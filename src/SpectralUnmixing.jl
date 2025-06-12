#  Copyright 2022 California Institute of Technology
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
# Author: Philip G. Brodrick, philip.g.brodrick@jpl.nasa.gov
module SpectralUnmixing

using ArchGDAL
using EllipsisNotation
using DelimitedFiles
using Logging
using Statistics
using Distributed
using Printf
using LinearAlgebra
using Combinatorics
using Random

include("Datasets.jl")
include("EndmemberLibrary.jl")
include("Solvers.jl")
include("Plotting.jl")

# Endmember Library Functions
export SpectralLibrary, load_data!, filter_by_class!, read_envi_wavelengths, interpolate_library_to_new_wavelengths!, remove_wavelength_region_inplace!, scale_library!
export reduce_endmembers_nmf!, reduce_endmembers_kmeans!, reduce_endmembers_pca!, brightness_normalize!, split_library
export prepare_combinations, prepare_options, scale_data

# Plotting Functions
export plot_mean_endmembers, plot_endmembers, plot_endmembers_individually

# Dataset functions
export initiate_output_datasets, set_band_names, write_results

# Unmixing and simulation functions
export unmix_line, unmix_pixel, simulate_pixel, unmix_and_write_line

"""
    wl_index(wavelengths::Vector{Float64}, target::Float64)

Return the index of the wavelength in a vector that is closest to the specified target
wavelength.
"""
function wl_index(wavelengths::Vector{Float64}, target)
    argmin(abs.(wavelengths .- target))
end

"""
    scale_data(refl::Matrix{Float64}, wavelengths::Vector{Float64}, criteria::String,
        bad_regions_wl=[[1300, 1500], [1800, 2000]])

Scale the reflectance data based on specified criteria.

# Arguments
- `refl::Matrix{Float64}`: Reflectance data either size (n,) or (m,n), where n is the
number of bands and m the number of pixels.
- `wavelengths::Vector{Float64}`: A vector of wavelengths of size (n,).
- `criteria::String`: Normalization options:
  - `"none"`: No scaling will be applied, and the original reflectance data will be
  returned.
  - `"brightness"`: The data will be normalized based on brightness, excluding specified
  bad wavelength regions [1300, 1500] and [1800, 2000].
  - A specific wavelength (as a string): The data will be normalized using the reflectance
  value at this wavelength.
- `bad_regions_wl`=[[1300, 1500], [1800, 2000]]: List of wavelength regions (in nm) to
ignore.

# Returns
- A matrix of scaled reflectance data with same dimensions as `refl`.
"""
function scale_data(refl::Matrix{Float64}, wavelengths::Vector{Float64}, criteria::String,
    bad_regions_wl=[[1300, 1500], [1800, 2000]])

    if criteria == "none"
        return refl
    elseif criteria == "brightness"
        good_bands = convert(Vector{Bool}, ones(length(wavelengths)))
        for br in bad_regions_wl
            good_bands[wl_index(wavelengths, br[1]):wl_index(wavelengths, br[2])] .= false
        end
        if length(size(refl)) == 2
            norm = sqrt.(mean(refl[:, good_bands] .^ 2, dims=2))
        else
            norm = sqrt.(mean(refl[good_bands] .^ 2))
        end
    else
        try
            target_wl = parse(Float64, criteria)
            if length(size(refl)) == 2
                norm = refl[:, wl_index(wavelengths, target_wl)] ./ 0.5
            else
                norm = refl[wl_index(wavelengths, target_wl)] ./ 0.5
            end
        catch e
            throw(ArgumentError(string("normalization must be [none, brightness, or a
            specific wavelength].  Provided:", criteria)))
        end
    end

    return refl ./ norm
end

"""
    get_sma_permutation(class_idx, num_endmembers::Vector{Int64},
                        combination_type::String, library_length::Int64)

Generate a permutation of endmember indices for SMA.

# Arguments
- `class_idx`: A collection (such as a vector of vectors) of size n where each element
is a set of indices corresponding to endmembers belonging to one class.
- `num_endmembers::Vector{Int64}`: Desired number of endmembers to select. The first
element specifies the number of endmembers to permute; if it is `-1`, all endmembers in
library will be included.
- `combination_type::String`: Permutation type options:
  - `"class-even"`: Randomly select one endmember from each class until the required number
  of endmembers is selected.
  - Any other types: Randomly selects endmembers from the entire library.
- `library_length::Int64`: Size of the entire endmember library.

# Returns
- A vector of integers representing the permuted endmember indices.
"""
function get_sma_permutation(class_idx, num_endmembers::Vector{Int64},
    combination_type::String, library_length::Int64)

    if num_endmembers[1] != -1
        if combination_type == "class-even"

            perm_class_idx = []
            for class_subset in class_idx
                push!(perm_class_idx, Random.shuffle(class_subset))
            end

            perm = []
            selector = 1
            while selector <= num_endmembers[1]
                _p = mod(selector, length(perm_class_idx)) + 1
                push!(perm, perm_class_idx[_p][1])
                deleteat!(perm_class_idx[_p], 1)

                if length(perm_class_idx[_p]) == 0
                    deleteat!(perm_class_idx, _p)
                end
                selector += 1
            end

        else
            perm = randperm(library_length)[1:num_endmembers[1]]
        end
    else
        perm = convert(Vector{Int64}, 1:library_length)
    end

    return perm
end

"""
    results_from_mc(results::Matrix{Float64}, cost::Vector{Float64}, mode::String)

Process Monte Carlo simulations and return variance and result based on `mode`.

# Arguments
- `results::Matrix{Float64}`: Outcomes of MC simulations, where each row represents a
simulation and each column corresponds to a different variable.
- `cost::Vector{Float64}`: A vector of Float64 values representing the cost associated
with each simulation, used to determine the best result when applicable.
- `mode::String`: How to process the results. It can be:
  - Contain `"best"`: Returns the result corresponding to the simulation with the lowest
  cost.
  - Any other value: Computes and returns the mean of the results across all simulations.

# Returns
- A tuple containing:
  - `output::Vector{Float64}`: MC simulation results
  - `output_var::Vector{Float64}`: Standard deviation of the results for each variable, or
  `nothing` if there is only one simulation.
"""
function results_from_mc(results::Matrix{Float64}, cost::Vector{Float64}, mode::String)

    if size(results)[1] == 1
        output_var = nothing
    else
        output_var = std(results, dims=1)[1, :]
    end

    if occursin("best", mode)
        best_idx = argmin(cost)
        output = results[best_idx, :]
    else
        output = mean(results, dims=1)[1, :]
    end

    return output, output_var
end

"""
    unmix_pixel(library::SpectralLibrary, img_dat_input::Array{Float64}, unc_dat,
                class_idx, options, mode::String, n_mc::Int64,
                num_endmembers::Vector{Int64}, normalization::String, optimization::String,
                max_combinations::Int64, combination_type::String)

Unmix a pixel's spectral data using a given spectral library with options for SMA or MESMA
and various Monte Carlo and optimization approaches.

# Arguments
- `library::SpectralLibrary`: Endmember library for unmixing.
- `img_dat_input::Array{Float64}`: Spectral data of the pixel.
- `unc_dat`: An optional array of uncertainty in the pixel spectral data.
- `class_idx`: A collection of indices indicating class memberships for different
endmembers in `library.spectra`. See also [`prepare_combinations`](@ref),
[`prepare_options`](@ref)
- `options`: A collection of potential endmember combinations for the unmixing process. See
also [`prepare_options`](@ref).
- `mode::String`: Determines the unmixing approach:
  - `"sma"`: Select endmembers randomly for MC unmixing, returns mean fractions.
  - `"sma-best"`: SMA and output best (lowest cost) MC fraction.
  - `"mesma"`: Evaluate different combinations of endmembers representing each class.
  - `"mesma-best"`: MESMA and output best (lowest cost) MC fraction.
- `n_mc::Int64`: Number of MC iterations.
- `num_endmembers::Vector{Int64}`: The number of endmembers to consider.
- `normalization::String`: The normalization method to apply to the spectral data. See
[`scale_data`](@ref).
- `optimization::String`: The optimization approach for unmixing: (e.g., `"bvls"`,
`"ldsqp"`, or `"inverse"`). Can optionally contain `"pinv"` or `"qr"` to specify the
inverse method, if applicable.
- `max_combinations::Int64`: Maximum number of combinations to consider.
- `combination_type::String`: The type of combination (e.g., `"class-even"`). See also
[`prepare_combinations`](@ref), [`prepare_options`](@ref), [`get_sma_permutation`](@ref).

# Returns
- A tuple containing:
  - `output_mixture::Vector{Float64}`: The estimated fraction of each class in the
  `library`, appended with the brightness.
  - `output_mixture_var::Vector{Float64}`: The variance of each class in the `library`,
  appended with the brightness variance.
  - `output_comp_frac::Vector{Float64}`: The estimated fraction of each endmember in the
  `library`, appended with the brightness.
  - `output_comp_frac_var::Vector{Float64}`: The variance of each endmember in the
  `library`, appended with the brightness variance.
"""
function unmix_pixel(library::SpectralLibrary, img_dat_input::Array{Float64}, unc_dat,
    class_idx, options, mode::String, n_mc::Int64, num_endmembers::Vector{Int64},
    normalization::String, optimization::String, max_combinations::Int64,
    combination_type::String)

    if length(size(img_dat_input)) == 1
        img_dat = reshape(img_dat_input, 1, length(img_dat_input))
    else
        img_dat = img_dat_input
    end

    mc_comp_frac = zeros(n_mc, size(library.spectra)[1] + 1)
    scores = zeros(n_mc)
    for mc in 1:n_mc #monte carlo loop
        Random.seed!(mc)

        d = img_dat
        if isnothing(unc_dat) == false
            d += (rand(size(d)...) .* 2 .- 1) .* unc_dat
        end

        if occursin("pinv", optimization)
            inverse_method = "pinv"
        elseif occursin("qr", optimization)
            inverse_method = "qr"
        else
            inverse_method = "default"
        end

        if mode == "sma" || mode == "sma-best"
            perm = get_sma_permutation(
                class_idx, num_endmembers, combination_type, size(library.spectra)[1]
            )
            G = library.spectra[perm, library.good_bands]

            G = scale_data(G, library.wavelengths[library.good_bands], normalization)'

            x0 = dolsq(G, d', method=inverse_method)

            x0 = x0[:]
            res = nothing
            if occursin("bvls", optimization)
                res, cost = bvls(
                    G, d[:], x0, zeros(size(x0)), ones(size(x0)), 1e-3, 100, 1,
                    inverse_method
                )
            elseif occursin("ldsqp", optimization)
                res, cost = opt_solve(G, d[:], x0, zeros(length(x0)), ones(length(x0)))
            elseif occursin("inverse", optimization)
                res = x0
                r = G * x0 - d[:]
                cost = dot(r, r)
            end
            mc_comp_frac[mc, perm] = res
            scores[mc] = cost

        elseif mode == "mesma" || mode == "mesma-best"
            solutions = []
            costs = zeros(size(options)[1]) .+ 1e12
            if max_combinations != -1 && length(options) > max_combinations
                perm = randperm(length(options))[1:max_combinations]
            else
                perm = convert(Vector{Int64}, 1:length(options))
            end

            for (_comb, comb) in enumerate(options[perm])
                comb = [c for c in comb]
                #G = hcat(library.spectra[comb,:], ones(size(library.spectra[comb,:])[1],1))
                G = scale_data(
                    library.spectra[comb, library.good_bands],
                    library.wavelengths[library.good_bands], normalization
                )'
                x0 = dolsq(G, d')
                x0 = x0[:]
                ls = nothing
                if optimization == "bvls"
                    ls, lc = bvls(
                        G, d[:], x0, zeros(size(x0)), ones(size(x0)), 1e-3, 10, 1,
                        inverse_method
                    )
                    costs[_comb] = lc
                elseif optimization == "ldsqp"
                    ls, lc = opt_solve(G, d[:], x0, 0, 1)
                    costs[_comb] = lc
                elseif optimization == "inverse"
                    ls = x0
                    r = G * x0 - d[:]
                    costs[_comb] = dot(r, r)
                end
                push!(solutions, ls)

            end
            best = argmin(costs)
            scores[mc] = best

            mc_comp_frac[mc, [ind for ind in options[perm][best]]] = solutions[best]

        else
            error("Invalid mode provided")
        end
    end

    # Calculate the sum of values (inverse of shade), and then normalize
    mc_comp_frac[mc_comp_frac.<0] .= 0
    mc_comp_frac[:, end] = sum(mc_comp_frac, dims=2)
    mc_comp_frac[:, 1:end-1] = mc_comp_frac[:, 1:end-1] ./ mc_comp_frac[:, end]

    # Aggregate results from per-library to per-unique-class
    mixture_results = zeros(size(mc_comp_frac)[1], length(library.class_valid_keys) + 1)
    for _i in 1:size(mc_comp_frac)[1]
        for (_class, cl) in enumerate(library.class_valid_keys)
            mixture_results[_i, _class] =
                sum(mc_comp_frac[_i, 1:end-1][cl.==library.classes])
        end
        mixture_results[_i, end] = mc_comp_frac[_i, end]
    end

    output_mixture, output_mixture_var = results_from_mc(mixture_results, scores, mode)
    output_comp_frac, output_comp_frac_var = results_from_mc(mc_comp_frac, scores, mode)

    return output_mixture, output_mixture_var, output_comp_frac, output_comp_frac_var
end

"""
    prepare_combinations(library::SpectralLibrary, combination_type::String)

Prepare the set of class indices for the specified combination type.

- Currently only supports `"class-even"`: returns the list of indices grouped by class.
"""
function prepare_combinations(library::SpectralLibrary, combination_type::String)
    class_idx = []
    if combination_type == "class-even"
        for uc in library.class_valid_keys
            push!(class_idx, (1:size(library.classes)[1])[library.classes.==uc])
        end
    end
    return class_idx
end

"""
    prepare_options(library::SpectralLibrary, combination_type::String,
                     num_endmembers::Vector{Int64}, class_idx)

Prepare combinations of endmembers based on the specified combination type.

- Combination type options:
    - `"class-even"``: Generate all combinations where one endmember is selected from each class.
    - `"all"`: Generate all possible combinations of `num_endmembers` spectra.
"""
function prepare_options(library::SpectralLibrary, combination_type::String,
    num_endmembers::Vector{Int64}, class_idx)

    # Prepare combinations if relevant
    options = []
    if combination_type == "class-even"
        options = collect(Iterators.product(class_idx...))[:]
    elseif combination_type == "all"
        for num in num_endmembers
            combo = [c for c in combinations(1:length(library.classes), num)]
            push!(options, combo...)
        end
    else
        error("Invalid combiation string")
    end

    return options
end

"""
    unmix_line(line::Int64, reflectance_file::String, mode::String,
                refl_nodata::Float64, refl_scale::Float64,
                normalization::String, library::SpectralLibrary,
                reflectance_uncertainty_file::String="", n_mc::Int64=1,
                combination_type::String="all",
                num_endmembers::Vector{Int64}=[2, 3],
                max_combinations::Int64=-1, optimization="bvls")

Unmix a specific line (row of pixels) of reflectance data. See also [`load_line`](@ref),
[`unmix_pixel`](@ref).

# Arguments
- `line::Int64`: Index of the line to unmix.
- `reflectance_file::String`: Path to the reflectance data.
- `mode::String`: The mode of unmixing to be used (e.g., "sma", "mesma-best"). See
[`unmix_pixel`](@ref)).
- `refl_nodata::Float64`: The no data value in the reflectance file.
- `refl_scale::Float64`: Scaling factor for the reflectance data.
- `normalization::String`: The normalization method to apply to the reflectance data. See
[`scale_data`](@ref).
- `library::SpectralLibrary`: Spectral library object containing endmembers for unmixing.
- `reflectance_uncertainty_file::String`: Optional path to the reflectance uncertainty file.
- `n_mc::Int64=1`: Number of Monte Carlo iterations to perform.
- `combination_type::String="all"`: The type of endmember combinations to prepare for
unmixing. See also [`prepare_combinations`](@ref), [`prepare_options`](@ref).
- `num_endmembers::Vector{Int64}=[2, 3]`: The number of endmembers to consider in
combinations.
- `max_combinations::Int64=-1`: Maximum number of combinations to consider, defaults to no
limit.
- `optimization::String`: The optimization method to use (default is "bvls"). See also
[`unmix_pixel`](@ref).

# Returns
- A tuple containing (see also [`unmix_pixel`](@ref)):
  - The line index.
  - Estimated class fractions for each pixel in the line.
  - A boolean mask of pixels with valid data.
  - Standard deviations of the mixture results (if applicable).
  - Complete fractions for each endmember in the unmixing of each pixel.

# Notes
- Initializes a random seed for reproducibility.
- Logs the execution time for unmixing a line.
- Returns `nothing` if input data is invalid or missing.
"""
function unmix_line(line::Int64, reflectance_file::String, mode::String,
    refl_nodata::Float64, refl_scale::Float64, normalization::String,
    library::SpectralLibrary, reflectance_uncertainty_file::String="", n_mc::Int64=1,
    combination_type::String="all", num_endmembers::Vector{Int64}=[2, 3],
    max_combinations::Int64=-1, optimization="bvls")

    Random.seed!(13)

    img_dat, unc_dat, good_data = load_line(
        reflectance_file, reflectance_uncertainty_file, line, library.good_bands,
        refl_nodata
    )
    if isnothing(img_dat)
        return line, nothing, good_data, nothing, nothing
    end

    mixture_results = fill(-9999.0, sum(good_data), size(library.class_valid_keys)[1] + 1)
    complete_fractions = zeros(size(img_dat)[1], size(library.spectra)[1] + 1)


    if n_mc > 1
        mixture_results_std =
            fill(-9999.0, sum(good_data), size(library.class_valid_keys)[1] + 1)
        complete_fractions_std = zeros(size(img_dat)[1], size(library.spectra)[1] + 1)
    else
        mixture_results_std = nothing
        complete_fractions_std = nothing
    end

    scale_data(img_dat, library.wavelengths[library.good_bands], normalization)
    img_dat = img_dat ./ refl_scale

    class_idx = prepare_combinations(library, combination_type)
    options = []
    # Prepare combinations if relevant
    if mode == "mesma" || mode == "mesma-best"
        options = prepare_options(library, combination_type, num_endmembers, class_idx)
    end

    start_time = time()

    # Solve for each pixel
    for _i in 1:size(img_dat)[1] # Pixel loop

        lid = img_dat[_i:_i, :]
        if isnothing(unc_dat)
            lud = nothing
        else
            lud = unc_dat[_i:_i, :]
        end

        loc_mixture_res, loc_mixture_var, loc_cf_res, loc_cf_var =
            unmix_pixel(library, lid,
                lud, class_idx, options, mode, n_mc, num_endmembers, normalization,
                optimization, max_combinations, combination_type
            )

        complete_fractions[_i, :] = loc_cf_res
        mixture_results[_i, :] = loc_mixture_res

        if n_mc > 1
            complete_fractions_std[_i, :] = loc_cf_var
            mixture_results_std[_i, :] = loc_mixture_var
        end

    end
    elapsed_time = time() - start_time
    if line != 1
        @info string("Line ", line, " run in ", round(elapsed_time, digits=4), " seconds")
    end
    return line, mixture_results, good_data, mixture_results_std, complete_fractions
end

"""
    unmix_and_write_line(line::Int64, reflectance_file::String, mode::String,
                          refl_nodata::Float64, refl_scale::Float64,
                          normalization::String, library::SpectralLibrary,
                          output_files::Vector{String}, write_complete_fractions::Bool,
                          reflectance_uncertainty_file::String="", n_mc::Int64=1,
                          combination_type::String="all",
                          num_endmembers::Vector{Int64}=[2, 3],
                          max_combinations::Int64=-1, optimization="bvls")

Unmix the specified `line` of `reflectance_file` and write the results to `output_files`.
Workhorse function for `unmix.jl`.

- See [`unmix_line`(@ref)] and [`write_line_results`](@ref) for arguments.
"""
function unmix_and_write_line(line::Int64, reflectance_file::String, mode::String,
    refl_nodata::Float64, refl_scale::Float64, normalization::String,
    library::SpectralLibrary, output_files::Vector{String}, write_complete_fractions::Bool,
    reflectance_uncertainty_file::String="", n_mc::Int64=1, combination_type::String="all",
    num_endmembers::Vector{Int64}=[2, 3], max_combinations::Int64=-1, optimization="bvls")

    line_results = unmix_line(
        line, reflectance_file, mode, refl_nodata, refl_scale, normalization, library,
        reflectance_uncertainty_file, n_mc, combination_type, num_endmembers,
        max_combinations, optimization
    )

    write_line_results(output_files, line_results, n_mc, write_complete_fractions)
end

"""
    simulate_pixel(library::SpectralLibrary, max_components::Int64,
                    combination_type::String, seed::Int64)

Simulate reflectance data for a pixel from a random normalized distribution of
endmembers in `library`.

# Arguments
- `library::SpectralLibrary`: Spectral library containing endmembers and class
information.
- `max_components::Int64`: Maximum number of endmembers to use in the simulation.
- `combination_type::String`: The type of combinations to consider when selecting
endmembers (e.g., "all" or "class-even").
- `seed::Int64`: Seed for the random number generator.

# Returns
- A tuple containing:
  - `simulated_rfl`: Simulated reflectance spectrum of the pixel.
  - `output_distribution`: A vector of endmember contributions for the pixel.
  - `output_distribution_classes`: A vector of contributions of each unique class.
"""
function simulate_pixel(library::SpectralLibrary, max_components::Int64,
    combination_type::String, seed::Int64)

    Random.seed!(seed)

    output_mixture = zeros(size(library.spectra)[2])
    output_mixture[:] .= NaN

    class_idx = prepare_combinations(library, combination_type)
    perm = get_sma_permutation(
        class_idx, [max_components], combination_type, size(library.spectra)[1]
    )

    G = library.spectra[perm, :]

    distribution = rand(max_components)
    distribution = distribution ./ sum(distribution)

    output_distribution = zeros(size(library.spectra)[1])
    output_distribution[perm] = distribution

    output_distribution_classes = zeros(size(library.class_valid_keys))
    for (_class, cl) in enumerate(library.class_valid_keys)
        output_distribution_classes[_class] = sum(output_distribution[cl.==library.classes])
    end

    simulated_rfl = G' * distribution
    return simulated_rfl, output_distribution, output_distribution_classes
end

end
