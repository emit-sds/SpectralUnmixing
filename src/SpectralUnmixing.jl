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

function wl_index(wavelengths::Vector{Float64}, target)
    argmin(abs.(wavelengths .- target))
end

function scale_data(refl::Matrix{Float64}, wavelengths::Vector{Float64}, criteria::String)

    if criteria == "none"
        return refl
    elseif criteria == "brightness"
        bad_regions_wl = [[1300,1500],[1800,2000]]
        good_bands = convert(Vector{Bool}, ones(length(wavelengths)))
        for br in bad_regions_wl
            good_bands[wl_index(wavelengths, br[1]):wl_index(wavelengths, br[2])] .= false
        end
        if length(size(refl)) == 2
            norm = sqrt.(mean(refl[:,good_bands].^2, dims=2))
        else
            norm = sqrt.(mean(refl[good_bands].^2))
        end
    else
        try
            target_wl = parse(Float64,criteria)
            if length(size(refl)) == 2
                norm = refl[:,wl_index(wavelengths, target_wl)] ./ 0.5
            else
                norm = refl[wl_index(wavelengths, target_wl)] ./ 0.5
            end
        catch e
            throw(ArgumentError(string("normalization must be [none, brightness, or a specific wavelength].  Provided:", criteria)))
        end
    end

    return refl ./ norm
end

function get_sma_permutation(class_idx, num_endmembers::Vector{Int64}, combination_type::String, library_length::Int64)
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
                deleteat!(perm_class_idx[_p],1)

                if length(perm_class_idx[_p]) == 0
                    deleteat!(perm_class_idx,_p)
                end
                selector += 1
            end

        else
            perm = randperm(library_length)[1:num_endmembers[1]]
        end
    else
        perm = convert(Vector{Int64},1:library_length)
    end

    return perm
end

function results_from_mc(results::Matrix{Float64}, cost::Vector{Float64}, mode::String)

    if size(results)[1] == 1
        output_var = nothing
    else
        output_var = std(results, dims=1)[1,:]
    end

    if occursin("best", mode)
        best_idx = argmin(cost)
        output = results[best_idx,:]
    else
        output = mean(results, dims=1)[1,:]
    end
    
    return output, output_var

end

function unmix_pixel(library::SpectralLibrary, img_dat_input::Array{Float64}, unc_dat, class_idx, options, mode::String, n_mc::Int64, num_endmembers::Vector{Int64}, normalization::String, optimization::String, max_combinations::Int64, combination_type::String)
    

    if length(size(img_dat_input)) == 1
        img_dat = reshape(img_dat_input, 1, length(img_dat_input))
    else
        img_dat = img_dat_input
    end

    mc_comp_frac = zeros(n_mc, size(library.spectra)[1]+1)
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
            perm = get_sma_permutation(class_idx, num_endmembers, combination_type, size(library.spectra)[1])
            G = library.spectra[perm, library.good_bands]

            G = scale_data(G, library.wavelengths[library.good_bands], normalization)'

            x0 = dolsq(G, d', method=inverse_method)

            x0 = x0[:]
            res = nothing
            if occursin("bvls", optimization)
                res, cost = bvls(G, d[:], x0, zeros(size(x0)), ones(size(x0)), 1e-3, 100, 1, inverse_method)
            elseif occursin("ldsqp", optimization)
                res, cost = opt_solve(G, d[:], x0, zeros(length(x0)), ones(length(x0)) )
            elseif occursin("inverse", optimization)
                res = x0
                r = G * x0 - d[:]
                cost = dot(r,r)
            end
            mc_comp_frac[mc, perm] = res
            scores[mc] = cost

        elseif mode == "mesma" || mode == "mesma-best"
            solutions = []
            costs = zeros(size(options)[1]).+1e12
            if max_combinations != -1 && length(options) > max_combinations
                perm = randperm(length(options))[1:max_combinations]
            else
                perm = convert(Vector{Int64},1:length(options))
            end
            
            for (_comb, comb) in enumerate(options[perm])
                comb = [c for c in comb]
                #G = hcat(library.spectra[comb,:], ones(size(library.spectra[comb,:])[1],1))
                G = scale_data(library.spectra[comb, library.good_bands], library.wavelengths[library.good_bands], normalization)'
                x0 = dolsq(G, d')
                x0 = x0[:]
                ls = nothing
                if optimization == "bvls"
                    ls, lc = bvls(G, d[:], x0, zeros(size(x0)), ones(size(x0)), 1e-3, 10, 1, inverse_method)
                    costs[_comb] = lc
                elseif optimization == "ldsqp"
                    ls, lc = opt_solve(G, d[:], x0, 0, 1)
                    costs[_comb] = lc
                elseif optimization == "inverse"
                    ls = x0
                    r = G * x0 - d[:]
                    costs[_comb] = dot(r,r)
                end
                push!(solutions,ls)

            end
            best = argmin(costs)
            scores[mc] = best

            mc_comp_frac[mc, [ind for ind in options[perm][best]]] = solutions[best]

        else
            error("Invalid mode provided")
        end
    end

    # Calculate the sum of values (inverse of shade), and then normalize
    mc_comp_frac[mc_comp_frac .< 0] .= 0
    mc_comp_frac[:,end] = sum(mc_comp_frac,dims=2)
    mc_comp_frac[:,1:end-1] = mc_comp_frac[:,1:end-1] ./ mc_comp_frac[:,end]

    # Aggregate results from per-library to per-unique-class
    mixture_results = zeros(size(mc_comp_frac)[1], length(library.class_valid_keys) + 1)
    for _i in 1:size(mc_comp_frac)[1]
        for (_class, cl) in enumerate(library.class_valid_keys)
            mixture_results[_i, _class] = sum(mc_comp_frac[_i,1:end-1][cl .== library.classes])
        end
        mixture_results[_i, end] = mc_comp_frac[_i,end]
    end

    output_mixture, output_mixture_var = results_from_mc(mixture_results, scores, mode)
    output_comp_frac, output_comp_frac_var = results_from_mc(mc_comp_frac, scores, mode)

    return output_mixture, output_mixture_var, output_comp_frac, output_comp_frac_var

end

function prepare_combinations(library::SpectralLibrary, combination_type::String)
    class_idx = []
    if combination_type == "class-even"
        for uc in library.class_valid_keys
            push!(class_idx, (1:size(library.classes)[1])[library.classes .== uc])
        end
    end
    return class_idx
end

function prepare_options(library::SpectralLibrary, combination_type::String, num_endmembers::Vector{Int64}, class_idx)

    # Prepare combinations if relevant
    options = []
    if combination_type == "class-even"
        options = collect(Iterators.product(class_idx...))[:]
    elseif combination_type == "all"
        for num in num_endmembers
            combo = [c for c in combinations(1:length(library.classes), num)]
            push!(options,combo...)
        end
    else
        error("Invalid combiation string")
    end

    return options
end



function unmix_line(line::Int64, reflectance_file::String, mode::String, refl_nodata::Float64,
                    refl_scale::Float64, normalization::String, library::SpectralLibrary,
                    reflectance_uncertainty_file::String = "", n_mc::Int64 = 1,
                    combination_type::String = "all", num_endmembers::Vector{Int64} = [2,3],
                    max_combinations::Int64 = -1, optimization="bvls")

    Random.seed!(13)
    
    img_dat, unc_dat, good_data = load_line(reflectance_file, reflectance_uncertainty_file, line, library.good_bands, refl_nodata)
    if isnothing(img_dat)
        return line, nothing, good_data, nothing, nothing
    end

    mixture_results = fill(-9999.0, sum(good_data), size(library.class_valid_keys)[1] + 1)
    complete_fractions = zeros(size(img_dat)[1], size(library.spectra)[1] + 1)
    
        
    if n_mc > 1
        mixture_results_std = fill(-9999.0, sum(good_data), size(library.class_valid_keys)[1] + 1)
        complete_fractions_std = zeros(size(img_dat)[1], size(library.spectra)[1] + 1)
    else
        mixture_results_std = nothing
        complete_fractions_std = nothing
    end

    scale_data(img_dat, library.wavelengths[library.good_bands], normalization)
    img_dat = img_dat ./ refl_scale

    class_idx = prepare_combinations(library, combination_type)
    # Prepare combinations if relevant
    if mode == "mesma" || mode == "mesma-best"
        options = prepare_options(library, combination_type, num_endmembers, class_idx)
    end

    start_time = time()

    # Solve for each pixel
    for _i in 1:size(img_dat)[1] # Pixel loop

        lid = img_dat[_i:_i,:]
        if isnothing(unc_dat) 
            lud = nothing 
        else 
            lud = unc_dat[_i:_i,:] 
        end

        loc_mixture_res, loc_mixture_var, loc_cf_res, loc_cf_var  = unmix_pixel(library, lid, 
            lud, class_idx, options, mode, n_mc, num_endmembers, normalization, optimization, 
            max_combinations, combination_type)

        complete_fractions[_i,:] = loc_cf_res
        mixture_results[_i,:] = loc_mixture_res

        if n_mc > 1
            complete_fractions_std[_i,:] = loc_cf_var
            mixture_results_std[_i,:] = loc_mixture_var
        end

    end
    elapsed_time = time() - start_time 
    if line != 1 
        @info string("Line " , line , " run in " , round(elapsed_time, digits=4) , " seconds")
    end
    return line, mixture_results, good_data, mixture_results_std, complete_fractions

end

function unmix_and_write_line(line::Int64, reflectance_file::String, mode::String, refl_nodata::Float64,
                    refl_scale::Float64, normalization::String, library::SpectralLibrary, output_files::Vector{String},
                    write_complete_fractions::Bool,
                    reflectance_uncertainty_file::String = "", n_mc::Int64 = 1,
                    combination_type::String = "all", num_endmembers::Vector{Int64} = [2,3],
                    max_combinations::Int64 = -1, optimization="bvls")
    
    line_results = unmix_line(line, reflectance_file, mode, refl_nodata, refl_scale, normalization, library, 
        reflectance_uncertainty_file, n_mc, combination_type, num_endmembers,
        max_combinations, optimization)

    write_line_results(output_files, line_results, n_mc, write_complete_fractions) 
    
end

function class_assign_fractions(complete_fractions, library::SpectralLibrary)
    # Aggregate results from per-library to per-unique-class
    if length(size(complete_fractions)) == 1
        cf = reshape(complete_fractions, (1, length(complete_fractions)))
    else
        cf = complete_fractions
    end
    mixture_results = zeros(size(cf)[1], length(library.class_valid_keys))

    for _i in 1:size(cf)[1]
        for (_class, cl) in enumerate(library.class_valid_keys)
            mixture_results[_i, _class] = sum(complete_fractions[_i,1:end-1][cl .== library.classes])
        end
    end
    
    if length(size(complete_fractions)) == 1
        return mixture_results[1,:]
    else
        return mixture_results
    end

end


function simulate_pixel(library::SpectralLibrary, max_components::Int64, combination_type::String, seed::Int64)

    Random.seed!(seed)
	
    output_mixture = zeros(size(library.spectra)[2])
    output_mixture[:] .= NaN

    class_idx = prepare_combinations(library, combination_type)
    perm = get_sma_permutation(class_idx, [max_components], combination_type, size(library.spectra)[1])

    G = library.spectra[perm,:]

    distribution = rand(max_components)
    distribution = distribution ./ sum(distribution)

    output_distribution = zeros(size(library.spectra)[1])
    output_distribution[perm] = distribution

    output_distribution_classes = zeros(size(library.class_valid_keys))
    for (_class, cl) in enumerate(library.class_valid_keys)
        output_distribution_classes[_class] = sum(output_distribution[cl .== library.classes])
    end

    simulated_rfl = G' * distribution
    return simulated_rfl, output_distribution, output_distribution_classes

end


end
