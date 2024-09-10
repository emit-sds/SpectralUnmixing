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

using DataFrames
using CSV
using Interpolations
using Logging
using Statistics
using Plots
using ModularIndices
using NMF
using MultivariateStats
using Clustering

function nanargmax(input::Array)
    x = copy(input)
    x[isnan.(x)] .= -Inf;
    return argmax(x);
end

function nanargmin(input::Array)
    x = copy(input)
    x[isnan.(x)] .= Inf;
    return argmin(x);
end

function read_envi_wavelengths(filename::String, nm::Bool=true)
    header_name = splitext(filename)[1] * ".hdr"
    header = readlines(header_name)
    found = false
    for line in header
        if occursin("wavelength = {", line) || occursin("wavelength= {", line) ||  occursin("wavelength={", line)
            header = line
            found = true
            break
        end
    end

    if !found
        @error "No wavelength found in " * header_name
        return nothing
    end

    wavelengths = [parse(Float64, strip(x)) for x in split(split(split(header, "{")[2], "}")[1],",")]

    if nm && all(wavelengths .< 25)
        wavelengths = wavelengths .* 1000
        @info "Converting wavelengths read from $filename to nm from microns."
    end
        
    return wavelengths
end

function get_good_bands_mask(wavelengths::Array{Float64}, wavelength_pairs)
    good_bands = ones(Bool, length(wavelengths))
    
    for wvp in wavelength_pairs
        wavelength_diff = wavelengths .- wvp[1]
        wavelength_diff[wavelength_diff .< 0] .= maximum(filter(!isnan, wavelength_diff))
        lower_index = nanargmin(wavelength_diff)

        wavelength_diff = wvp[2] .- wavelengths
        wavelength_diff[wavelength_diff .< 0] .= maximum(filter(!isnan, wavelength_diff))
        upper_index = nanargmin(wavelength_diff)
        good_bands[lower_index:upper_index] .= false
    end
    return good_bands
end

mutable struct SpectralLibrary
    file_name::String
    class_header_name::String
    spectral_starting_column::Int64
    truncate_end_columns::Int64
    class_valid_keys
    scale_factor::Float64
    wavelength_regions_ignore
    SpectralLibrary(file_name::String, class_header_name::String, spectral_starting_column::Int64 = 2, truncate_end_columns::Int64 = 0, class_valid_keys = nothing, 
                    scale_factor = 1.0, wavelength_regions_ignore= [0,440,1310,1490,1770,2050,2440,2880]) = 
                    new(file_name, class_header_name, spectral_starting_column, truncate_end_columns, class_valid_keys, scale_factor, wavelength_regions_ignore)

    spectra
    classes
    good_bands
    wavelengths
end

function load_data!(library::SpectralLibrary)

    df = DataFrame(CSV.File(library.file_name))

    paired_sets = []
    for i in 1:2:size(library.wavelength_regions_ignore)[1]
        push!(paired_sets, [library.wavelength_regions_ignore[i],library.wavelength_regions_ignore[i+1]])
    end
    library.wavelength_regions_ignore = paired_sets
    @info "Ignoring wavelength regions: " * string(library.wavelength_regions_ignore)

    try
        library.spectra = Matrix(df[:,library.spectral_starting_column:end-library.truncate_end_columns])
    catch e
        throw(ArgumentError("Could not read spectra from endmember library.  Try adjusting spectral_starting_column or truncate_end_columns, or reworking the library."))
    end

    try
        library.classes = convert(Array{AbstractString}, df[!,library.class_header_name])
    catch e
        throw(ArgumentError("Could not read classes from endmember library using class key "*library.class_header_name * 
                            ".  Try adjusting class_header_name, or reworking the library."))
    end

    try
        library.wavelengths = parse.(Float64, names(df[:,library.spectral_starting_column:end-library.truncate_end_columns]))
    catch e
        throw(ArgumentError("Could not read wavelengths from endmember library.  Try adjusting class_header_name, or reworking the library."))
    end


    wl_order = sortperm(library.wavelengths)
    library.spectra = library.spectra[:,wl_order]
    library.wavelengths = library.wavelengths[wl_order]

    if isnothing(library.class_valid_keys)
        library.class_valid_keys = unique(library.classes)
    end

    good_bands = get_good_bands_mask(library.wavelengths, library.wavelength_regions_ignore)
    library.good_bands = good_bands

    ignore_regions_nm = false
    for wri in library.wavelength_regions_ignore
        if wri[1] > 25 || wri[2] > 25
            ignore_regions_nm = true
        end
    end

    if ignore_regions_nm && all(library.wavelengths .< 25)
        library.wavelengths = library.wavelengths .* 1000
        @info "Unit mismatch between library wavelengths and wavelength regions to ignore detected.  Converting library wavelengths to nm from microns"
    end

    return library
end


function save_data(library::SpectralLibrary, output_filename::String; class_label_name::String = "Label")

    df = DataFrame()
    insertcols!(df, 1, class_label_name => library.classes)
    for (_wv, wv) in enumerate(library.wavelengths)
        insertcols!(df, 1 + _wv, string(wv) => library.spectra[:,_wv])
    end
    CSV.write(output_filename, df)

end

function filter_by_class!(library::SpectralLibrary)
    if isnothing(library.class_valid_keys)
        @info "No class valid keys provided, no filtering occuring"
        return
    end

    valid_classes = zeros(Bool, size(library.spectra)[1])
    for cla in library.class_valid_keys
        valid_classes[library.classes .== cla] .= true
    end

    library.spectra = library.spectra[valid_classes,:]
    library.classes = library.classes[valid_classes]
end

function remove_wavelength_region_inplace!(library::SpectralLibrary, set_as_nans::Bool=false)
    if set_as_nans
        library.spectra[:,.!library.good_bands] .= NaN
        library.wavelengths[.!library.good_bands] .= NaN
    else
        library.spectra = library.spectra[:, library.good_bands]
        library.wavelengths = library.wavelengths[library.good_bands]
        library.good_bands = library.good_bands[library.good_bands]
    end
end

function interpolate_library_to_new_wavelengths!(library::SpectralLibrary, new_wavelengths::Array{Float64})
    old_spectra = copy(library.spectra)

    library.spectra = zeros((size(library.spectra)[1], length(new_wavelengths)))
    for _s in 1:size(old_spectra)[1]
        fit = LinearInterpolation(library.wavelengths, old_spectra[_s,:], extrapolation_bc=Flat());
        library.spectra[_s,:] = fit(new_wavelengths)
    end
    library.wavelengths = new_wavelengths

    good_bands = get_good_bands_mask(library.wavelengths, library.wavelength_regions_ignore)
    library.good_bands = good_bands
end

function scale_library!(library::SpectralLibrary, scaling_factor=nothing)
    if isnothing(scaling_factor)
        library.spectra =  library.spectra ./ library.scale_factor
    else
        library.spectra /= scaling_factor
    end
end

function brightness_normalize!(library::SpectralLibrary)
    library.spectra = library.spectra ./ sqrt.(mean(library.spectra[:,library.good_bands].^2, dims=2))
end

function split_library(library::SpectralLibrary, split_fraction::Float64)
    perm = randperm(size(library.spectra)[1])

    split_1 = perm[1:Int(round(split_fraction * length(perm)))]
    split_2 = perm[Int(round(split_fraction * length(perm))):end]

    output_library_1 = deepcopy(library)
    output_library_2 = deepcopy(library)

    output_library_1.spectra = output_library_1.spectra[split_1,:]
    output_library_2.spectra = output_library_2.spectra[split_2,:]

    output_library_1.classes = output_library_1.classes[split_1]
    output_library_2.classes = output_library_2.classes[split_2]

    return output_library_1, output_library_2

end

function reduce_endmembers_nmf!(library::SpectralLibrary, max_endmembers_per_class::Int64)

    reduced_library = []
    new_classes = copy(library.classes)[1:max_endmembers_per_class*size(library.class_valid_keys)[1]]
    for (_c, cla) in enumerate(library.class_valid_keys)
        library_subset = library.spectra[library.classes .== cla, library.good_bands]

        out_spectra = zeros(max_endmembers_per_class, size(library.spectra)[2])
        out_spectra[:,:] .= NaN

        r = nnmf(library_subset, max_endmembers_per_class; maxiter=500, tol=1.0e-2)
        out_spectra[:,library.good_bands] = r.H

        push!(reduced_library,out_spectra)
        new_classes[(_c-1) * max_endmembers_per_class + 1:_c * max_endmembers_per_class] .= cla
    end
    library.spectra = cat(reduced_library...,dims=1)
    library.classes = new_classes
end

function reduce_endmembers_pca!(library::SpectralLibrary, max_endmembers_per_class::Int64)

    reduced_library = []
    new_classes = copy(library.classes)[1:max_endmembers_per_class*size(library.class_valid_keys)[1]]
    for (_c, cla) in enumerate(library.class_valid_keys)
        library_subset = library.spectra[library.classes .== cla, library.good_bands]

        out_spectra = zeros(max_endmembers_per_class, size(library.spectra)[2])
        out_spectra[:,:] .= NaN

        model = fit(PCA, library_subset', maxoutdim=max_endmembers_per_class, pratio=1)
        out_spectra[:,library.good_bands] = projection(model)'

        push!(reduced_library,out_spectra)
        new_classes[(_c-1) * max_endmembers_per_class + 1:_c * max_endmembers_per_class] .= cla
    end
    library.spectra = cat(reduced_library...,dims=1)
    library.classes = new_classes
end

function reduce_endmembers_kmeans!(library::SpectralLibrary, max_endmembers_per_class::Int64)

    reduced_library = []
    new_classes = copy(library.classes)[1:max_endmembers_per_class*size(library.class_valid_keys)[1]]
    for (_c, cla) in enumerate(library.class_valid_keys)
        library_subset = library.spectra[library.classes .== cla, library.good_bands]

        out_spectra = zeros(max_endmembers_per_class, size(library.spectra)[2])
        out_spectra[:,:] .= NaN

        model = kmeans(library_subset', max_endmembers_per_class, maxiter=1000)
        out_spectra[:,library.good_bands] = model.centers'

        push!(reduced_library,out_spectra)
        new_classes[(_c-1) * max_endmembers_per_class + 1:_c * max_endmembers_per_class] .= cla
    end
    library.spectra = cat(reduced_library...,dims=1)
    library.classes = new_classes
end
