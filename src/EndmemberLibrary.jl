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

"""
    nanargmin(input::Array)

Return the index of the minimum value in the input array, ignoring NaN values.
"""
function nanargmax(input::Array)
    x = copy(input)
    x[isnan.(x)] .= -Inf
    return argmax(x)
end

"""
    nanargmin(input::Array)

Return the index of the maximum value in the input array, ignoring NaN values.
"""
function nanargmin(input::Array)
    x = copy(input)
    x[isnan.(x)] .= Inf
    return argmin(x)
end

"""
    read_envi_wavelengths(filename::String, nm::Bool=true)

Read wavelength data from an ENVI header file with option to convert from microns to
nanometers.

- If no wavelengths are found, returns `nothing`.

# Notes
- The function constructs the expected header filename by replacing the extension of the
given filename with `.hdr`.
- This function logs an error if no wavelength data is found in the header file.
"""
function read_envi_wavelengths(filename::String, nm::Bool=true)
    header_name = splitext(filename)[1] * ".hdr"
    header = readlines(header_name)
    found = false
    for line in header
        if occursin("wavelength = {", line) ||
           occursin("wavelength= {", line) ||
           occursin("wavelength={", line)
            header = line
            found = true
            break
        end
    end

    if !found
        @error "No wavelength found in " * header_name
        return nothing
    end

    wavelengths = [
        parse(Float64, strip(x)) for x in split(split(split(header, "{")[2], "}")[1], ",")
    ]

    if nm && all(wavelengths .< 25)
        wavelengths = wavelengths .* 1000
        @info "Converting wavelengths read from $filename to nm from microns."
    end

    return wavelengths
end

"""
    get_good_bands_mask(wavelengths::Array{Float64}, wavelength_pairs)

Return a boolean mask for `wavelengths` *outside* the specified wavelength ranges given
by `wavelength_pairs`.
"""
function get_good_bands_mask(wavelengths::Array{Float64}, wavelength_pairs)
    good_bands = ones(Bool, length(wavelengths))

    for wvp in wavelength_pairs
        wavelength_diff = wavelengths .- wvp[1]
        wavelength_diff[wavelength_diff.<0] .= maximum(filter(!isnan, wavelength_diff))
        lower_index = nanargmin(wavelength_diff)

        wavelength_diff = wvp[2] .- wavelengths
        wavelength_diff[wavelength_diff.<0] .= maximum(filter(!isnan, wavelength_diff))
        upper_index = nanargmin(wavelength_diff)
        good_bands[lower_index:upper_index] .= false
    end

    # conditions = [((r[1] .<= wavelengths) .& (wavelengths .<= r[2])) for r in wavelength_pairs]
    # bad_bands = reduce(.|, conditions)
    # good_bands = Bool.(1 .- bad_bands)

    return good_bands
end

"""
    mutable struct SpectralLibrary

A structure for managing spectral data from a file, supporting class labels, scaling, and
wavelength filtering.
"""
mutable struct SpectralLibrary
    file_name::String
    class_header_name::String
    spectral_starting_column::Int64
    truncate_end_columns::Int64
    class_valid_keys
    scale_factor::Float64
    wavelength_regions_ignore
    @doc """
        SpectralLibrary(file_name::String,
            class_header_name::String,
            spectral_starting_column::Int64=2,
            truncate_end_columns::Int64=0,
            class_valid_keys=nothing,
            scale_factor=1.0,
            wavelength_regions_ignore=[0, 440, 1310, 1490, 1770, 2050, 2440, 2880])

    Constructor for the `SpectralLibrary` struct.

    # Arguments
    - `file_name::String`: Endmember Library data file name, CSV-like.
    - `class_header_name::String`: Column for class labels in data file.
    - `spectral_starting_column::Int64=2`: Column for start of spectral data
    - `truncate_end_columns::Int64=0`: Number of columns to ignore from end of data.
    - `class_valid_keys=nothing`: List of valid class labels.
    - `scale_factor::Float64=1.0`: Spectral scaling factor.
    - `wavelength_regions_ignore=[0,440,1310,1490,1770,2050,2440,2880]`: List of wavelength
    regions to ignore as beginning/end pairs in list.
    """
    SpectralLibrary(file_name::String,
        class_header_name::String,
        spectral_starting_column::Int64=2,
        truncate_end_columns::Int64=0,
        class_valid_keys=nothing,
        scale_factor=1.0,
        wavelength_regions_ignore=[0, 440, 1310, 1490, 1770, 2050, 2440, 2880]) =
        new(file_name, class_header_name, spectral_starting_column,
            truncate_end_columns, class_valid_keys, scale_factor, wavelength_regions_ignore)
    spectra
    classes
    good_bands
    wavelengths
end

"""
    load_data!(library::SpectralLibrary)

Load and preprocess spectral data associated with the `SpectralLibrary` object.

- Load spectral data from a CSV-like file `library.file_name`.
- Process wavelength ignore regions into paired sets `library.wavelength_regions_ignore`.
- Load `library.spectra`, `library.classes`, and `library.class_valid_keys`.
- Sort and load `library.wavelengths` as nanometers and rearrange spectral data accordingly.
- Generate mask `library.good_bands` based on `library.wavelength_regions_ignore`.

# Returns
- The updated `SpectralLibrary` instance with loaded and processed data.
"""
function load_data!(library::SpectralLibrary)

    df = DataFrame(CSV.File(library.file_name))

    paired_sets = []
    for i in 1:2:size(library.wavelength_regions_ignore)[1]
        push!(paired_sets, [library.wavelength_regions_ignore[i],
            library.wavelength_regions_ignore[i+1]])
    end
    library.wavelength_regions_ignore = paired_sets
    @info "Ignoring wavelength regions: " * string(library.wavelength_regions_ignore)

    try
        library.spectra = Matrix(
            df[:, library.spectral_starting_column:end-library.truncate_end_columns])
    catch e
        throw(ArgumentError("Could not read spectra from endmember library.  Try adjusting
        spectral_starting_column or truncate_end_columns, or reworking the library."))
    end

    try
        library.classes = convert(Array{AbstractString}, df[!, library.class_header_name])
    catch e
        throw(ArgumentError("Could not read classes from endmember library using class key"
                            * library.class_header_name *
                            ".  Try adjusting class_header_name, or reworking the library."
        ))
    end

    try
        library.wavelengths = parse.(Float64, names(
            df[:, library.spectral_starting_column:end-library.truncate_end_columns])
        )
    catch e
        throw(ArgumentError("Could not read wavelengths from endmember library.  Try
        adjusting class_header_name, or reworking the library."))
    end


    wl_order = sortperm(library.wavelengths)
    library.spectra = library.spectra[:, wl_order]
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
        @info "Unit mismatch between library wavelengths and wavelength regions to ignore
        detected.  Converting library wavelengths to nm from microns."
    end

    return library
end

"""
    save_data(library::SpectralLibrary, output_filename::String;
    class_label_name::String="Label")

Wrte the spectra and classes from a `SpectralLibrary` to a CSV file.

- The output CSV will have the class labels as the first column under `class_label_name`,
followed by columns for each wavelength with the spectra values filling the corresponding
rows.
"""
function save_data(library::SpectralLibrary, output_filename::String;
    class_label_name::String="Label")

    df = DataFrame()
    insertcols!(df, 1, class_label_name => library.classes)
    for (_wv, wv) in enumerate(library.wavelengths)
        insertcols!(df, 1 + _wv, string(wv) => library.spectra[:, _wv])
    end
    CSV.write(output_filename, df)
end

"""
    filter_by_class!(library::SpectralLibrary)

Filter spectra in the `SpectralLibrary` based on `library.class_valid_keys`.

- Update `library.spectra` and `library.classes` to include only those matching `library.
class_valid_keys`.

# Notes
- This function modifies the `library` in place.
- If `library.class_valid_keys` is `nothing`, logs a message indicating that no filtering
will occur, and exits without making any changes.
"""
function filter_by_class!(library::SpectralLibrary)
    if isnothing(library.class_valid_keys)
        @info "No class valid keys provided, no filtering occuring"
        return
    end

    valid_classes = zeros(Bool, size(library.spectra)[1])
    for cla in library.class_valid_keys
        valid_classes[library.classes.==cla] .= true
    end

    library.spectra = library.spectra[valid_classes, :]
    library.classes = library.classes[valid_classes]
end

"""
    remove_wavelength_region_inplace!(library::SpectralLibrary, set_as_nans::Bool=false)

Remove wavelength regions from `library` outside `library.good_bands` either by setting
them to NaN or by filtering them out.

# Notes
- This function modifies the `library` in place.
"""
function remove_wavelength_region_inplace!(library::SpectralLibrary,
    set_as_nans::Bool=false)
    if set_as_nans
        library.spectra[:, .!library.good_bands] .= NaN
        library.wavelengths[.!library.good_bands] .= NaN
    else
        library.spectra = library.spectra[:, library.good_bands]
        library.wavelengths = library.wavelengths[library.good_bands]
        library.good_bands = library.good_bands[library.good_bands]
    end
end

"""
    interpolate_library_to_new_wavelengths!(library::SpectralLibrary, new_wavelengths::Array
    {Float64})

Interpolate `library` spectra to new wavelengths.

- Perform linear interpolation, applying a flat extrapolation boundary condition.
- Update `library.wavelengths`, `library.good_bands`, and `library.spectra` matrix to the
resampled spectra.

# Notes
- This function modifies the `library` in place
"""
function interpolate_library_to_new_wavelengths!(library::SpectralLibrary,
    new_wavelengths::Array{Float64})
    old_spectra = copy(library.spectra)

    library.spectra = zeros((size(library.spectra)[1], length(new_wavelengths)))
    for _s in 1:size(old_spectra)[1]
        fit = LinearInterpolation(
            library.wavelengths, old_spectra[_s, :], extrapolation_bc=Flat()
        )
        library.spectra[_s, :] = fit(new_wavelengths)
    end
    library.wavelengths = new_wavelengths

    good_bands = get_good_bands_mask(library.wavelengths, library.wavelength_regions_ignore)
    library.good_bands = good_bands
end

"""
    scale_library!(library::SpectralLibrary, scaling_factor=nothing)

Scale the spectral data in the `SpectralLibrary` by the library's defined scale factor
(default) or a specified scaling factor.
"""
function scale_library!(library::SpectralLibrary, scaling_factor=nothing)
    if isnothing(scaling_factor)
        library.spectra = library.spectra ./ library.scale_factor
    else
        library.spectra /= scaling_factor
    end
end

"""
    brightness_normalize!(library::SpectralLibrary)

Normalize `library.spectra` based on the RMS brightness of `library.good_bands`.

# Notes
- This function modifies the `library` in place.
- Only considers `library.good_bands` in the brightness calculation.
"""
function brightness_normalize!(library::SpectralLibrary)
    library.spectra = library.spectra ./
                      sqrt.(mean(library.spectra[:, library.good_bands] .^ 2, dims=2))
end

"""
    split_library(library::SpectralLibrary, split_fraction::Float64)

Split a `SpectralLibrary` into two new libraries based on a specified fraction of the
total spectra.

# Returns
- `Tuple{SpectralLibrary, SpectralLibrary}` with each output library a random,
mutually exclusive subset of the original library containing `split_fraction` and `1 -
split_fraction` of the spectra, respectively.

# Notes
- The split is random; consecutive calls with the same `library` may yield different
results.
"""
function split_library(library::SpectralLibrary, split_fraction::Float64)

    if !(0 < split_fraction < 1)
        throw(ArgumentError("split_fraction must be between 0 and 1 (exclusive)."))
    end
    perm = randperm(size(library.spectra)[1])

    split_1 = perm[1:Int(round(split_fraction * length(perm)))]
    split_2 = perm[Int(round(split_fraction * length(perm))):end]

    output_library_1 = deepcopy(library)
    output_library_2 = deepcopy(library)

    output_library_1.spectra = output_library_1.spectra[split_1, :]
    output_library_2.spectra = output_library_2.spectra[split_2, :]

    output_library_1.classes = output_library_1.classes[split_1]
    output_library_2.classes = output_library_2.classes[split_2]

    return output_library_1, output_library_2

end

"""
    reduce_endmembers_nmf!(library::SpectralLibrary, max_endmembers_per_class::Int64)

Reduce the number of endmembers in the spectral library using Non-negative Matrix
Factorization (NMF).

- For each class in `library.class_valid_keys`, apply NMF to the spectra subset to identify
the specified maximum number of endmembers.
- Update `library.classes` and `library.spectra` to contain the reduced set of endmembers.

# Notes
- This function modifies the `library` in place.
- Only `library.good_bands` will be considered in the NMF computation.
- The NMF uses a maximum of 500 iterations with a tolerance of 1.0e-2 for convergence.
"""
function reduce_endmembers_nmf!(library::SpectralLibrary, max_endmembers_per_class::Int64)

    reduced_library = []
    new_classes =
        copy(library.classes)[1:max_endmembers_per_class*size(library.class_valid_keys)[1]]
    for (_c, cla) in enumerate(library.class_valid_keys)
        library_subset = library.spectra[library.classes.==cla, library.good_bands]

        out_spectra = zeros(max_endmembers_per_class, size(library.spectra)[2])
        out_spectra[:, :] .= NaN

        r = nnmf(library_subset, max_endmembers_per_class; maxiter=500, tol=1.0e-2)
        out_spectra[:, library.good_bands] = r.H

        push!(reduced_library, out_spectra)
        new_classes[(_c-1)*max_endmembers_per_class+1:_c*max_endmembers_per_class] .= cla
    end
    library.spectra = cat(reduced_library..., dims=1)
    library.classes = new_classes
end

"""
    reduce_endmembers_pca!(library::SpectralLibrary, max_endmembers_per_class::Int64)

Reduce the number of endmembers in the spectral library using Principal Component Analysis
(PCA).

- For each class in `library.class_valid_keys`, apply PCA to the spectra subset to identify
the specified maximum number of endmembers.
- Update `library.classes` and `library.spectra` to contain the reduced set of endmembers.

# Notes
- This function modifies the `library` in place.
- Only `library.good_bands` will be considered in the PCA computation.
- The PCA is truncated to `max_endmembers_per_class` PCs per endmember class.
"""
function reduce_endmembers_pca!(library::SpectralLibrary, max_endmembers_per_class::Int64)

    reduced_library = []
    new_classes =
        copy(library.classes)[1:max_endmembers_per_class*size(library.class_valid_keys)[1]]
    for (_c, cla) in enumerate(library.class_valid_keys)
        library_subset = library.spectra[library.classes.==cla, library.good_bands]

        out_spectra = zeros(max_endmembers_per_class, size(library.spectra)[2])
        out_spectra[:, :] .= NaN

        model = fit(PCA, library_subset', maxoutdim=max_endmembers_per_class, pratio=1)
        out_spectra[:, library.good_bands] = projection(model)'

        push!(reduced_library, out_spectra)
        new_classes[(_c-1)*max_endmembers_per_class+1:_c*max_endmembers_per_class] .= cla
    end
    library.spectra = cat(reduced_library..., dims=1)
    library.classes = new_classes
end

"""
    reduce_endmembers_kmeans!(library::SpectralLibrary, max_endmembers_per_class::Int64)

Reduce the number of endmembers in the spectral library using K-means clustering.

- For each class in `library.class_valid_keys`, apply K-means clustering to the spectra
subset to identify the specified maximum number of endmembers from cluster centers.
- Update `library.classes` and `library.spectra` to contain the reduced set of endmembers.

# Notes
- This function modifies the `library` in place.
- Only `library.good_bands` will be considered in the K-means computation.
- The K-means algorithm is run with a maximum of 1000 iterations.
"""
function reduce_endmembers_kmeans!(library::SpectralLibrary,
    max_endmembers_per_class::Int64)

    reduced_library = []
    new_classes =
        copy(library.classes)[1:max_endmembers_per_class*size(library.class_valid_keys)[1]]
    for (_c, cla) in enumerate(library.class_valid_keys)
        library_subset = library.spectra[library.classes.==cla, library.good_bands]

        out_spectra = zeros(max_endmembers_per_class, size(library.spectra)[2])
        out_spectra[:, :] .= NaN

        model = kmeans(library_subset', max_endmembers_per_class, maxiter=1000)
        out_spectra[:, library.good_bands] = model.centers'

        push!(reduced_library, out_spectra)
        new_classes[(_c-1)*max_endmembers_per_class+1:_c*max_endmembers_per_class] .= cla
    end
    library.spectra = cat(reduced_library..., dims=1)
    library.classes = new_classes
end
