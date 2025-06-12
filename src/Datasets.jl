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

using ArchGDAL
using GDAL

"""
    set_band_names(filename::String, band_names::Vector)

Set the names of the bands in a raster dataset specified by the given filename.

# Arguments
- `filename::String`: Path to the raster file. The file must be in a format supported by
GDAL. The data file is updated in place.
- `band_names::Vector`: A vector of strings containing the new names for each band in the
raster dataset. The length of this vector must match the number of bands in the raster.
"""
function set_band_names(filename::String, band_names::Vector)
    GC.gc()
    local_ods = GDAL.gdalopen(filename, GDAL.GA_Update)
    for _b in 1:(length(band_names))
        GDAL.gdalsetdescription(GDAL.gdalgetrasterband(local_ods, _b), band_names[_b])
    end
    local_ods = nothing
end

"""
    initiate_output_datasets(output_files::Vector{String}, x_len::Int64, y_len::Int64,
                             output_bands::Vector, reference_dataset::ArchGDAL.IDataset)

Initialize multiple output raster datasets in ENVI format with specified dimensions and
band information using a reference dataset for projection and geotransformation.

# Arguments
- `output_files::Vector{String}`: File paths of the output datasets. Each entry corresponds
to a separate output raster.
- `x_len::Int64`: Width of the output datasets.
- `y_len::Int64`: Height of the output datasets.
- `output_bands::Vector`: Number of bands for each output dataset. The length of this
vector must match the length of `output_files`.
- `reference_dataset::ArchGDAL.IDataset`: Reference dataset object from which projection
and geotransformation information will be copied.

# Notes
- All bands in the output datasets are initialized with the no-data value of `-9999`.
"""
function initiate_output_datasets(output_files::Vector{String}, x_len::Int64, y_len::Int64,
    output_bands::Vector, reference_dataset::ArchGDAL.IDataset)
    GC.gc()
    for _o in 1:length(output_files)
        outDataset = ArchGDAL.create(
            output_files[_o], driver=ArchGDAL.getdriver("ENVI"), width=x_len,
            height=y_len, nbands=output_bands[_o], dtype=Float64
        )
        ArchGDAL.setproj!(outDataset, ArchGDAL.getproj(reference_dataset))
        @info "Output Image Size (x,y,b): $x_len, $y_len, $output_bands.
        Creating output fractional cover dataset."
        try
            ArchGDAL.setgeotransform!(
                outDataset, ArchGDAL.getgeotransform(reference_dataset)
            )
        catch e
            println("No geotransorm available, proceeding without")
        end

        for _b in 1:length(output_bands)
            ArchGDAL.setnodatavalue!(ArchGDAL.getband(outDataset, _b), -9999)
        end

        output = zeros(x_len, y_len, output_bands[_o]) .- 9999
        ArchGDAL.write!(
            outDataset, output, [1:size(output)[end];], 0, 0, size(output)[1],
            size(output)[2]
        )

        outDataset = nothing
    end
end

"""
    write_results(output_files::Vector{String}, output_bands::Vector, x_len::Int64,
                  y_len::Int64, results, args)

Write processed results to specified output ENVI raster files, including primary outputs,
uncertainty estimates, and complete fractions, based on the provided parameters.

# Arguments
- `output_files::Vector{String}`: File paths of the output datasets to be created or
updated. Each entry corresponds to a separate output raster.
- `output_bands::Vector`: Number of bands for each output dataset. The length of this
vector must match the number of output files.
- `x_len::Int64`: Width of the output datasets.
- `y_len::Int64`: Height of the output datasets.
- `results`: A collection of results, where each result is expected to be a tuple or
similar structure. Specific elements contain the data to be written to the output datasets.
- `args`: An object or structure containing additional parameters that dictate the
writing behavior. Available options are:
    - `args.n_mc`: If greater than 1, an uncertainty output is also written to the second
    output file.
    - `args.write_complete_fractions`: If set to 1, write complete fractions to the
    subsequent output file.

# Notes
- The primary output is initialized with a no-data value of `-9999` and written to the
  first output file specified in `output_files`.
- Each output array is permuted to match the desired dimensions before writing.
"""
function write_results(output_files::Vector{String}, output_bands::Vector, x_len::Int64,
    y_len::Int64, results, args)

    GC.gc()

    @info "Write primary output"
    output = zeros(y_len, x_len, output_bands[1]) .- 9999
    for res in results
        if isnothing(res[2]) == false
            output[res[1], res[3], :] = res[2]
        end
    end
    output = permutedims(output, (2, 1, 3))
    outDataset = ArchGDAL.read(output_files[1], flags=1, alloweddrivers=["ENVI"])
    ArchGDAL.write!(
        outDataset, output, [1:size(output)[end];], 0, 0, size(output)[1], size(output)[2]
    )
    outDataset = nothing

    ods_idx = 2
    if args.n_mc > 1
        @info "Write uncertainty output"
        output = zeros(y_len, x_len, output_bands[ods_idx]) .- 9999
        for res in results
            if isnothing(res[4]) == false
                output[res[1], res[3], :] = res[4]
            end
        end

        output = permutedims(output, (2, 1, 3))
        outDataset = ArchGDAL.read(output_files[ods_idx], flags=1, alloweddrivers=["ENVI"])
        ArchGDAL.write!(
            outDataset, output, [1:size(output)[end];], 0, 0, size(output)[1],
            size(output)[2]
        )
        outDataset = nothing
        ods_idx += 1
    end

    if args.write_complete_fractions == 1
        @info "Write complete fractions output"
        output = zeros(y_len, x_len, output_bands[ods_idx]) .- 9999
        for res in results
            if isnothing(res[5]) == false
                output[res[1], res[3], :] = res[5]
            end
        end

        output = permutedims(output, (2, 1, 3))
        outDataset = ArchGDAL.read(output_files[ods_idx], flags=1, alloweddrivers=["ENVI"])
        ArchGDAL.write!(
            outDataset, output, [1:size(output)[end];], 0, 0, size(output)[1],
            size(output)[2]
        )
        outDataset = nothing
        ods_idx += 1
    end
end

"""
    write_line_results(output_files::Vector{String}, results, n_mc::Int64,
                       write_complete_fractions::Bool)

Write line-based [`unmix_line`](@ref) results to specified output ENVI raster files.

# Arguments
- `output_files::Vector{String}`: File paths of the output datasets to be created or
updated. The first entry is for pixel mixture fractions, subsequent entries
are for uncertainty and complete fractions, in that order, if applicable.
- `results`: Results to be written, as produced by [`unmix_line`](@ref).
- `n_mc::Int64`: If greater than one, uncertainty outputs will be written to the
subsequent output file.
- `write_complete_fractions::Bool`: Indicates whether to write complete fractions
to the subsequent output file.
"""
function write_line_results(output_files::Vector{String}, results, n_mc::Int64,
    write_complete_fractions::Bool)

    GC.gc()

    if isnothing(results[2]) == false

        output = zeros(size(results[3])[1], 1, size(results[2])[2]) .- 9999
        output[results[3], 1, :] = results[2]

        outDataset = ArchGDAL.read(output_files[1], flags=1, alloweddrivers=["ENVI"])
        ArchGDAL.write!(
            outDataset, output, [1:size(output)[end];], 0, results[1] - 1,
            size(output)[1], 1
        )
        outDataset = nothing
    end
    ods_idx = 2

    if n_mc > 1
        if isnothing(results[4]) == false
            output = zeros(size(results[3])[1], 1, size(results[4])[2]) .- 9999
            output[results[3], 1, :] = results[4]

            outDataset = ArchGDAL.read(
                output_files[ods_idx], flags=1, alloweddrivers=["ENVI"]
            )
            ArchGDAL.write!(
                outDataset, output, [1:size(output)[end];], 0, results[1] - 1,
                size(output)[1], 1
            )
            outDataset = nothing
        end
        ods_idx += 1
    end

    if write_complete_fractions == 1
        if isnothing(results[5]) == false
            output = zeros(size(results[3])[1], 1, size(results[5])[2]) .- 9999
            output[results[3], 1, :] = results[5]

            outDataset = ArchGDAL.read(
                output_files[ods_idx], flags=1, alloweddrivers=["ENVI"]
            )
            ArchGDAL.write!(
                outDataset, output, [1:size(output)[end];], 0, results[1] - 1,
                size(output)[1], 1
            )
            outDataset = nothing
        end
        ods_idx += 1
    end
end

"""
    load_line(reflectance_file::String, reflectance_uncertainty_file::String,
               line::Int64, good_bands::Array{Bool}, refl_nodata::Float64)

Load a specific line (row) of reflectance data and its associated uncertainty from raster
files, filtering based on specified good bands and no-data values.

# Arguments
- `reflectance_file::String`: Path to the ENVI reflectance raster file.
- `reflectance_uncertainty_file::String`: Path to the uncertainty raster file corresponding
to the reflectance data. An empty string indicates no uncertainty data.
- `line::Int64`: Line (row) of the raster data to be loaded.
- `good_bands::Array{Bool}`: Boolean array indicating which bands of the data to keep and
process.
- `refl_nodata::Float64`: The no-data value used in the reflectance data.

# Returns
- A tuple containing:
  - `img_dat::Array{Float64}`: The filtered reflectance data for the specified line and
  good bands. If no valid data exists, this will be `nothing`.
  - `unc_dat::Array{Float64}`: The uncertainty data, filtered similarly. If no uncertainty
  data is loaded or valid, this will be `nothing`.
  - `good_data`: A boolean array indicating which pixels contain valid reflectance data.
"""
function load_line(reflectance_file::String, reflectance_uncertainty_file::String,
    line::Int64, good_bands::Array{Bool}, refl_nodata::Float64)

    reflectance_dataset = ArchGDAL.read(reflectance_file, alloweddrivers=["ENVI"])
    img_dat = convert(
        Array{Float64},
        ArchGDAL.readraster(reflectance_file, alloweddrivers=["ENVI"])[:, line, :]
    )
    img_dat = img_dat[:, good_bands]
    good_data = .!all(img_dat .== refl_nodata, dims=2)[:, 1]
    img_dat = img_dat[good_data, :]

    if sum(good_data) >= 1
        if reflectance_uncertainty_file != ""
            unc_dat = convert(
                Array{Float64},
                ArchGDAL.readraster(
                    reflectance_uncertainty_file, alloweddrivers=["ENVI"]
                )[:, line, :]
            )
            unc_dat = unc_dat[:, good_bands]
            unc_dat = unc_dat[good_data, :]
        else
            unc_dat = nothing
        end
    else
        return nothing, nothing, good_data
    end

    return img_dat, unc_dat, good_data
end
