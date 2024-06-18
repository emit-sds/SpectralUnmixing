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


function set_band_names(filename::String, band_names::Vector)
    GC.gc()
    local_ods = GDAL.gdalopen(filename, GDAL.GA_Update)
    for _b in 1:(length(band_names))
        GDAL.gdalsetdescription(GDAL.gdalgetrasterband(local_ods, _b), band_names[_b])
    end
    local_ods = nothing
end

function initiate_output_datasets(output_files::Vector{String}, x_len::Int64, y_len::Int64, output_bands::Vector, reference_dataset::ArchGDAL.IDataset)
    GC.gc()
    for _o in 1:length(output_files)
        outDataset = ArchGDAL.create(output_files[_o], driver=ArchGDAL.getdriver("ENVI"), width=x_len,
                                     height=y_len, nbands=output_bands[_o], dtype=Float64)
        ArchGDAL.setproj!(outDataset, ArchGDAL.getproj(reference_dataset))
        @info "Output Image Size (x,y,b): $x_len, $y_len, $output_bands.  Creating output fractional cover dataset."
        try
            ArchGDAL.setgeotransform!(outDataset, ArchGDAL.getgeotransform(reference_dataset))
        catch e
            println("No geotransorm avaialble, proceeding without")
        end

        for _b in 1:length(output_bands)
            ArchGDAL.setnodatavalue!(ArchGDAL.getband(outDataset,_b), -9999)
        end

        output = zeros(x_len, y_len, output_bands[_o]) .- 9999
        ArchGDAL.write!(outDataset, output, [1:size(output)[end];], 0, 0, size(output)[1], size(output)[2])

        outDataset = nothing

    end
end


function write_results(output_files::Vector{String}, output_bands::Vector, x_len::Int64, y_len::Int64, results, args)
    GC.gc()

    @info "Write primary output"
    output = zeros(y_len, x_len, output_bands[1]) .- 9999
    for res in results
        if isnothing(res[2]) == false
            output[res[1],res[3], :] = res[2]
        end
    end
    output = permutedims( output, (2,1,3))
    outDataset = ArchGDAL.read(output_files[1], flags=1, alloweddrivers =["ENVI"])
    ArchGDAL.write!(outDataset, output, [1:size(output)[end];], 0, 0, size(output)[1], size(output)[2])
    outDataset = nothing

    ods_idx = 2
    if args.n_mc > 1
        @info "Write uncertainty output"
        output = zeros(y_len, x_len, output_bands[ods_idx]) .- 9999
        for res in results
            if isnothing(res[4]) == false
                output[res[1],res[3], :] = res[4]
            end
        end

        output = permutedims( output, (2,1,3))
        outDataset = ArchGDAL.read(output_files[ods_idx], flags=1, alloweddrivers =["ENVI"])
        ArchGDAL.write!(outDataset, output, [1:size(output)[end];], 0, 0, size(output)[1], size(output)[2])
        outDataset = nothing
        ods_idx += 1
    end

    if args.write_complete_fractions == 1
        @info "Write complete fractions output"
        output = zeros(y_len, x_len, output_bands[ods_idx]) .- 9999
        for res in results
            if isnothing(res[5]) == false
                output[res[1],res[3], :] = res[5]
            end
        end

        output = permutedims( output, (2,1,3))
        outDataset = ArchGDAL.read(output_files[ods_idx], flags=1, alloweddrivers =["ENVI"])
        ArchGDAL.write!(outDataset, output, [1:size(output)[end];], 0, 0, size(output)[1], size(output)[2])
        outDataset = nothing
        ods_idx += 1
    end

end

function write_line_results(output_files::Vector{String}, results, n_mc::Int64, write_complete_fractions::Bool)
    GC.gc()
    
    if isnothing(results[2]) == false
        
        output = zeros(size(results[3])[1], 1, size(results[2])[2]) .- 9999
        output[results[3], 1, :] = results[2]

        outDataset = ArchGDAL.read(output_files[1], flags=1, alloweddrivers =["ENVI"])
        ArchGDAL.write!(outDataset, output, [1:size(output)[end];], 0, results[1]-1, size(output)[1], 1)
        outDataset = nothing
    end
    ods_idx = 2

    if n_mc > 1
        if isnothing(results[4]) == false
            output = zeros(size(results[3])[1], 1, size(results[4])[2]) .- 9999
            output[results[3], 1, :] = results[4]
            
            outDataset = ArchGDAL.read(output_files[ods_idx], flags=1, alloweddrivers =["ENVI"])
            ArchGDAL.write!(outDataset, output, [1:size(output)[end];], 0, results[1]-1, size(output)[1], 1)
            outDataset = nothing
        end
        ods_idx += 1
    end

    if write_complete_fractions == 1
        if isnothing(results[5]) == false
            output = zeros(size(results[3])[1], 1, size(results[5])[2]) .- 9999
            output[results[3], 1, :] = results[5]
            
            outDataset = ArchGDAL.read(output_files[ods_idx], flags=1, alloweddrivers =["ENVI"])
            ArchGDAL.write!(outDataset, output, [1:size(output)[end];], 0, results[1]-1, size(output)[1], 1)
            outDataset = nothing
        end
        ods_idx += 1
    end

end

function load_line(reflectance_file::String, reflectance_uncertainty_file::String, line::Int64,
                       good_bands::Array{Bool}, refl_nodata::Float64)
        
        reflectance_dataset = ArchGDAL.read(reflectance_file, alloweddrivers =["ENVI"])
        img_dat = convert(Array{Float64},ArchGDAL.readraster(reflectance_file, alloweddrivers =["ENVI"])[:,line,:])
        img_dat = img_dat[:, good_bands]
        good_data = .!all(img_dat .== refl_nodata, dims=2)[:,1]
        img_dat = img_dat[good_data,:]

        if sum(good_data) >= 1
            if reflectance_uncertainty_file != ""
                unc_dat = convert(Array{Float64},ArchGDAL.readraster(reflectance_uncertainty_file, alloweddrivers =["ENVI"])[:,line,:])
                unc_dat = unc_dat[:, good_bands]
                unc_dat = unc_dat[good_data,:]
            else
                unc_dat = nothing
            end
        else
            return nothing, nothing, good_data
        end

        return img_dat, unc_dat, good_data
end
