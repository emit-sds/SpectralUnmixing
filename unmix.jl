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
using ArgParse2
using DelimitedFiles
using Logging
using Distributed
using Printf

@everywhere using SpectralUnmixing

function main()

    parser = ArgumentParser(prog = "Spectral Unmixer",
                            description = "Unmix spectra in one of many ways")

    add_argument!(parser, "reflectance_file", type = String, help = "Input reflectance file")
    add_argument!(parser, "endmember_file", type = String, help = "Endmember reflectance deck")
    add_argument!(parser, "endmember_class_header", type = String, help = "header of column to use in endmember_file")
    add_argument!(parser, "output_file_base", type = String, help = "Output file base name")
    add_argument!(parser, "--spectral_starting_column", type = Int64, default = 2, help = "Column of library file that spectral information starts on")
    add_argument!(parser, "--truncate_end_columns", type = Int64, default = 0, help = "Columns to remove from spectral library at end of file")
    add_argument!(parser, "--reflectance_uncertainty_file", type = String, default = "", help = "Channelized uncertainty for reflectance input image")
    add_argument!(parser, "--n_mc", type = Int64, default = 1, help = "number of monte carlo runs to use, requires reflectance uncertainty file")
    add_argument!(parser, "--mode", type = String, default = "sma", help = "operating mode.  Options = [sma, mesma, plots]")
    add_argument!(parser, "--refl_nodata", type = Float64, default = -9999.0, help = "nodata value expected in input reflectance data")
    add_argument!(parser, "--refl_scale", type = Float64, default = 1.0, help = "scale image data (divide it by) this amount")
    add_argument!(parser, "--normalization", type = String, default = "none", help = "flag to indicate the scaling type. Options = [none, brightness, specific wavelength]")
    add_argument!(parser, "--combination_type", type = String, default = "class-even", help = "style of combinations.  Options = [all, class-even]")
    add_argument!(parser, "--max_combinations", type = Int64, default = -1, help = "set the maximum number of enmember combinations (relevant only to mesma)")
    add_argument!(parser, "--num_endmembers", type = Int64, default = [3], nargs="+", help = "set the maximum number of enmember to use")
    add_argument!(parser, "--write_complete_fractions", type=Bool, default = 1, help = "flag to indicate if per-endmember fractions should be written out")
    add_argument!(parser, "--write_spectral_residual", type=Bool, default = 1, help = "flag to indicate if spectral residuals should be written out")
    add_argument!(parser, "--optimizer", type=String, default = "bvls", help = "Choice of core optimization.  Options = [inverse, bvls, ldsqp]")
    add_argument!(parser, "--start_line", type=Int64, default = 1, help = "line of image to start on")
    add_argument!(parser, "--end_line", type=Int64, default = -1, help = "line of image to stop on (-1 does the full image)")
    add_argument!(parser, "--endmember_classes", type=String, default = [""], help = "Class names to use from the endmember library.  By default, all will be used")
    add_argument!(parser, "--wavelength_ignore_regions", nargs="+", type=Float64, default = [0,440,1310,1490,1770,2050,2440,2880], help = "Wavelength pairs of regions to ignore in the spectrum.")
    add_argument!(parser, "--log_file", type = String, default = nothing, help = "log file to write to")
    args = parse_args(parser)

    if isnothing(args.log_file)
        logger = Logging.SimpleLogger()
    else
        logger = Logging.SimpleLogger(open(args.log_file, "w+"))
    end

    if size(args.wavelength_ignore_regions)[1] % 2 != 0
        throw(ArgumentError("wavelength_ignore_regions must be an even number of values"))
    end
    Logging.global_logger(logger)

    valid_keys = nothing
    if args.endmember_classes[1] != ""
        valid_keys = args.endmember_classes
    end

    @info string("Unmixing was processed on: ", gethostname())
    @info string("Reflectance file processed: ", args.reflectance_file)
    @info string("Arguments: $(args)")
    
    library_scale_factor=1.0
    endmember_library = SpectralLibrary(args.endmember_file, args.endmember_class_header, args.spectral_starting_column, args.truncate_end_columns, valid_keys, library_scale_factor, args.wavelength_ignore_regions)
    load_data!(endmember_library)
    filter_by_class!(endmember_library)

    refl_file_wl = read_envi_wavelengths(args.reflectance_file)
    interpolate_library_to_new_wavelengths!(endmember_library, refl_file_wl)

    remove_wavelength_region_inplace!(endmember_library, true)

    reflectance_dataset = ArchGDAL.read(args.reflectance_file, alloweddrivers =["ENVI"])
    x_len = Int64(ArchGDAL.width(reflectance_dataset))
    y_len = Int64(ArchGDAL.height(reflectance_dataset))
    
    if args.start_line > y_len || args.start_line < 1;
       throw(ArgumentError(string("start_line must be less than length of scene, and greater than 1.  Provided: ", args.start_line, ", Scene length: ", y_len)))
    end

    if args.end_line == -1; 
        end_line = y_len; 
    else
        end_line = args.end_line
    end

    if end_line > y_len || end_line < args.start_line;
       throw(ArgumentError(string("end_line must be less than length of scene (or -1 for full scene), and greater than start_line.  Provided: ", args.end_line, ", start_line: ", args.start_line)))
    end
    @info string("Running from lines: ", args.start_line, " - ", end_line)

    if args.reflectance_uncertainty_file != ""
        reflectance_uncertainty_dataset = ArchGDAL.read(args.reflectance_uncertainty_file, alloweddrivers =["ENVI"])
        if ArchGDAL.width(reflectance_uncertainty_dataset) != x_len error("Reflectance_uncertainty_file size mismatch") end
        if ArchGDAL.height(reflectance_uncertainty_dataset) != y_len error("Reflectance_uncertainty_file size mismatch") end
        reflectance_uncertainty_dataset = nothing
    end


    if args.mode == "plots"
        plot_mean_endmembers(endmember_library, output_name=string(args.output_file_base, "_mean_endmembers.png"))
        plot_endmembers(endmember_library, output_name=string(args.output_file_base, "_endmembers.png"))
        plot_endmembers_individually(endmember_library, output_name=string(args.output_file_base, "_endmembers_individually.png"))
        exit()
    end


    n_classes = length(unique(endmember_library.classes))
    output_bands = [n_classes + 1]
    output_files = [string(args.output_file_base , "_fractional_cover")]
    
    if args.n_mc > 1
        push!(output_bands, n_classes + 1)
        push!(output_files, string(args.output_file_base , "_fractional_cover_uncertainty"))
    end

    if args.write_complete_fractions == 1
        push!(output_bands, size(endmember_library.spectra)[1] + 1)
        push!(output_files,string(args.output_file_base , "_complete_fractions") )
    end

    # if True, write a file with spectral residuals
    if args.write_spectral_residual == 1
            push!(output_bands, size(endmember_library.spectra)[1] + 1)
            push!(output_files,string(args.output_file_base , "_spectral_resudial") )
        end

    output_band_names = copy(endmember_library.class_valid_keys)
    println(output_band_names)
    push!(output_band_names, "Brightness")
    println(output_band_names)

    initiate_output_datasets(output_files, x_len, y_len, output_bands, reflectance_dataset)
    set_band_names(output_files[1], output_band_names)
    if args.n_mc > 1
        set_band_names(output_files[2], output_band_names)
    end
   
    # if True, write a file with spectral residuals
    if args.write_spectral_residual == 1
        z_len = [size(refl_file_wl)[1]]
        initiate_output_datasets([output_files[3]], x_len, y_len, z_len, reflectance_dataset)
        band_names = read_envi_wavelengths(args.reflectance_file)
        set_band_names(output_files[3], string.(band_names))
    end 

    @info string("Unmix output files: ", output_files)
    @info string("total number of workers available: ", nworkers())  
    results = @time pmap(line->unmix_and_write_line(line,args.reflectance_file, args.mode, args.refl_nodata,
               args.refl_scale, args.normalization, endmember_library, output_files, args.write_complete_fractions,
               args.write_spectral_residual, args.reflectance_uncertainty_file, args.n_mc,
               args.combination_type, args.num_endmembers, args.max_combinations, args.optimizer), args.start_line:end_line)
    
end


main()
