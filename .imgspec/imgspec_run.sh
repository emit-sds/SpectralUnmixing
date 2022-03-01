#!/bin/bash

# Description:
#
# The top-level run script to execute unmix.jl on ImgSPEC. This script accepts the inputs described in the
# algorithm_config.yaml file and pre-processes them as needed to pass into SpectralUnmixing/unmix.jl.  This script
# is currently compatible with AVIRIS Classic, AVIRIS-NG, and PRISMA data.
#
# Inputs:
#
# $1: URL of endmember library
# $2: Column name from endmember library that describes the library classes (default is "class")
#
# In addition to the positional arguments, this script expects a downloaded reflectance granule to be present in a
# folder called "input".


imgspec_dir=$( cd "$(dirname "$0")" ; pwd -P )
specun_dir=$(dirname ${imgspec_dir})

input="input"
mkdir -p output

# Spectral Unmixing paths
unmix_exe="$specun_dir/unmix.jl"
endmember_library_path="$input/endmember_library.csv"

# Get reflectance path from input granule
echo "Looking for input granule gzip and extracting to get reflectance path..."
rfl_path=$(python ${imgspec_dir}/get_paths_from_granules.py -p rfl)
echo "Found input reflectance file: $rfl_path"

# Process positional args to get EcoSIS CSV files
echo "Getting endmember library..."
curl --retry 10 --output $endmember_library_path $1

# Get output base from reflectance path
rfl_name=$(basename $rfl_path)
output_base="output"
if [[ rfl_name == f* ]]; then
    output_base=$(echo rfl_name | cut -c1-16)
elif [[ rfl_name == ang* ]]; then
    output_base=$(echo rfl_name | cut -c1-18)
elif [[ rfl_name == PRS* ]]; then
    output_base=$(echo rfl_name | cut -c1-38)
fi

cmd="julia $unmix_exe $rfl_path $endmember_library_path $2 $output_base"
echo "Executing command: $cmd"
$cmd
