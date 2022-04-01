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
# $3: Number of cores to use (default is "1")
# $4: Nodata value expected in input reflectance data (default is "None" - will use the repo's default)
# $5: Scale image data (divide it by) this amount (default is "None" - will use the repo's default)
# $6: Flag to indicate the scaling type. Options = [none, brightness, specific wavelength] (default is "None" - will
#     use the repo's default)
#
# In addition to the positional arguments, this script expects a downloaded reflectance granule to be present in a
# folder called "input".

# Define some useful directories
imgspec_dir=$( cd "$(dirname "$0")" ; pwd -P )
specun_dir=$(dirname ${imgspec_dir})
input="input"
mkdir -p output

# Activate conda environment
source activate spectral-unmixing

# Export JULIA_PROJECT again in case it doesn't carry over from install.sh
export JULIA_PROJECT=$specun_dir

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
if [[ $rfl_name == f* ]]; then
    output_base=$(echo $rfl_name | cut -c1-16)
elif [[ $rfl_name == ang* ]]; then
    output_base=$(echo $rfl_name | cut -c1-18)
elif [[ $rfl_name == PRS* ]]; then
    output_base=$(echo $rfl_name | cut -c1-38)
fi
output_base_path="output/$output_base"
echo "Output base path: $output_base_path"

# Build command and execute
cmd="julia"

# Add number of cores
if [[ $3 != "1" ]]; then
    cmd="$cmd -p $3"
fi

# Add the required args (and assume mode is sma for now)
cmd="$cmd $unmix_exe $rfl_path $endmember_library_path $2 $output_base_path --mode=sma"

# Add the optional args
if [[ $4 != "None" ]]; then
    cmd="$cmd --refl_nodata=$4"
fi
if [[ $5 != "None" ]]; then
    cmd="$cmd --refl_scale=$5"
fi
if [[ $6 != "None" ]]; then
    cmd="$cmd --normalization=$6"
fi

echo "Executing command: $cmd"
$cmd
