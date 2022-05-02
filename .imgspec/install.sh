#!/bin/bash

# Install script to 1) install Julia via conda, and 2) install the Julia dependencies for this project

set -x

# Define some useful directories
imgspec_dir=$( cd "$(dirname "$0")" ; pwd -P )
specun_dir=$(dirname ${imgspec_dir})

# Install Julia and then install Julia dependencies
conda create -n spectral-unmixing -y -c conda-forge julia=1.7 python=3.9
source activate spectral-unmixing
pushd $specun_dir
julia -e 'using Pkg; Pkg.activate("."); Pkg.add(path="https://github.com/kmsquire/ArgParse2.jl"); Pkg.instantiate()'
export JULIA_PROJECT=$specun_dir
popd
