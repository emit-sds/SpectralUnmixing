#!/bin/bash

# Install script to 1) install Julia via conda, and 2) install the Julia dependencies for this project

set -x

# Define some useful directories
imgspec_dir=$( cd "$(dirname "$0")" ; pwd -P )
specun_dir=$(dirname ${imgspec_dir})
cur_dir=$(pwd -P)

# Install Julia and then install Julia dependencies
conda install -y -c conda-forge julia=1.7
cd $specun_dir
julia -e 'using Pkg; Pkg.activate("."); Pkg.add(path="https://github.com/kmsquire/ArgParse2.jl"); Pkg.instantiate()'
export JULIA_PROJECT=$specun_dir

# Return to original directory
cd $cur_dir
