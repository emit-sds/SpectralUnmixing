<h1 align="center">
SpectralUnmixing
</h1>
A general, fast, flexible, and including spectral unmixing package.  Oriented towards VSWIR imaging spectroscopy data but applicable for different sensor types.  Includes options for different treatments of endmember library assemblages, including MESMA and bootstrapping (aka Monte Carlo) strategies.

## Installation
This package has not been registered yet, but will be soon.  In the interim, after cloning and navigating into the repository, it can be installed from the Julia REPL, by running

```
julia -e 'using Pkg; Pkg.activate("."); Pkg.precompile()'
export JULIA_PROJECT=${PWD}
```

## Using the script
Currently the package supports reading and writing ENVI raster data.

Basic:

```
julia unmix.jl REFLECTANCE_IMAGE ENDMEMBER_LIBRARY ENDMEMBER_COLUMN OUTPUT_BASE --mode sma
```


Parallel implementation (with 10 cores):

```
julia -p 10 unmix.jl REFLECTANCE_IMAGE ENDMEMBER_LIBRARY ENDMEMBER_COLUMN OUTPUT_BASE --mode sma
```

Bootstrapping uncertainty:

```
julia -p 10 unmix.jl REFLECTANCE_IMAGE ENDMEMBER_LIBRARY ENDMEMBER_COLUMN OUTPUT_BASE --mode sma --n_mc 50
```

Normalization:

```
julia -p 10 unmix.jl REFLECTANCE_IMAGE ENDMEMBER_LIBRARY ENDMEMBER_COLUMN OUTPUT_BASE --mode sma --n_mc 50 --normalization brightness
```

Preset maximum number of endmembers used for unmixing:

```
julia -p 10 unmix.jl REFLECTANCE_IMAGE ENDMEMBER_LIBRARY ENDMEMBER_COLUMN OUTPUT_BASE --mode sma --n_mc 50 --normalization brightness --num_endmembers 10
```

## EMIT-style Runs
Following [Ochoa et al.](https://d197for5662m48.cloudfront.net/documents/publicationstatus/232672/preprint_pdf/973acea360e10b97752976bf19e5c071.pdf), to run SpectralUnmixing in the same manner as EMIT, use:

```
julia -p 64 unmix.jl REFLECTANCE_IMAGE ENDMEMBER_LIBRARY ENDMEMBER_COLUMN OUTPUT_BASE --mode sma-best --normalization brightness --num_endmember 30 --n_mc 20 --reflectance_uncertainty_file REFLECTANCE_UNCERTAINTY_IMAGE
```
