<h1 align="center">
<br>
<a href="https://github.com/emit-sds/SpectralUnmixing.jl"><img src="docs/src/assets/logo.svg" alt="SpectralUnmixing.jl" width="200"></a>
<br>
SpectralUnmixing.jl
<br>
</h1>


[![version](https://github.com/emit-sds/SpectralUnmixing/actions/workflows/unit-tests.yml/badge.svg)](https://github.com/emit-sds/SpectralUnmixing/actions/workflows/unit-tests.yml/)
[![](https://img.shields.io/github/license/emit-sds/SpectralUnmixing)](https://github.com/emit-sds/SpectralUnmixing/blob/master/LICENSE)
[![](https://img.shields.io/badge/docs-latest-blue)](https://emit-sds.github.io/SpectralUnmixing.jl)

A general, fast, flexible, and including spectral unmixing package.  Oriented towards VSWIR imaging spectroscopy data but applicable for different sensor types.  Includes options for different treatments of endmember library assemblages, including MESMA and bootstrapping (aka monte carlo) strategies.


## Installation
This package is registered and can be added using the Julia package manager.
1. Install [Julia](https://julialang.org/install/)
2. Install the package from the Julia Pkg REPL:
```
julia
julia> ]
pkg> add SpectralUnmixing
```
or equivalently,
```
julia 'using Pkg; Pkg.add("SpectralUnmixing")'
```
Remember to activate the appropriate [environment](https://pkgdocs.julialang.org/v1/) via the Pkg REPL or by using the `--project` flag or setting the `JULIA_PROJECT` environment variable.
3. Install the CLI script
```
julia> using SpectralUnmixing
julia> CLI.install()
```

If you would like to install a local version of the repository, first clone a local copy and navigate into the base SpectralUnmixing directory.  Then run:
```
julia --project='.' -e 'using Pkg; Pkg.activate(".");'
export JULIA_PROJECT=${PWD}
```

## Using the CLI script
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
