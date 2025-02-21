<h1 align="center">
SpectralUnmixing
</h1>

[![version](https://github.com/emit-sds/SpectralUnmixing/actions/workflows/unit-tests.yml/badge.svg)](https://github.com/emit-sds/SpectralUnmixing/actions/workflows/unit-tests.yml/)
[![](https://img.shields.io/github/license/emit-sds/SpectralUnmixing)](https://github.com/emit-sds/SpectralUnmixing/blob/master/LICENSE)
[![](https://img.shields.io/badge/docs-latest-blue)](https://emit-sds.github.io/SpectralUnmixing/dev/)



A general, fast, flexible, and including spectral unmixing package.  Oriented towards VSWIR imaging spectroscopy data but applicable for different sensor types.  Includes options for different treatments of endmember library assemblages, including MESMA and bootstrapping (aka monte carlo) strategies.

## Installation
This package is registered and may be added using:
```
julia 'using Pkg; Pkg.add("SpectralUnmixing")'
```
Remember to use the --project flag or to set the JULIA_PROJECT environment variable to activate the appropriate environment.

If you would like to install a local version of the repository, first pull a local copy and navigate into the base SpectralUnmixing directory.  Then run:

```
julia --project='.' -e 'using Pkg; Pkg.activate(".");'
export JULIA_PROJECT=${PWD}
```

## Using the script
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

