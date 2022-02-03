<h1 align="center">
SpectralUnmixing
</h1>
A general, fast, flexible, and including spectral unmixing package.  Oriented towards VSWIR imaging spectroscopy data but applicable for different sensor types.  Includes options for different treatments of endmember library assemblages, including MESMA and bootstrapping (aka monte carlo) strategies.

## Installation
This package has not been registered yet, but will be soon.  In the interim, after cloning and navigating into the repository, it can be installed from the Julia REPL, by running

```
julia
]activate .
]add https://github.com/kmsquire/ArgParse2.jl
Pkg.instantiate()
exit()
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

