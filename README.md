# spectral-unmixing
A general, fast, flexible, and including spectral unmixing package.  Oriented towards VSWIR imaging spectroscopy data but applicable for different sensor types.  Includes options for different treatments of endmember library assemblages, including MESMA and bootstrapping (aka monte carlo) strategies.

This package has not been registered yet, but will be soon.  In the interim, it can be installed directly from this repo. From the Julia REPL or in a Jupyter cell, run

]activate .  # Optionally activate a project
]add https://github.com/emit-sds/SpectralUnmixing
Or, equivalently, from the command line with

julia --project=. -e 'import Pkg; Pkg.add("https://github.com/emit-sds/SpectralUnmixing")'
(Remove --project=. if not working in a project.)

Base usage:

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

