# spectral-unmixing
A general, fast, flexible, and including spectral unmixing package.  Oriented towards VSWIR imaging spectroscopy data but applicable for different sensor types.  Includes options for different treatments of endmember library assemblages, including MESMA and bootstrapping (aka monte carlo) strategies.

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

