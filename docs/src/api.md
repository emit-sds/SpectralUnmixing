# API Documentation
```@meta
CurrentModule = SpectralUnmixing
```

## Endmember Library Functions
```@docs
SpectralLibrary
SpectralLibrary(file_name::String, class_header_name::String)
load_data!
filter_by_class!
read_envi_wavelengths
interpolate_library_to_new_wavelengths!
remove_wavelength_region_inplace!
scale_library!
reduce_endmembers_nmf!
reduce_endmembers_kmeans!
reduce_endmembers_pca!
brightness_normalize!
split_library
prepare_combinations
prepare_options
scale_data
get_good_bands_mask
save_data
```

## Unmixing and Simulation Functions
```@docs
unmix_pixel
simulate_pixel
unmix_line
unmix_and_write_line
get_sma_permutation
results_from_mc
```

## Dataset Functions
```@docs
initiate_output_datasets
set_band_names
write_results
write_line_results
load_line
```

## Plotting Functions
```@docs
plot_mean_endmembers
plot_endmembers
plot_endmembers_individually
```

## Solver Functions
```@docs
opt_solve
dolsq
bvls
compute_kkt_optimality
```

# Utility Functions
```@docs
wl_index
nanargmax
nanargmin
```

# CLI
```@docs
CLI.install
```