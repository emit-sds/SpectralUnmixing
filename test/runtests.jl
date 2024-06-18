
using Test
    
using SpectralUnmixing

datafile = "data/basic_endmember_library.csv"
classname = "Class"

@info "Basic Loading"
lib = load_data!(SpectralLibrary(datafile, classname))
@test isnothing(lib) == false

@time @test isnothing(lib.spectra) == false
@time @test length(size(lib.spectra)) == 2

@info "Filter By Class"
tmplib = deepcopy(lib)
filter_by_class!(tmplib)
@time @test length(lib.classes) == length(tmplib.classes)

tmplib =  load_data!(SpectralLibrary(datafile, classname, 2, 0, [unique(lib.classes)[1]]))
filter_by_class!(tmplib)
@time @test length(lib.classes) != length(tmplib.classes)
@time @test length(unique(tmplib.classes)) == 1
tmplib=nothing

@info "Remove wavelength region inplace"

tmplib =  load_data!(SpectralLibrary(datafile, classname))

remove_wavelength_region_inplace!(tmplib)
@time @test size(lib.spectra) != size(tmplib.spectra)

@info "Wavelength interpolation"
tmplib = deepcopy(lib)
new_wl = 400 .+ Vector{Float64}(1:200)*10.
interpolate_library_to_new_wavelengths!(tmplib, new_wl)
@time @test size(tmplib.spectra)[2] == length(new_wl)
@time @test size(lib.spectra) != size(tmplib.spectra)


@info "Scale Library"
tmplib = deepcopy(lib)
@time scale_library!(tmplib)
@test lib.spectra == tmplib.spectra

@time scale_library!(tmplib, 2)
@test lib.spectra != tmplib.spectra

tmplib = deepcopy(lib)
@time brightness_normalize!(tmplib)
@test lib.spectra != tmplib.spectra


@info "Simulate Mixture - class-even type"
remove_wavelength_region_inplace!(lib)
class_idx = []
for uc in lib.class_valid_keys
    push!(class_idx, (1:size(lib.classes)[1])[lib.classes .== uc])
end
num_components = 3
seed = 13
@time simulated_rfl, true_complete_fractions, true_mixture_fractions = simulate_pixel(lib, num_components, "class-even", seed)
@test sum(isnothing.(simulated_rfl)) + sum(isnan.(simulated_rfl)) == 0
println(lib.spectra' * true_complete_fractions)
println(simulated_rfl[:])
@test sum(lib.spectra' * true_complete_fractions .≈ simulated_rfl[:]) == length(simulated_rfl[:])


for mode in ["sma", "sma-best", "mesma", "mesma-best"]
    n_mc = 2
    mode = "sma-best"
    num_endmembers=[3]
    normalization="brightness"
    optimization="bvls"

    max_combinations=10000
    combination_type="class_even"

    unmixing_library = lib

    options = []

    class_idx = []
    if combination_type == "class-even"
        for uc in library.class_valid_keys
            push!(class_idx, (1:size(unmixing_library.classes)[1])[unmixing_library.classes .== uc])
        end
    end

    @info "Unmix Pixel - Mode: " * mode
    @time mr, mv, cfr, cfv = unmix_pixel(unmixing_library, simulated_rfl, nothing, class_idx, options, mode, n_mc, 
            num_endmembers, normalization, optimization, max_combinations, combination_type)
    @test sum(mr[1:end-1]) ≈ 1
    @test size(mr) == size(mv)
    @test sum(cfr[1:end-1]) ≈ 1
    @test size(cfr) == size(cfv)
end

