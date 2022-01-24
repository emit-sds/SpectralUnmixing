
using Test
    
include("../src/endmember_library.jl")

datafile = "data/basic_endmember_library.csv"
classname = "Class"

@info "Basic Loading"
lib = load_data!(SpectralLibrary(datafile, classname, 18))
@test isnothing(lib) == false

@time @test isnothing(lib.spectra) == false
@time @test length(size(lib.spectra)) == 2

@info "Filter By Class"
tmplib = deepcopy(lib)
filter_by_class!(tmplib)
@time @test length(lib.classes) == length(tmplib.classes)

tmplib =  load_data!(SpectralLibrary(datafile, classname, 18, [unique(lib.classes)[1]]))
filter_by_class!(tmplib)
@time @test length(lib.classes) != length(tmplib.classes)
@time @test length(unique(tmplib.classes)) == 1
tmplib=nothing

@info "Remove wavelength region inplace"

tmplib =  load_data!(SpectralLibrary(datafile, classname, 18))

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
scale_library!(tmplib)
@time @test lib.spectra == tmplib.spectra

scale_library!(tmplib, 2)
@time @test lib.spectra != tmplib.spectra

tmplib = deepcopy(lib)
brightness_normalize!(tmplib)
@time @test lib.spectra != tmplib.spectra



