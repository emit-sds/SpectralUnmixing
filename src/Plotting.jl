
using Plots
using Statistics


function plot_mean_endmembers(endmember_library; output_name = "")
    plots = []
    for (_u, u) in enumerate(endmember_library.class_valid_keys)
        mean_spectra = mean(endmember_library.spectra[endmember_library.classes .== u,:],dims=1)[:]
        plot!(endmember_library.wavelengths, mean_spectra, label=u, xlim=[300,3200])
    end
    xlabel!("Wavelength [nm]")
    ylabel!("Reflectance")
    xticks!([500, 1000, 1500, 2000, 2500, 3000])
    plot!(dpi=200)
    if output_name != ""
        savefig(output_name)
    end
    plot!()
end

function plot_endmembers(endmember_library::SpectralLibrary; output_name::String = "")

    for (_u, u) in enumerate(endmember_library.class_valid_keys)
        if _u == 1
            plot(endmember_library.wavelengths, endmember_library.spectra[endmember_library.classes .== u,:]', lab="", xlim=[300,3200], color=palette(:tab10)[Mod(_u)],dpi=200)
        else
            plot!(endmember_library.wavelengths, endmember_library.spectra[endmember_library.classes .== u,:]', lab="",xlim=[300,3200], color=palette(:tab10)[Mod(_u)])
        end
    end
    xlabel!("Wavelenth [nm]")
    ylabel!("Reflectance")
    xticks!([500, 1000, 1500, 2000, 2500, 3000])
    for (_u, u) in enumerate(endmember_library.class_valid_keys)
        plot!([1:2],[0,0.3], color=palette(:tab10)[Mod(_u)], labels=u)
    end
    if output_name != ""
        savefig(output_name)
    end
    plot!()
end

function plot_endmembers_individually(endmember_library::SpectralLibrary; output_name::String = "")
    plots = []
    spectra = endmember_library.spectra
    classes = endmember_library.classes
    for (_u, u) in enumerate(endmember_library.class_valid_keys)
        sp = spectra[classes .== u,:]
        sp[broadcast(isnan,sp)] .= 0
        toplot = spectra[classes .== u,:]

        if size(toplot)[2] == sum(endmember_library.good_bands)
            push!(plots, plot(endmember_library.wavelengths[endmember_library.good_bands], toplot', title=u, xlabel="Wavelength [nm]", ylabel="Reflectance"))
        else 
            push!(plots, plot(endmember_library.wavelengths, toplot', title=u, xlabel="Wavelength [nm]", ylabel="Reflectance"))
        end
        xticks!([500, 1000, 1500, 2000, 2500])
    end
    plot(plots...,size=(1000,600),dpi=200)
    plot!(legend=:outertopright)
    if output_name != ""
        savefig(output_name)
    end
    plot!()
end