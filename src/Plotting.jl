#  Copyright 2022 California Institute of Technology
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
# Author: Philip G. Brodrick, philip.g.brodrick@jpl.nasa.gov

using Plots
using Statistics

"""
    plot_mean_endmembers(endmember_library::SpectralLibrary; output_name::String="")

Plot the mean spectra of endmembers of each class from a spectral library.

- Saves file to `output_name` if provided, otherwise defaults to no saving.
"""
function plot_mean_endmembers(endmember_library::SpectralLibrary; output_name::String="")
    plots = []
    for (_u, u) in enumerate(endmember_library.class_valid_keys)
        mean_spectra =
            mean(endmember_library.spectra[endmember_library.classes.==u, :], dims=1)[:]
        plot!(endmember_library.wavelengths, mean_spectra, label=u, xlim=[300, 3200])
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

"""
    plot_endmembers(endmember_library::SpectralLibrary; output_name::String="")

Plot all endmember spectra from a spectral library.

- Endmember spectra are colored and labeled by class.
- Saves file to `output_name` if provided, otherwise defaults to no saving.
"""
function plot_endmembers(endmember_library::SpectralLibrary; output_name::String="")
    for (_u, u) in enumerate(endmember_library.class_valid_keys)
        if _u == 1
            plot(endmember_library.wavelengths,
                endmember_library.spectra[endmember_library.classes.==u, :]',
                lab="",
                xlim=[300, 3200],
                color=palette(:tab10)[Mod(_u)],
                dpi=200
            )
        else
            plot!(endmember_library.wavelengths,
                endmember_library.spectra[endmember_library.classes.==u, :]',
                lab="",
                xlim=[300, 3200],
                color=palette(:tab10)[Mod(_u)]
            )
        end
    end
    xlabel!("Wavelenth [nm]")
    ylabel!("Reflectance")
    xticks!([500, 1000, 1500, 2000, 2500, 3000])
    for (_u, u) in enumerate(endmember_library.class_valid_keys)
        plot!([1:2], [0, 0.3], color=palette(:tab10)[Mod(_u)], labels=u)
    end
    if output_name != ""
        savefig(output_name)
    end
    plot!()
end

"""
    plot_endmembers_individually(endmember_library::SpectralLibrary;
    output_name::String="", legend_on::Bool=false)

Plot all endmember spectra individually from a spectral library, grouped by class.

- Create seperate plots for each class, legends are disabled by default.
- Saves file to `output_name` if provided, otherwise defaults to no saving.
"""
function plot_endmembers_individually(endmember_library::SpectralLibrary;
    output_name::String="", legend_on::Bool=false)

    plots = []
    spectra = endmember_library.spectra
    classes = endmember_library.classes
    for (_u, u) in enumerate(endmember_library.class_valid_keys)
        sp = spectra[classes.==u, :]
        sp[broadcast(isnan, sp)] .= 0
        toplot = spectra[classes.==u, :]

        if size(toplot)[2] == sum(endmember_library.good_bands)
            push!(plots,
                plot(endmember_library.wavelengths[endmember_library.good_bands],
                    toplot',
                    title=u,
                    xlabel="Wavelength [nm]",
                    ylabel="Reflectance"
                )
            )
        else
            push!(plots,
                plot(endmember_library.wavelengths,
                    toplot',
                    title=u,
                    xlabel="Wavelength [nm]",
                    ylabel="Reflectance"
                )
            )
        end
        xticks!([500, 1000, 1500, 2000, 2500])
    end
    plot(plots..., size=(1000, 600), dpi=200)
    if legend_on
        plot!(legend=:outertopright)
    else
        plot!(legend=false)
    end
    if output_name != ""
        savefig(output_name)
    end
    plot!()
end