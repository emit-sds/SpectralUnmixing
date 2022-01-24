


@everywhere begin
    using ArchGDAL
    using EllipsisNotation
    using DelimitedFiles
    using Logging
    using Statistics
    using PyCall
    using Distributed
    using Printf
    using LinearAlgebra
    using Combinatorics
    using Random
    include("solvers.jl")
    include("endmember_library.jl")
    include("datasets.jl")

    function wl_index(wavelengths::Array{Float64}, target)
        argmin(abs.(wavelengths .- target))
    end

    function scale_data(refl::Array{Float64}, wavelengths::Array{Float64}, criteria::String)

        if criteria == "none"
            return refl
        elseif criteria == "brightness"
            bad_regions_wl = [[1300,1500],[1800,2000]]
            good_bands = convert(Array{Bool}, ones(length(wavelengths)))
            for br in bad_regions_wl
                good_bands[wl_index(wavelengths, br[1]):wl_index(wavelengths, br[2])] .= false
            end
            norm = sqrt.(mean(refl[:,good_bands].^2, dims=2))
        else
            try
                target_wl = parse(Float64,criteria)
                norm = refl[:,wl_index(wavelengths, target_wl)] ./ 0.5
            catch e
                throw(ArgumentError(string("normalization must be [none, brightness, or a specific wavelength].  Provided:", criteria)))
            end
        end

        return refl ./ norm
    end

    function unmix_line(line::Int64, reflectance_file::String, mode::String, refl_nodata::Float64,
                        refl_scale::Float64, normalization::String, library::SpectralLibrary,
                        reflectance_uncertainty_file::String = "", n_mc::Int64 = 1,
                        combination_type::String = "all", num_endmembers::Vector{Int64} = [2,3],
                        max_combinations::Int64 = -1, optimization="bvls")

        Random.seed!(13)
        println(line)
        img_dat, unc_dat, good_data = load_line(reflectance_file, reflectance_uncertainty_file, line, library.good_bands, refl_nodata)
        mixture_results = fill(-9999.0, sum(good_data), size(library.class_valid_keys)[1] + 1)
        if n_mc > 1
            mixture_results_std = fill(-9999.0, sum(good_data), size(library.class_valid_keys)[1] + 1)
        else
            mixture_results_std = nothing
        end

        if isnothing(img_dat)
            return line, nothing, good_data, nothing, nothing
        end
        scale_data(img_dat, library.wavelengths[library.good_bands], normalization)
        img_dat = img_dat ./ refl_scale

        if combination_type == "class-even"
            class_idx = []
            for uc in library.class_valid_keys
                push!(class_idx, (1:size(library.classes)[1])[library.classes .== uc])
            end
        end

        # Prepare combinations if relevant
        if mode == "mesma"
            if combination_type == "class-even"
                options = collect(Iterators.product(class_idx...))[:]
            elseif combination_type == "all"
                options = []
                for num in num_endmembers
                    combo = [c for c in combinations(1:length(library.classes), num)]
                    push!(options,combo...)
                end
            else
                error("Invalid combiation string")
            end
        end


        # Solve complete fraction set (based on full library deck)
        complete_fractions = zeros(size(img_dat)[1], size(library.spectra)[1] + 1)
        complete_fractions_std = zeros(size(img_dat)[1], size(library.spectra)[1] + 1)
        for _i in 1:size(img_dat)[1] # Pixel loop

            mc_comp_frac = unmix_pixel(library, img_dat[_i:_i,:], unc_dat[_i:_i,:], class_idx, options,
                                n_mc, num_endmembers, normalization, optimization, max_combinations, combination_type)

            complete_fractions[_i,:] = mean(mc_comp_frac,dims=1)
            complete_fractions_std[_i,:] = std(mc_comp_frac,dims=1)

            # Aggregate results from per-library to per-unique-class
            for (_class, cl) in enumerate(library.class_valid_keys)
                mixture_results[_i, _class] = sum(complete_fractions[_i,1:end-1][cl .== library.classes])
            end
            mixture_results[_i, end] = complete_fractions[_i,end]

            #Aggregate uncertainty if relevant
            if n_mc > 1
                for (_class, cl) in enumerate(library.class_valid_keys)
                    mixture_results_std[_i, _class] = std(sum(mc_comp_frac[:,1:end-1][:,cl .== library.classes], dims=2))
                end
                mixture_results_std[_i, end] = std(mc_comp_frac[:,end])
            end

        end

        return line, mixture_results, good_data, mixture_results_std, complete_fractions


    end



    function unmix_pixel(library, img_dat, unc_dat, class_idx, options, n_mc, num_endmembers, normalization, optimization, max_combinations, combination_type)

        mc_comp_frac = zeros(n_mc, size(library.spectra)[1]+1)
        for mc in 1:n_mc #monte carlo loop
            Random.seed!(mc)

            d = img_dat
            if isnothing(unc_dat) == false
                d += (rand(size(d)) .* 2 .- 1) .* unc_dat
            end


            if mode == "sma"
                if num_endmembers[1] != -1
                    if combination_type == "class-even"

                        perm_class_idx = []
                        for class_subset in class_idx
                            push!(perm_class_idx, Random.shuffle(class_subset))
                        end

                        perm = []
                        selector = 1
                        while selector <= num_endmembers[1]
                            _p = mod(selector, length(perm_class_idx)) + 1
                            push!(perm, perm_class_idx[_p][1])
                            deleteat!(perm_class_idx[_p],1)

                            if length(perm_class_idx[_p]) == 0
                                deleteat!(perm_class_idx,_p)
                            end
                            selector += 1
                        end

                    else
                        perm = randperm(size(library.spectra)[1])[1:num_endmembers[1]]
                    end

                    G = library.spectra[perm,:]
                else
                    perm = convert(Vector{Int64},1:size(library.spectra)[1])
                    G = library.spectra
                end

                G = scale_data(G, library.wavelengths[library.good_bands], normalization)'

                x0 = dolsq(G, d')
                x0 = x0[:]
                res = nothing
                if optimization == "bvls"
                    res, cost = bvls(G, d[:], x0, zeros(size(x0)), ones(size(x0)), 1e-3, 100, 1)
                elseif optimization == "ldsqp"
                    res, cost = opt_solve(G, d[:], x0, 0, 1 )
                elseif optimization == "inverse"
                    res = x0
                end
                mc_comp_frac[mc, perm] = res

            elseif occursin("mesma", mode)
                solutions = []
                costs = zeros(size(options)[1]).+1e12

                if max_combinations != -1 && length(options) > max_combinations
                    perm = randperm(length(options))[1:max_combinations]
                else
                    perm = 1:length(options)
                end

                for (_comb, comb) in enumerate(options[perm])
                    comb = [c for c in comb]
                    #G = hcat(library.spectra[comb,:], ones(size(library.spectra[comb,:])[1],1))
                    G = scale_data(library.spectra[comb,:], library.wavelengths[library.good_bands], normalization)'

                    x0 = dolsq(G, d')
                    ls = nothing
                    if optimization == "bvls"
                        ls, lc = bvls(G, d[:], x0, zeros(size(x0)), ones(size(x0)), 1e-3, 10, 1)
                        costs[_comb] = lc
                    elseif optimization == "ldsqp"
                        ls, lc = opt_solve(G, d[:], x0, 0, 1)
                        costs[_comb] = lc
                    elseif optimization == "inverse"
                        ls = x0
                        r = G * x0 - d[:]
                        costs[_comb] = dot(r,r)
                    end
                    push!(solutions,ls)

                end
                best = argmin(costs)

                mc_comp_frac[mc, [ind for ind in options[perm][best]]] = solutions[best]
            else
                error("Invalid mode provided")
            end
        end

        # Calculate the sum of values (inverse of shade), and then normalize
        mc_comp_frac[mc_comp_frac .< 0] .= 0
        mc_comp_frac[:,end] = sum(mc_comp_frac,dims=2)
        mc_comp_frac[:,1:end-1] = mc_comp_frac[:,1:end-1] ./ mc_comp_frac[:,end]

        return mc_comp_frac

    end

end
