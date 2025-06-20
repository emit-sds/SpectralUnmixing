module CLI

export install

"""
    install(; dest_dir::String=joinpath(DEPOT_PATH[1], "bin"),
        sym_name::String="unmix.jl")

Create a symlink to the `unmix.jl` CLI script.

# Arguments
- `dest_dir::String`: The directory where the symlink will be created (default: `~/.julia/bin`).
- `sym_name::String`: The name of the symlink to be created (default: `unmix.jl`).
"""
function install(; dest_dir::String=joinpath(DEPOT_PATH[1], "bin"),
    sym_name::String="unmix.jl")

    pkgpath = dirname(Base.find_package("SpectralUnmixing"))
    src = joinpath(pkgpath, "..", "unmix.jl") |> normpath
    dest = joinpath(dest_dir, sym_name)

    mkpath(dest_dir)
    if isfile(dest) || islink(dest)
        rm(dest)
    end
    symlink(src, dest)
    println("Symlink created at $dest")
end

end # module