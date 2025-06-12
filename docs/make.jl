using Documenter, SpectralUnmixing

# copy readme to src/index.md
readme_title = """
<h1 align="center">
<br>
<a href="https://github.com/emit-sds/SpectralUnmixing.jl"><img src="docs/src/assets/logo.svg" alt="SpectralUnmixing.jl" width="200"></a>
<br>
SpectralUnmixing.jl
<br>
</h1>
"""
index_title = """
# SpectralUnmixing.jl
"""
readme_str = read(joinpath(@__DIR__, "..", "README.md"), String)
write(
    joinpath(@__DIR__, "src", "index.md"),
    replace(readme_str, readme_title => index_title),
)

makedocs(
    sitename="SpectralUnmixing.jl Documentation",
    modules=[SpectralUnmixing],
    pages=[
        "Introduction" => "index.md",
        "Examples" => "examples.md",
        "Contributing" => "contributing.md",
        "API" => "api.md",
    ]
)

deploydocs(
    repo="github.com/emit-sds/SpectralUnmixing.jl.git",
    devbranch="dev",
)