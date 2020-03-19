using Documenter, DiskArrayTools

makedocs(
    modules = [DiskArrayTools],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Fabian Gans",
    sitename = "DiskArrayTools.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/meggart/DiskArrayTools.jl.git",
    push_preview = true
)
