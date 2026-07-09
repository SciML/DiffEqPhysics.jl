using Documenter, DiffEqPhysics

makedocs(
    sitename = "DiffEqPhysics.jl",
    authors = "Chris Rackauckas",
    modules = [DiffEqPhysics],
    clean = true,
    doctest = false,
    linkcheck = false,
    warnonly = false,
    format = Documenter.HTML(
        canonical = "https://docs.sciml.ai/DiffEqPhysics/stable/"
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md"
    ]
)

deploydocs(
    repo = "github.com/SciML/DiffEqPhysics.jl.git";
    push_preview = true
)
