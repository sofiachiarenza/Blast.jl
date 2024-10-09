using Documenter
using Blast

ENV["GKSwstype"] = "100"

push!(LOAD_PATH,"../src/")

makedocs(
    modules = [Blast],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true",
    sidebar_sitename=true),
    sitename = "Blast.jl",
    authors  = "Sofia Chiarenza",
    pages = [
        "Home" => "index.md"
        "The algorithm" => "alg.md"
        "API" => "api.md"
    ]
)

deploydocs(
    repo = "github.com/sofiachiarenza/Blast.jl.git",
    devbranch = "develop"
)