using Documenter, QGDipoles, RecipesBase

format = Documenter.HTML(
    assets = ["assets/favicon.ico"],
    collapselevel = 2,
    prettyurls = get(ENV, "CI", nothing) == "true",
)

makedocs(
    authors = "Matthew N. Crowe",
    modules = [QGDipoles],
    sitename = "QGDipoles.jl",
    format = format,
    pages = Any[
        "Home"=>"index.md",
        "Installation"=>"Installation.md",
        "Methodology"=>"Methodology.md",
        "Examples"=>"Examples.md",
        "List of Functions"=>"Functions.md",
    ],
)

deploydocs(repo = "github.com/mncrowe/QGDipoles.jl.git")
