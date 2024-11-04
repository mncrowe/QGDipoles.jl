using Documenter, QGDipoles

makedocs(
  authors = "Matthew N. Crowe",
  modules = [QGDipoles],
 sitename = "QGDipoles.jl",
   format = Documenter.HTML(assets = ["assets/favicon.ico"]),
    pages = Any[
                "Home" => "index.md",
		"Examples" => "Examples.md",
		"List of Functions" => "Functions.md"
		]
	)

deploydocs(
    repo = "github.com/mncrowe/QGDipoles.jl.git",
)
