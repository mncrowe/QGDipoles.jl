using Documenter, QGDipoles

makedocs(
  authors = "Matthew N. Crowe",
 sitename = "QGDipoles.jl",
    pages = Any[
                "Home" => "index.md",
		"List of Functions" => "Functions.md"
		]
	)

deploydocs(
    repo = "github.com/mncrowe/QGDipoles.jl.git",
)
