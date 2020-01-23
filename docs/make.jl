using Documenter, OhMyQSIM

makedocs(sitename="OhMyQSIM",
         pages=Any[
         "Getting Started"    =>    "index.md",
         "Algorithms"         =>    "algorithms.md"
         ])

deploydocs(
    repo = "github.com/nmoran/OhMyQSIM.jl.git",
)
