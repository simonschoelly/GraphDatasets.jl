using GraphDatasets
using Documenter

makedocs(;
    modules=[GraphDatasets],
    authors="Simon Schoelly <sischoel@gmail.com> and contributors",
    repo="https://github.com/simonschoelly/GraphDatasets.jl/blob/{commit}{path}#L{line}",
    sitename="GraphDatasets.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://simonschoelly.github.io/GraphDatasets.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/simonschoelly/GraphDatasets.jl",
)
