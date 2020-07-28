using Porteau
using Documenter

makedocs(;
    modules=[Porteau],
    authors="Deltares",
    repo="https://github.com/openearth/Porteau.jl/blob/{commit}{path}#L{line}",
    sitename="Porteau.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://openearth.github.io/Porteau.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/openearth/Porteau.jl",
)
