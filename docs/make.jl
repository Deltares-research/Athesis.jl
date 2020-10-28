using Athesis
using Documenter

makedocs(;
    modules=[Athesis],
    authors="mjr-deltares <martijn.russcher@deltares.nl> and contributors",
    repo="https://github.com/mjr-deltares/Athesis.jl/blob/{commit}{path}#L{line}",
    sitename="Athesis.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mjr-deltares.github.io/Athesis.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mjr-deltares/Athesis.jl",
)
