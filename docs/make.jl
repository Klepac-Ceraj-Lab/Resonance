using Resonance
using Documenter

DocMeta.setdocmeta!(Resonance, :DocTestSetup, :(using Resonance); recursive=true)

makedocs(;
    modules=[Resonance],
    authors="Kevin Bonham, PhD <kbonham@wellesley.edu> and contributors",
    repo="https://github.com/kescobo/Resonance.jl/blob/{commit}{path}#{line}",
    sitename="Resonance.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kescobo.github.io/Resonance.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kescobo/Resonance.jl",
    devbranch="main",
)
