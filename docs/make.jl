using Resonance
using Documenter

DocMeta.setdocmeta!(Resonance, :DocTestSetup, :(using Resonance); recursive=true)

makedocs(;
    modules=[Resonance],
    authors="Kevin Bonham, PhD <kbonham@wellesley.edu> and contributors",
    repo="https://github.com/Klepac-Ceraj-lab/Resonance/blob/{commit}{path}#{line}",
    sitename="Resonance",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Klepac-Ceraj-lab.github.io/Resonance",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Setup" => "setup.md",
        "Figures" => [
            "Figure 1" => "Figures/figure1.md",
            "Figure 2" => "Figures/figure1.md",
            "Figure 3" => "Figures/figure1.md",
        ]
    ],
)

deploydocs(;
    repo="github.com/Klepac-Ceraj-lab/Resonance",
    devbranch="main",
)
