# Figure3 - RandomForest prediction of cogScores from DEM, Taxonomic and Functional Profiles

```julia
using Chain
using CSV
using Statistics
using Random
using DataFrames
using StableRNGs
using MLJ
using CairoMakie
using ColorSchemes
using DecisionTree
using JLD2
using Resonance
using Distributions
using CategoricalArrays
ml_rng = StableRNG(0)
```

## Loading the pretrained models

```julia
## 1. Metadata
mdata = Resonance.load(Metadata())
seqs = subset(mdata, "sample"=> ByRow(!ismissing)) 
transform!(seqs, "sample"=> ByRow(String)=> "sample")
seqs.edfloat = map(x-> ismissing(x) ? missing : Float64(levelcode(x)), seqs.education)

## 2. Taxonomic Profiles
taxa = Resonance.load(TaxonomicProfiles(); timepoint_metadata = seqs) # this can take a bit
species = filter(f-> taxrank(f) == :species, taxa)

RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree
# concurrent cogScore regression from taxonomic profiles
regression_currentCogScores_00to06mo_onlydemo = JLD2.load(
    modelfiles("regression_currentCogScores_00to06mo_onlydemo.jld"),
    "regression_currentCogScores_00to06mo_onlydemo"
);
regression_currentCogScores_00to06mo_onlytaxa = JLD2.load(
    modelfiles("regression_currentCogScores_00to06mo_onlytaxa.jld"),
    "regression_currentCogScores_00to06mo_onlytaxa"
);
regression_currentCogScores_00to06mo_demoplustaxa = JLD2.load(
    modelfiles("regression_currentCogScores_00to06mo_demoplustaxa.jld"),
    "regression_currentCogScores_00to06mo_demoplustaxa"
);
regression_currentCogScores_00to06mo_onlyecs = JLD2.load(
    modelfiles("regression_currentCogScores_00to06mo_onlyecs.jld"),
    "regression_currentCogScores_00to06mo_onlyecs"
);
regression_currentCogScores_00to06mo_demoplusecs = JLD2.load(
    modelfiles("regression_currentCogScores_00to06mo_demoplusecs.jld"),
    "regression_currentCogScores_00to06mo_demoplusecs"
);
regression_currentCogScores_18to120mo_onlydemo = JLD2.load(
    modelfiles("regression_currentCogScores_18to120mo_onlydemo.jld"),
    "regression_currentCogScores_18to120mo_onlydemo"
);
regression_currentCogScores_18to120mo_onlytaxa = JLD2.load(
    modelfiles("regression_currentCogScores_18to120mo_onlytaxa.jld"),
    "regression_currentCogScores_18to120mo_onlytaxa"
);
regression_currentCogScores_18to120mo_demoplustaxa = JLD2.load(
    modelfiles("regression_currentCogScores_18to120mo_demoplustaxa.jld"),
    "regression_currentCogScores_18to120mo_demoplustaxa"
);
regression_currentCogScores_18to120mo_onlyecs = JLD2.load(
    modelfiles("regression_currentCogScores_18to120mo_onlyecs.jld"),
    "regression_currentCogScores_18to120mo_onlyecs"
);
regression_currentCogScores_18to120mo_demoplusecs = JLD2.load(
    modelfiles("regression_currentCogScores_18to120mo_demoplusecs.jld"),
    "regression_currentCogScores_18to120mo_demoplusecs"
);
regression_futureCogScores_demoplustaxa = JLD2.load(
    modelfiles("regression_futureCogScores_demoplustaxa.jld"),
    "regression_futureCogScores_demoplustaxa"
);
regression_futureCogScores_onlytaxa = JLD2.load(
    modelfiles("regression_futureCogScores_onlytaxa.jld"),
    "regression_futureCogScores_onlytaxa"
);
```

## Initializing Figure 3

```julia
figure = Figure(resolution = (1920, 1536));

ABCD_subfig = GridLayout(figure[1,1:2]; alignmode=Outside())
AB_subfig = GridLayout(ABCD_subfig[1,1])
CD_subfig = GridLayout(ABCD_subfig[1,2])
E_subfig = GridLayout(figure[2,1], alignmode=Outside())
F_subfig = GridLayout(figure[2,2], alignmode=Outside())

# colsize!(figure.layout, 1, Relative(0.3))
# colsize!(figure.layout, 2, Relative(0.7))
# 
# rowsize!(figure.layout, 1, Relative(0.6))
# rowsize!(figure.layout, 2, Relative(0.2))
# rowsize!(figure.layout, 3, Relative(0.2))
```

### Plot panels A and B - Comparative Importances

```julia
axA = Axis(
    AB_subfig[1, 1];
    xlabel = "-log(p) for LM coefficients",
    ylabel = "relative importance",
    title = "over 18mo",
)
t10 = ColorSchemes.tableau_10
let plot_colorset = [t10[1], t10[3], t10[7], (:white, 0.)]
# plot_comparative_lmvsrf_scatterplots!(axA, regression_currentCogScores_00to06mo_onlytaxa, tablefiles("figure2", "lms_species_00to06.csv"); plot_colorset = plot_colorset, strokewidth=1)
# plot_comparative_lmvsrf_scatterplots!(axB, regression_currentCogScores_18to120mo_onlytaxa, tablefiles("figure2", "lms_species_18to120.csv"); plot_colorset = plot_colorset, strokewidth=1)
    plot_comparative_lmvsrf_scatterplots!(axA, regression_currentCogScores_18to120mo_onlytaxa,
        tablefiles("figure2/lms_species_18to120.csv");
        plot_colorset, strokewidth=1
    ) 

    Legend(
        AB_subfig[1, 2],
        [
            MarkerElement(; marker=:circle, color=plot_colorset[1], strokewidth=1),
            MarkerElement(; marker=:circle, color=plot_colorset[2], strokewidth=1),
            MarkerElement(; marker=:circle, color=plot_colorset[3], strokewidth=1),
            MarkerElement(; marker=:circle, color=plot_colorset[4], strokewidth=1),
        ],
        [
            "> 60% ranked\nimportance",
            "q < 0.2 in LM",
            "Both",
            "None"
        ];
        tellheight = false,
        tellwidth = true,
    )
end

axB = GridLayout(AB_subfig[2,1])

let plot_colorset = [t10[1], t10[6], t10[5], (:white, 0.)]
    Resonance.plot_comparative_rfvsrf_scatterplots!(axB,
        regression_currentCogScores_18to120mo_onlytaxa,
        regression_futureCogScores_onlytaxa;
        plot_colorset, strokewidth=1
    ) 

    Legend(
        AB_subfig[2, 2],
        [
            MarkerElement(; marker=:circle, color=plot_colorset[1], strokewidth=1),
            MarkerElement(; marker=:circle, color=plot_colorset[2], strokewidth=1),
            MarkerElement(; marker=:circle, color=plot_colorset[3], strokewidth=1),
            MarkerElement(; marker=:circle, color=plot_colorset[4], strokewidth=1),
        ],
        [
            "over 18mo",
            "future",
            "Both",
            "None"
        ];
        tellheight = false,
        tellwidth = true,
    )
end

```

```julia
nbars_toplot = 15

joined_importances_00to06 = compute_joined_importances(
    regression_currentCogScores_00to06mo_onlytaxa,
    regression_currentCogScores_00to06mo_demoplustaxa;
    imp_fun = weighted_hpimportances
)
joined_importances_18to120 = compute_joined_importances(
    regression_currentCogScores_18to120mo_onlytaxa,
    regression_currentCogScores_18to120mo_demoplustaxa;
    imp_fun = weighted_hpimportances
)

joined_importances_future = compute_joined_importances(
    regression_futureCogScores_onlytaxa,
    regression_futureCogScores_demoplustaxa;
    imp_fun = weighted_hpimportances
)


joined_importances_future = compute_joined_importances(
    regression_futureCogScores_onlytaxa,
    regression_futureCogScores_demoplustaxa;
    imp_fun = weighted_hpimportances
)

axC = Axis(
    CD_subfig[1, 1];
    xlabel = "relative importance",
    yticks = (reverse(collect(1:nbars_toplot)), format_species_labels(joined_importances_18to120.variable[1:nbars_toplot])),
    ylabel = "Feature",
    title = "18 to 120 months",
    # yticklabelrotation= -pi/2
)

plot_comparativedemo_importance_barplots!(axC, joined_importances_18to120; n_rows = nbars_toplot)

axD = Axis(
    CD_subfig[1, 2];
    xlabel = "relative importance",
    yticks = (reverse(collect(1:nbars_toplot)), format_species_labels(joined_importances_future.variable[1:nbars_toplot])),
    ylabel = "Feature",
    title = "Future",
    # yticklabelrotation= -pi/2
)

plot_comparativedemo_importance_barplots!(axD, joined_importances_future; n_rows = nbars_toplot)

Legend(
    CD_subfig[1, 1], [MarkerElement(; marker=:rect, color=:gray), MarkerElement(; marker=:star8, color=:red)], ["Microbiome alone", "Microbiome + demographics"];
    tellheight = false,
    tellwidth = false,
    margin = (10, 10, 10, 10),
    halign = :right, valign = :bottom, orientation = :vertical
)
```

### Plot panel E - Deep dives on taxa
```julia
e1 = GridLayout(E_subfig[1,1])
e2 = GridLayout(E_subfig[1,2])
e3 = GridLayout(E_subfig[2,1])
e4 = GridLayout(E_subfig[2,2])
plot_taxon_deepdive!(e1, species, :filter_18to120, "Blautia_wexlerae";)
plot_taxon_deepdive!(e2, species, :filter_18to120, "Bifidobacterium_pseudocatenulatum";)
plot_taxon_deepdive!(e3, species, :filter_18to120, "Faecalibacterium_prausnitzii";)
plot_taxon_deepdive!(e4, species, :filter_18to120, "Eubacterium_eligens";)


f1 = GridLayout(F_subfig[1,1])
f2 = GridLayout(F_subfig[1,2])
f3 = GridLayout(F_subfig[2,1])
f4 = GridLayout(F_subfig[2,2])
plot_future_deepdive!(f1, mdata, species, "Klebsiella_pneumoniae";)
plot_future_deepdive!(f2, mdata, species, "Klebsiella_variicola";)
plot_future_deepdive!(f3, mdata, species, "Ruminococcus_gnavus";)
plot_future_deepdive!(f4, mdata, species, "Bacteroides_ovatus";)
```

### Panel labels
```julia
Label(AB_subfig[1, 1, TopLeft()], "A", fontsize = 26,font = :bold, padding = (0, 5, 5, 0), halign = :right)
Label(AB_subfig[2, 1, TopLeft()], "B", fontsize = 26,font = :bold, padding = (0, 5, 5, 0), halign = :right)
Label(CD_subfig[1, 1, TopLeft()], "C", fontsize = 26,font = :bold, padding = (0, 5, 5, 0), halign = :right)
Label(CD_subfig[1, 2, TopLeft()], "D", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
Label(E_subfig[1, 1, TopLeft()], "E", fontsize = 26,font = :bold, padding = (0, 5, 5, 0), halign = :right)
Label(F_subfig[1, 1, TopLeft()], "F", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
# 
# Label(E_subfig[1:2, 1, Left()], "0 to 6 months", fontsize = 20, font = :bold, padding = (0, 80, 0, 0), halign = :center, valign = :center, rotation = pi/2)
# Label(F_subfig[1:2, 1, Left()], "18 to 120 months", fontsize = 20, font = :bold, padding = (0, 80, 0, 0), halign = :center, valign = :center, rotation = pi/2)
#figure
# Label(AB_subfig[1, 1, TopLeft()], "A", fontsize = 26,font = :bold, padding = (0, 5, 5, 0), halign = :right)
# Label(AB_subfig[2, 1, TopLeft()], "B", fontsize = 26,font = :bold, padding = (0, 5, 5, 0), halign = :right)
# Label(CD_subfig[1, 1, TopLeft()], "C", fontsize = 26,font = :bold, padding = (0, 5, 5, 0), halign = :right)
# Label(CD_subfig[1, 2, TopLeft()], "D", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
# Label(E_subfig[1, 1, TopLeft()], "E", fontsize = 26,font = :bold, padding = (0, 5, 5, 0), halign = :right)
# Label(E_subfig[1, 5, TopLeft()], "F", fontsize = 26,font = :bold, padding = (0, 5, 5, 0), halign = :right)
# Label(F_subfig[1, 5, TopLeft()], "G", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
# 
# Label(E_subfig[1:2, 1, Left()], "0 to 6 months", fontsize = 20, font = :bold, padding = (0, 80, 0, 0), halign = :center, valign = :center, rotation = pi/2)
# Label(F_subfig[1:2, 1, Left()], "18 to 120 months", fontsize = 20, font = :bold, padding = (0, 80, 0, 0), halign = :center, valign = :center, rotation = pi/2)
colsize!(ABCD_subfig, 1, Relative(1/3))
figure
```

```julia
save(figurefiles("Figure3.png"), figure)
save(figurefiles("Figure3.svg"), figure)
save("manuscript/assets/Figure3.png", figure)
```

# Supplement to Figure 3
```julia
## 1. Building the supplementary-style tables
supptblA = singlemodel_importances_suppltable(regression_currentCogScores_00to06mo_onlytaxa)
supptblB = singlemodel_importances_suppltable(regression_currentCogScores_18to120mo_onlytaxa)

## 2. Building the Figure
figure = Figure(resolution = (1200, 1200));

this_barcolor = :lightblue
this_curvecolor = :orange

plot_importances_pareto!(figure[1,1], supptblA, "Pareto Plot - 0 to 6 months"; barcolor = this_barcolor, curvecolor = this_curvecolor)
plot_importances_pareto!(figure[2,1], supptblB[1:70, :], "Pareto Plot - 18 to 120 months"; barcolor = this_barcolor, curvecolor = this_curvecolor)

Legend(
    figure[3, 1],
    [
        PolyElement(; color = this_barcolor, strokewidth = 1, strokecolor = :black),
        LineElement(; color = this_curvecolor, linewidth = 5)
    ] , [
        "Individual relative Importance",
        "Cumulative relative importance"

    ],
    tellheight = true,
    tellwidth = false,
    margin = (10, 10, 10, 10),
    orientation = :horizontal
)

## Saving
save(figurefiles("Supp_Figure4.svg"), figure)
save("manuscript/assets/Supp_Figure4.png", figure)
```
