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
ml_rng = StableRNG(0)
```

## Loading the pretrained models

```julia
RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree
# concurrent cogScore reression from taxonomic profiles
regression_currentCogScores_00to06mo_onlydemo = JLD2.load(modelfiles("manuscript", "regression_currentCogScores_00to06mo_onlydemo.jld"))["regression_currentCogScores_00to06mo_onlydemo"]
regression_currentCogScores_00to06mo_onlytaxa = JLD2.load(modelfiles("manuscript", "regression_currentCogScores_00to06mo_onlytaxa.jld"))["regression_currentCogScores_00to06mo_onlytaxa"]
regression_currentCogScores_00to06mo_demoplustaxa = JLD2.load(modelfiles("manuscript", "regression_currentCogScores_00to06mo_demoplustaxa.jld"))["regression_currentCogScores_00to06mo_demoplustaxa"]
regression_currentCogScores_00to06mo_onlyecs = JLD2.load(modelfiles("manuscript", "regression_currentCogScores_00to06mo_onlyecs.jld"))["regression_currentCogScores_00to06mo_onlyecs"]
regression_currentCogScores_00to06mo_demoplusecs = JLD2.load(modelfiles("manuscript", "regression_currentCogScores_00to06mo_demoplusecs.jld"))["regression_currentCogScores_00to06mo_demoplusecs"]
regression_currentCogScores_18to120mo_onlydemo = JLD2.load(modelfiles("manuscript", "regression_currentCogScores_18to120mo_onlydemo.jld"))["regression_currentCogScores_18to120mo_onlydemo"]
regression_currentCogScores_18to120mo_onlytaxa = JLD2.load(modelfiles("manuscript", "regression_currentCogScores_18to120mo_onlytaxa.jld"))["regression_currentCogScores_18to120mo_onlytaxa"]
regression_currentCogScores_18to120mo_demoplustaxa = JLD2.load(modelfiles("manuscript", "regression_currentCogScores_18to120mo_demoplustaxa.jld"))["regression_currentCogScores_18to120mo_demoplustaxa"]
regression_currentCogScores_18to120mo_onlyecs = JLD2.load(modelfiles("manuscript", "regression_currentCogScores_18to120mo_onlyecs.jld"))["regression_currentCogScores_18to120mo_onlyecs"]
regression_currentCogScores_18to120mo_demoplusecs = JLD2.load(modelfiles("manuscript", "regression_currentCogScores_18to120mo_demoplusecs.jld"))["regression_currentCogScores_18to120mo_demoplusecs"]
```

## Initializing Figure 3

```julia
figure = Figure(resolution = (1920, 1536))

AB_subfig = GridLayout(figure[1,1], alignmode=Outside())
CD_subfig = GridLayout(figure[1,2], alignmode=Outside())
E_subfig = GridLayout(figure[2,1:2], alignmode=Outside())
F_subfig = GridLayout(figure[3,1:2], alignmode=Outside())

colsize!(figure.layout, 1, Relative(0.3))
colsize!(figure.layout, 2, Relative(0.7))

rowsize!(figure.layout, 1, Relative(0.6))
rowsize!(figure.layout, 2, Relative(0.2))
rowsize!(figure.layout, 3, Relative(0.2))
```

### Plot panels A and B - Comparative Importances

```julia
axA = Axis(
    AB_subfig[1, 1];
    xlabel = "-log(p) for LM coefficients",
    ylabel = "Fitness-weighted importance from  RFs",
    title = "00 to 06 months",
)

axB = Axis(
    AB_subfig[2, 1];
    xlabel = "-log(p) for LM coefficients",
    ylabel = "Fitness-weighted importance from  RFs",
    title = "18 to 120 months",
)

plot_colorset = [:gray, ColorSchemes.tableau_10[3], ColorSchemes.tableau_10[1], ColorSchemes.tableau_10[7]]
plot_comparative_lmvsrf_scatterplots!(axA, regression_currentCogScores_00to06mo_onlytaxa, "/brewster/kevin/scratch/derived/tables/lms_species_00to06.csv"; plot_colorset = plot_colorset)
plot_comparative_lmvsrf_scatterplots!(axB, regression_currentCogScores_18to120mo_onlytaxa, "/brewster/kevin/scratch/derived/tables/lms_species_18to120.csv"; plot_colorset = plot_colorset)

Legend(
    AB_subfig[3, 1],
    [
        MarkerElement(; marker=:circle, color=plot_colorset[3]),
        MarkerElement(; marker=:circle, color=plot_colorset[2]),
        MarkerElement(; marker=:circle, color=plot_colorset[4]),
        MarkerElement(; marker=:circle, color=plot_colorset[1]),
    ],
    [
        ">60% ranked importance",
        "q < 0.2 in LM",
        "Both",
        "None"
    ],
    "Input composition",
    tellheight = true,
    tellwidth = false,
    nbanks = 2,
    orientation = :horizontal
)
```

```julia
nbars_toplot = 25

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

axC = Axis(
    CD_subfig[1, 1];
    xlabel = "Importance",
    yticks = (reverse(collect(1:nbars_toplot)), replace.(joined_importances_00to06.variable[1:nbars_toplot], "_"=>" ")),
    yticklabelfont = "TeX Gyre Heros Makie Italic",
    ylabel = "Feature",
    title = "00 to 06 months",
    # yticklabelrotation= -pi/2
)

axD = Axis(
    CD_subfig[1, 2];
    xlabel = "Importance",
    yticks = (reverse(collect(1:nbars_toplot)), replace.(joined_importances_18to120.variable[1:nbars_toplot], "_"=>" ")),
    yticklabelfont = "TeX Gyre Heros Makie Italic",
    ylabel = "Feature",
    title = "18 to 120 months",
    # yticklabelrotation= -pi/2
)

plot_comparativedemo_importance_barplots!(axC, joined_importances_00to06; n_rows = nbars_toplot)
plot_comparativedemo_importance_barplots!(axD, joined_importances_18to120; n_rows = nbars_toplot)

Legend(
    CD_subfig[1, 2], [MarkerElement(; marker=:rect, color=:gray), MarkerElement(; marker=:star8, color=:red)], ["Microbiome alone", "Microbiome + demographics"], "Input composition",
    tellheight = false,
    tellwidth = false,
    margin = (10, 10, 10, 10),
    halign = :right, valign = :bottom, orientation = :vertical
)
```

### Plot panel E - Deep dives on taxa
```julia
mdata = Resonance.load(Metadata())
spec = Resonance.load(TaxonomicProfiles(); timepoint_metadata = mdata)

plot_taxon_deepdive!(E_subfig, 1, spec, :filter_00to06, "Blautia_wexlerae";)
plot_taxon_deepdive!(E_subfig, 2, spec, :filter_00to06, "Gordonibacter_pamelaeae";)
plot_taxon_deepdive!(E_subfig, 3, spec, :filter_00to06, "Bifidobacterium_longum";)
plot_taxon_deepdive!(E_subfig, 4, spec, :filter_00to06, "Ruminococcus_gnavus";)
plot_taxon_deepdive!(E_subfig, 5, spec, :filter_00to06, "Eggerthella_lenta";)
plot_taxon_deepdive!(E_subfig, 6, spec, :filter_00to06, "Erysipelatoclostridium_ramosum";)

plot_taxon_deepdive!(F_subfig, 1, spec, :filter_18to120, "Blautia_wexlerae";)
plot_taxon_deepdive!(F_subfig, 2, spec, :filter_18to120, "Gordonibacter_pamelaeae";)
plot_taxon_deepdive!(F_subfig, 3, spec, :filter_18to120, "Bifidobacterium_longum";)
plot_taxon_deepdive!(F_subfig, 4, spec, :filter_18to120, "Ruminococcus_gnavus";)
plot_taxon_deepdive!(F_subfig, 5, spec, :filter_18to120, "Faecalibacterium_prausnitzii";)
plot_taxon_deepdive!(F_subfig, 6, spec, :filter_18to120, "Alistipes_finegoldii";)
```

### Panel labels
```julia
Label(AB_subfig[1, 1, TopLeft()], "A", textsize = 26,font = :bold, padding = (0, 5, 5, 0), halign = :right)
Label(AB_subfig[2, 1, TopLeft()], "B", textsize = 26,font = :bold, padding = (0, 5, 5, 0), halign = :right)
Label(CD_subfig[1, 1, TopLeft()], "C", textsize = 26,font = :bold, padding = (0, 5, 5, 0), halign = :right)
Label(CD_subfig[1, 2, TopLeft()], "D", textsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
Label(E_subfig[1, 1, TopLeft()], "E", textsize = 26,font = :bold, padding = (0, 5, 5, 0), halign = :right)
Label(E_subfig[1, 5, TopLeft()], "F", textsize = 26,font = :bold, padding = (0, 5, 5, 0), halign = :right)
Label(F_subfig[1, 5, TopLeft()], "G", textsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)

Label(E_subfig[1:2, 1, Left()], "00 to 06 months", textsize = 20, font = :bold, padding = (0, 80, 0, 0), halign = :center, valign = :center, rotation = pi/2)
Label(F_subfig[1:2, 1, Left()], "18 to 120 months", textsize = 20, font = :bold, padding = (0, 80, 0, 0), halign = :center, valign = :center, rotation = pi/2)

save("manuscript/assets/Figure3.png", figure)
```

## Table 3 generation - [TODO: Move to `tables.md`]
```julia

conf_interval(vv, critical_z=1.96) = critical_z * std(vv) / sqrt(length(vv))

prod_summary_table = vcat(
    combine(
        groupby(regression_currentCogScores_00to06mo_onlydemo.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to06_onlydemo") => :model
    ),
    combine(
        groupby(regression_currentCogScores_00to06mo_onlytaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to06_onlytaxa") => :model
    ),
    combine(
        groupby(regression_currentCogScores_00to06mo_demoplustaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to06_demoplustaxa") => :model
    ),
    combine(
        groupby(regression_currentCogScores_00to06mo_onlyecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to06_onlyecs") => :model
    ),
    combine(
        groupby(regression_currentCogScores_00to06mo_demoplusecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to06_demoplusecs") => :model
    ),
    combine(
        groupby(regression_currentCogScores_18to120mo_onlydemo.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "18to120_onlydemo") => :model
    ),
    combine(
        groupby(regression_currentCogScores_18to120mo_onlytaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "18to120_onlytaxa") => :model
    ),
    combine(
        groupby(regression_currentCogScores_18to120mo_demoplustaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "18to120_demoplustaxa") => :model
    ),
    combine(
        groupby(regression_currentCogScores_18to120mo_onlyecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "18to120_onlyecs") => :model
    ),
    combine(
        groupby(regression_currentCogScores_18to120mo_demoplusecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "18to120_demoplusecs") => :model
    )
)

prod_summary_table = select(prod_summary_table, [:model, :Test_RMSE_mean, :Test_RMSE_CI, :Test_Cor_mean, :Test_Cor_CI])
```