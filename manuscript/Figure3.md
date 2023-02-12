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
using DecisionTree
using JLD2
using Resonance
using Distributions
ml_rng = StableRNG(0)
```

## Defining relevant functions

### Calculating fitness-weighted importances
```julia
function calculate_fitness(train_cor::Vector{Float64}, test_cor::Vector{Float64})

    positive_train_cor = map( x -> maximum([x, 0.0]), train_cor)
    positive_test_cor = map( x -> maximum([x, 0.0]), test_cor)

    return positive_train_cor .* positive_test_cor

end

function weightedimportances(m; change_hashnames=false, hashnamestable::Union{Nothing, DataFrame} = nothing)

    merits = m.merits
    importances = m.importances
    fitnesses = calculate_fitness(merits.Train_Cor, merits.Test_Cor)
    weighted_importances = map(x -> sum(x .* fitnesses)/sum(fitnesses .> 0.0), collect(eachrow(Matrix(importances[:, 2:end]))))
    weighted_importances_df = sort( DataFrame(:variable => importances[:,1], :weightedImportance => weighted_importances), :weightedImportance; rev=true)
    sum(fitnesses .> 0.0) / nrow(merits)
    contributing_rows = merits[fitnesses .> 0.0, :]
    sum(contributing_rows.Train_Cor .> contributing_rows.Test_Cor)
    sum(contributing_rows.Train_Cor .< contributing_rows.Test_Cor)

    if change_hashnames
        newnames = leftjoin(weighted_importances_df, hashnamestable, on = :variable => :hashname).longname
        weighted_importances_df.variable = newnames
    end

    return weighted_importances_df
end
```

### Comparative Importance plots
```julia
function compute_joined_importances(modelwithoutdemo, modelwithdemo; imp_fun = hpimportances)

    importanceswithoutdemo = imp_fun(modelwithoutdemo)
    importanceswithdemo = imp_fun(modelwithdemo)
    rename!(importanceswithoutdemo, "weightedImportance" => "ImportanceWithoutDemo")
    rename!(importanceswithdemo, "weightedImportance" => "ImportanceWithDemo")

    joined_importances = dropmissing(outerjoin(importanceswithdemo, importanceswithoutdemo, on = [ :variable => :variable ]))
    sort!(joined_importances, :ImportanceWithoutDemo, rev = true)
    insertcols!(joined_importances, 1, :rank => collect(1:nrow(joined_importances)))

    return joined_importances

end

function plot_comparativedemo_importance_barplots!(plot_axis::Axis, joined_importances::DataFrame; n_rows = 40)

    # Plot barplot
    barplot!(plot_axis, reverse(collect(1:n_rows)), joined_importances.ImportanceWithoutDemo[1:n_rows], color = "medium turquoise", direction = :x)
    scatter!(plot_axis, joined_importances.ImportanceWithDemo[1:n_rows], reverse(collect(1:n_rows)), marker = :star8, color = "purple", direction = :x)
    
    return plot_axis

end
```

### Comparative Linear/RF importances scatterplots

```julia
function plot_comparative_lmvsrf_scatterplots!(plot_axis::Axis, rf_model::ProbeData, lm_path::String;)

    rf_model_importances = weightedimportances(rf_model; change_hashnames=false)
    lm_coefs = CSV.read(lm_path, DataFrame)
    lm_coefs = select(lm_coefs, [:feature, :coef, :pvalue])
    lm_coefs.coef = sqrt.(lm_coefs.coef .^ 2)
    plot_comaprative_df = outerjoin(rf_model_importances, lm_coefs, on = [ :variable => :feature ])

    # Plot scatterplot
    scatter!(plot_axis, log.(plot_comaprative_df.pvalue).*(-1), plot_comaprative_df.weightedImportance, color = "blue")

    return plot_axis
end
```

### Taxon deepdive quantile plots
```julia
function plot_taxon_deepdive!(figure_layout::GridLayout, figure_col::Int, taxonomic_profile::CommunityProfile, filter_row::Symbol, taxon_to_dive::String;)

    filtered_spec = taxonomic_profile[:, get(spec, filter_row)]
    (lq, uq) = quantile(get(filtered_spec, :cogScore), [0.25, 0.75])
    lqidx = findall(x-> x <= lq, get(filtered_spec, :cogScore))
    mqidx = findall(x-> lq < x <= uq, get(filtered_spec, :cogScore))
    uqidx = findall(x-> x > uq, get(filtered_spec, :cogScore))

    taxab = vec(abundances(filtered_spec[Regex(taxon_to_dive), :]))
    g1 = filter(!=(0), taxab[lqidx])
    g2 = filter(!=(0), taxab[mqidx])
    g3 = filter(!=(0), taxab[uqidx])

    g1p = count(!=(0), taxab[lqidx]) / length(lqidx)
    g2p = count(!=(0), taxab[mqidx]) / length(mqidx)
    g3p = count(!=(0), taxab[uqidx]) / length(uqidx)

    ax1 = Axis(figure_layout[1,figure_col]; ylabel="relative abundance\n$taxon_to_dive")
    ax2 = Axis(figure_layout[2,figure_col]; xticks = (1:3, ["lower", "mid", "upper"]), xlabel="quartile", ylabel = "Intra-q\nPrevalence")
    
    Random.seed!(0)
    scatter!(ax1, 1 .+ rand(Normal(0, 0.1), length(g1)), g1)
    Random.seed!(0)
    scatter!(ax1, 2 .+ rand(Normal(0, 0.1), length(g2)), g2)
    Random.seed!(0)
    scatter!(ax1, 3 .+ rand(Normal(0, 0.1), length(g3)), g3)
    barplot!(ax2, [1,2,3], [g1p, g2p, g3p]; color=Makie.wong_colors()[1:3])
    
    hidexdecorations!(ax1; grid=false)
    rowsize!(figure_layout, 2, Relative(1/3))

    return figure_layout
end
```

## Loading the pretrained models

```julia
RandomForestClassifier = MLJ.@load RandomForestClassifier pkg=DecisionTree
RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree
# concurrent cogScore reression from taxonomic profiles
JLD2.@load "models/2023-01-24/regression_currentCogScores_00to06mo_onlydemo.jld"
JLD2.@load "models/2023-01-24/regression_currentCogScores_00to06mo_onlytaxa.jld"
JLD2.@load "models/2023-01-23/regression_currentCogScores_00to06mo_demoplustaxa.jld"
JLD2.@load "models/2023-01-31/regression_currentCogScores_00to06mo_onlyecs.jld"
JLD2.@load "models/2023-01-31/regression_currentCogScores_00to06mo_demoplusecs.jld"
JLD2.@load "models/2023-01-24/regression_currentCogScores_18to120mo_onlydemo.jld"
JLD2.@load "models/2023-01-24/regression_currentCogScores_18to120mo_onlytaxa.jld"
JLD2.@load "models/2023-01-24/regression_currentCogScores_18to120mo_demoplustaxa.jld"
JLD2.@load "models/2023-01-31/regression_currentCogScores_18to120mo_onlyecs.jld"
JLD2.@load "models/2023-01-31/regression_currentCogScores_18to120mo_demoplusecs.jld"
```

## Initializing Figure 3

```julia
figure = Figure(resolution = (1920, 1536))

AB_subfig = GridLayout(figure[1,1])
CD_subfig = GridLayout(figure[1,2])
E_subfig = GridLayout(figure[2,1:2])
F_subfig = GridLayout(figure[3,1:2])


colsize!(figure.layout, 1, Relative(0.7))
colsize!(figure.layout, 2, Relative(0.3))

rowsize!(figure.layout, 1, Relative(0.6))
rowsize!(figure.layout, 2, Relative(0.2))
rowsize!(figure.layout, 3, Relative(0.2))
```

### Plot panels A and B - Comparative Importances

```julia
nbars_toplot = 30

joined_importances_00to06 = compute_joined_importances(
    regression_currentCogScores_00to06mo_onlytaxa,
    regression_currentCogScores_00to06mo_demoplustaxa;
    imp_fun = weightedimportances
)
joined_importances_18to120 = compute_joined_importances(
    regression_currentCogScores_18to120mo_onlytaxa,
    regression_currentCogScores_18to120mo_demoplustaxa;
    imp_fun = weightedimportances
)

axA = Axis(
    AB_subfig[1, 1];
    xlabel = "Importance",
    yticks = (reverse(collect(1:nbars_toplot)), replace.(joined_importances_00to06.variable[1:nbars_toplot], "_"=>" ")),
    yticklabelfont = "TeX Gyre Heros Makie Italic",
    ylabel = "Feature",
    title = "00 to 06 months",
    # yticklabelrotation= -pi/2
)

axB = Axis(
    AB_subfig[1, 2];
    xlabel = "Importance",
    yticks = (reverse(collect(1:nbars_toplot)), replace.(joined_importances_18to120.variable[1:nbars_toplot], "_"=>" ")),
    yticklabelfont = "TeX Gyre Heros Makie Italic",
    ylabel = "Feature",
    title = "18 to 120 months",
    # yticklabelrotation= -pi/2
)

plot_comparativedemo_importance_barplots!(axA, joined_importances_00to06; n_rows = nbars_toplot)
plot_comparativedemo_importance_barplots!(axB, joined_importances_18to120; n_rows = nbars_toplot)

Legend(AB_subfig[2, 1:2], [MarkerElement(; marker=:rect, color=:aqua), MarkerElement(; marker=:star8, color=:purple)], ["Microbiome alone", "Microbiome + demographics"], "Input composition", orientation=:horizontal)
```

```julia
axC = Axis(
    CD_subfig[1, 1];
    xlabel = "-log(p) for LM coefficients",
    ylabel = "Fitness-weighted importance from  RFs",
    title = "00 to 06 months",
)

axD = Axis(
    CD_subfig[2, 1];
    xlabel = "-log(p) for LM coefficients",
    ylabel = "Fitness-weighted importance from  RFs",
    title = "18 to 120 months",
)

plot_comparative_lmvsrf_scatterplots!(axC, regression_currentCogScores_00to06mo_onlytaxa, "/brewster/kevin/scratch/derived/tables/lms_species_00to06.csv";)
plot_comparative_lmvsrf_scatterplots!(axD, regression_currentCogScores_18to120mo_onlytaxa, "/brewster/kevin/scratch/derived/tables/lms_species_18to120.csv";)
```

### Plot panel E - Deep dives on taxa
```julia
mdata = Resonance.load(Metadata())
spec = Resonance.load(TaxonomicProfiles(); timepoint_metadata = mdata)

plot_taxon_deepdive!(E_subfig, 1, spec, :filter_00to06, "Blautia_wexlerae";)
plot_taxon_deepdive!(E_subfig, 2, spec, :filter_00to06, "Erysipelatoclostridium_ramosum";)
plot_taxon_deepdive!(E_subfig, 3, spec, :filter_00to06, "Bifidobacterium_longum";)
plot_taxon_deepdive!(E_subfig, 4, spec, :filter_00to06, "Ruminococcus_gnavus";)
plot_taxon_deepdive!(E_subfig, 5, spec, :filter_00to06, "Gordonibacter_pamelaeae";)
plot_taxon_deepdive!(E_subfig, 6, spec, :filter_00to06, "Eggerthella_lenta";)

plot_taxon_deepdive!(F_subfig, 1, spec, :filter_18to120, "Alistipes_finegoldii";)
plot_taxon_deepdive!(F_subfig, 2, spec, :filter_18to120, "Faecalibacterium_prausnitzii";)
plot_taxon_deepdive!(F_subfig, 3, spec, :filter_18to120, "Ruminococcus_gnavus";)
plot_taxon_deepdive!(F_subfig, 4, spec, :filter_18to120, "Bifidobacterium_longum";)
plot_taxon_deepdive!(F_subfig, 5, spec, :filter_18to120, "Erysipelatoclostridium_ramosum";)
plot_taxon_deepdive!(F_subfig, 6, spec, :filter_18to120, "Eubacterium_eligens";)
```