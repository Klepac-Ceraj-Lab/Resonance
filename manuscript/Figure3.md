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

## Defining relevant functions

### Calculating fitness-weighted importances
```julia
function calculate_fitness(train_cor::Vector{Float64}, test_cor::Vector{Float64})

    positive_train_cor = map( x -> maximum([x, 0.0]), train_cor)
    positive_test_cor = map( x -> maximum([x, 0.0]), test_cor)

    return positive_train_cor .* positive_test_cor

end

function weighted_hpimportances(m, hp = 1; change_hashnames=false, hashnamestable::Union{Nothing, DataFrame} = nothing)
    merits = m.merits
    importances = m.importances
    fitnesses = calculate_fitness(merits.Train_Cor, merits.Test_Cor)
    fitness_subset_idx = findall(merits.Hyperpar_Idx .== hp)

    mean_importances = map(x -> sum(x .* fitnesses[fitness_subset_idx])/sum(fitnesses[fitness_subset_idx] .> 0.0), collect(eachrow(Matrix(importances[:, fitness_subset_idx .+ 1]))))
    mean_importances_df = sort( DataFrame(:variable => importances[:,1], :weightedImportance => mean_importances), :weightedImportance; rev=true)

    if change_hashnames
        newnames = leftjoin(mean_importances_df, hashnamestable, on = :variable => :hashname).longname
        mean_importances_df.variable = newnames
    end

    return mean_importances_df
end
```

### Comparative Importance plots
```julia
function compute_joined_importances(modelwithoutdemo, modelwithdemo; imp_fun = weighted_hpimportances)

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
    barplot!(plot_axis, reverse(collect(1:n_rows)), joined_importances.ImportanceWithoutDemo[1:n_rows], color = "gray", direction = :x)
    scatter!(plot_axis, joined_importances.ImportanceWithDemo[1:n_rows], reverse(collect(1:n_rows)), marker = :star8, color = "red", direction = :x)
    
    return plot_axis

end
```

### Comparative Linear/RF importances scatterplots

```julia

function attribute_colors(lm_qvalue, rf_cumImp, plot_colorset = [:dodgerblue, :orange, :green, :purple])
    if (lm_qvalue & rf_cumImp)
        return plot_colorset[4]
    elseif (!lm_qvalue & rf_cumImp)
        return plot_colorset[3]
    elseif (lm_qvalue & !rf_cumImp)
        return plot_colorset[2]
    else
        return plot_colorset[1]
    end
end

function plot_comparative_lmvsrf_scatterplots!(
    plot_axis::Axis,
    rf_model::ProbeData,
    lm_path::String;
    exclude_from_importances= [ "ageMonths" ],
    cumulative_importance_threshold = 0.6,
    plot_colorset = [:dodgerblue, :orange, :green, :purple])

    ## 1. calculate the importances and the cumulative sum, excluding the relevant variables
    rf_model_importances = weighted_hpimportances(rf_model; change_hashnames=false)
    rf_model_importances.relativeWeightedImportance = rf_model_importances.weightedImportance ./ sum(rf_model_importances.weightedImportance[.!(rf_model_importances.variable .∈ Ref(exclude_from_importances))])
    rf_model_importances[rf_model_importances.variable .∈ Ref(exclude_from_importances), :relativeWeightedImportance] .= 0.0
    rf_model_importances.cumulativeWeightedImportance = cumsum(rf_model_importances.relativeWeightedImportance)
    
    ## 2. get the data for linear models, from `Figure2-calculations.jl`
    lm_coefs = CSV.read(lm_path, DataFrame)
    lm_coefs = select(lm_coefs, [:feature, :coef, :pvalue, :qvalue])

    ## 3. Join the data
    plot_comparative_df = dropmissing(outerjoin(rf_model_importances, lm_coefs, on = [ :variable => :feature ]))
    
    ## 4. Attribute color to each point
    is_significative_lm = map(q-> q < 0.2 ? true : false, plot_comparative_df.qvalue)
    is_over_threshold = map(ci-> ci <= cumulative_importance_threshold ? true : false, plot_comparative_df.cumulativeWeightedImportance)
    point_colors = map((a, b) -> attribute_colors(a, b, plot_colorset), is_significative_lm, is_over_threshold)

    # @show DataFrame(:variable => plot_comparative_df.variable, :color => point_colors) # For Debug

    # 5. Plot scatterplot
    scatter!(
        plot_axis,
        log.(plot_comparative_df.pvalue).*(-1), plot_comparative_df.weightedImportance,
        color = point_colors
    )

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

    ax1 = Axis(
        figure_layout[1,figure_col];
        ylabel="relative abundance",#\n$taxon_to_dive",
        title = replace(taxon_to_dive, "_"=>" "),
        titlefont = "TeX Gyre Heros Makie Italic")
    ax2 = Axis(figure_layout[2,figure_col]; xticks = (1:3, ["lower", "mid", "upper"]), xlabel="quartile", ylabel = "Intra-q\nPrevalence")
    ylims!(ax2, [0, 1.0])

    boxplot!(ax1, repeat([ 1.0 ], length(g1)), g1, color = (Makie.wong_colors()[1], 0.5), show_notch = false, show_outliers=false)
    boxplot!(ax1, repeat([ 2.0 ], length(g2)), g2, color = (Makie.wong_colors()[2], 0.5), show_notch = false, show_outliers=false)
    boxplot!(ax1, repeat([ 3.0 ], length(g3)), g3, color = (Makie.wong_colors()[3], 0.5), show_notch = false, show_outliers=false)    
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
JLD2.@load "models/2023-02-15/regression_currentCogScores_00to06mo_onlydemo.jld"
JLD2.@load "models/2023-02-15/regression_currentCogScores_00to06mo_onlytaxa.jld"
JLD2.@load "models/2023-02-15/regression_currentCogScores_00to06mo_demoplustaxa.jld"
JLD2.@load "models/2023-02-15/regression_currentCogScores_00to06mo_onlyecs.jld"
JLD2.@load "models/2023-02-15/regression_currentCogScores_00to06mo_demoplusecs.jld"
JLD2.@load "models/2023-02-15/regression_currentCogScores_18to120mo_onlydemo.jld"
JLD2.@load "models/2023-02-15/regression_currentCogScores_18to120mo_onlytaxa.jld"
JLD2.@load "models/2023-02-15/regression_currentCogScores_18to120mo_demoplustaxa.jld"
JLD2.@load "models/2023-02-15/regression_currentCogScores_18to120mo_onlyecs.jld"
JLD2.@load "models/2023-02-15/regression_currentCogScores_18to120mo_demoplusecs.jld"
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