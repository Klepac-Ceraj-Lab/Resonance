using Resonance
using CairoMakie
using Statistics
using MLJ
using DecisionTree
using JLD2
using StableRNGs
using GLM
using ProgressMeter
using CategoricalArrays

#####
# Loading Data
#####

RandomForestClassifier = MLJ.@load RandomForestClassifier pkg=DecisionTree
RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree
# cogScore classification from taxonomic profiles
JLD2.@load "models/classification_currentCogScores_00to06_fromtaxa_results.jld"
JLD2.@load "models/classification_currentCogScores_06to12_fromtaxa_results.jld"
JLD2.@load "models/classification_currentCogScores_12to18_fromtaxa_results.jld"
JLD2.@load "models/classification_currentCogScores_18to24_fromtaxa_results.jld"
JLD2.@load "models/classification_futureCogScores_allselected_fromtaxa_results.jld"
# cogScore regression from taxonomic profiles
JLD2.@load "models/regression_currentCogScores_00to06_fromtaxa_results.jld"
JLD2.@load "models/regression_currentCogScores_06to12_fromtaxa_results.jld"
JLD2.@load "models/regression_currentCogScores_12to18_fromtaxa_results.jld"
JLD2.@load "models/regression_currentCogScores_18to24_fromtaxa_results.jld"
JLD2.@load "models/regression_futureCogScores_allselected_fromtaxa_results.jld"

ensemble_classification_models = UnivariatePredictorEnsemble(
    [ "Concurrent, 0-6mo", "Concurrent, 6-12mo", "Concurrent, 12-18mo", "Concurrent, 18-24mo", "Future, all samples" ],
    [ :concurrent_00to06, :concurrent_06to12, :concurrent_12to18, :concurrent_18to24, :future_allsamples ],
    [
        classification_currentCogScores_00to06_fromtaxa_results,
        classification_currentCogScores_06to12_fromtaxa_results,
        classification_currentCogScores_12to18_fromtaxa_results,
        classification_currentCogScores_18to24_fromtaxa_results,
        classification_futureCogScores_allselected_fromtaxa_results
    ]
)

ensemble_regression_models = UnivariatePredictorEnsemble(
    [ "Concurrent, 0-6mo", "Concurrent, 6-12mo", "Concurrent, 12-18mo", "Concurrent, 18-24mo", "Future, all samples" ],
    [ :concurrent_00to06, :concurrent_06to12, :concurrent_12to18, :concurrent_18to24, :future_allsamples ],
    [
        regression_currentCogScores_00to06_fromtaxa_results,
        regression_currentCogScores_06to12_fromtaxa_results,
        regression_currentCogScores_12to18_fromtaxa_results,
        regression_currentCogScores_18to24_fromtaxa_results,
        regression_futureCogScores_allselected_fromtaxa_results
    ]
)

#####
# Computation
#####

## Classification

classification_individual_abundances = multimodel_individual_abundances(ensemble_classification_models)
classification_individual_prevalences = multimodel_individual_prevalences(ensemble_classification_models)
classification_individual_correlations = multimodel_individual_correlations(ensemble_classification_models)
classification_aggregate_importances = get_multimodel_aggregate_summaryimportances(ensemble_classification_models)
classification_aggregate_importances.Order = 1:nrow(classification_aggregate_importances)

plot_aggregate_df = DataFrames.outerjoin(classification_aggregate_importances, classification_individual_correlations, on = :Variable, makeunique = true)
plot_aggregate_df = DataFrames.outerjoin(plot_aggregate_df, classification_individual_prevalences, on = :Variable, makeunique = true)
plot_aggregate_df = DataFrames.outerjoin(plot_aggregate_df, classification_individual_abundances, on = :Variable, makeunique = true)
sort!(plot_aggregate_df, :Order)

#####
# Plotting
#####

# Logistic Regressions

logistic_regression_fig = Figure(resolution=(2400, 2000))

multimodel_logistic_regression!(
    logistic_regression_fig,
    ensemble_classification_models,
    [ "Gordonibacter_pamelaeae", "Bifidobacterium_breve", "Ruminococcus_gnavus", "Clostridium_innocuum", "Bifidobacterium_longum" ];
)

save("figures/logistic_regressions.png", logistic_regression_fig)

# Heatmaps

ens = ensemble_classification_models
n_plot = 50
heatmap_figure = Figure(resolution = (1800, 1800) )

## Correlation

correlation_colormap = cgrad(
    [:red, :white, :blue],
    [0.0, 0.5, 1.0]
)

ax_correlation = Axis(
    heatmap_figure[1,1];
    ylabel = "Most important Taxa",
    yticks = (collect(1:n_plot), plot_aggregate_df.Variable[1:n_plot]),
    xlabel = "Model",
    xticks = (collect(1:5), ens.screen_names),
    title = "Average Point-biserial correlation through the $(length(ens.predictors)) models",
    xticklabelrotation = pi/2,
    yreversed = true
)

corr_mat = permutedims(Matrix(plot_aggregate_df[1:n_plot, 4:8]))

correlation_hm = heatmap!(ax_correlation, corr_mat; colormap=correlation_colormap, colorrange=(-0.8, 0.8), highclip = :blue, lowclip = :red)

for ci in CartesianIndices(corr_mat)
    text!(ax_correlation,
    string(round(corr_mat[ci], digits=2)),
    ; position=(ci[1],ci[2]),
    align=(:center, :center)
    )
end

## Prevalence

prevalence_colormap = :viridis

ax_prevalence = Axis(
    heatmap_figure[1,2];
    #ylabel = "Most important Taxa",
    yticks = (collect(1:n_plot), plot_aggregate_df.Variable[1:n_plot]),
    yticksvisible = false,
    yticklabelsvisible = false,
    xlabel = "Model",
    xticks = (collect(1:5), ens.screen_names),
    title = "Average feature Prevalence through the $(length(ens.predictors)) models",
    xticklabelrotation = pi/2,
    yreversed = true
)

prevalence_mat = permutedims(Matrix(plot_aggregate_df[1:n_plot, 9:13]))

prevalence_hm = heatmap!(ax_prevalence, prevalence_mat; colormap=prevalence_colormap)

for ci in CartesianIndices(prevalence_mat)
    text!(ax_prevalence,
    string(round(prevalence_mat[ci], digits=4)),
    ; position=(ci[1],ci[2]),
    align=(:center, :center)
    )
end

## Prevalence

abundance_colormap = :viridis

ax_abundance = Axis(
    heatmap_figure[1,3];
    #ylabel = "Most important Taxa",
    yticks = (collect(1:n_plot), plot_aggregate_df.Variable[1:n_plot]),
    yticksvisible = false,
    yticklabelsvisible = false,
    xlabel = "Model",
    xticks = (collect(1:5), ens.screen_names),
    title = "Average feature Abundance through the $(length(ens.predictors)) models",
    xticklabelrotation = pi/2,
    yreversed = true
)

abundance_mat = permutedims(Matrix(plot_aggregate_df[1:n_plot, 14:18]))

abundance_hm = heatmap!(ax_abundance, abundance_mat; colormap=abundance_colormap)

for ci in CartesianIndices(abundance_mat)
    text!(ax_abundance,
    string(round(abundance_mat[ci], digits=4)),
    ; position=(ci[1],ci[2]),
    align=(:center, :center)
    )
end

Colorbar(heatmap_figure[2, 1], correlation_hm, label = "Point-Biserial correlation", vertical = false)
Colorbar(heatmap_figure[2, 2], prevalence_hm, label = "Prevalence", vertical = false)
Colorbar(heatmap_figure[2, 3], abundance_hm, label = "Abundance", vertical = false)

heatmap_figure

save("figures/correlation_prevalence_abundance_classification_heatmap.png", heatmap_figure)