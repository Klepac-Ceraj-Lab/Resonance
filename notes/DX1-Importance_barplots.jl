using Resonance
using CairoMakie
using ColorSchemes
using Statistics
using MLJ
using DecisionTree
using JLD2
using StableRNGs
using GLM

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

#####
# Individual model importances averaged over all splits
#####

individual_models = [
    classification_currentCogScores_00to06_fromtaxa_results,
    classification_currentCogScores_06to12_fromtaxa_results,
    classification_currentCogScores_12to18_fromtaxa_results,
    classification_currentCogScores_18to24_fromtaxa_results,
    classification_futureCogScores_allselected_fromtaxa_results,
    regression_currentCogScores_00to06_fromtaxa_results,
    regression_currentCogScores_06to12_fromtaxa_results,
    regression_currentCogScores_12to18_fromtaxa_results,
    regression_currentCogScores_18to24_fromtaxa_results,
    regression_futureCogScores_allselected_fromtaxa_results
]

individual_plot_titles = [
    "Classification model, concurrent cogScores, samples 0 to 6 months (n = 73)",
    "Classification model, concurrent cogScores, samples 6 to 12 months (n = 61)",
    "Classification model, concurrent cogScores, samples 12 to 18 months (n = 39)",
    "Classification model, concurrent cogScores, samples 18 to 24 months (n = 50)",
    "Classification model, future cogScores, all qualifying samples\n(sample < 12mo, closest later assesment < 24mo)",
    "Regression model, concurrent cogScores, samples 0 to 6 months (n = 73)",
    "Regression model, concurrent cogScores, samples 6 to 12 months (n = 61)",
    "Regression model, concurrent cogScores, samples 12 to 18 months (n = 39)",
    "Regression model, concurrent cogScores, samples 18 to 24 months (n = 50)",
    "Regression model, future cogScores, all qualifying samples\n(sample < 12mo, closest later assesment < 24mo)"
]

individual_output_paths = [
    "figures/importance_plots/Importances_classification_currentCogScores_00to06_fromtaxa.png",
    "figures/importance_plots/Importances_classification_currentCogScores_06to12_fromtaxa.png",
    "figures/importance_plots/Importances_classification_currentCogScores_12to18_fromtaxa.png",
    "figures/importance_plots/Importances_classification_currentCogScores_18to24_fromtaxa.png",
    "figures/importance_plots/Importances_classification_futureCogScores_allselected_fromtaxa.png",
    "figures/importance_plots/Importances_regression_currentCogScores_00to06_fromtaxa.png",
    "figures/importance_plots/Importances_regression_currentCogScores_06to12_fromtaxa.png",
    "figures/importance_plots/Importances_regression_currentCogScores_12to18_fromtaxa.png",
    "figures/importance_plots/Importances_regression_currentCogScores_18to24_fromtaxa.png",
    "figures/importance_plots/Importances_regression_futureCogScores_allselected_fromtaxa.png"
]

for i in eachindex(individual_models)

    figure = Figure( ; resolution = (1000, 800) )

    singlemodel_avgimportance_barplot!(
        figure,
        individual_models[i],
        (1,1), individual_plot_titles[i];
        n = 30
    )

    save(individual_output_paths[i], figure)

end

#####
# Ensembled average importances/topN
#####

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

## Plot summary figure for Classification models
figure = Figure( ; resolution = (1000, 1000) )
multimodel_avgimportance_barplot!(
    figure,
    ensemble_classification_models,
    (1,1), "Summary of variable importances throughout all classification models"
)
save("figures/importance_plots/Summarybarplot_classification_fromtaxa.png", figure)

## Plot summary figure for Classification models
figure = Figure( ; resolution = (1000, 1000) )
multimodel_avgimportance_barplot!(
    figure,
    ensemble_classification_models,
    (1,1), "Summary of variable importances throughout all regressions models"
)
save("figures/importance_plots/Summarybarplot_regression_fromtaxa.png", figure)