
using Resonance
using CairoMakie
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

# cogScore classification from taxonomic profiles
report_merits(classification_currentCogScores_00to06_fromtaxa_results)
report_merits(classification_currentCogScores_06to12_fromtaxa_results)
report_merits(classification_currentCogScores_12to18_fromtaxa_results)
report_merits(classification_currentCogScores_18to24_fromtaxa_results)
report_merits(classification_futureCogScores_allselected_fromtaxa_results)
# cogScore regression from taxonomic profiles
report_merits(regression_currentCogScores_00to06_fromtaxa_results)
report_merits(regression_currentCogScores_06to12_fromtaxa_results)
report_merits(regression_currentCogScores_12to18_fromtaxa_results)
report_merits(regression_currentCogScores_18to24_fromtaxa_results)
report_merits(regression_futureCogScores_allselected_fromtaxa_results)

#####
# Plotting singlemodel importances
#####

get_singlemodel_singlesplit_importance(classification_currentCogScores_00to06_fromtaxa_results)
get_singlemodel_singlesplit_importance(regression_currentCogScores_00to06_fromtaxa_results)

get_singlemodel_allsplits_importances(classification_currentCogScores_00to06_fromtaxa_results)
get_singlemodel_allsplits_importances(regression_currentCogScores_00to06_fromtaxa_results)

get_singlemodel_summary_importances(classification_currentCogScores_00to06_fromtaxa_results)
get_singlemodel_summary_importances(regression_currentCogScores_00to06_fromtaxa_results)

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

get_multimodel_individual_summaryimportances(ensemble_classification_models)
get_multimodel_individual_summaryimportances(ensemble_regression_models)

get_multimodel_aggregate_summaryimportances(ensemble_classification_models)
get_multimodel_aggregate_summaryimportances(ensemble_regression_models)

#####
# Top N Importances
#####

get_singlemodel_binarytopn_importances(classification_currentCogScores_00to06_fromtaxa_results)
get_singlemodel_binarytopn_importances(regression_currentCogScores_00to06_fromtaxa_results)

get_multimodel_individual_binarytopns(ensemble_classification_models)
get_multimodel_individual_binarytopns(ensemble_regression_models)

get_multimodel_aggregate_binarytopns(ensemble_classification_models)
get_multimodel_aggregate_binarytopns(ensemble_regression_models)

