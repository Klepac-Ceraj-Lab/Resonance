######
# This notebook was used to convert old Dictionaries to the newer data structures
######

using Resonance
using CairoMakie
using Statistics
using MLJ
using DecisionTree
using JLD2
using StableRNGs
using GLM
RandomForestClassifier = MLJ.@load RandomForestClassifier pkg=DecisionTree
RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree

##################
# CLASSIFICATION #
##################

###
# classification_currentCogScores_00to06_fromtaxa_results
JLD2.@load "/home/guilherme/Documents/models/classification_currentCogScores_00to06_fromtaxa_results.jld"

classification_currentCogScores_00to06_fromtaxa_results = UnivariateRandomForestClassifier(
    "classification_currentCogScores_00to06_fromtaxa_results",
    classification_currentCogScores_00to06_fromtaxa_results[:inputs_outputs],
    classification_currentCogScores_00to06_fromtaxa_results[:n_trials],
    classification_currentCogScores_00to06_fromtaxa_results[:dataset_partitions],
    classification_currentCogScores_00to06_fromtaxa_results[:models],
    classification_currentCogScores_00to06_fromtaxa_results[:selected_trial],
    classification_currentCogScores_00to06_fromtaxa_results[:train_accuracies],
    classification_currentCogScores_00to06_fromtaxa_results[:test_accuracies],
)

JLD2.@save "models/classification_currentCogScores_00to06_fromtaxa_results.jld" classification_currentCogScores_00to06_fromtaxa_results

###
# classification_currentCogScores_06to12_fromtaxa_results
JLD2.@load "/home/guilherme/Documents/models/classification_currentCogScores_06to12_fromtaxa_results.jld"

classification_currentCogScores_06to12_fromtaxa_results = UnivariateRandomForestClassifier(
    "classification_currentCogScores_06to12_fromtaxa_results",
    classification_currentCogScores_06to12_fromtaxa_results[:inputs_outputs],
    classification_currentCogScores_06to12_fromtaxa_results[:n_trials],
    classification_currentCogScores_06to12_fromtaxa_results[:dataset_partitions],
    classification_currentCogScores_06to12_fromtaxa_results[:models],
    classification_currentCogScores_06to12_fromtaxa_results[:selected_trial],
    classification_currentCogScores_06to12_fromtaxa_results[:train_accuracies],
    classification_currentCogScores_06to12_fromtaxa_results[:test_accuracies],
)

JLD2.@save "models/classification_currentCogScores_06to12_fromtaxa_results.jld" classification_currentCogScores_06to12_fromtaxa_results

###
# classification_currentCogScores_12to18_fromtaxa_results
JLD2.@load "/home/guilherme/Documents/models/classification_currentCogScores_12to18_fromtaxa_results.jld"

classification_currentCogScores_12to18_fromtaxa_results = UnivariateRandomForestClassifier(
    "classification_currentCogScores_12to18_fromtaxa_results",
    classification_currentCogScores_12to18_fromtaxa_results[:inputs_outputs],
    classification_currentCogScores_12to18_fromtaxa_results[:n_trials],
    classification_currentCogScores_12to18_fromtaxa_results[:dataset_partitions],
    classification_currentCogScores_12to18_fromtaxa_results[:models],
    classification_currentCogScores_12to18_fromtaxa_results[:selected_trial],
    classification_currentCogScores_12to18_fromtaxa_results[:train_accuracies],
    classification_currentCogScores_12to18_fromtaxa_results[:test_accuracies],
)

JLD2.@save "models/classification_currentCogScores_12to18_fromtaxa_results.jld" classification_currentCogScores_12to18_fromtaxa_results

###
# classification_currentCogScores_18to24_fromtaxa_results
JLD2.@load "/home/guilherme/Documents/models/classification_currentCogScores_18to24_fromtaxa_results.jld"

classification_currentCogScores_18to24_fromtaxa_results = UnivariateRandomForestClassifier(
    "classification_currentCogScores_18to24_fromtaxa_results",
    classification_currentCogScores_18to24_fromtaxa_results[:inputs_outputs],
    classification_currentCogScores_18to24_fromtaxa_results[:n_trials],
    classification_currentCogScores_18to24_fromtaxa_results[:dataset_partitions],
    classification_currentCogScores_18to24_fromtaxa_results[:models],
    classification_currentCogScores_18to24_fromtaxa_results[:selected_trial],
    classification_currentCogScores_18to24_fromtaxa_results[:train_accuracies],
    classification_currentCogScores_18to24_fromtaxa_results[:test_accuracies],
)

JLD2.@save "models/classification_currentCogScores_18to24_fromtaxa_results.jld" classification_currentCogScores_18to24_fromtaxa_results

###
# classification_futureCogScores_allselected_fromtaxa_results
JLD2.@load "/home/guilherme/Documents/models/classification_futureCogScores_allselected_fromtaxa_results.jld"

classification_futureCogScores_allselected_fromtaxa_results = UnivariateRandomForestClassifier(
    "classification_futureCogScores_allselected_fromtaxa_results",
    classification_futureCogScores_allselected_fromtaxa_results[:inputs_outputs],
    classification_futureCogScores_allselected_fromtaxa_results[:n_trials],
    classification_futureCogScores_allselected_fromtaxa_results[:dataset_partitions],
    classification_futureCogScores_allselected_fromtaxa_results[:models],
    classification_futureCogScores_allselected_fromtaxa_results[:selected_trial],
    classification_futureCogScores_allselected_fromtaxa_results[:train_accuracies],
    classification_futureCogScores_allselected_fromtaxa_results[:test_accuracies],
)

JLD2.@save "models/classification_futureCogScores_allselected_fromtaxa_results.jld" classification_futureCogScores_allselected_fromtaxa_results

##############
# REGRESSION #
##############

###
# regression_currentCogScores_00to06_fromtaxa_results
JLD2.@load "/home/guilherme/Documents/models/regression_currentCogScores_00to06_fromtaxa_results.jld"

regression_currentCogScores_00to06_fromtaxa_results = UnivariateRandomForestRegressor(
    "regression_currentCogScores_00to06_fromtaxa_results",
    regression_currentCogScores_00to06_fromtaxa_results[:inputs_outputs],
    regression_currentCogScores_00to06_fromtaxa_results[:n_trials],
    regression_currentCogScores_00to06_fromtaxa_results[:dataset_partitions],
    regression_currentCogScores_00to06_fromtaxa_results[:models],
    regression_currentCogScores_00to06_fromtaxa_results[:slope_corrections],
    regression_currentCogScores_00to06_fromtaxa_results[:selected_trial],
    regression_currentCogScores_00to06_fromtaxa_results[:train_maes],
    regression_currentCogScores_00to06_fromtaxa_results[:test_maes],
    regression_currentCogScores_00to06_fromtaxa_results[:train_mapes],
    regression_currentCogScores_00to06_fromtaxa_results[:test_mapes],
    regression_currentCogScores_00to06_fromtaxa_results[:train_correlations],
    regression_currentCogScores_00to06_fromtaxa_results[:test_correlations],
)

JLD2.@save "models/regression_currentCogScores_00to06_fromtaxa_results.jld" regression_currentCogScores_00to06_fromtaxa_results

###
# regression_currentCogScores_06to12_fromtaxa_results
JLD2.@load "/home/guilherme/Documents/models/regression_currentCogScores_06to12_fromtaxa_results.jld"

regression_currentCogScores_06to12_fromtaxa_results = UnivariateRandomForestRegressor(
    "regression_currentCogScores_06to12_fromtaxa_results",
    regression_currentCogScores_06to12_fromtaxa_results[:inputs_outputs],
    regression_currentCogScores_06to12_fromtaxa_results[:n_trials],
    regression_currentCogScores_06to12_fromtaxa_results[:dataset_partitions],
    regression_currentCogScores_06to12_fromtaxa_results[:models],
    regression_currentCogScores_06to12_fromtaxa_results[:slope_corrections],
    regression_currentCogScores_06to12_fromtaxa_results[:selected_trial],
    regression_currentCogScores_06to12_fromtaxa_results[:train_maes],
    regression_currentCogScores_06to12_fromtaxa_results[:test_maes],
    regression_currentCogScores_06to12_fromtaxa_results[:train_mapes],
    regression_currentCogScores_06to12_fromtaxa_results[:test_mapes],
    regression_currentCogScores_06to12_fromtaxa_results[:train_correlations],
    regression_currentCogScores_06to12_fromtaxa_results[:test_correlations],
)

JLD2.@save "models/regression_currentCogScores_06to12_fromtaxa_results.jld" regression_currentCogScores_06to12_fromtaxa_results

###
# regression_currentCogScores_12to18_fromtaxa_results
JLD2.@load "/home/guilherme/Documents/models/regression_currentCogScores_12to18_fromtaxa_results.jld"

regression_currentCogScores_12to18_fromtaxa_results = UnivariateRandomForestRegressor(
    "regression_currentCogScores_12to18_fromtaxa_results",
    regression_currentCogScores_12to18_fromtaxa_results[:inputs_outputs],
    regression_currentCogScores_12to18_fromtaxa_results[:n_trials],
    regression_currentCogScores_12to18_fromtaxa_results[:dataset_partitions],
    regression_currentCogScores_12to18_fromtaxa_results[:models],
    regression_currentCogScores_12to18_fromtaxa_results[:slope_corrections],
    regression_currentCogScores_12to18_fromtaxa_results[:selected_trial],
    regression_currentCogScores_12to18_fromtaxa_results[:train_maes],
    regression_currentCogScores_12to18_fromtaxa_results[:test_maes],
    regression_currentCogScores_12to18_fromtaxa_results[:train_mapes],
    regression_currentCogScores_12to18_fromtaxa_results[:test_mapes],
    regression_currentCogScores_12to18_fromtaxa_results[:train_correlations],
    regression_currentCogScores_12to18_fromtaxa_results[:test_correlations],
)

JLD2.@save "models/regression_currentCogScores_12to18_fromtaxa_results.jld" regression_currentCogScores_12to18_fromtaxa_results

###
# regression_currentCogScores_18to24_fromtaxa_results
JLD2.@load "/home/guilherme/Documents/models/regression_currentCogScores_18to24_fromtaxa_results.jld"

regression_currentCogScores_18to24_fromtaxa_results = UnivariateRandomForestRegressor(
    "regression_currentCogScores_18to24_fromtaxa_results",
    regression_currentCogScores_18to24_fromtaxa_results[:inputs_outputs],
    regression_currentCogScores_18to24_fromtaxa_results[:n_trials],
    regression_currentCogScores_18to24_fromtaxa_results[:dataset_partitions],
    regression_currentCogScores_18to24_fromtaxa_results[:models],
    regression_currentCogScores_18to24_fromtaxa_results[:slope_corrections],
    regression_currentCogScores_18to24_fromtaxa_results[:selected_trial],
    regression_currentCogScores_18to24_fromtaxa_results[:train_maes],
    regression_currentCogScores_18to24_fromtaxa_results[:test_maes],
    regression_currentCogScores_18to24_fromtaxa_results[:train_mapes],
    regression_currentCogScores_18to24_fromtaxa_results[:test_mapes],
    regression_currentCogScores_18to24_fromtaxa_results[:train_correlations],
    regression_currentCogScores_18to24_fromtaxa_results[:test_correlations],
)

JLD2.@save "models/regression_currentCogScores_18to24_fromtaxa_results.jld" regression_currentCogScores_18to24_fromtaxa_results

###
# regression_futureCogScores_allselected_fromtaxa_results
JLD2.@load "/home/guilherme/Documents/models/regression_futureCogScores_allselected_fromtaxa_results.jld"

regression_futureCogScores_allselected_fromtaxa_results = UnivariateRandomForestRegressor(
    "regression_futureCogScores_allselected_fromtaxa_results",
    regression_futureCogScores_allselected_fromtaxa_results[:inputs_outputs],
    regression_futureCogScores_allselected_fromtaxa_results[:n_trials],
    regression_futureCogScores_allselected_fromtaxa_results[:dataset_partitions],
    regression_futureCogScores_allselected_fromtaxa_results[:models],
    regression_futureCogScores_allselected_fromtaxa_results[:slope_corrections],
    regression_futureCogScores_allselected_fromtaxa_results[:selected_trial],
    regression_futureCogScores_allselected_fromtaxa_results[:train_maes],
    regression_futureCogScores_allselected_fromtaxa_results[:test_maes],
    regression_futureCogScores_allselected_fromtaxa_results[:train_mapes],
    regression_futureCogScores_allselected_fromtaxa_results[:test_mapes],
    regression_futureCogScores_allselected_fromtaxa_results[:train_correlations],
    regression_futureCogScores_allselected_fromtaxa_results[:test_correlations],
)

JLD2.@save "models/regression_futureCogScores_allselected_fromtaxa_results.jld" regression_futureCogScores_allselected_fromtaxa_results