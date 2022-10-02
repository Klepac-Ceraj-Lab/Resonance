#####
# Notebook D12 - Regression of current CogScores above/below average from current functional data
#####

using Chain
using CSV
using Statistics
using Random
using DataFrames
using StableRNGs
using ProgressMeter
using MLJ
using CairoMakie
using DecisionTree
using JLD2
using Resonance
ml_rng = StableRNG(0)

#####
# Loading Data
#####

mdata = @chain Resonance.load(Metadata()) begin
    select( [ :subject, :timepoint, :ageMonths, :cogScore, :omni ] )
    dropmissing( [ :omni ] )
end

ecs = @chain Resonance.load(ECProfiles(); timepoint_metadata = mdata) begin
    filter(!hastaxon, _ )
    filter(f -> !in(name(f), ("UNMAPPED", "UNGROUPED")), _ )
end

ecs_df = DataFrame([eachcol(collect(ecs.abundances'))...], map( x -> string(x), ecs.features), copycols=true)
ecs_df = ecs_df[ !, sortperm(names(ecs_df)) ]
rename!(ecs_df, [ names(ecs_df) .=> replace.(names(ecs_df), "s__" => "")]...)
insertcols!(ecs_df, 1, :sample => collect(keys(ecs.sidx))[collect(values(ecs.sidx))])

cogscore_functions_df = leftjoin(mdata, ecs_df, on = :omni => :sample, matchmissing=:error) |>  y -> rename!(y, :omni => :sample) |>  y -> sort(y, [ :subject, :timepoint ]);

#####
# Training and saving models
#####

RandomForestRegressor= MLJ.@load RandomForestRegressor pkg=DecisionTree

tuning_space = (
    maxnodes_range = collect(1:1:15) ,
    nodesize_range = collect(1:1:20),
    sampsize_range = [0.5, 0.6, 0.7, 0.8],
    mtry_range = collect(5:10:200),
    ntrees_range = [100, 300, 500, 700]
)

## 0 to 6 months
regression_currentCogScores_00to06_fromfunctions_results = train_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_00to06_fromfunctions_results",
    cogscore_functions_df,
    x -> dropmissing(filter_age_bracket(x, 0.0, 6.0)),
    6:554,
    :cogScore;
    n_splits = 10,
    tuning_space = tuning_space,
    train_rng = ml_rng
)
JLD2.@save "models/regression_currentCogScores_00to06_fromfunctions_results.jld" regression_currentCogScores_00to06_fromfunctions_results

## 6 to 12 months
regression_currentCogScores_06to12_fromfunctions_results = train_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_06to12_fromfunctions_results",
    cogscore_functions_df,
    x -> dropmissing(filter_age_bracket(x, 6.0, 12.0)),
    6:554,
    :cogScore;
    n_splits = 10,
    tuning_space = tuning_space,
    train_rng = ml_rng
)
JLD2.@save "models/regression_currentCogScores_06to12_fromfunctions_results.jld" regression_currentCogScores_06to12_fromfunctions_results

## 12 to 18 months
regression_currentCogScores_12to18_fromfunctions_results = train_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_12to18_fromfunctions_results",
    cogscore_functions_df,
    x -> dropmissing(filter_age_bracket(x, 12.0, 18.0)),
    6:554,
    :cogScore;
    n_splits = 10,
    tuning_space = tuning_space,
    train_rng = ml_rng
)
JLD2.@save "models/regression_currentCogScores_12to18_fromfunctions_results.jld" regression_currentCogScores_12to18_fromfunctions_results

## 18 to 24 months
regression_currentCogScores_18to24_fromfunctions_results = train_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_18to24_fromfunctions_results",
    cogscore_functions_df,
    x -> dropmissing(filter_age_bracket(x, 18.0, 24.0)),
    6:554,
    :cogScore;
    n_splits = 10,
    tuning_space = tuning_space,
    train_rng = ml_rng
)

JLD2.@save "models/regression_currentCogScores_18to24_fromfunctions_results.jld" regression_currentCogScores_18to24_fromfunctions_results