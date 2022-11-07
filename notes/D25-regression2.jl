#####
# Notebook D02 - Regression of current CogScores percentiles
#####

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
ml_rng = StableRNG(0)

#####
# Loading Data
#####

## 1. Subject-timepoint metadata, cogScore, SES (MaternalEd) and sample ref#

mdata_df = @chain Resonance.load(Metadata()) begin
    select( [ :subject, :timepoint, :ageMonths, :sex, :education, :cogScorePercentile, :omni, :read_depth] )
    dropmissing( [ :ageMonths, :sex, :education, :omni ] )
end

mdata_df.education = coerce(int.(skipmissing(mdata_df.education), type = Int), OrderedFactor)
mdata_df.sex = coerce(int.(skipmissing(mdata_df.sex), type = Int), OrderedFactor)

## 2. Taxonomic Profiles

taxa = @chain Resonance.load(TaxonomicProfiles(); timepoint_metadata = mdata_df) begin
    filter( f-> taxrank(f) == :species, _ )
end

taxa_df = DataFrame([eachcol(collect(taxa.abundances'))...], map( x -> string(x), taxa.features), copycols=true)
taxa_df = taxa_df[ !, sortperm(names(taxa_df)) ]
rename!(taxa_df, [ names(taxa_df) .=> replace.(names(taxa_df), "s__" => "")]...)
taxa_df.Collinsella_massiliensis = myxor.(taxa_df.Collinsella_massiliensis, taxa_df."[Collinsella]_massiliensis")
select!(taxa_df, Not("[Collinsella]_massiliensis"))
insertcols!(taxa_df, 1, :sample => collect(keys(taxa.sidx))[collect(values(taxa.sidx))])

mdata_taxa_df = leftjoin(mdata_df[:, 1:end-1], taxa_df, on = :omni => :sample, matchmissing=:error) |>  y -> rename!(y, :omni => :sample) |>  y -> sort(y, [ :subject, :timepoint ]);

## 3. Functional Profiles

ecs = @chain Resonance.load(ECProfiles(); timepoint_metadata = mdata_df) begin
    filter(!hastaxon, _ )
    filter(f -> !occursin(name(f), "UNMAPPED"), _ ) # "UNGROUPED" will be kept for this analysis
end

ecs_df = DataFrame([eachcol(collect(ecs.abundances'))...], map( x -> string(x), ecs.features), copycols=true)
ecs_df = ecs_df[ !, sortperm(names(ecs_df)) ]
rename!(ecs_df, [ names(ecs_df) .=> replace.(names(ecs_df), "s__" => "")]...)

insertcols!(ecs_df, 1, :sample => collect(keys(ecs.sidx))[collect(values(ecs.sidx))])

mdata_ecs_df = leftjoin(mdata_df, ecs_df, on = :omni => :sample, matchmissing=:error) |>  y -> rename!(y, :omni => :sample) |>  y -> sort(y, [ :subject, :timepoint ]);

##### IDEA 2022-10-26

for coll in names(mdata_ecs_df)[10:end]
    mdata_ecs_df[:, coll] .= mdata_ecs_df[:, coll] .* 1e6 ./ mdata_ecs_df[:, :read_depth]
end

select!(mdata_ecs_df, Not(:read_depth))

#####
# Training models
#####

RandomForestRegressor= MLJ.@load RandomForestRegressor pkg=DecisionTree

# # Production Tuning Grid
# onlydemo_tuning_space = (
#     maxnodes_range = collect(1:1:15),
#     nodesize_range = collect(1:1:20),
#     sampsize_range = [0.5, 0.6, 0.7, 0.8],
#     mtry_range = [ 1 ],
#     ntrees_range = [100, 300, 500, 700]
#     )

# taxa_tuning_space = (
#     maxnodes_range = collect(1:2:9),
#     nodesize_range = collect(1:2:15),
#     sampsize_range = [0.5, 0.6, 0.7],
#     mtry_range = collect(5:5:100),
#     ntrees_range = [100, 300, 500]
#     )

# ecs_tuning_space = (
#     maxnodes_range = collect(1:2:9),
#     nodesize_range = collect(1:2:15),
#     sampsize_range = [0.5, 0.6, 0.7],
#     mtry_range = collect(25:25:500),
#     ntrees_range = [100, 300, 500]
#     )

# Local test tuning grid
onlydemo_tuning_space = (
    maxnodes_range = [2, 4],
    nodesize_range = [2, 3],
    sampsize_range = [0.5, 0.6],
    mtry_range = [ 1 ],
    ntrees_range = [10, 30]
)

taxa_tuning_space = (
    maxnodes_range = [2, 4],
    nodesize_range = [2, 3],
    sampsize_range = [0.5, 0.6],
    mtry_range = [ 3, 5, 8, 10, 15, 20, 30, 40 ],
    ntrees_range = [10, 20, 30, 40]
)

ecs_tuning_space = (
    maxnodes_range = [2, 4, 6, 8, 10],
    nodesize_range = [2, 3, 4, 6],
    sampsize_range = [0.5, 0.6],
    mtry_range = [ 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500 ],
    ntrees_range = [10, 20, 30, 40, 50, 60, 70, 80]
)

#####
# 00 to 06 months
#####

## 1. Only SES

regression_currentCogScores_06to12mo_onlydemo = train_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_06to12mo_onlydemo",
    mdata_taxa_df[:, 1:8],
    x -> unique(dropmissing(filter_age_bracket(x, 6.0, 12.0)), :subject),
    [3,4,5],
    :cogScorePercentile;
    n_splits = 4,
    tuning_space = onlydemo_tuning_space,
    train_rng = ml_rng
)

report_merits(regression_currentCogScores_06to12mo_onlydemo)

JLD2.@save "models/regression_currentCogScores_06to12mo_onlydemo.jld" regression_currentCogScores_06to12mo_onlydemo

## 2. Only taxonomic profiles

regression_currentCogScores_06to12mo_onlytaxa = train_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_06to12mo_onlytaxa",
    mdata_taxa_df,
    x -> unique(dropmissing(filter_age_bracket(x, 6.0, 12.0)), :subject),
    8:556,
    :cogScorePercentile;
    n_splits = 4,
    tuning_space = taxa_tuning_space,
    train_rng = ml_rng
)

report_merits(regression_currentCogScores_06to12mo_onlytaxa)

JLD2.@save "models/regression_currentCogScores_06to12mo_onlytaxa.jld" regression_currentCogScores_06to12mo_onlytaxa

## 3. SES + taxonomic profiles

regression_currentCogScores_06to12mo_demoplustaxa = train_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_06to12mo_demoplustaxa",
    mdata_taxa_df,
    x -> unique(dropmissing(filter_age_bracket(x, 6.0, 12.0)), :subject),
    [3, 4, 5, collect(8:556)...],
    :cogScorePercentile;
    n_splits = 4,
    tuning_space = taxa_tuning_space,
    train_rng = ml_rng
)

report_merits(regression_currentCogScores_06to12mo_demoplustaxa)

JLD2.@save "models/regression_currentCogScores_06to12mo_demoplustaxa.jld" regression_currentCogScores_06to12mo_demoplustaxa

## 4. Only functional profiles

regression_currentCogScores_06to12mo_onlyecs = train_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_06to12mo_onlyecs",
    mdata_ecs_df,
    x -> unique(dropmissing(filter_age_bracket(x, 6.0, 12.0)), :subject),
    8:2440,
    :cogScorePercentile;
    n_splits = 4,
    tuning_space = ecs_tuning_space,
    train_rng = ml_rng
)

report_merits(regression_currentCogScores_06to12mo_onlyecs)

JLD2.@save "models/regression_currentCogScores_06to12mo_onlyecs.jld" regression_currentCogScores_06to12mo_onlyecs

## 5. SES + taxonomic profiles

regression_currentCogScores_06to12mo_demoplusecs = train_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_06to12mo_demoplusecs",
    mdata_ecs_df,
    x -> unique(dropmissing(filter_age_bracket(x, 6.0, 12.0)), :subject),
    [3, 4, 5, collect(8:2440)...],
    :cogScorePercentile;
    n_splits = 4,
    tuning_space = ecs_tuning_space,
    train_rng = ml_rng
)

report_merits(regression_currentCogScores_06to12mo_demoplusecs)

JLD2.@save "models/regression_currentCogScores_06to12mo_demoplusecs.jld" regression_currentCogScores_06to12mo_demoplusecs

#####
# 18 to 120 months
#####

## 6. Only SES

regression_currentCogScores_12to18mo_onlydemo = train_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_12to18mo_onlydemo",
    mdata_taxa_df[:, 1:8],
    x -> unique(dropmissing(filter_age_bracket(x, 12.0, 18.0)), :subject),
    [3,4,5],
    :cogScorePercentile;
    n_splits = 4,
    tuning_space = onlydemo_tuning_space,
    train_rng = ml_rng
)

report_merits(regression_currentCogScores_12to18mo_onlydemo)

JLD2.@save "models/regression_currentCogScores_12to18mo_onlydemo.jld" regression_currentCogScores_12to18mo_onlydemo

## 7. Only taxonomic profiles

regression_currentCogScores_12to18mo_onlytaxa = train_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_12to18mo_onlytaxa",
    mdata_taxa_df,
    x -> unique(dropmissing(filter_age_bracket(x, 12.0, 18.0)), :subject),
    8:556,
    :cogScorePercentile;
    n_splits = 4,
    tuning_space = taxa_tuning_space,
    train_rng = ml_rng
)

report_merits(regression_currentCogScores_12to18mo_onlytaxa)

JLD2.@save "models/regression_currentCogScores_12to18mo_onlytaxa.jld" regression_currentCogScores_12to18mo_onlytaxa

## 8. SES + taxonomic profiles

regression_currentCogScores_12to18mo_demoplustaxa = train_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_12to18mo_demoplustaxa",
    mdata_taxa_df,
    x -> unique(dropmissing(filter_age_bracket(x, 12.0, 18.0)), :subject),
    [3, 4, 5, collect(8:556)...],
    :cogScorePercentile;
    n_splits = 4,
    tuning_space = taxa_tuning_space,
    train_rng = ml_rng
)

report_merits(regression_currentCogScores_12to18mo_demoplustaxa)

JLD2.@save "models/regression_currentCogScores_12to18mo_demoplustaxa.jld" regression_currentCogScores_12to18mo_demoplustaxa

## 9. Only functional profiles

regression_currentCogScores_12to18mo_onlyecs = train_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_12to18mo_onlyecs",
    mdata_ecs_df,
    x -> unique(dropmissing(filter_age_bracket(x, 12.0, 18.0)), :subject),
    8:2440,
    :cogScorePercentile;
    n_splits = 4,
    tuning_space = ecs_tuning_space,
    train_rng = ml_rng
)

report_merits(regression_currentCogScores_12to18mo_onlyecs)

JLD2.@save "models/regression_currentCogScores_12to18mo_onlyecs.jld" regression_currentCogScores_12to18mo_onlyecs

## 10. SES + taxonomic profiles

regression_currentCogScores_12to18mo_demoplusecs = train_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_12to18mo_demoplusecs",
    mdata_ecs_df,
    x -> unique(dropmissing(filter_age_bracket(x, 12.0, 18.0)), :subject),
    [3, 4, 5, collect(8:2440)...],
    :cogScorePercentile;
    n_splits = 4,
    tuning_space = ecs_tuning_space,
    train_rng = ml_rng
)

report_merits(regression_currentCogScores_12to18mo_demoplusecs)

JLD2.@save "models/regression_currentCogScores_12to18mo_demoplusecs.jld" regression_currentCogScores_12to18mo_demoplusecs

# regression_results = [
#     regression_currentCogScores_06to12mo_onlydemo,
#     regression_currentCogScores_06to12mo_onlytaxa,
#     regression_currentCogScores_06to12mo_demoplustaxa,
#     regression_currentCogScores_06to12mo_onlyecs,
#     regression_currentCogScores_06to12mo_demoplusecs,
#     regression_currentCogScores_12to18mo_onlydemo,
#     regression_currentCogScores_12to18mo_onlytaxa,
#     regression_currentCogScores_12to18mo_demoplustaxa,
#     regression_currentCogScores_12to18mo_onlyecs,
#     regression_currentCogScores_12to18mo_demoplusecs
# ]

# selected_models = [ findmin(report_merits(m)[1:10, :Test_RMSE]) for m in regression_results ]
# selected_models = [ findmax(report_merits(m)[1:10, :Test_COR]) for m in regression_results ]
# selected_maes = [ el[1] for el in selected_models ]
# selected_indexes = [ el[2] for el in selected_models ]
# selected_correlations = [ report_merits(m)[i, :Test_COR] for (m,i) in zip(regression_results, selected_indexes) ]