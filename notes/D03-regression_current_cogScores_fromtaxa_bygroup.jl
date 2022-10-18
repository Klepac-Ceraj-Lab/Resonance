#####
# Notebook D03 - Regression of current CogScores above/below average from current taxonomic data
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

mdata = @chain Resonance.load(Metadata()) begin
    select( [ :subject, :timepoint, :ageMonths, :cogScore, :omni ] )
    dropmissing( [ :omni ] )
end

taxa = @chain Resonance.load(TaxonomicProfiles(); timepoint_metadata = mdata) begin
    filter( f-> taxrank(f) == :species, _ )
end

taxa_df = DataFrame([eachcol(collect(taxa.abundances'))...], map( x -> string(x), taxa.features), copycols=true)
taxa_df = taxa_df[ !, sortperm(names(taxa_df)) ]
rename!(taxa_df, [ names(taxa_df) .=> replace.(names(taxa_df), "s__" => "")]...)
taxa_df.Collinsella_massiliensis = myxor.(taxa_df.Collinsella_massiliensis, taxa_df."[Collinsella]_massiliensis")
select!(taxa_df, Not("[Collinsella]_massiliensis"))
insertcols!(taxa_df, 1, :sample => collect(keys(taxa.sidx))[collect(values(taxa.sidx))])

cogscore_taxa_df = leftjoin(mdata, taxa_df, on = :omni => :sample, matchmissing=:error) |>  y -> rename!(y, :omni => :sample) |>  y -> sort(y, [ :subject, :timepoint ]);

#####
# Training and saving models
#####

RandomForestRegressor= MLJ.@load RandomForestRegressor pkg=DecisionTree

tuning_space = (
        maxnodes_range = collect(1:1:15) ,
        nodesize_range = collect(1:1:20),
        sampsize_range = [0.5, 0.6, 0.7, 0.8],
        mtry_range = collect(5:5:100),
        ntrees_range = [100, 300, 500, 700]
    )

## 0 to 6 months
regression_currentCogScores_00to06_fromtaxa_results = train_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_00to06_fromtaxa_results",
    cogscore_taxa_df,
    x -> dropmissing(filter_age_bracket(x, 0.0, 6.0)),
    6:554,
    :cogScore;
    n_splits = 10,
    tuning_space = tuning_space,
    train_rng = ml_rng
)
JLD2.@save "models/regression_currentCogScores_00to06_fromtaxa_results.jld" regression_currentCogScores_00to06_fromtaxa_results

## 6 to 12 months
regression_currentCogScores_06to12_fromtaxa_results = train_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_06to12_fromtaxa_results",
    cogscore_taxa_df,
    x -> dropmissing(filter_age_bracket(x, 6.0, 12.0)),
    6:554,
    :cogScore;
    n_splits = 10,
    tuning_space = tuning_space,
    train_rng = ml_rng
)
JLD2.@save "models/regression_currentCogScores_06to12_fromtaxa_results.jld" regression_currentCogScores_06to12_fromtaxa_results

## 12 to 18 months
regression_currentCogScores_12to18_fromtaxa_results = train_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_12to18_fromtaxa_results",
    cogscore_taxa_df,
    x -> dropmissing(filter_age_bracket(x, 12.0, 18.0)),
    6:554,
    :cogScore;
    n_splits = 10,
    tuning_space = tuning_space,
    train_rng = ml_rng
)
JLD2.@save "models/regression_currentCogScores_12to18_fromtaxa_results.jld" regression_currentCogScores_12to18_fromtaxa_results

## 18 to 24 months
regression_currentCogScores_18to24_fromtaxa_results = train_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_18to24_fromtaxa_results",
    cogscore_taxa_df,
    x -> dropmissing(filter_age_bracket(x, 18.0, 24.0)),
    6:554,
    :cogScore;
    n_splits = 10,
    tuning_space = tuning_space,
    train_rng = ml_rng
)

JLD2.@save "models/regression_currentCogScores_18to24_fromtaxa_results.jld" regression_currentCogScores_18to24_fromtaxa_results