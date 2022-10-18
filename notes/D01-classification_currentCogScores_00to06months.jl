#####
# Notebook D01 - Classification of binary current CogScores above/below 50th percentile for 00 to 06 months
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
    select( [ :subject, :timepoint, :ageMonths, :cogScore, :education, :omni ] )
    dropmissing( [ :education, :omni ] )
end

insertcols!(mdata_df, 6, :educationInt => coerce(int.(skipmissing(mdata_df.education), type = Int), OrderedFactor))

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

mdata_taxa_df = leftjoin(mdata_df, taxa_df, on = :omni => :sample, matchmissing=:error) |>  y -> rename!(y, :omni => :sample) |>  y -> sort(y, [ :subject, :timepoint ]);

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

#####
# Training models
#####

RandomForestClassifier= MLJ.@load RandomForestClassifier pkg=DecisionTree

onlyses_tuning_space = (
        maxnodes_range = collect(1:1:10) ,
        nodesize_range = collect(1:2:10),
        sampsize_range = [0.5, 0.6, 0.7, 0.8],
        mtry_range = 1,
        ntrees_range = [100, 300, 500, 700]
    )

regular_tuning_space = (
        maxnodes_range = collect(1:5:15) ,
        nodesize_range = collect(1:5:20),
        sampsize_range = [0.5, 0.6],
        mtry_range = collect(5:15:100),
        ntrees_range = [300, 500]
    )

## 1. Only SES

insertcols!(mdata_df, 7, :educationInte => coerce(int.(skipmissing(mdata_df.education), type = Int), OrderedFactor)) # To tackle the problem with AbstractVetor as Input

classification_currentCogScores_00to16mo_onlyses = train_randomforest(
    Resonance.Classification(),
    "classification_currentCogScores_00to16mo_onlyses",
    mdata_df,
    x -> dropmissing(filter_age_bracket(x, 0.0, 6.0)),
    meanclass,
    [6, 7],
    :cogScore;
    n_splits = 10,
    tuning_space = ses_tuning_space,
    train_rng = ml_rng
)

report_merits(classification_currentCogScores_00to16mo_onlyses)

JLD2.@save "models/classification_currentCogScores_00to16mo_onlyses.jld" classification_currentCogScores_00to16mo_onlyses

## 2. Only taxonomic profiles

classification_currentCogScores_00to16mo_onlytaxa = train_randomforest(
    Resonance.Classification(),
    "classification_currentCogScores_00to16mo_onlytaxa",
    mdata_taxa_df,
    x -> dropmissing(filter_age_bracket(x, 6.0, 12.0)),
    meanclass,
    6:554,
    :cogScore;
    n_splits = 10,
    tuning_space = regular_tuning_space,
    train_rng = ml_rng
)
JLD2.@save "models/classification_currentCogScores_00to16mo_onlytaxa.jld" classification_currentCogScores_00to16mo_onlytaxa

## 3. SES + taxonomic profiles

classification_currentCogScores_12to18_fromtaxa_results = train_randomforest(
    Resonance.Classification(),
    "classification_currentCogScores_12to18_fromtaxa_results",
    cogscore_taxa_df,
    x -> dropmissing(filter_age_bracket(x, 12.0, 18.0)),
    meanclass,
    6:554,
    :cogScore;
    n_splits = 10,
    tuning_space = tuning_space,
    train_rng = ml_rng
)
JLD2.@save "models/classification_currentCogScores_12to18_fromtaxa_results.jld" classification_currentCogScores_12to18_fromtaxa_results

## 18 to 24 months
classification_currentCogScores_18to24_fromtaxa_results = train_randomforest(
    Resonance.Classification(),
    "classification_currentCogScores_18to24_fromtaxa_results",
    cogscore_taxa_df,
    x -> dropmissing(filter_age_bracket(x, 18.0, 24.0)),
    meanclass,
    6:554,
    :cogScore;
    n_splits = 10,
    tuning_space = tuning_space,
    train_rng = ml_rng
)

JLD2.@save "models/classification_currentCogScores_18to24_fromtaxa_results.jld" classification_currentCogScores_18to24_fromtaxa_results