#####
# Notebook D01 - Classification of binary current CogScores above/below 50th percentile
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
    select( [ :subject, :timepoint, :ageMonths, :cogScore, :cogScorePercentile, :education, :omni ] )
    dropmissing( [ :education, :omni ] )
end

insertcols!(mdata_df, 8, :educationInt => coerce(int.(skipmissing(mdata_df.education), type = Int), OrderedFactor))

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
    maxnodes_range = collect(1:1:15) ,
    nodesize_range = collect(1:1:20),
    sampsize_range = [0.5, 0.6, 0.7, 0.8],
    mtry_range = 1,
    ntrees_range = [100, 300, 500, 700]
    )

taxa_tuning_space = (
    maxnodes_range = collect(1:1:15) ,
    nodesize_range = collect(1:1:20),
    sampsize_range = [0.5, 0.6, 0.7, 0.8],
    mtry_range = collect(5:5:100),
    ntrees_range = [100, 300, 500, 700]
    )

ecs_tuning_space = (
    maxnodes_range = [1, 2, 3, 4, 5, 7, 9, 11, 13, 15],
    nodesize_range = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20],
    sampsize_range = [0.5, 0.6, 0.7, 0.8],
    mtry_range = collect(5:20:500),
    ntrees_range = [100, 300, 500, 700]
    )

upperhalf_percentile(x::Vector{T} where T <: Real) = coerce(x .>= 0.50, OrderedFactor)

insertcols!(mdata_df, 8, :educationInte => coerce(int.(skipmissing(mdata_df.education), type = Int), OrderedFactor)) # To tackle the problem with AbstractVetor as Input

# #####
# # 00 to 06 months
# #####

# ## 1. Only SES

# classification_currentCogScores_00to06mo_onlyses = train_randomforest(
#     Resonance.Classification(),
#     "classification_currentCogScores_00to06mo_onlyses",
#     mdata_df,
#     x -> dropmissing(filter_age_bracket(x, 0.0, 6.0)),
#     upperhalf_percentile, 
#     [8,9], # To tackle the problem with AbstractVetor as Input
#     :cogScorePercentile;
#     n_splits = 10,
#     tuning_space = onlyses_tuning_space,
#     train_rng = ml_rng
# )

# report_merits(classification_currentCogScores_00to06mo_onlyses)

# JLD2.@save "models/classification_currentCogScores_00to06mo_onlyses.jld" classification_currentCogScores_00to06mo_onlyses

# ## 2. Only taxonomic profiles

# classification_currentCogScores_00to06mo_onlytaxa = train_randomforest(
#     Resonance.Classification(),
#     "classification_currentCogScores_00to06mo_onlytaxa",
#     mdata_taxa_df,
#     x -> dropmissing(filter_age_bracket(x, 0.0, 6.0)),
#     upperhalf_percentile,
#     9:557,
#     :cogScorePercentile;
#     n_splits = 10,
#     tuning_space = taxa_tuning_space,
#     train_rng = ml_rng
# )

# report_merits(classification_currentCogScores_00to06mo_onlytaxa)

# JLD2.@save "models/classification_currentCogScores_00to06mo_onlytaxa.jld" classification_currentCogScores_00to06mo_onlytaxa

# ## 3. SES + taxonomic profiles

# classification_currentCogScores_00to06mo_sesplustaxa = train_randomforest(
#     Resonance.Classification(),
#     "classification_currentCogScores_00to06mo_sesplustaxa",
#     mdata_taxa_df,
#     x -> dropmissing(filter_age_bracket(x, 0.0, 6.0)),
#     upperhalf_percentile,
#     8:557,
#     :cogScorePercentile;
#     n_splits = 10,
#     tuning_space = taxa_tuning_space,
#     train_rng = ml_rng
# )

# report_merits(classification_currentCogScores_00to06mo_sesplustaxa)

# JLD2.@save "models/classification_currentCogScores_00to06mo_sesplustaxa.jld" classification_currentCogScores_00to06mo_sesplustaxa

# ## 4. Only functional profiles

# classification_currentCogScores_00to06mo_onlyecs = train_randomforest(
#     Resonance.Classification(),
#     "classification_currentCogScores_00to06mo_onlyecs",
#     mdata_ecs_df,
#     x -> dropmissing(filter_age_bracket(x, 0.0, 6.0)),
#     upperhalf_percentile,
#     9:2441,
#     :cogScorePercentile;
#     n_splits = 10,
#     tuning_space = ecs_tuning_space,
#     train_rng = ml_rng
# )

# report_merits(classification_currentCogScores_00to06mo_onlytaxa)

# JLD2.@save "models/classification_currentCogScores_00to06mo_onlyecs.jld" classification_currentCogScores_00to06mo_onlyecs

# ## 5. SES + taxonomic profiles

# classification_currentCogScores_00to06mo_sesplusecs = train_randomforest(
#     Resonance.Classification(),
#     "classification_currentCogScores_00to06mo_sesplusecs",
#     mdata_ecs_df,
#     x -> dropmissing(filter_age_bracket(x, 0.0, 6.0)),
#     upperhalf_percentile,
#     8:2441,
#     :cogScorePercentile;
#     n_splits = 10,
#     tuning_space = ecs_tuning_space,
#     train_rng = ml_rng
# )

# report_merits(classification_currentCogScores_00to06mo_sesplusecs)

# JLD2.@save "models/classification_currentCogScores_00to06mo_sesplusecs.jld" classification_currentCogScores_00to06mo_sesplusecs

#####
# 18 to 120 months
#####

## 6. Only SES

classification_currentCogScores_18to120mo_onlyses = train_randomforest(
    Resonance.Classification(),
    "classification_currentCogScores_18to120mo_onlyses",
    mdata_df,
    x -> dropmissing(filter_age_bracket(x, 18.0, 120.0)),
    upperhalf_percentile, 
    [8,9], # To tackle the problem with AbstractVetor as Input
    :cogScorePercentile;
    n_splits = 10,
    tuning_space = onlyses_tuning_space,
    train_rng = ml_rng
)

report_merits(classification_currentCogScores_18to120mo_onlyses)

JLD2.@save "models/classification_currentCogScores_18to120mo_onlyses.jld" classification_currentCogScores_18to120mo_onlyses

## 7. Only taxonomic profiles

classification_currentCogScores_18to120mo_onlytaxa = train_randomforest(
    Resonance.Classification(),
    "classification_currentCogScores_18to120mo_onlytaxa",
    mdata_taxa_df,
    x -> dropmissing(filter_age_bracket(x, 18.0, 120.0)),
    upperhalf_percentile,
    9:557,
    :cogScorePercentile;
    n_splits = 10,
    tuning_space = taxa_tuning_space,
    train_rng = ml_rng
)

report_merits(classification_currentCogScores_18to120mo_onlytaxa)

JLD2.@save "models/classification_currentCogScores_18to120mo_onlytaxa.jld" classification_currentCogScores_18to120mo_onlytaxa

## 8. SES + taxonomic profiles

classification_currentCogScores_18to120mo_sesplustaxa = train_randomforest(
    Resonance.Classification(),
    "classification_currentCogScores_18to120mo_sesplustaxa",
    mdata_taxa_df,
    x -> dropmissing(filter_age_bracket(x, 18.0, 120.0)),
    upperhalf_percentile,
    8:557,
    :cogScorePercentile;
    n_splits = 10,
    tuning_space = taxa_tuning_space,
    train_rng = ml_rng
)

report_merits(classification_currentCogScores_18to120mo_sesplustaxa)

JLD2.@save "models/classification_currentCogScores_18to120mo_sesplustaxa.jld" classification_currentCogScores_18to120mo_sesplustaxa

## 9. Only functional profiles

classification_currentCogScores_18to120mo_onlyecs = train_randomforest(
    Resonance.Classification(),
    "classification_currentCogScores_18to120mo_onlyecs",
    mdata_ecs_df,
    x -> dropmissing(filter_age_bracket(x, 18.0, 120.0)),
    upperhalf_percentile,
    9:2441,
    :cogScorePercentile;
    n_splits = 10,
    tuning_space = ecs_tuning_space,
    train_rng = ml_rng
)

report_merits(classification_currentCogScores_18to120mo_onlytaxa)

JLD2.@save "models/classification_currentCogScores_18to120mo_onlyecs.jld" classification_currentCogScores_18to120mo_onlyecs

## 10. SES + taxonomic profiles

classification_currentCogScores_18to120mo_sesplusecs = train_randomforest(
    Resonance.Classification(),
    "classification_currentCogScores_18to120mo_sesplusecs",
    mdata_ecs_df,
    x -> dropmissing(filter_age_bracket(x, 18.0, 120.0)),
    upperhalf_percentile,
    8:2441,
    :cogScorePercentile;
    n_splits = 10,
    tuning_space = ecs_tuning_space,
    train_rng = ml_rng
)

report_merits(classification_currentCogScores_18to120mo_sesplusecs)

JLD2.@save "models/classification_currentCogScores_18to120mo_sesplusecs.jld" classification_currentCogScores_18to120mo_sesplusecs