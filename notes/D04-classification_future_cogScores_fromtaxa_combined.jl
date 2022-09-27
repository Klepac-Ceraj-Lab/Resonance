#####
# Notebook D04 - Classification of binary future CogScores above/below average from current taxonomic data
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
end # no dropmissing(:omni) because subjects/timepoints without sample but with assessment must be considered on the prediction routine

taxa = @chain Resonance.load(TaxonomicProfiles(); timepoint_metadata = mdata) begin
    filter( f-> taxrank(f) == :species, _ )
end

taxa_df = DataFrame([eachcol(collect(taxa.abundances'))...], map( x -> string(x), taxa.features), copycols=true)
taxa_df = taxa_df[ !, sortperm(names(taxa_df)) ]
rename!(taxa_df, [ names(taxa_df) .=> replace.(names(taxa_df), "s__" => "")]...)
taxa_df.Collinsella_massiliensis = myxor.(taxa_df.Collinsella_massiliensis, taxa_df."[Collinsella]_massiliensis")
select!(taxa_df, Not("[Collinsella]_massiliensis"))
insertcols!(taxa_df, 1, :sample => collect(keys(taxa.sidx))[collect(values(taxa.sidx))])

cogscore_taxa_df = leftjoin(mdata, taxa_df, on = :omni => :sample, matchmissing=:equal) |>  y -> rename!(y, :omni => :sample) |>  y -> sort(y, [ :subject, :timepoint ]);

#####
# Training and saving models
#####

RandomForestClassifier= MLJ.@load RandomForestClassifier pkg=DecisionTree

tuning_space = (
        maxnodes_range = collect(1:1:15) ,
        nodesize_range = collect(1:1:20),
        sampsize_range = [0.5, 0.6, 0.7, 0.8],
        mtry_range = collect(5:5:100),
        ntrees_range = [100, 300, 500, 700]
    )

classification_futureCogScores_allselected_fromtaxa_results = train_randomforest(
    Resonance.Classification(),
    "classification_futureCogScores_allselected_fromtaxa_results",
    cogscore_taxa_df,
    x -> prepare_future_prediction_df(x, Symbol.(names(taxa_df)), [ :cogScore ], 12.0, 24.0),
    meanclass,
    10:558,
    :futureCogScore;
    n_splits = 10,
    tuning_space = tuning_space,
    train_rng = ml_rng
)
JLD2.@save "models/classification_futureCogScores_allselected_fromtaxa_results.jld" classification_futureCogScores_allselected_fromtaxa_results