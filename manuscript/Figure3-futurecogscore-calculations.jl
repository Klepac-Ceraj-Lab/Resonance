##########
# Confuture cognitive assessment Score regression notebook
##########

#####
# 0. Loading libraries
#####

# using Distributed
# addprocs(8)  # 8 procs.

# @everywhere begin
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
    using CategoricalArrays
    using Resonance
    using GLM
    # using ShapML        
    ml_rng = StableRNG(0)
# end

function build_future_row(df)
    this_subject = df.subject[1]
    this_sex = df.sex[1]
    this_education = df.education[1]
    this_present_timepoint = minimum(df.timepoint)
    this_present_ageMonths = minimum(df.ageMonths)
    this_future_timepoint = maximum(df.timepoint)
    this_future_ageMonths = maximum(df.ageMonths)
    this_present_sample = df.sample[df.filter_future_omni]
    this_future_cog = df.cogScore[df.filter_future_cog]

    DataFrame(
        :subject => this_subject,
        :present_timepoint => this_present_timepoint,
        :future_timepoint => this_future_timepoint,
        :sample => this_present_sample,
        :future_cogScore => this_future_cog,
        :sex => this_sex,
        :education => this_education,
        :present_ageMonths => this_present_ageMonths,
        :future_ageMonths => this_future_ageMonths
    )
end

#####
# Loading Data
#####

## 1. Metadata
mdata = Resonance.load(Metadata())
future_mdata = vcat( [ build_future_row(df) for df in groupby(mdata[ mdata.filter_future_omni .| mdata.filter_future_cog, : ], :subject) ]... )
seqs = subset(future_mdata, "sample"=> ByRow(!ismissing)) 
transform!(seqs, "sample"=> ByRow(String)=> "sample")
seqs.edfloat = map(x-> ismissing(x) ? missing : Float64(levelcode(x)), seqs.education)

## 2. Taxonomic Profiles
taxa = Resonance.load(TaxonomicProfiles(); timepoint_metadata = seqs) # this can take a bit
species = filter(f-> taxrank(f) == :species, taxa)
mdata_taxa_df = sort(Resonance.comm2wide(species), [ :subject ]);

## 3. Functional Profiles
ecs = Resonance.load(ECProfiles(); timepoint_metadata = seqs)
mdata_ecs_df = sort(Resonance.comm2wide(ecs), [ :subject ]);

### 3.1. Hashing EC names with fixed-length alphabetic hashes due to length and special characters; storing text files for remapping.
oldnames = names(mdata_ecs_df)[20:end]
for oldname in names(mdata_ecs_df)[20:end]
    rename!(mdata_ecs_df, oldname => randstring(['A':'Z'; 'a':'z'], 12))
end
newnames = names(mdata_ecs_df)[20:end]
open(scratchfiles("future_longnames.txt"), "w") do io
    for i in oldnames
        println(io, i)
    end
end
open(scratchfiles("future_hashnames.txt"), "w") do io
    for i in newnames
        println(io, i)
    end
end

#####
# Applying prevalence filters
#####
taxa_prevalence_threshold_lower = 0.10
taxa_prevalence_threshold_upper = 1.00
ecs_prevalence_threshold_lower = 0.25
ecs_prevalence_threshold_upper = 0.85

## 00 to 06 months

filtered_mdata_taxa_df_future = filter_prevalences(
    dropmissing(mdata_taxa_df),
    :future_cogScore,
    [:subject, :sex, :education, :present_ageMonths, :future_ageMonths],
    [:sample, :read_depth, :present_timepoint, :future_timepoint, :edfloat ],
    lbound = taxa_prevalence_threshold_lower,
    ubound = taxa_prevalence_threshold_upper
)

filtered_mdata_ecs_df_future = filter_prevalences(
    dropmissing(mdata_ecs_df),
    :future_cogScore,
    [:subject, :sex, :education, :present_ageMonths, :future_ageMonths],
    [:sample, :read_depth, :present_timepoint, :future_timepoint, :edfloat ],
    lbound = ecs_prevalence_threshold_lower,
    ubound = ecs_prevalence_threshold_upper
)
# ## Exports

isdir(tablefiles("figure3")) || mkdir(tablefiles("figure3"))

CSV.write(tablefiles("figure3", "future_taxa.csv"), filtered_mdata_taxa_df_future) # @Hugemiler what are these used for?
CSV.write(tablefiles("figure3", "future_ecs.csv"), filtered_mdata_ecs_df_future)

#####
# Training models
#####

RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree

onlydemo_tuning_space = (
    maxnodes_range = [ -1 ],
    nodesize_range = [ 5 ],
    sampsize_range =  [ 0.7 ],
    mtry_range = [ -1 ],
    ntrees_range = [ 500 ]
)

taxa_tuning_space = (
    maxnodes_range = [ -1 ],
    nodesize_range = [ 5 ],
    sampsize_range =  [ 0.7 ],
    mtry_range = [ -1 ],
    ntrees_range = [ 500 ]
)

ecs_tuning_space = (
    maxnodes_range = [ -1 ],
    nodesize_range = [ 5 ],
    sampsize_range =  [ 0.7 ],
    mtry_range = [ -1 ],
    ntrees_range = [ 500 ]
)

#####
# Future
#####

## 1. Only SES
regression_futureCogScores_onlydemo = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_futureCogScores_onlydemo",
    filtered_mdata_taxa_df_future,
    identity,
    [3, 4, 5, 6],
    :future_cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = onlydemo_tuning_space,
    train_rng=ml_rng
)
    
jldsave(modelfiles("regression_futureCogScores_onlydemo.jld"); regression_futureCogScores_onlydemo)

## 2. Only taxonomic profiles
regression_futureCogScores_onlytaxa = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_futureCogScores_onlytaxa",
    filtered_mdata_taxa_df_future,
    identity,
    collect(5:ncol(filtered_mdata_taxa_df_future)), # 4 to include age
    :future_cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = taxa_tuning_space,
    train_rng=ml_rng
)
    
jldsave(modelfiles("regression_futureCogScores_onlytaxa.jld"); regression_futureCogScores_onlytaxa)

## 3. SES + taxonomic profiles
regression_futureCogScores_demoplustaxa = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_futureCogScores_demoplustaxa",
    filtered_mdata_taxa_df_future,
    identity,
    collect(3:ncol(filtered_mdata_taxa_df_future)),
    :future_cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = taxa_tuning_space,
    train_rng=ml_rng
)

jldsave(modelfiles("regression_futureCogScores_demoplustaxa.jld"); regression_futureCogScores_demoplustaxa)

## 4. Only functional profiles
regression_futureCogScores_onlyecs = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_futureCogScores_onlyecs",
    filtered_mdata_ecs_df_future,
    identity,
    collect(5:ncol(filtered_mdata_ecs_df_future)),
    :future_cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = ecs_tuning_space,
    train_rng=ml_rng
)

jldsave(modelfiles("regression_futureCogScores_onlyecs.jld"); regression_futureCogScores_onlyecs)

## 5. SES + functional profiles
regression_futureCogScores_demoplusecs = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_futureCogScores_demoplusecs",
    filtered_mdata_ecs_df_future,
    identity,
    collect(3:ncol(filtered_mdata_ecs_df_future)),
    :future_cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = ecs_tuning_space,
    train_rng=ml_rng
)

jldsave(modelfiles("regression_futureCogScores_demoplusecs.jld"); regression_futureCogScores_demoplusecs)