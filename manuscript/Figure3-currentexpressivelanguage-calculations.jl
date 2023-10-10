##########
# Concurrent cognitive assessment Score regression notebook
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

#####
# Loading Data
#####

## 1. Metadata
mdata = Resonance.load(Metadata())
rename!(mdata, "Mullen::mullen_ExpressivelanguageT" => "ExpressiveLanguage")
seqs = subset(mdata, "sample"=> ByRow(!ismissing)) 
transform!(seqs, "sample"=> ByRow(String)=> "sample")
seqs.edfloat = map(x-> ismissing(x) ? missing : Float64(levelcode(x)), seqs.education)

## 2. Taxonomic Profiles
taxa = Resonance.load(TaxonomicProfiles(); timepoint_metadata = seqs) # this can take a bit
species = filter(f-> taxrank(f) == :species, taxa)
mdata_taxa_df = sort(Resonance.comm2wide(species), [ :subject, :timepoint ]);

## 3. Functional Profiles
ecs = Resonance.load(ECProfiles(); timepoint_metadata = seqs)
mdata_ecs_df = sort(Resonance.comm2wide(ecs), [ :subject, :timepoint ]);

### 3.1. Hashing EC names with fixed-length alphabetic hashes due to length and special characters; storing text files for remapping.
oldnames = names(mdata_ecs_df)[22:end]
for oldname in names(mdata_ecs_df)[22:end]
    rename!(mdata_ecs_df, oldname => randstring(['A':'Z'; 'a':'z'], 12))
end
newnames = names(mdata_ecs_df)[22:end]
open(scratchfiles("concurrentEL_longnames.txt"), "w") do io
    for i in oldnames
        println(io, i)
    end
end
open(scratchfiles("concurrentEL_hashnames.txt"), "w") do io
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

filtered_mdata_taxa_df_00to06 = filter_prevalences(
    dropmissing(mdata_taxa_df[mdata_taxa_df.filter_00to06, :], Not(["sample", "read_depth", "edfloat", "race", "filter_00to120", "filter_00to06", "filter_18to120", "filter_future_omni", "filter_future_cog", "omni", "cogScore", "Mullen::mullen_FineMotorT", "Mullen::mullen_GrossMotorT", "Mullen::mullen_ReceptiveLanguageT", "Mullen::mullen_VisualReceptionT"])),
    :ExpressiveLanguage,
    [:subject, :timepoint, :sex, :education, :ageMonths],
    Symbol.(["sample", "read_depth", "edfloat", "race", "filter_00to120", "filter_00to06", "filter_18to120", "filter_future_omni", "filter_future_cog", "omni", "cogScore", "Mullen::mullen_FineMotorT", "Mullen::mullen_GrossMotorT", "Mullen::mullen_ReceptiveLanguageT", "Mullen::mullen_VisualReceptionT"]);
    lbound = taxa_prevalence_threshold_lower,
    ubound = taxa_prevalence_threshold_upper
)

filtered_mdata_ecs_df_00to06 = filter_prevalences(
    dropmissing(mdata_ecs_df[mdata_ecs_df.filter_00to06, :], Not(["sample", "read_depth", "edfloat", "race", "filter_00to120", "filter_00to06", "filter_18to120", "filter_future_omni", "filter_future_cog", "omni", "cogScore", "Mullen::mullen_FineMotorT", "Mullen::mullen_GrossMotorT", "Mullen::mullen_ReceptiveLanguageT", "Mullen::mullen_VisualReceptionT"])),
    :ExpressiveLanguage,
    [:subject, :timepoint, :sex, :education, :ageMonths],
    Symbol.(["sample", "read_depth", "edfloat", "race", "filter_00to120", "filter_00to06", "filter_18to120", "filter_future_omni", "filter_future_cog", "omni", "cogScore", "Mullen::mullen_FineMotorT", "Mullen::mullen_GrossMotorT", "Mullen::mullen_ReceptiveLanguageT", "Mullen::mullen_VisualReceptionT"]);
    lbound = ecs_prevalence_threshold_lower,
    ubound = ecs_prevalence_threshold_upper
)

filtered_mdata_taxa_df_18to120 = filter_prevalences(
    dropmissing(mdata_taxa_df[mdata_taxa_df.filter_18to120, :], Not(["sample", "read_depth", "edfloat", "race", "filter_00to120", "filter_00to06", "filter_18to120", "filter_future_omni", "filter_future_cog", "omni", "cogScore", "Mullen::mullen_FineMotorT", "Mullen::mullen_GrossMotorT", "Mullen::mullen_ReceptiveLanguageT", "Mullen::mullen_VisualReceptionT"])),
    :ExpressiveLanguage,
    [:subject, :timepoint, :sex, :education, :ageMonths],
    Symbol.(["sample", "read_depth", "edfloat", "race", "filter_00to120", "filter_00to06", "filter_18to120", "filter_future_omni", "filter_future_cog", "omni", "cogScore", "Mullen::mullen_FineMotorT", "Mullen::mullen_GrossMotorT", "Mullen::mullen_ReceptiveLanguageT", "Mullen::mullen_VisualReceptionT"]);
    lbound = taxa_prevalence_threshold_lower,
    ubound = taxa_prevalence_threshold_upper
)

filtered_mdata_ecs_df_18to120 = filter_prevalences(
    dropmissing(mdata_ecs_df[mdata_ecs_df.filter_18to120, :], Not(["sample", "read_depth", "edfloat", "race", "filter_00to120", "filter_00to06", "filter_18to120", "filter_future_omni", "filter_future_cog", "omni", "cogScore", "Mullen::mullen_FineMotorT", "Mullen::mullen_GrossMotorT", "Mullen::mullen_ReceptiveLanguageT", "Mullen::mullen_VisualReceptionT"])),
    :ExpressiveLanguage,
    [:subject, :timepoint, :sex, :education, :ageMonths],
    Symbol.(["sample", "read_depth", "edfloat", "race", "filter_00to120", "filter_00to06", "filter_18to120", "filter_future_omni", "filter_future_cog", "omni", "cogScore", "Mullen::mullen_FineMotorT", "Mullen::mullen_GrossMotorT", "Mullen::mullen_ReceptiveLanguageT", "Mullen::mullen_VisualReceptionT"]);
    lbound = ecs_prevalence_threshold_lower,
    ubound = ecs_prevalence_threshold_upper
)

filtered_mdata_taxa_df_00to120 = filter_prevalences(
    dropmissing(mdata_taxa_df[mdata_taxa_df.filter_00to120, :], Not(["sample", "read_depth", "edfloat", "race", "filter_00to120", "filter_00to06", "filter_18to120", "filter_future_omni", "filter_future_cog", "omni", "cogScore", "Mullen::mullen_FineMotorT", "Mullen::mullen_GrossMotorT", "Mullen::mullen_ReceptiveLanguageT", "Mullen::mullen_VisualReceptionT"])),
    :ExpressiveLanguage,
    [:subject, :timepoint, :sex, :education, :ageMonths],
    Symbol.(["sample", "read_depth", "edfloat", "race", "filter_00to120", "filter_00to06", "filter_18to120", "filter_future_omni", "filter_future_cog", "omni", "cogScore", "Mullen::mullen_FineMotorT", "Mullen::mullen_GrossMotorT", "Mullen::mullen_ReceptiveLanguageT", "Mullen::mullen_VisualReceptionT"]);
    lbound = taxa_prevalence_threshold_lower,
    ubound = taxa_prevalence_threshold_upper
)

filtered_mdata_ecs_df_00to120 = filter_prevalences(
    dropmissing(mdata_ecs_df[mdata_ecs_df.filter_00to120, :], Not(["sample", "read_depth", "edfloat", "race", "filter_00to120", "filter_00to06", "filter_18to120", "filter_future_omni", "filter_future_cog", "omni", "cogScore", "Mullen::mullen_FineMotorT", "Mullen::mullen_GrossMotorT", "Mullen::mullen_ReceptiveLanguageT", "Mullen::mullen_VisualReceptionT"])),
    :ExpressiveLanguage,
    [:subject, :timepoint, :sex, :education, :ageMonths],
    Symbol.(["sample", "read_depth", "edfloat", "race", "filter_00to120", "filter_00to06", "filter_18to120", "filter_future_omni", "filter_future_cog", "omni", "cogScore", "Mullen::mullen_FineMotorT", "Mullen::mullen_GrossMotorT", "Mullen::mullen_ReceptiveLanguageT", "Mullen::mullen_VisualReceptionT"]);
    lbound = ecs_prevalence_threshold_lower,
    ubound = ecs_prevalence_threshold_upper
)

transform!(filtered_mdata_taxa_df_00to06, :ExpressiveLanguage => (x -> Float64.(x)) => :ExpressiveLanguage; renamecols = false)
transform!(filtered_mdata_ecs_df_00to06, :ExpressiveLanguage => (x -> Float64.(x)) => :ExpressiveLanguage; renamecols = false)
transform!(filtered_mdata_taxa_df_18to120, :ExpressiveLanguage => (x -> Float64.(x)) => :ExpressiveLanguage; renamecols = false)
transform!(filtered_mdata_ecs_df_18to120, :ExpressiveLanguage => (x -> Float64.(x)) => :ExpressiveLanguage; renamecols = false)
transform!(filtered_mdata_taxa_df_00to120, :ExpressiveLanguage => (x -> Float64.(x)) => :ExpressiveLanguage; renamecols = false)
transform!(filtered_mdata_ecs_df_00to120, :ExpressiveLanguage => (x -> Float64.(x)) => :ExpressiveLanguage; renamecols = false)

# ## Exports

isdir(tablefiles("figure3")) || mkdir(tablefiles("figure3"))

CSV.write(tablefiles("figure3", "expressivelanguage_00to06_taxa.csv"), filtered_mdata_taxa_df_00to06) # @Hugemiler what are these used for?
CSV.write(tablefiles("figure3", "expressivelanguage_00to06_ecs.csv"), filtered_mdata_ecs_df_00to06)
CSV.write(tablefiles("figure3", "expressivelanguage_18to120_taxa.csv"), filtered_mdata_taxa_df_18to120)
CSV.write(tablefiles("figure3", "expressivelanguage_18to120_ecs.csv"), filtered_mdata_ecs_df_18to120)
CSV.write(tablefiles("figure3", "expressivelanguage_00to120_taxa.csv"), filtered_mdata_taxa_df_00to120)
CSV.write(tablefiles("figure3", "expressivelanguage_00to120_ecs.csv"), filtered_mdata_ecs_df_00to120)

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
# 00 to 06 months
#####

## 1. Only SES
regression_currentExpressiveLanguages_00to06mo_onlydemo = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentExpressiveLanguages_00to06mo_onlydemo",
    filtered_mdata_taxa_df_00to06,
    identity,
    [4, 5, 6],
    :ExpressiveLanguage;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = onlydemo_tuning_space,
    train_rng=ml_rng
)
    
jldsave(modelfiles("regression_currentExpressiveLanguages_00to06mo_onlydemo.jld"); regression_currentExpressiveLanguages_00to06mo_onlydemo)

## 2. Only taxonomic profiles
regression_currentExpressiveLanguages_00to06mo_onlytaxa = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentExpressiveLanguages_00to06mo_onlytaxa",
    filtered_mdata_taxa_df_00to06,
    identity,
    collect(6:ncol(filtered_mdata_taxa_df_00to06)), # 4 to include age
    :ExpressiveLanguage;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = taxa_tuning_space,
    train_rng=ml_rng
)
    
jldsave(modelfiles("regression_currentExpressiveLanguages_00to06mo_onlytaxa.jld"); regression_currentExpressiveLanguages_00to06mo_onlytaxa)

## 3. SES + taxonomic profiles
regression_currentExpressiveLanguages_00to06mo_demoplustaxa = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentExpressiveLanguages_00to06mo_demoplustaxa",
    filtered_mdata_taxa_df_00to06,
    identity,
    collect(4:ncol(filtered_mdata_taxa_df_00to06)),
    :ExpressiveLanguage;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = taxa_tuning_space,
    train_rng=ml_rng
)

jldsave(modelfiles("regression_currentExpressiveLanguages_00to06mo_demoplustaxa.jld"); regression_currentExpressiveLanguages_00to06mo_demoplustaxa)

## 4. Only functional profiles
regression_currentExpressiveLanguages_00to06mo_onlyecs = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentExpressiveLanguages_00to06mo_onlyecs",
    filtered_mdata_ecs_df_00to06,
    identity,
    collect(6:ncol(filtered_mdata_ecs_df_00to06)),
    :ExpressiveLanguage;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = ecs_tuning_space,
    train_rng=ml_rng
)

jldsave(modelfiles("regression_currentExpressiveLanguages_00to06mo_onlyecs.jld"); regression_currentExpressiveLanguages_00to06mo_onlyecs)

## 5. SES + functional profiles
regression_currentExpressiveLanguages_00to06mo_demoplusecs = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentExpressiveLanguages_00to06mo_demoplusecs",
    filtered_mdata_ecs_df_00to06,
    identity,
    collect(4:ncol(filtered_mdata_ecs_df_00to06)),
    :ExpressiveLanguage;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = ecs_tuning_space,
    train_rng=ml_rng
)

jldsave(modelfiles("regression_currentExpressiveLanguages_00to06mo_demoplusecs.jld"); regression_currentExpressiveLanguages_00to06mo_demoplusecs)

#####
# 18 to 120 months
#####

# ## 6. Only SES
regression_currentExpressiveLanguages_18to120mo_onlydemo = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentExpressiveLanguages_18to120mo_onlydemo",
    filtered_mdata_taxa_df_18to120,
    identity,
    [4, 5, 6],
    :ExpressiveLanguage;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = onlydemo_tuning_space,
    train_rng=ml_rng
)
    
jldsave(modelfiles("regression_currentExpressiveLanguages_18to120mo_onlydemo.jld"); regression_currentExpressiveLanguages_18to120mo_onlydemo)

## 7. Only taxonomic profiles
regression_currentExpressiveLanguages_18to120mo_onlytaxa = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentExpressiveLanguages_18to120mo_onlytaxa",
    filtered_mdata_taxa_df_18to120,
    identity,
    collect(6:ncol(filtered_mdata_taxa_df_18to120)),
    :ExpressiveLanguage;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = taxa_tuning_space,
    train_rng=ml_rng
)
    
jldsave(modelfiles("regression_currentExpressiveLanguages_18to120mo_onlytaxa.jld"); regression_currentExpressiveLanguages_18to120mo_onlytaxa)

## 8. SES + taxonomic profiles
regression_currentExpressiveLanguages_18to120mo_demoplustaxa = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentExpressiveLanguages_18to120mo_demoplustaxa",
    filtered_mdata_taxa_df_18to120,
    identity,
    collect(4:ncol(filtered_mdata_taxa_df_18to120)),
    :ExpressiveLanguage;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = taxa_tuning_space,
    train_rng=ml_rng
)

jldsave(modelfiles("regression_currentExpressiveLanguages_18to120mo_demoplustaxa.jld"); regression_currentExpressiveLanguages_18to120mo_demoplustaxa)

## 9. Only functional profiles
regression_currentExpressiveLanguages_18to120mo_onlyecs = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentExpressiveLanguages_18to120mo_onlyecs",
    filtered_mdata_ecs_df_18to120,
    identity,
    collect(6:ncol(filtered_mdata_ecs_df_18to120)),
    :ExpressiveLanguage;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = ecs_tuning_space,
    train_rng=ml_rng
)

jldsave(modelfiles("regression_currentExpressiveLanguages_18to120mo_onlyecs.jld"); regression_currentExpressiveLanguages_18to120mo_onlyecs)

## 10. SES + functional profiles
regression_currentExpressiveLanguages_18to120mo_demoplusecs = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentExpressiveLanguages_18to120mo_demoplusecs",
    filtered_mdata_ecs_df_18to120,
    identity,
    collect(4:ncol(filtered_mdata_ecs_df_18to120)),
    :ExpressiveLanguage;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = ecs_tuning_space,
    train_rng=ml_rng
)

jldsave(modelfiles("regression_currentExpressiveLanguages_18to120mo_demoplusecs.jld"); regression_currentExpressiveLanguages_18to120mo_demoplusecs)

#####
# 18 to 120 months
#####

## 6. Only SES
regression_currentExpressiveLanguages_00to120mo_onlydemo = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentExpressiveLanguages_00to120mo_onlydemo",
    filtered_mdata_taxa_df_00to120,
    identity,
    [4, 5, 6],
    :ExpressiveLanguage;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = onlydemo_tuning_space,
    train_rng=ml_rng
)
    
jldsave(modelfiles("regression_currentExpressiveLanguages_00to120mo_onlydemo.jld"); regression_currentExpressiveLanguages_00to120mo_onlydemo)

## 7. Only taxonomic profiles
regression_currentExpressiveLanguages_00to120mo_onlytaxa = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentExpressiveLanguages_00to120mo_onlytaxa",
    filtered_mdata_taxa_df_00to120,
    identity,
    collect(6:ncol(filtered_mdata_taxa_df_00to120)),
    :ExpressiveLanguage;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = taxa_tuning_space,
    train_rng=ml_rng
)
    
jldsave(modelfiles("regression_currentExpressiveLanguages_00to120mo_onlytaxa.jld"); regression_currentExpressiveLanguages_00to120mo_onlytaxa)

## 8. SES + taxonomic profiles
regression_currentExpressiveLanguages_00to120mo_demoplustaxa = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentExpressiveLanguages_00to120mo_demoplustaxa",
    filtered_mdata_taxa_df_00to120,
    identity,
    collect(4:ncol(filtered_mdata_taxa_df_00to120)),
    :ExpressiveLanguage;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = taxa_tuning_space,
    train_rng=ml_rng
)

jldsave(modelfiles("regression_currentExpressiveLanguages_00to120mo_demoplustaxa.jld"); regression_currentExpressiveLanguages_00to120mo_demoplustaxa)

## 9. Only functional profiles
regression_currentExpressiveLanguages_00to120mo_onlyecs = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentExpressiveLanguages_00to120mo_onlyecs",
    filtered_mdata_ecs_df_00to120,
    identity,
    collect(6:ncol(filtered_mdata_ecs_df_00to120)),
    :ExpressiveLanguage;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = ecs_tuning_space,
    train_rng=ml_rng
)

jldsave(modelfiles("regression_currentExpressiveLanguages_00to120mo_onlyecs.jld"); regression_currentExpressiveLanguages_00to120mo_onlyecs)

## 10. SES + functional profiles
regression_currentExpressiveLanguages_00to120mo_demoplusecs = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentExpressiveLanguages_00to120mo_demoplusecs",
    filtered_mdata_ecs_df_00to120,
    identity,
    collect(4:ncol(filtered_mdata_ecs_df_00to120)),
    :ExpressiveLanguage;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = ecs_tuning_space,
    train_rng=ml_rng
)

jldsave(modelfiles("regression_currentExpressiveLanguages_00to120mo_demoplusecs.jld"); regression_currentExpressiveLanguages_00to120mo_demoplusecs)