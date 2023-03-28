##########
# Concurrent cognitive assessment Score regression notebook
##########

#####
# 0. Loading libraries
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

## 1. Taxonomic Profiles
taxa = Resonance.load(TaxonomicProfiles())
mdata_taxa_df = sort(Resonance.comm2wide(taxa), [ :subject, :timepoint ]);

## 2. Functional Profiles
ecs = Resonance.load(ECProfiles())
mdata_ecs_df = sort(Resonance.comm2wide(ecs), [ :subject, :timepoint ]);

### 2.1. Hashing EC names with fixed-length alphabetic hashes due to length and special characters; storing text files for remapping.
oldnames = names(mdata_ecs_df)[15:2414]
for oldname in names(mdata_ecs_df)[15:2414]
    rename!(mdata_ecs_df, oldname => randstring(['A':'Z'; 'a':'z'], 12))
end
newnames = names(mdata_ecs_df)[15:2414]
open("longnames.txt", "w") do io
    for i in oldnames
        println(io, i)
    end
end
open("hashnames.txt", "w") do io
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
    dropmissing(mdata_taxa_df[mdata_taxa_df.filter_00to06, :], Not([:race, :read_depth])),
    :cogScore,
    [:subject, :timepoint, :sex, :education, :ageMonths],
    [:sample, :sample_base, :race, :date, :read_depth, :filter_00to120, :filter_00to06, :filter_18to120];
    lbound = taxa_prevalence_threshold_lower,
    ubound = taxa_prevalence_threshold_upper
)

filtered_mdata_ecs_df_00to06 = filter_prevalences(
    dropmissing(mdata_ecs_df[mdata_ecs_df.filter_00to06, :], Not([:race, :read_depth])),
    :cogScore,
    [:subject, :timepoint, :sex, :education, :ageMonths],
    [:sample, :sample_base, :race, :date, :read_depth, :filter_00to120, :filter_00to06, :filter_18to120];
    lbound = ecs_prevalence_threshold_lower,
    ubound = ecs_prevalence_threshold_upper
)

filtered_mdata_taxa_df_18to120 = filter_prevalences(
    dropmissing(mdata_taxa_df[mdata_taxa_df.filter_18to120, :], Not([:race, :read_depth])),
    :cogScore,
    [:subject, :timepoint, :sex, :education, :ageMonths],
    [:sample, :sample_base, :race, :date, :read_depth, :filter_00to120, :filter_00to06, :filter_18to120];
    lbound = taxa_prevalence_threshold_lower,
    ubound = taxa_prevalence_threshold_upper
)

filtered_mdata_ecs_df_18to120 = filter_prevalences(
    dropmissing(mdata_ecs_df[mdata_ecs_df.filter_18to120, :], Not([:race, :read_depth])),
    :cogScore,
    [:subject, :timepoint, :sex, :education, :ageMonths],
    [:sample, :sample_base, :race, :date, :read_depth, :filter_00to120, :filter_00to06, :filter_18to120];
    lbound = ecs_prevalence_threshold_lower,
    ubound = ecs_prevalence_threshold_upper
)

# ## Exports
CSV.write("00to06_taxa.csv", filtered_mdata_taxa_df_00to06)
CSV.write("00to06_ecs.csv", filtered_mdata_ecs_df_00to06)
CSV.write("18to120_taxa.csv", filtered_mdata_taxa_df_18to120)
CSV.write("18to120_ecs.csv", filtered_mdata_ecs_df_18to120)

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
regression_currentCogScores_00to06mo_onlydemo = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_00to06mo_onlydemo",
    filtered_mdata_taxa_df_00to06,
    identity,
    [4, 5, 6],
    :cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = onlydemo_tuning_space,
    train_rng=ml_rng
)
    
jldsave(modelfiles("regression_currentCogScores_00to06mo_onlydemo.jld"), regression_currentCogScores_00to06mo_onlydemo)

## 2. Only taxonomic profiles
regression_currentCogScores_00to06mo_onlytaxa = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_00to06mo_onlytaxa",
    filtered_mdata_taxa_df_00to06,
    identity,
    collect(6:ncol(filtered_mdata_taxa_df_00to06)), # 4 to include age
    :cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = taxa_tuning_space,
    train_rng=ml_rng
)
    
jldsave(modelfiles("regression_currentCogScores_00to06mo_onlytaxa.jld"), regression_currentCogScores_00to06mo_onlytaxa)

## 3. SES + taxonomic profiles
regression_currentCogScores_00to06mo_demoplustaxa = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_00to06mo_demoplustaxa",
    filtered_mdata_taxa_df_00to06,
    identity,
    collect(4:ncol(filtered_mdata_taxa_df_00to06)),
    :cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = taxa_tuning_space,
    train_rng=ml_rng
)

jldsave(modelfiles("regression_currentCogScores_00to06mo_demoplustaxa.jld"), regression_currentCogScores_00to06mo_demoplustaxa)

## 4. Only functional profiles
regression_currentCogScores_00to06mo_onlyecs = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_00to06mo_onlyecs",
    filtered_mdata_ecs_df_00to06,
    identity,
    collect(4:ncol(filtered_mdata_ecs_df_00to06)),
    :cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = ecs_tuning_space,
    train_rng=ml_rng
)

jldsave(modelfiles("regression_currentCogScores_00to06mo_onlyecs.jld"), regression_currentCogScores_00to06mo_onlyecs)

## 5. SES + functional profiles
regression_currentCogScores_00to06mo_demoplusecs = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_00to06mo_demoplusecs",
    filtered_mdata_ecs_df_00to06,
    identity,
    collect(4:ncol(filtered_mdata_ecs_df_00to06)),
    :cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = ecs_tuning_space,
    train_rng=ml_rng
)

jldsave(modelfiles("regression_currentCogScores_00to06mo_demoplusecs.jld"), regression_currentCogScores_00to06mo_demoplusecs)

#####
# 18 to 120 months
#####

# ## 6. Only SES
regression_currentCogScores_18to120mo_onlydemo = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_18to120mo_onlydemo",
    filtered_mdata_taxa_df_18to120,
    identity,
    [4, 5, 6],
    :cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = onlydemo_tuning_space,
    train_rng=ml_rng
)
    
jldsave(modelfiles("regression_currentCogScores_18to120mo_onlydemo.jld"), regression_currentCogScores_18to120mo_onlydemo)

## 7. Only taxonomic profiles
regression_currentCogScores_18to120mo_onlytaxa = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_18to120mo_onlytaxa",
    filtered_mdata_taxa_df_18to120,
    identity,
    collect(6:ncol(filtered_mdata_taxa_df_18to120)),
    :cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = taxa_tuning_space,
    train_rng=ml_rng
)
    
jldsave(modelfiles("regression_currentCogScores_18to120mo_onlytaxa.jld"), regression_currentCogScores_18to120mo_onlytaxa)

## 8. SES + taxonomic profiles
regression_currentCogScores_18to120mo_demoplustaxa = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_18to120mo_demoplustaxa",
    filtered_mdata_taxa_df_18to120,
    identity,
    collect(4:ncol(filtered_mdata_taxa_df_18to120)),
    :cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = taxa_tuning_space,
    train_rng=ml_rng
)

jldsave(modelfiles("regression_currentCogScores_18to120mo_demoplustaxa.jld"), regression_currentCogScores_18to120mo_demoplustaxa)

## 9. Only functional profiles
regression_currentCogScores_18to120mo_onlyecs = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_18to120mo_onlyecs",
    filtered_mdata_ecs_df_18to120,
    identity,
    collect(6:ncol(filtered_mdata_ecs_df_18to120)),
    :cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = ecs_tuning_space,
    train_rng=ml_rng
)

jldsave(modelfiles("regression_currentCogScores_18to120mo_onlyecs.jld"), regression_currentCogScores_18to120mo_onlyecs)

## 10. SES + functional profiles
regression_currentCogScores_18to120mo_demoplusecs = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_18to120mo_demoplusecs",
    filtered_mdata_ecs_df_18to120,
    identity,
    collect(4:ncol(filtered_mdata_ecs_df_18to120)),
    :cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = ecs_tuning_space,
    train_rng=ml_rng
)

jldsave(modelfiles("regression_currentCogScores_18to120mo_demoplusecs.jld"), regression_currentCogScores_18to120mo_demoplusecs)