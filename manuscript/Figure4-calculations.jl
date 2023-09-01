#####
# Notebook DN01 - Regression of binary current CogScores above/below 50th percentile -- NULL DISTRIBUTION
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
using CategoricalArrays
ml_rng = StableRNG(0)

isdir(tablefiles("figure4")) || mkdir(tablefiles("figure4"))

#####
# Loading Data
#####

## 1. Metadata
mdata = Resonance.load(Metadata())
seqs = subset(mdata, "sample"=> ByRow(!ismissing)) 
transform!(seqs, "sample"=> ByRow(String)=> "sample")
seqs.edfloat = map(x-> ismissing(x) ? missing : Float64(levelcode(x)), seqs.education)

## 2. Taxonomic Profiles
taxa = Resonance.load(TaxonomicProfiles(); timepoint_metadata = seqs) # this can take a bit
species = filter(f-> taxrank(f) == :species, taxa)
mdata_taxa_df = sort(Resonance.comm2wide(species), [ :subject, :timepoint ]);

## 2. Brain data
brain = Resonance.load(Neuroimaging())
mdata_brain_df = @chain Resonance.comm2wide(brain; feature_func = string) begin
    sort([ :subject, :timepoint ]);
    select!(Not([
        "hassample",
        "subject",
        "timepoint",
        "ageMonths",
        "sex",
        "education",
        "cogScore",
        "race",
        "omni",
        "Mullen::mullen_ExpressivelanguageT",
        "Mullen::mullen_FineMotorT",
        "Mullen::mullen_GrossMotorT",
        "Mullen::mullen_ReceptiveLanguageT",
        "Mullen::mullen_VisualReceptionT"]
    ))
end
# TBV = mdata_brain_df."White-matter" .+ mdata_brain_df."Gray-matter" # This is not necessary anymore, because the columns are already TBV normalized on the input.

for col in names(mdata_brain_df)[2:end]
    mdata_brain_df[!, col] = 100.0 * mdata_brain_df[:, col]
end

## Filter taxonomic profile rows to match brain samples from 18 to 20 months
mdata_taxa_brain_df = innerjoin(mdata_taxa_df, select!(mdata_brain_df, Not([:filter_00to120, :filter_00to06, :filter_18to120, :filter_future_omni, :filter_future_cog])); on = :sample)

#####
# Applying prevalence filters
#####
taxa_prevalence_threshold_lower = 0.10
taxa_prevalence_threshold_upper = 1.00

filtered_mdata_taxa_brain_df_18to120 = filter_prevalences(
    dropmissing(mdata_taxa_brain_df[mdata_taxa_brain_df.filter_18to120, :], Not(["sample", "read_depth", "edfloat", "race", "filter_00to120", "filter_00to06", "filter_18to120", "filter_future_omni", "filter_future_cog", "omni", "Mullen::mullen_ExpressivelanguageT", "Mullen::mullen_FineMotorT", "Mullen::mullen_GrossMotorT", "Mullen::mullen_ReceptiveLanguageT", "Mullen::mullen_VisualReceptionT"])),
    :cogScore,
    [:subject, :timepoint, :sex, :education, :ageMonths],
    Symbol.(["sample", "read_depth", "edfloat", "race", "filter_00to120", "filter_00to06", "filter_18to120", "filter_future_omni", "filter_future_cog", "omni", "Mullen::mullen_ExpressivelanguageT", "Mullen::mullen_FineMotorT", "Mullen::mullen_GrossMotorT", "Mullen::mullen_ReceptiveLanguageT", "Mullen::mullen_VisualReceptionT"]);
    lbound = taxa_prevalence_threshold_lower,
    ubound = taxa_prevalence_threshold_upper
)[:, 2:end]

CSV.write(tablefiles("figure4", "18to120_brain_taxa.csv"), filtered_mdata_taxa_brain_df_18to120)

inputs_taxa = filtered_mdata_taxa_brain_df_18to120[:, 1:166] # Ends on "Ruminococcus_sp_CAG_177" 

#####
# Training models
#####

RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree

taxa_tuning_space = (
    maxnodes_range = [ -1 ],
    nodesize_range = [ 5 ],
    sampsize_range =  [ 0.7 ],
    mtry_range = [ -1 ],
    ntrees_range = [ 500 ]
)

#####
# 18 to 120 months
#####

# Summary of brain segment aggregations:
# Removed "left-lateral-ventricle", "right-lateral-ventricle", "left-inferior-lateral-ventricle", "right-inferior-lateral-ventricle", because vessels are and not necessarily interesting;
# Removed "left-cerebellum-exterior", "right-cerebellum-exterior"; not useful and bad predictions.
# Removed "left-ventral-DC", "right-ventral-DC". See http://neuromorphometrics.com/Seg/html/segmentation/ventral%20diencephalon.html
# Combined "superior-temporal", "middle-temporal", "transverse-temporal" and "inferior-temporal" into "temporal"
filtered_mdata_taxa_brain_df_18to120."left-temporal" = filtered_mdata_taxa_brain_df_18to120."left-superior-temporal" .+ filtered_mdata_taxa_brain_df_18to120."left-middle-temporal" .+ filtered_mdata_taxa_brain_df_18to120."left-transverse-temporal" .+ filtered_mdata_taxa_brain_df_18to120."left-inferior-temporal"
filtered_mdata_taxa_brain_df_18to120."right-temporal" = filtered_mdata_taxa_brain_df_18to120."right-superior-temporal" .+ filtered_mdata_taxa_brain_df_18to120."right-middle-temporal" .+ filtered_mdata_taxa_brain_df_18to120."right-transverse-temporal" .+ filtered_mdata_taxa_brain_df_18to120."right-inferior-temporal"
# Combined "lateral-orbitofrontal" and "medial-orbitofrontal" into "orbitofrontal"
filtered_mdata_taxa_brain_df_18to120."left-orbitofrontal" = filtered_mdata_taxa_brain_df_18to120."left-lateral-orbitofrontal" .+ filtered_mdata_taxa_brain_df_18to120."left-medial-orbitofrontal"
filtered_mdata_taxa_brain_df_18to120."right-orbitofrontal" = filtered_mdata_taxa_brain_df_18to120."right-lateral-orbitofrontal" .+ filtered_mdata_taxa_brain_df_18to120."right-medial-orbitofrontal"
# Combined "inferior-parietal" and "superior-parietal" into "parietal"
filtered_mdata_taxa_brain_df_18to120."left-parietal" = filtered_mdata_taxa_brain_df_18to120."left-inferior-parietal" .+ filtered_mdata_taxa_brain_df_18to120."left-superior-parietal"
filtered_mdata_taxa_brain_df_18to120."right-parietal" = filtered_mdata_taxa_brain_df_18to120."right-inferior-parietal" .+ filtered_mdata_taxa_brain_df_18to120."right-superior-parietal"
# Combined "caudal-middle-frontal" and "rostral-middle-frontal" into "middle-frontal"
filtered_mdata_taxa_brain_df_18to120."left-middle-frontal" = filtered_mdata_taxa_brain_df_18to120."left-caudal-middle-frontal" .+ filtered_mdata_taxa_brain_df_18to120."left-rostral-middle-frontal"
filtered_mdata_taxa_brain_df_18to120."right-middle-frontal" = filtered_mdata_taxa_brain_df_18to120."right-caudal-middle-frontal" .+ filtered_mdata_taxa_brain_df_18to120."right-rostral-middle-frontal"
# Combined "caudal-anterior-cingulate" and "ostral-anterior-cingulate" into "anterior-cingulate"
filtered_mdata_taxa_brain_df_18to120."left-anterior-cingulate" = filtered_mdata_taxa_brain_df_18to120."left-caudal-anterior-cingulate" .+ filtered_mdata_taxa_brain_df_18to120."left-rostral-anterior-cingulate"
filtered_mdata_taxa_brain_df_18to120."right-anterior-cingulate" = filtered_mdata_taxa_brain_df_18to120."right-caudal-anterior-cingulate" .+ filtered_mdata_taxa_brain_df_18to120."right-rostral-anterior-cingulate"

include("manuscript/Figure4-definitions.jl")

# ## 7. Only taxonomic profiles

brain_models = Dict{String, ProbeData}()

for brain_segment in ordered_brain_segments_list

    tmp_df = deepcopy(inputs_taxa)
    prediction_df = insertcols!(tmp_df, 1, :target => filtered_mdata_taxa_brain_df_18to120[:, brain_segment])
    model_name = "regression_"*string(brain_segment)*"_18to120mo_onlytaxa"

    push!(
        brain_models,
        brain_segment => probe_prod_randomforest(
            Resonance.Regression(),
            model_name,
            prediction_df,
            identity,
            collect(4:ncol(prediction_df)),
            :target;
            n_folds = 3,
            n_replicas = 100,
            n_rngs = 10,
            tuning_space = taxa_tuning_space,
            train_rng=ml_rng
        )
    )

    println("finished $(string(brain_segment))!")

end
    
jldsave(modelfiles("brain_models.jld"); brain_models)