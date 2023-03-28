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
ml_rng = StableRNG(0)

#####
# Loading Data
#####

## 1. Taxonomic Profiles
taxa = Resonance.load(TaxonomicProfiles())
mdata_taxa_df = sort(Resonance.comm2wide(taxa), [ :subject, :timepoint ]);

## 2. Brain data
brain = Resonance.load(Neuroimaging())
mdata_brain_df = sort(Resonance.comm2wide(brain; feature_func = string), [ :subject, :timepoint ]);
mdata_brain_df = dropmissing(mdata_brain_df[mdata_brain_df.filter_18to120, :])
TBV = mdata_brain_df."White-matter" .+ mdata_brain_df."Gray-matter"

for col in names(mdata_brain_df)[16:end-2]
    mdata_brain_df[!, col] = 100 * mdata_brain_df[:, col] ./ TBV
end

## Filter taxonomic profile rows to match brain samples from 18 to 20 months

taxa_row_indexes = [ findfirst(mdata_taxa_df.sample .== el) for el in mdata_brain_df.sample ]
mdata_taxa_df = mdata_taxa_df[taxa_row_indexes, :]

#####
# Applying prevalence filters
#####
taxa_prevalence_threshold_lower = 0.10
taxa_prevalence_threshold_upper = 1.00

filtered_mdata_taxa_df_18to120 = filter_prevalences(
    dropmissing(mdata_taxa_df[mdata_taxa_df.filter_18to120, :], Not([:race, :read_depth])),
    :cogScore,
    [:subject, :timepoint, :sex, :education, :ageMonths],
    [:sample, :sample_base, :race, :date, :read_depth, :filter_00to120, :filter_00to06, :filter_18to120];
    lbound = taxa_prevalence_threshold_lower,
    ubound = taxa_prevalence_threshold_upper
)[:, 2:end]

CSV.write("18to120_brain_taxa.csv", hcat(filtered_mdata_taxa_df_18to120, mdata_brain_df[:, 16:end-2]))

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
mdata_brain_df."left-temporal" = mdata_brain_df."left-superior-temporal" .+ mdata_brain_df."left-middle-temporal" .+ mdata_brain_df."left-transverse-temporal" .+ mdata_brain_df."left-inferior-temporal"
mdata_brain_df."right-temporal" = mdata_brain_df."right-superior-temporal" .+ mdata_brain_df."right-middle-temporal" .+ mdata_brain_df."right-transverse-temporal" .+ mdata_brain_df."right-inferior-temporal"
# Combined "lateral-orbitofrontal" and "medial-orbitofrontal" into "orbitofrontal"
mdata_brain_df."left-orbitofrontal" = mdata_brain_df."left-lateral-orbitofrontal" .+ mdata_brain_df."left-medial-orbitofrontal"
mdata_brain_df."right-orbitofrontal" = mdata_brain_df."right-lateral-orbitofrontal" .+ mdata_brain_df."right-medial-orbitofrontal"
# Combined "inferior-parietal" and "superior-parietal" into "parietal"
mdata_brain_df."left-parietal" = mdata_brain_df."left-inferior-parietal" .+ mdata_brain_df."left-superior-parietal"
mdata_brain_df."right-parietal" = mdata_brain_df."right-inferior-parietal" .+ mdata_brain_df."right-superior-parietal"
# Combined "caudal-middle-frontal" and "rostral-middle-frontal" into "middle-frontal"
mdata_brain_df."left-middle-frontal" = mdata_brain_df."left-caudal-middle-frontal" .+ mdata_brain_df."left-rostral-middle-frontal"
mdata_brain_df."right-middle-frontal" = mdata_brain_df."right-caudal-middle-frontal" .+ mdata_brain_df."right-rostral-middle-frontal"
# Combined "caudal-anterior-cingulate" and "ostral-anterior-cingulate" into "anterior-cingulate"
mdata_brain_df."left-anterior-cingulate" = mdata_brain_df."left-caudal-anterior-cingulate" .+ mdata_brain_df."left-rostral-anterior-cingulate"
mdata_brain_df."right-anterior-cingulate" = mdata_brain_df."right-caudal-anterior-cingulate" .+ mdata_brain_df."right-rostral-anterior-cingulate"

brain_segments_list = [
    "left-temporal", "right-temporal",
    "left-orbitofrontal", "right-orbitofrontal",
    "left-parietal", "right-parietal",
    "left-middle-frontal", "right-middle-frontal",
    "left-anterior-cingulate", "right-anterior-cingulate",
    "left-lateral-occipital", "right-lateral-occipital",
    "left-cerebellum-white-matter", "right-cerebellum-white-matter",
    "left-thalamus-proper", "right-thalamus-proper",
    "left-caudate", "right-caudate",
    "left-putamen", "right-putamen",
    "left-pallidum", "right-pallidum",
    "left-hippocampus", "right-hippocampus",
    "left-amygdala", "right-amygdala",
    "left-accumbens-area", "right-accumbens-area",
    "left-basal-forebrain", "right-basal-forebrain", 
    "left-cuneus", "right-cuneus", 
    "left-entorhinal", "right-entorhinal",
    "left-fusiform", "right-fusiform",
    "left-isthmus-cingulate", "right-isthmus-cingulate",
    "left-lingual", "right-lingual",
    "left-parahippocampal", "right-parahippocampal",
    "left-paracentral", "right-paracentral",
    "left-pars-opercularis", "right-pars-opercularis",
    "left-pars-orbitalis", "right-pars-orbitalis",
    "left-pars-triangularis", "right-pars-triangularis",
    "left-pericalcarine", "right-pericalcarine",
    "left-postcentral", "right-postcentral",
    "left-posterior-cingulate", "right-posterior-cingulate",
    "left-precentral", "right-precentral",
    "left-precuneus", "right-precuneus",
    "left-superior-frontal", "right-superior-frontal",
    "left-supramarginal", "right-supramarginal",
    "left-insula", "right-insula",
    "cerebellar-vermal-lobules-I-V", "cerebellar-vermal-lobules-VI-VII", "cerebellar-vermal-lobules-VIII-X",
    "Brain-stem", "CSF"
]

# ## 7. Only taxonomic profiles

brain_models = Dict{String, ProbeData}()

for brain_segment in brain_segments_list

    tmp_df = deepcopy(filtered_mdata_taxa_df_18to120)
    prediction_df = insertcols!(tmp_df, 1, :target => mdata_brain_df[:, brain_segment])
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
    
JLD2.@save "models/brain_models.jld" brain_models

@show [ mean(brain_models[el].merits.Test_Cor) for el in brain_segments_list ]

mean_brain_merits = sort(
    reduce(
        vcat,
        [
            DataFrame(
                :variable => i,
                :Train_RMSE => mean(j.merits.Train_RMSE),
                :Test_RMSE => mean(j.merits.Test_RMSE),
                :Train_MAPE => mean(j.merits.Train_MAPE),
                :Test_MAPE => mean(j.merits.Test_MAPE),
                :Train_Cor => mean(j.merits.Train_Cor),
                :Test_Cor => mean(j.merits.Test_Cor)
            ) for (i,j) in brain_models
        ]
    ), :Test_Cor, rev=true
)

CSV.write("brain_segmented_merits_SuppTable.csv", mean_brain_merits)

combine(mean_brain_merits,
    [ :Test_Cor => mean => :Test_Cor_mean,
    :Test_Cor => std => :Test_Cor_std,
    :Test_Cor => maximum => :Test_Cor_maximum,
    :Test_Cor => (x -> quantile(x, 0.75)) => :Test_Cor_upperquartile,
    :Test_Cor => (x -> quantile(x, 0.50)) => :Test_Cor_median,
    :Test_Cor => (x -> quantile(x, 0.25)) => :Test_Cor_lowerquartile,
    :Test_Cor => minimum => :Test_Cor_minimum ]
    )


combine(mean_brain_merits,
    [
        :Test_MAPE => mean => :Test_MAPE_mean,
        :Test_MAPE => std => :Test_MAPE_std,
        :Test_MAPE => maximum => :Test_MAPE_maximum,
        :Test_MAPE => (x -> quantile(x, 0.75)) => :Test_MAPE_upperquartile,
        :Test_MAPE => (x -> quantile(x, 0.50)) => :Test_MAPE_median,
        :Test_MAPE => (x -> quantile(x, 0.25)) => :Test_MAPE_lowerquartile,
        :Test_MAPE => minimum => :Test_MAPE_minimum ]
)