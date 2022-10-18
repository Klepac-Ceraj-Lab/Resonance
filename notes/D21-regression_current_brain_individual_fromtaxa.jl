#####
# Notebook D21 - Independent regression of current Brain Segment from taxonomic data
#####

using Resonance
using Chain
using CSV
using Statistics
using Random
using DataFrames
using StableRNGs
using ProgressMeter
using MLJ
using DecisionTree
using GLM
using JLD2
ml_rng = StableRNG(0)

#####
# Loading Data
#####

mdata = @chain Resonance.load(Metadata()) begin
    #subset(:has_segmentation)
    select( [ :subject, :timepoint, :ageMonths, :omni ] )
    dropmissing(:omni)
end

taxa = Resonance.load(Resonance.TaxonomicProfiles(); timepoint_metadata = mdata)
taxa_df = DataFrame([eachcol(collect(taxa.abundances'))...], map( x -> string(x), taxa.features), copycols=true)
taxa_df = taxa_df[:, [ el.rank for el in taxa.features ] .== :species ]
taxa_df.s__Collinsella_massiliensis = myxor.(taxa_df.s__Collinsella_massiliensis, taxa_df."s__[Collinsella]_massiliensis")
select!(taxa_df, Not("s__[Collinsella]_massiliensis"))
insertcols!(taxa_df, 1, :sample => collect(keys(taxa.sidx))[collect(values(taxa.sidx))])

brain = Resonance.load(Resonance.Neuroimaging())
brain_df = DataFrame([eachcol(collect(brain.abundances'))...], map( x -> string(x), brain.features), copycols=true)
insertcols!(brain_df, 1, :sample => collect(keys(brain.sidx))[collect(values(brain.sidx))])

subject_taxa_df = @chain mdata begin
    leftjoin(taxa_df, on = :omni => :sample, matchmissing=:notequal)
    leftjoin(brain_df, on = :omni => :sample, matchmissing=:notequal)
    # rename!(:omni => :sample)
    select!( Not(:omni) )
    sort!( [ :subject, :timepoint ])
    dropmissing()
end ## Columns 1:3 are metadata, columns 4:552 are taxonomic profile, columns 553:651 are brain data

subject_taxa_df[:, 553:650] .= 100.0 .* subject_taxa_df[:, 553:650]

#####
# Joining certain columns
#####

joined_segments_df = @chain subject_taxa_df begin
    # 1. Removing ventricles and other vessels (too small)
    select(Not(["right-inferior-lateral-ventricle", "left-inferior-lateral-ventricle"]))
    select(Not(["3rd-ventricle", "4th-ventricle"]))
    select(Not("left-vessel")) # Curiously, we do not have a right vessel on the dataframe
end

#####
# Training Models
#####

original_prediction_segments = [
    # "left-lateral-ventricle", "right-lateral-ventricle",
    # "left-inferior-lateral-ventricle", "right-inferior-lateral-ventricle",
    "left-cerebellum-exterior", "right-cerebellum-exterior",
    "left-cerebellum-white-matter", "right-cerebellum-white-matter",
    "left-thalamus-proper", "right-thalamus-proper",
    "left-caudate", "right-caudate",
    "left-putamen", "right-putamen",
    "left-pallidum", "right-pallidum",
    "left-hippocampus", "right-hippocampus",
    "left-amygdala", "right-amygdala",
    "left-accumbens-area", "right-accumbens-area",
    "left-ventral-DC", "right-ventral-DC",
    # "left-vessel",
    "left-basal-forebrain", "right-basal-forebrain",
    "cerebellar-vermal-lobules-I-V", "cerebellar-vermal-lobules-VI-VII", "cerebellar-vermal-lobules-VIII-X",
    "left-caudal-anterior-cingulate", "right-caudal-anterior-cingulate",
    "left-caudal-middle-frontal", "right-caudal-middle-frontal",
    "left-cuneus", "right-cuneus", 
    "left-entorhinal", "right-entorhinal",
    "left-fusiform", "right-fusiform",
    "left-inferior-parietal", "right-inferior-parietal",
    "left-inferior-temporal", "right-inferior-temporal",
    "left-isthmus-cingulate", "right-isthmus-cingulate",
    "left-lateral-occipital", "right-lateral-occipital",
    "left-lateral-orbitofrontal", "right-lateral-orbitofrontal",
    "left-lingual", "right-lingual", 
    "left-medial-orbitofrontal", "right-medial-orbitofrontal",
    "left-middle-temporal", "right-middle-temporal",
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
    "left-rostral-anterior-cingulate", "right-rostral-anterior-cingulate",
    "left-rostral-middle-frontal", "right-rostral-middle-frontal",
    "left-superior-frontal", "right-superior-frontal",
    "left-superior-parietal", "right-superior-parietal",
    "left-superior-temporal", "right-superior-temporal",
    "left-supramarginal", "right-supramarginal",
    "left-transverse-temporal", "right-transverse-temporal",
    "left-insula", "right-insula",
    # "3rd-ventricle", "4th-ventricle",
    "Brain-stem"
]

tuning_space = (
    maxnodes_range = collect(1:3:15),
    nodesize_range = collect(3:3:15),
    sampsize_range = [0.6, 0.7, 0.8],
    mtry_range = collect(10:15:100),
    ntrees_range = [100, 300, 500]
    )

RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree

brain_results = Vector{UnivariateRandomForestRegressor}()

for (coll) in original_prediction_segments
    println("Column $(coll)")

    regression_results = train_randomforest(
        Regression(),
        coll,
        subject_taxa_df,
        identity,
        4:552,
        n_splits = 10,
        coll;
        tuning_space = tuning_space,
        split_proportion=0.75,
        train_rng=ml_rng
        )

    push!(brain_results, regression_results)
end

JLD2.@save "models/brain_univariate_regression.jld" brain_results