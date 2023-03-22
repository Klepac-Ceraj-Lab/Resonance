#####
# Table 4. Summary statistics for the neuroanatomy prediction benchmarks
#####

using Resonance
using CategoricalArrays
using Statistics
using Chain
using JLD2
using MLJ
using DataFrames.PrettyTables

isdir(tablefiles()) || mkpath(tablefiles())

## 1. Loading the model files

RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree
# concurrent brain volumes regression from taxonomic profiles
JLD2.@load "models/2023-02-15/brain_models.jld"

ordered_brain_segments_list = [
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

## 2. Calculating the Figures of Merit
mean_brain_merits = reduce(
    vcat,
    [ DataFrame(:variable => ordered_brain_segments_list[i], :Train_MAPE => mean(brain_models[ordered_brain_segments_list[i]].merits.Train_MAPE), :Test_MAPE => mean(brain_models[ordered_brain_segments_list[i]].merits.Test_MAPE), :Train_Cor => mean(brain_models[ordered_brain_segments_list[i]].merits.Train_Cor), :Test_Cor => mean(brain_models[ordered_brain_segments_list[i]].merits.Test_Cor)) for i in eachindex(ordered_brain_segments_list) ]
)

## 3. Building the actual table

table4 = DataFrame(
    :statisic => [ "Mean", "Standard deviation", "Maximum", "75th percentile", "50th percentile", "25th percentile", "Minimum" ],
    :test_MAPE => [ round(f(mean_brain_merits.Test_MAPE);digits=3) for f in [ mean, std, maximum, x -> quantile(x, 0.75), x -> quantile(x, 0.50), x -> quantile(x, 0.25), minimum ] ],
    :test_Cor => [ round(f(mean_brain_merits.Test_Cor);digits=3) for f in [ mean, std, maximum, x -> quantile(x, 0.75), x -> quantile(x, 0.50), x -> quantile(x, 0.25), minimum ] ]
)

CSV.write(tablefiles("Table4.csv"), table4)

pretty_table(table4;
    header = [ "Statistic", "Mean absolute proportional error (MAPE)", "Correlation coefficient (R)" ],
    backend = Val(:latex)
)