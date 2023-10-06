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
brain_models = JLD2.load(modelfiles("brain_models.jld"))["brain_models"]
include("manuscript/Figure4-definitions.jl")

## 2. Calculating the Figures of Merit
mean_brain_merits = reduce(
    vcat,
    [ DataFrame(:variable => ordered_brain_segments_list[i], :Train_MAPE => mean(brain_models[ordered_brain_segments_list[i]].merits.Train_MAPE), :Test_MAPE => mean(brain_models[ordered_brain_segments_list[i]].merits.Test_MAPE), :Train_Cor => mean(brain_models[ordered_brain_segments_list[i]].merits.Train_Cor), :Test_Cor => mean(brain_models[ordered_brain_segments_list[i]].merits.Test_Cor)) for i in eachindex(ordered_brain_segments_list) ]
)
mean_brain_merits = subset(mean_brain_merits, :variable => x -> x .∈ Ref(interesting_segments) )

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