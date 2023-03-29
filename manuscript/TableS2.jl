#####
# Supplementary Table 2. Summary statistics for the neuroanatomy prediction benchmarks
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
brain_models = JLD2.load(modelfiles("manuscript", "brain_models.jld"))["brain_models"]
include("manuscript/Figure4-definitions.jl")

## 2. Calculating the Figures of Merit
tableS2 = mapreduce( vcat, eachindex(ordered_brain_segments_list)) do i
    DataFrame(:variable => ordered_brain_segments_list[i],
              :Test_MAPE => round(mean(brain_models[ordered_brain_segments_list[i]].merits.Test_MAPE); digits = 4), 
              :Test_Cor => round(mean(brain_models[ordered_brain_segments_list[i]].merits.Test_Cor); digits = 4))
end
tableS2.variable = uppercasefirst.(map(x -> replace(x, "-"=>" "), tableS2.variable))

CSV.write(joinpath("manuscript", "assets", "TableS2.csv"), tableS2)

pretty_table(tableS2;
    header = [ "Segment", "Mean absolute proportional error (MAPE)", "Correlation coefficient (R)" ],
    backend = Val(:latex)
)