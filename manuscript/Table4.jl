#####
# Table 3. Benchmark metrics for the cognitive assessment score prediction models. Confidence intervals are calculated from the distribution of metrics from repeated CV at a confidence level of 95%
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
RandomForestClassifier = MLJ.@load RandomForestClassifier pkg=DecisionTree
RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree
# concurrent cogScore reression from taxonomic profiles
regression_futureCogScores_onlydemo = JLD2.load(modelfiles("regression_futureCogScores_onlydemo.jld"))["regression_futureCogScores_onlydemo"]
regression_futureCogScores_onlytaxa = JLD2.load(modelfiles("regression_futureCogScores_onlytaxa.jld"))["regression_futureCogScores_onlytaxa"]
regression_futureCogScores_demoplustaxa = JLD2.load(modelfiles("regression_futureCogScores_demoplustaxa.jld"))["regression_futureCogScores_demoplustaxa"]
regression_futureCogScores_onlyecs = JLD2.load(modelfiles("regression_futureCogScores_onlyecs.jld"))["regression_futureCogScores_onlyecs"]
regression_futureCogScores_demoplusecs = JLD2.load(modelfiles("regression_futureCogScores_demoplusecs.jld"))["regression_futureCogScores_demoplusecs"]
## 2. Calculating the Figures of Merit

conf_interval(vv, critical_z=1.96) = critical_z * std(vv) / sqrt(length(vv))

prod_summary_table = vcat(
    DataFrames.combine(
        groupby(regression_futureCogScores_onlydemo.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "future_onlydemo") => :model
    ),
    DataFrames.combine(
        groupby(regression_futureCogScores_onlytaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "future_onlytaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_futureCogScores_demoplustaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "future_demoplustaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_futureCogScores_onlyecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "future_onlyecs") => :model
    ),
    DataFrames.combine(
        groupby(regression_futureCogScores_demoplusecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "future_demoplusecs") => :model
    )
)

## 3. Building the actual table

table4 = DataFrame(
    :subject_ages => repeat([ "future" ], 5),
    :microbial_feature => [ "-", "taxa", "taxa", "genes", "genes" ],
    :demo_provided => [ "+", "-", "+", "-", "+" ],
    :test_rmse => [ string(i)*" ± "*string(j) for (i,j) in zip(prod_summary_table.Test_RMSE_mean, prod_summary_table.Test_RMSE_CI) ],
    :test_cor => [ string(i)*" ± "*string(j) for (i,j) in zip(prod_summary_table.Test_Cor_mean, prod_summary_table.Test_Cor_CI) ]
)

CSV.write(tablefiles("Table4.csv"), table4)

pretty_table(table4;
    header = [ "Subject Ages (months)", "Microbial feature", "Demo.", "Test set Root-mean-square error (± C.I.)",  "Test set correlation (± C.I.)" ],
    backend = Val(:latex)
)