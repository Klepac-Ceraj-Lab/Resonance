#####
# Table 5. Benchmark metrics for the cognitive assessment Mullen subscores
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
regression_currentExpressiveLanguages_18to120mo_onlydemo = JLD2.load(modelfiles("regression_currentExpressiveLanguages_18to120mo_onlydemo.jld"))["regression_currentExpressiveLanguages_18to120mo_onlydemo"]
regression_currentExpressiveLanguages_18to120mo_onlytaxa = JLD2.load(modelfiles("regression_currentExpressiveLanguages_18to120mo_onlytaxa.jld"))["regression_currentExpressiveLanguages_18to120mo_onlytaxa"]
regression_currentExpressiveLanguages_18to120mo_demoplustaxa = JLD2.load(modelfiles("regression_currentExpressiveLanguages_18to120mo_demoplustaxa.jld"))["regression_currentExpressiveLanguages_18to120mo_demoplustaxa"]
regression_currentExpressiveLanguages_18to120mo_onlyecs = JLD2.load(modelfiles("regression_currentExpressiveLanguages_18to120mo_onlyecs.jld"))["regression_currentExpressiveLanguages_18to120mo_onlyecs"]
regression_currentExpressiveLanguages_18to120mo_demoplusecs = JLD2.load(modelfiles("regression_currentExpressiveLanguages_18to120mo_demoplusecs.jld"))["regression_currentExpressiveLanguages_18to120mo_demoplusecs"]
regression_currentGrossMotors_18to120mo_onlydemo = JLD2.load(modelfiles("regression_currentGrossMotors_18to120mo_onlydemo.jld"))["regression_currentGrossMotors_18to120mo_onlydemo"]
regression_currentGrossMotors_18to120mo_onlytaxa = JLD2.load(modelfiles("regression_currentGrossMotors_18to120mo_onlytaxa.jld"))["regression_currentGrossMotors_18to120mo_onlytaxa"]
regression_currentGrossMotors_18to120mo_demoplustaxa = JLD2.load(modelfiles("regression_currentGrossMotors_18to120mo_demoplustaxa.jld"))["regression_currentGrossMotors_18to120mo_demoplustaxa"]
regression_currentGrossMotors_18to120mo_onlyecs = JLD2.load(modelfiles("regression_currentGrossMotors_18to120mo_onlyecs.jld"))["regression_currentGrossMotors_18to120mo_onlyecs"]
regression_currentGrossMotors_18to120mo_demoplusecs = JLD2.load(modelfiles("regression_currentGrossMotors_18to120mo_demoplusecs.jld"))["regression_currentGrossMotors_18to120mo_demoplusecs"]
regression_currentVisualReceptions_18to120mo_onlydemo = JLD2.load(modelfiles("regression_currentVisualReceptions_18to120mo_onlydemo.jld"))["regression_currentVisualReceptions_18to120mo_onlydemo"]
regression_currentVisualReceptions_18to120mo_onlytaxa = JLD2.load(modelfiles("regression_currentVisualReceptions_18to120mo_onlytaxa.jld"))["regression_currentVisualReceptions_18to120mo_onlytaxa"]
regression_currentVisualReceptions_18to120mo_demoplustaxa = JLD2.load(modelfiles("regression_currentVisualReceptions_18to120mo_demoplustaxa.jld"))["regression_currentVisualReceptions_18to120mo_demoplustaxa"]
regression_currentVisualReceptions_18to120mo_onlyecs = JLD2.load(modelfiles("regression_currentVisualReceptions_18to120mo_onlyecs.jld"))["regression_currentVisualReceptions_18to120mo_onlyecs"]
regression_currentVisualReceptions_18to120mo_demoplusecs = JLD2.load(modelfiles("regression_currentVisualReceptions_18to120mo_demoplusecs.jld"))["regression_currentVisualReceptions_18to120mo_demoplusecs"]

## 2. Calculating the Figures of Merit

conf_interval(vv, critical_z=1.96) = critical_z * std(vv) / sqrt(length(vv))

prod_summary_table = vcat(
    DataFrames.combine(
        groupby(regression_currentExpressiveLanguages_18to120mo_onlydemo.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 4)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "EL_onlydemo") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentExpressiveLanguages_18to120mo_onlytaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 4)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "EL_onlytaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentExpressiveLanguages_18to120mo_demoplustaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 4)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "EL_demoplustaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentExpressiveLanguages_18to120mo_onlyecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 4)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "EL_onlyecs") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentExpressiveLanguages_18to120mo_demoplusecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 4)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "EL_demoplusecs") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentGrossMotors_18to120mo_onlydemo.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 4)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 5)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "GM_onlydemo") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentGrossMotors_18to120mo_onlytaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 4)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 5)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "GM_onlytaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentGrossMotors_18to120mo_demoplustaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 4)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 5)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "GM_demoplustaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentGrossMotors_18to120mo_onlyecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 4)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 5)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "GM_onlyecs") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentGrossMotors_18to120mo_demoplusecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 4)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 5)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "GM_demoplusecs") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentVisualReceptions_18to120mo_onlydemo.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 4)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 5)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "VR_onlydemo") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentVisualReceptions_18to120mo_onlytaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 4)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 5)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "VR_onlytaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentVisualReceptions_18to120mo_demoplustaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 4)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 5)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "VR_demoplustaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentVisualReceptions_18to120mo_onlyecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 4)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 5)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "VR_onlyecs") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentVisualReceptions_18to120mo_demoplusecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 4)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 5)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "VR_demoplusecs") => :model
    )
)

## 3. Building the actual table

table5 = DataFrame(
    :target_subscale => vcat( repeat([ "Expressive Language" ], 5), repeat( [ "Gross Motor" ], 5), repeat( [ "Visual Reception" ], 5)),
    :microbial_feature => repeat([ "-", "taxa", "taxa", "genes", "genes" ], 3),
    :demo_provided => repeat([ "+", "-", "+", "-", "+" ], 3),
    :test_rmse => [ string(i)*" ± "*string(j) for (i,j) in zip(prod_summary_table.Test_RMSE_mean, prod_summary_table.Test_RMSE_CI) ],
    :test_cor => [ string(i)*" ± "*string(j) for (i,j) in zip(prod_summary_table.Test_Cor_mean, prod_summary_table.Test_Cor_CI) ]
)

CSV.write(tablefiles("Table5.csv"), table5)

pretty_table(table5;
    header = [ "Target Mullen scale", "Microbial feature", "Demo.", "Test set Root-mean-square error (± C.I.)", "Test set correlation (± C.I.)" ],
    backend = Val(:latex)
)