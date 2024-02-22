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
regression_currentExpressiveLanguages_00to06mo_onlydemo = JLD2.load(modelfiles("regression_currentExpressiveLanguages_00to06mo_onlydemo.jld"))["regression_currentExpressiveLanguages_00to06mo_onlydemo"]
regression_currentExpressiveLanguages_00to06mo_onlytaxa = JLD2.load(modelfiles("regression_currentExpressiveLanguages_00to06mo_onlytaxa.jld"))["regression_currentExpressiveLanguages_00to06mo_onlytaxa"]
regression_currentExpressiveLanguages_00to06mo_demoplustaxa = JLD2.load(modelfiles("regression_currentExpressiveLanguages_00to06mo_demoplustaxa.jld"))["regression_currentExpressiveLanguages_00to06mo_demoplustaxa"]
regression_currentExpressiveLanguages_00to06mo_onlyecs = JLD2.load(modelfiles("regression_currentExpressiveLanguages_00to06mo_onlyecs.jld"))["regression_currentExpressiveLanguages_00to06mo_onlyecs"]
regression_currentExpressiveLanguages_00to06mo_demoplusecs = JLD2.load(modelfiles("regression_currentExpressiveLanguages_00to06mo_demoplusecs.jld"))["regression_currentExpressiveLanguages_00to06mo_demoplusecs"]
regression_currentExpressiveLanguages_18to120mo_onlydemo = JLD2.load(modelfiles("regression_currentExpressiveLanguages_18to120mo_onlydemo.jld"))["regression_currentExpressiveLanguages_18to120mo_onlydemo"]
regression_currentExpressiveLanguages_18to120mo_onlytaxa = JLD2.load(modelfiles("regression_currentExpressiveLanguages_18to120mo_onlytaxa.jld"))["regression_currentExpressiveLanguages_18to120mo_onlytaxa"]
regression_currentExpressiveLanguages_18to120mo_demoplustaxa = JLD2.load(modelfiles("regression_currentExpressiveLanguages_18to120mo_demoplustaxa.jld"))["regression_currentExpressiveLanguages_18to120mo_demoplustaxa"]
regression_currentExpressiveLanguages_18to120mo_onlyecs = JLD2.load(modelfiles("regression_currentExpressiveLanguages_18to120mo_onlyecs.jld"))["regression_currentExpressiveLanguages_18to120mo_onlyecs"]
regression_currentExpressiveLanguages_18to120mo_demoplusecs = JLD2.load(modelfiles("regression_currentExpressiveLanguages_18to120mo_demoplusecs.jld"))["regression_currentExpressiveLanguages_18to120mo_demoplusecs"]
regression_currentExpressiveLanguages_00to120mo_onlydemo = JLD2.load(modelfiles("regression_currentExpressiveLanguages_00to120mo_onlydemo.jld"))["regression_currentExpressiveLanguages_00to120mo_onlydemo"]
regression_currentExpressiveLanguages_00to120mo_onlytaxa = JLD2.load(modelfiles("regression_currentExpressiveLanguages_00to120mo_onlytaxa.jld"))["regression_currentExpressiveLanguages_00to120mo_onlytaxa"]
regression_currentExpressiveLanguages_00to120mo_demoplustaxa = JLD2.load(modelfiles("regression_currentExpressiveLanguages_00to120mo_demoplustaxa.jld"))["regression_currentExpressiveLanguages_00to120mo_demoplustaxa"]
regression_currentExpressiveLanguages_00to120mo_onlyecs = JLD2.load(modelfiles("regression_currentExpressiveLanguages_00to120mo_onlyecs.jld"))["regression_currentExpressiveLanguages_00to120mo_onlyecs"]
regression_currentExpressiveLanguages_00to120mo_demoplusecs = JLD2.load(modelfiles("regression_currentExpressiveLanguages_00to120mo_demoplusecs.jld"))["regression_currentExpressiveLanguages_00to120mo_demoplusecs"]
regression_currentGrossMotors_00to06mo_onlydemo = JLD2.load(modelfiles("regression_currentGrossMotors_00to06mo_onlydemo.jld"))["regression_currentGrossMotors_00to06mo_onlydemo"]
regression_currentGrossMotors_00to06mo_onlytaxa = JLD2.load(modelfiles("regression_currentGrossMotors_00to06mo_onlytaxa.jld"))["regression_currentGrossMotors_00to06mo_onlytaxa"]
regression_currentGrossMotors_00to06mo_demoplustaxa = JLD2.load(modelfiles("regression_currentGrossMotors_00to06mo_demoplustaxa.jld"))["regression_currentGrossMotors_00to06mo_demoplustaxa"]
regression_currentGrossMotors_00to06mo_onlyecs = JLD2.load(modelfiles("regression_currentGrossMotors_00to06mo_onlyecs.jld"))["regression_currentGrossMotors_00to06mo_onlyecs"]
regression_currentGrossMotors_00to06mo_demoplusecs = JLD2.load(modelfiles("regression_currentGrossMotors_00to06mo_demoplusecs.jld"))["regression_currentGrossMotors_00to06mo_demoplusecs"]
regression_currentGrossMotors_18to120mo_onlydemo = JLD2.load(modelfiles("regression_currentGrossMotors_18to120mo_onlydemo.jld"))["regression_currentGrossMotors_18to120mo_onlydemo"]
regression_currentGrossMotors_18to120mo_onlytaxa = JLD2.load(modelfiles("regression_currentGrossMotors_18to120mo_onlytaxa.jld"))["regression_currentGrossMotors_18to120mo_onlytaxa"]
regression_currentGrossMotors_18to120mo_demoplustaxa = JLD2.load(modelfiles("regression_currentGrossMotors_18to120mo_demoplustaxa.jld"))["regression_currentGrossMotors_18to120mo_demoplustaxa"]
regression_currentGrossMotors_18to120mo_onlyecs = JLD2.load(modelfiles("regression_currentGrossMotors_18to120mo_onlyecs.jld"))["regression_currentGrossMotors_18to120mo_onlyecs"]
regression_currentGrossMotors_18to120mo_demoplusecs = JLD2.load(modelfiles("regression_currentGrossMotors_18to120mo_demoplusecs.jld"))["regression_currentGrossMotors_18to120mo_demoplusecs"]
regression_currentGrossMotors_00to120mo_onlydemo = JLD2.load(modelfiles("regression_currentGrossMotors_00to120mo_onlydemo.jld"))["regression_currentGrossMotors_00to120mo_onlydemo"]
regression_currentGrossMotors_00to120mo_onlytaxa = JLD2.load(modelfiles("regression_currentGrossMotors_00to120mo_onlytaxa.jld"))["regression_currentGrossMotors_00to120mo_onlytaxa"]
regression_currentGrossMotors_00to120mo_demoplustaxa = JLD2.load(modelfiles("regression_currentGrossMotors_00to120mo_demoplustaxa.jld"))["regression_currentGrossMotors_00to120mo_demoplustaxa"]
regression_currentGrossMotors_00to120mo_onlyecs = JLD2.load(modelfiles("regression_currentGrossMotors_00to120mo_onlyecs.jld"))["regression_currentGrossMotors_00to120mo_onlyecs"]
regression_currentGrossMotors_00to120mo_demoplusecs = JLD2.load(modelfiles("regression_currentGrossMotors_00to120mo_demoplusecs.jld"))["regression_currentGrossMotors_00to120mo_demoplusecs"]
regression_currentVisualReceptions_00to06mo_onlydemo = JLD2.load(modelfiles("regression_currentVisualReceptions_00to06mo_onlydemo.jld"))["regression_currentVisualReceptions_00to06mo_onlydemo"]
regression_currentVisualReceptions_00to06mo_onlytaxa = JLD2.load(modelfiles("regression_currentVisualReceptions_00to06mo_onlytaxa.jld"))["regression_currentVisualReceptions_00to06mo_onlytaxa"]
regression_currentVisualReceptions_00to06mo_demoplustaxa = JLD2.load(modelfiles("regression_currentVisualReceptions_00to06mo_demoplustaxa.jld"))["regression_currentVisualReceptions_00to06mo_demoplustaxa"]
regression_currentVisualReceptions_00to06mo_onlyecs = JLD2.load(modelfiles("regression_currentVisualReceptions_00to06mo_onlyecs.jld"))["regression_currentVisualReceptions_00to06mo_onlyecs"]
regression_currentVisualReceptions_00to06mo_demoplusecs = JLD2.load(modelfiles("regression_currentVisualReceptions_00to06mo_demoplusecs.jld"))["regression_currentVisualReceptions_00to06mo_demoplusecs"]
regression_currentVisualReceptions_18to120mo_onlydemo = JLD2.load(modelfiles("regression_currentVisualReceptions_18to120mo_onlydemo.jld"))["regression_currentVisualReceptions_18to120mo_onlydemo"]
regression_currentVisualReceptions_18to120mo_onlytaxa = JLD2.load(modelfiles("regression_currentVisualReceptions_18to120mo_onlytaxa.jld"))["regression_currentVisualReceptions_18to120mo_onlytaxa"]
regression_currentVisualReceptions_18to120mo_demoplustaxa = JLD2.load(modelfiles("regression_currentVisualReceptions_18to120mo_demoplustaxa.jld"))["regression_currentVisualReceptions_18to120mo_demoplustaxa"]
regression_currentVisualReceptions_18to120mo_onlyecs = JLD2.load(modelfiles("regression_currentVisualReceptions_18to120mo_onlyecs.jld"))["regression_currentVisualReceptions_18to120mo_onlyecs"]
regression_currentVisualReceptions_18to120mo_demoplusecs = JLD2.load(modelfiles("regression_currentVisualReceptions_18to120mo_demoplusecs.jld"))["regression_currentVisualReceptions_18to120mo_demoplusecs"]
regression_currentVisualReceptions_00to120mo_onlydemo = JLD2.load(modelfiles("regression_currentVisualReceptions_00to120mo_onlydemo.jld"))["regression_currentVisualReceptions_00to120mo_onlydemo"]
regression_currentVisualReceptions_00to120mo_onlytaxa = JLD2.load(modelfiles("regression_currentVisualReceptions_00to120mo_onlytaxa.jld"))["regression_currentVisualReceptions_00to120mo_onlytaxa"]
regression_currentVisualReceptions_00to120mo_demoplustaxa = JLD2.load(modelfiles("regression_currentVisualReceptions_00to120mo_demoplustaxa.jld"))["regression_currentVisualReceptions_00to120mo_demoplustaxa"]
regression_currentVisualReceptions_00to120mo_onlyecs = JLD2.load(modelfiles("regression_currentVisualReceptions_00to120mo_onlyecs.jld"))["regression_currentVisualReceptions_00to120mo_onlyecs"]
regression_currentVisualReceptions_00to120mo_demoplusecs = JLD2.load(modelfiles("regression_currentVisualReceptions_00to120mo_demoplusecs.jld"))["regression_currentVisualReceptions_00to120mo_demoplusecs"]

## 2. Calculating the Figures of Merit

conf_interval(vv, critical_z=1.96) = critical_z * std(vv) / sqrt(length(vv))

prod_summary_table = vcat(
    DataFrames.combine(
        groupby(regression_currentExpressiveLanguages_00to06mo_onlydemo.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to06_onlydemo") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentExpressiveLanguages_00to06mo_onlytaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to06_onlytaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentExpressiveLanguages_00to06mo_demoplustaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to06_demoplustaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentExpressiveLanguages_00to06mo_onlyecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to06_onlyecs") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentExpressiveLanguages_00to06mo_demoplusecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to06_demoplusecs") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentExpressiveLanguages_18to120mo_onlydemo.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "18to120_onlydemo") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentExpressiveLanguages_18to120mo_onlytaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "18to120_onlytaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentExpressiveLanguages_18to120mo_demoplustaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "18to120_demoplustaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentExpressiveLanguages_18to120mo_onlyecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "18to120_onlyecs") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentExpressiveLanguages_18to120mo_demoplusecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "18to120_demoplusecs") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentExpressiveLanguages_00to120mo_onlydemo.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to120_onlydemo") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentExpressiveLanguages_00to120mo_onlytaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to120_onlytaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentExpressiveLanguages_00to120mo_demoplustaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to120_demoplustaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentExpressiveLanguages_00to120mo_onlyecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to120_onlyecs") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentExpressiveLanguages_00to120mo_demoplusecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to120_demoplusecs") => :model
    ), ### END OF ExpressiveLanguage
    DataFrames.combine(
        groupby(regression_currentGrossMotors_00to06mo_onlydemo.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to06_onlydemo") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentGrossMotors_00to06mo_onlytaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to06_onlytaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentGrossMotors_00to06mo_demoplustaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to06_demoplustaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentGrossMotors_00to06mo_onlyecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to06_onlyecs") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentGrossMotors_00to06mo_demoplusecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to06_demoplusecs") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentGrossMotors_18to120mo_onlydemo.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "18to120_onlydemo") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentGrossMotors_18to120mo_onlytaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "18to120_onlytaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentGrossMotors_18to120mo_demoplustaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "18to120_demoplustaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentGrossMotors_18to120mo_onlyecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "18to120_onlyecs") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentGrossMotors_18to120mo_demoplusecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "18to120_demoplusecs") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentGrossMotors_00to120mo_onlydemo.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to120_onlydemo") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentGrossMotors_00to120mo_onlytaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to120_onlytaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentGrossMotors_00to120mo_demoplustaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to120_demoplustaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentGrossMotors_00to120mo_onlyecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to120_onlyecs") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentGrossMotors_00to120mo_demoplusecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to120_demoplusecs") => :model
    ), ### END OF GrossMotor
    DataFrames.combine(
        groupby(regression_currentVisualReceptions_00to06mo_onlydemo.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to06_onlydemo") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentVisualReceptions_00to06mo_onlytaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to06_onlytaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentVisualReceptions_00to06mo_demoplustaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to06_demoplustaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentVisualReceptions_00to06mo_onlyecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to06_onlyecs") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentVisualReceptions_00to06mo_demoplusecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 2)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 2)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to06_demoplusecs") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentVisualReceptions_18to120mo_onlydemo.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "18to120_onlydemo") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentVisualReceptions_18to120mo_onlytaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "18to120_onlytaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentVisualReceptions_18to120mo_demoplustaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "18to120_demoplustaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentVisualReceptions_18to120mo_onlyecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "18to120_onlyecs") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentVisualReceptions_18to120mo_demoplusecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "18to120_demoplusecs") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentVisualReceptions_00to120mo_onlydemo.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to120_onlydemo") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentVisualReceptions_00to120mo_onlytaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to120_onlytaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentVisualReceptions_00to120mo_demoplustaxa.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to120_demoplustaxa") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentVisualReceptions_00to120mo_onlyecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to120_onlyecs") => :model
    ),
    DataFrames.combine(
        groupby(regression_currentVisualReceptions_00to120mo_demoplusecs.merits, :Hyperpar_Idx),
        :Test_RMSE => (x -> round(mean(x); digits = 2)) => :Test_RMSE_mean,
        :Test_RMSE => (x -> round(conf_interval(x); digits = 2)) => :Test_RMSE_CI,
        :Test_Cor => (x -> round(mean(x); digits = 3)) => :Test_Cor_mean,
        :Test_Cor => (x -> round(conf_interval(x); digits = 3)) => :Test_Cor_CI,
        :Hyperpar_Idx=> (x-> "00to120_demoplusecs") => :model
    ) ### END OF VisualReception
)

## 3. Building the actual table

table3 = DataFrame(
    :subscore => vcat( repeat([ "Expressive Language" ], 15), repeat( [ "Gross Motor" ], 15), repeat( [ "Visual Reception" ], 15) ),
    :subject_ages => repeat( vcat( repeat([ "0 to 6" ], 5), repeat( [ "18 to 120" ], 5), repeat( [ "00 to 120" ], 5)), 3),
    :microbial_feature => repeat([ "-", "taxa", "taxa", "genes", "genes" ], 9),
    :demo_provided => repeat([ "+", "-", "+", "-", "+" ], 9),
    :test_rmse => [ string(i)*" ± "*string(j) for (i,j) in zip(prod_summary_table.Test_RMSE_mean, prod_summary_table.Test_RMSE_CI) ],
    :test_cor => [ string(i)*" ± "*string(j) for (i,j) in zip(prod_summary_table.Test_Cor_mean, prod_summary_table.Test_Cor_CI) ]
)

CSV.write(tablefiles("Table3_subscores.csv"), table3)

pretty_table(table3;
    header = [  "Subscale", "Subject Ages (months)", "Microbial feature", "Demo.", "Test set correlation (± C.I.)", "Test set Root-mean-square error (± C.I.)" ],
    backend = Val(:latex)
)