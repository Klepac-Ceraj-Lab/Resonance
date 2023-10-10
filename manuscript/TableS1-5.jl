#####
# Supplementary Table S1 thourgh S5. List of most important taxa for Random Forest models, MICROBIOME ONLY, ordered by relative weighted importance rank
#####
using Resonance
using CategoricalArrays
using Statistics
using Chain
using JLD2
using MLJ
using DataFrames.PrettyTables

isdir(tablefiles()) || mkpath(tablefiles())
RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree

## 1. Loading the model files

### 1.1. Concurrent composite scores
regression_currentCogScores_00to06mo_onlydemo = JLD2.load(modelfiles("regression_currentCogScores_00to06mo_onlydemo.jld"))["regression_currentCogScores_00to06mo_onlydemo"]
regression_currentCogScores_00to06mo_onlytaxa = JLD2.load(modelfiles("regression_currentCogScores_00to06mo_onlytaxa.jld"))["regression_currentCogScores_00to06mo_onlytaxa"]
regression_currentCogScores_00to06mo_demoplustaxa = JLD2.load(modelfiles("regression_currentCogScores_00to06mo_demoplustaxa.jld"))["regression_currentCogScores_00to06mo_demoplustaxa"]
regression_currentCogScores_00to06mo_onlyecs = JLD2.load(modelfiles("regression_currentCogScores_00to06mo_onlyecs.jld"))["regression_currentCogScores_00to06mo_onlyecs"]
regression_currentCogScores_00to06mo_demoplusecs = JLD2.load(modelfiles("regression_currentCogScores_00to06mo_demoplusecs.jld"))["regression_currentCogScores_00to06mo_demoplusecs"]
regression_currentCogScores_18to120mo_onlydemo = JLD2.load(modelfiles("regression_currentCogScores_18to120mo_onlydemo.jld"))["regression_currentCogScores_18to120mo_onlydemo"]
regression_currentCogScores_18to120mo_onlytaxa = JLD2.load(modelfiles("regression_currentCogScores_18to120mo_onlytaxa.jld"))["regression_currentCogScores_18to120mo_onlytaxa"]
regression_currentCogScores_18to120mo_demoplustaxa = JLD2.load(modelfiles("regression_currentCogScores_18to120mo_demoplustaxa.jld"))["regression_currentCogScores_18to120mo_demoplustaxa"]
regression_currentCogScores_18to120mo_onlyecs = JLD2.load(modelfiles("regression_currentCogScores_18to120mo_onlyecs.jld"))["regression_currentCogScores_18to120mo_onlyecs"]
regression_currentCogScores_18to120mo_demoplusecs = JLD2.load(modelfiles("regression_currentCogScores_18to120mo_demoplusecs.jld"))["regression_currentCogScores_18to120mo_demoplusecs"]
regression_currentCogScores_00to120mo_onlydemo = JLD2.load(modelfiles("regression_currentCogScores_00to120mo_onlydemo.jld"))["regression_currentCogScores_00to120mo_onlydemo"]
regression_currentCogScores_00to120mo_onlytaxa = JLD2.load(modelfiles("regression_currentCogScores_00to120mo_onlytaxa.jld"))["regression_currentCogScores_00to120mo_onlytaxa"]
regression_currentCogScores_00to120mo_demoplustaxa = JLD2.load(modelfiles("regression_currentCogScores_00to120mo_demoplustaxa.jld"))["regression_currentCogScores_00to120mo_demoplustaxa"]
regression_currentCogScores_00to120mo_onlyecs = JLD2.load(modelfiles("regression_currentCogScores_00to120mo_onlyecs.jld"))["regression_currentCogScores_00to120mo_onlyecs"]
regression_currentCogScores_00to120mo_demoplusecs = JLD2.load(modelfiles("regression_currentCogScores_00to120mo_demoplusecs.jld"))["regression_currentCogScores_00to120mo_demoplusecs"]

### 1.2. Concurrent Expressive Language scores
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

### 1.2. Concurrent Gross Motor scores
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

### 1.2. Concurrent Visual Reception scores
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

### 1.2. Future composite scores
regression_futureCogScores_onlydemo = JLD2.load(modelfiles("regression_futureCogScores_onlydemo.jld"))["regression_futureCogScores_onlydemo"]
regression_futureCogScores_onlytaxa = JLD2.load(modelfiles("regression_futureCogScores_onlytaxa.jld"))["regression_futureCogScores_onlytaxa"]
regression_futureCogScores_demoplustaxa = JLD2.load(modelfiles("regression_futureCogScores_demoplustaxa.jld"))["regression_futureCogScores_demoplustaxa"]
regression_futureCogScores_onlyecs = JLD2.load(modelfiles("regression_futureCogScores_onlyecs.jld"))["regression_futureCogScores_onlyecs"]
regression_futureCogScores_demoplusecs = JLD2.load(modelfiles("regression_futureCogScores_demoplusecs.jld"))["regression_futureCogScores_demoplusecs"]

## 2. Defining Suppl Table Functions

models_to_report = [
    regression_currentCogScores_18to120mo_onlytaxa,
    regression_currentExpressiveLanguages_18to120mo_onlytaxa,
    regression_currentGrossMotors_18to120mo_onlytaxa,
    regression_currentVisualReceptions_18to120mo_onlytaxa,
    regression_futureCogScores_onlytaxa
]

## 4. LaTeX outputs
spec_hl = LatexHighlighter((val, i, j) -> j == 2, ["textit"])
spec_ft = (val, i, j) -> j == 2 ? replace(val, "_"=> " ") : val

for m in models_to_report

    supptbl = singlemodel_importances_suppltable(m)
    prettyformat_suppltable!(supptbl)
    CSV.write(tablefiles("figure3_merits_importances", "importances_"*m.name*".csv"), supptbl)

    @show latextable = pretty_table(supptbl; highlighters=spec_hl, formatters=spec_ft,
        header = [ "Rank", "Variable", "Average fitness-weighted importance", "Relative fitness-weighted importance", "Rank-cumulative relative importance " ],
        backend = Val(:latex) )
    # write(tablefiles("figure3_merits_importances", "latex_importances_"*m.name*".txt"), latextable) # @Kevin I don't know how to save that.

end

# 4. LaTeX outputs
spec_hl = LatexHighlighter((val, i, j) -> j == 2, ["textit"])
spec_ft = (val, i, j) -> j == 2 ? replace(val, "_"=> " ") : val

for m in filter(f-> contains(f, r"TableS\d+\.csv"), readdir("manuscript/assets"; join=true))
    @info m
    supptbl = CSV.read(m, DataFrame)

    pretty_table(supptbl; highlighters=spec_hl, formatters=spec_ft,
        show_subheader=false,
        backend = Val(:latex) )
    # write(tablefiles("figure3_merits_importances", "latex_importances_"*m.name*".txt"), latextable) # @Kevin I don't know how to save that.

end