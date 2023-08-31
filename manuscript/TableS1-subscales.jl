#####
# Supplementary Table S1A-B. List of most important taxa for Random Forest models on the (A) 00-06 and (B) 18-120 months age bracket, MICROBIOME ONLY, ordered by relative weighted importance rank
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
## 2. Defining Suppl Table Functions
function singlemodel_importances_suppltable(rf_model)
    rf_model_importances = weighted_hpimportances(rf_model; normalize_importances=false, change_hashnames=false)
    rf_model_importances.relativeWeightedImportance = rf_model_importances.weightedImportance ./ sum(skipmissing(rf_model_importances.weightedImportance))
    rf_model_importances.cumulativeWeightedImportance = cumsum(rf_model_importances.relativeWeightedImportance)
    return rf_model_importances
end

function prettyformat_suppltable!(suptbl)
    suptbl = @chain suptbl begin
        insertcols!(1, :rank => 1:nrow(suptbl))
        transform!(:weightedImportance => (x -> string.(round.(x; digits = 5))) => :weightedImportance, renamecols = false)
        transform!(:relativeWeightedImportance => (x -> string.(round.(100*x; digits = 2)).*" %"), renamecols = false)
        transform!(:cumulativeWeightedImportance => (x -> string.(round.(100*x; digits = 2)).*" %"), renamecols = false)
    end    
    return suptbl
end

models_to_report = [
    regression_currentExpressiveLanguages_00to06mo_onlydemo,
    regression_currentExpressiveLanguages_00to06mo_onlytaxa,
    regression_currentExpressiveLanguages_00to06mo_demoplustaxa,
    regression_currentExpressiveLanguages_00to06mo_onlyecs,
    regression_currentExpressiveLanguages_00to06mo_demoplusecs,
    regression_currentExpressiveLanguages_18to120mo_onlydemo,
    regression_currentExpressiveLanguages_18to120mo_onlytaxa,
    regression_currentExpressiveLanguages_18to120mo_demoplustaxa,
    regression_currentExpressiveLanguages_18to120mo_onlyecs,
    regression_currentExpressiveLanguages_18to120mo_demoplusecs,
    regression_currentExpressiveLanguages_00to120mo_onlydemo,
    regression_currentExpressiveLanguages_00to120mo_onlytaxa,
    regression_currentExpressiveLanguages_00to120mo_demoplustaxa,
    regression_currentExpressiveLanguages_00to120mo_onlyecs,
    regression_currentExpressiveLanguages_00to120mo_demoplusecs,
    regression_currentGrossMotors_00to06mo_onlydemo,
    regression_currentGrossMotors_00to06mo_onlytaxa,
    regression_currentGrossMotors_00to06mo_demoplustaxa,
    regression_currentGrossMotors_00to06mo_onlyecs,
    regression_currentGrossMotors_00to06mo_demoplusecs,
    regression_currentGrossMotors_18to120mo_onlydemo,
    regression_currentGrossMotors_18to120mo_onlytaxa,
    regression_currentGrossMotors_18to120mo_demoplustaxa,
    regression_currentGrossMotors_18to120mo_onlyecs,
    regression_currentGrossMotors_18to120mo_demoplusecs,
    regression_currentGrossMotors_00to120mo_onlydemo,
    regression_currentGrossMotors_00to120mo_onlytaxa,
    regression_currentGrossMotors_00to120mo_demoplustaxa,
    regression_currentGrossMotors_00to120mo_onlyecs,
    regression_currentGrossMotors_00to120mo_demoplusecs,    
    regression_currentVisualReceptions_00to06mo_onlydemo,
    regression_currentVisualReceptions_00to06mo_onlytaxa,
    regression_currentVisualReceptions_00to06mo_demoplustaxa,
    regression_currentVisualReceptions_00to06mo_onlyecs,
    regression_currentVisualReceptions_00to06mo_demoplusecs,
    regression_currentVisualReceptions_18to120mo_onlydemo,
    regression_currentVisualReceptions_18to120mo_onlytaxa,
    regression_currentVisualReceptions_18to120mo_demoplustaxa,
    regression_currentVisualReceptions_18to120mo_onlyecs,
    regression_currentVisualReceptions_18to120mo_demoplusecs,
    regression_currentVisualReceptions_00to120mo_onlydemo,
    regression_currentVisualReceptions_00to120mo_onlytaxa,
    regression_currentVisualReceptions_00to120mo_demoplustaxa,
    regression_currentVisualReceptions_00to120mo_onlyecs,
    regression_currentVisualReceptions_00to120mo_demoplusecs,    
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