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

## 1. Loading the model files
RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree
# concurrent cogScore reression from taxonomic profiles
regression_currentCogScores_00to06mo_onlytaxa = JLD2.load(modelfiles("2023-02-15", "regression_currentCogScores_00to06mo_onlytaxa.jld"))["regression_currentCogScores_00to06mo_onlytaxa"]
regression_currentCogScores_18to120mo_onlytaxa = JLD2.load(modelfiles("2023-02-15", "regression_currentCogScores_18to120mo_onlytaxa.jld"))["regression_currentCogScores_18to120mo_onlytaxa"]

## 2. Relevant functions

function calculate_fitness(train_cor::Vector{Float64}, test_cor::Vector{Float64})
    positive_train_cor = map( x -> maximum([x, 0.0]), train_cor)
    positive_test_cor = map( x -> maximum([x, 0.0]), test_cor)
    return positive_train_cor .* positive_test_cor
end

function weighted_hpimportances(m, hp = 1; change_hashnames=false, hashnamestable::Union{Nothing, DataFrame} = nothing)
    merits = m.merits
    importances = m.importances
    fitnesses = calculate_fitness(merits.Train_Cor, merits.Test_Cor)
    fitness_subset_idx = findall(merits.Hyperpar_Idx .== hp)

    mean_importances = map(x -> sum(x .* fitnesses[fitness_subset_idx])/sum(fitnesses[fitness_subset_idx] .> 0.0), collect(eachrow(Matrix(importances[:, fitness_subset_idx .+ 1]))))
    mean_importances_df = sort( DataFrame(:variable => importances[:,1], :weightedImportance => mean_importances), :weightedImportance; rev=true)

    if change_hashnames
        newnames = leftjoin(mean_importances_df, hashnamestable, on = :variable => :hashname).longname
        mean_importances_df.variable = newnames
    end

    return mean_importances_df
end

function suppl_table(rf_model)
    rf_model_importances = weighted_hpimportances(rf_model; change_hashnames=false)
    rf_model_importances.relativeWeightedImportance = rf_model_importances.weightedImportance ./ sum(skipmissing(rf_model_importances.weightedImportance))
    rf_model_importances.cumulativeWeightedImportance = cumsum(rf_model_importances.relativeWeightedImportance)
    return rf_model_importances
end

function tblformat_suppl_table!(suptbl)
    suptbl = @chain suptbl begin
        insertcols!(1, :rank => 1:nrow(suptbl))
        transform!(:relativeWeightedImportance => (x -> string.(round.(100*x; digits = 2)).*" %"), renamecols = false)
        transform!(:cumulativeWeightedImportance => (x -> string.(round.(100*x; digits = 2)).*" %"), renamecols = false)
    end    
    return suptbl
end

## 3. Building the tables

supptblA = suppl_table(regression_currentCogScores_00to06mo_onlytaxa)
tblformat_suppl_table!(supptblA)
CSV.write(joinpath("manuscript", "assets", "TableS1A.csv"), supptblA)

supptblB = suppl_table(regression_currentCogScores_18to120mo_onlytaxa)
tblformat_suppl_table!(supptblB)
CSV.write(joinpath("manuscript", "assets", "TableS1B.csv"), supptblB)

## 4. LaTeX outputs

spec_hl = LatexHighlighter((val, i, j) -> j == 2, ["textit"])
spec_ft = (val, i, j) -> j == 2 ? replace(val, "_"=> " ") : val

pretty_table(supptblA; highlighters=spec_hl, formatters=spec_ft,
    header = [ "Rank", "Variable", "Average fitness-weighted importance", "Relative fitness-weighted importance", "Rank-cumulative relative importance " ],
    backend = Val(:latex)
)

pretty_table(supptblB; highlighters=spec_hl, formatters=spec_ft,
    header = [ "Rank", "Variable", "Average fitness-weighted importance", "Relative fitness-weighted importance", "Rank-cumulative relative importance " ],
    backend = Val(:latex)
)