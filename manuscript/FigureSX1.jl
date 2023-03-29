#####
# Supplementary Figure: Pareto plot for taxa importances from FIgure 3
#####
using Resonance
using CategoricalArrays
using Statistics
using Chain
using JLD2
using MLJ
using CairoMakie

isdir(tablefiles()) || mkpath(tablefiles())

## 1. Loading the model files
RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree
regression_currentCogScores_00to06mo_onlytaxa = JLD2.load(modelfiles("manuscript", "regression_currentCogScores_00to06mo_onlytaxa.jld"))["regression_currentCogScores_00to06mo_onlytaxa"]
regression_currentCogScores_18to120mo_onlytaxa = JLD2.load(modelfiles("manuscript", "regression_currentCogScores_18to120mo_onlytaxa.jld"))["regression_currentCogScores_18to120mo_onlytaxa"]

## 2. Defining Suppl Table Functions
function singlemodel_importances_suppltable(rf_model)
    rf_model_importances = weighted_hpimportances(rf_model; normalize_importances=false, change_hashnames=false)
    rf_model_importances.relativeWeightedImportance = rf_model_importances.weightedImportance ./ sum(skipmissing(rf_model_importances.weightedImportance))
    rf_model_importances.cumulativeWeightedImportance = cumsum(rf_model_importances.relativeWeightedImportance)
    return rf_model_importances
end

## 3. Building the tables
supptblA = singlemodel_importances_suppltable(regression_currentCogScores_00to06mo_onlytaxa)
supptblB = singlemodel_importances_suppltable(regression_currentCogScores_18to120mo_onlytaxa)

## 4. Building the Figure
figure = Figure(resolution = (1200, 1200))
barcolor = :lightblue
curvecolor = :orange

ax1L = Axis(
    figure[1,1],
    xlabel = "Feature",
    xticks = (collect(1:nrow(supptblA)), replace.(supptblA.variable, "_"=>" ")),
    xticklabelsize=8,
    xticklabelrotation= pi/4,
    xticklabelfont="TeX Gyre Heros Makie Italic",
    title = "0 to 6 months",
    ylabel = "Individual relative importance")
ax1R = Axis(
    figure[1,1],
    xlabel = "Feature",
    xticks = (collect(1:nrow(supptblA))),
    ylabel = "Cumulative relative importance",
    yticks = 0:10:100,
    yaxisposition = :right)

ylims!(ax1L, [0.0, 10.0])
ylims!(ax1R, [0.0, 100.0])
hidexdecorations!(ax1R)
# hidespines!(ax1R)
linkxaxes!(ax1L, ax1R)

barplot!(
    ax1L,
    collect(1:nrow(supptblA)), (100 .* supptblA.relativeWeightedImportance),
    color = barcolor, strokecolor = :black, strokewidth = 1)
scatter!(
    ax1R,
    collect(1:nrow(supptblA)), (100 .* supptblA.cumulativeWeightedImportance),
    color = curvecolor, markersize = 15)
lines!(
    ax1R,
    collect(1:nrow(supptblA)), (100 .* supptblA.cumulativeWeightedImportance),
    color = curvecolor, linewidth = 5)

ax1L = Axis(
    figure[1,1],
    xlabel = "Feature",
    xticks = (collect(1:nrow(supptblA)), replace.(supptblA.variable, "_"=>" ")),
    xticklabelsize=8,
    xticklabelrotation= pi/4,
    xticklabelfont="TeX Gyre Heros Makie Italic",
    title = "0 to 6 months",
    ylabel = "Individual relative importance")
ax1R = Axis(
    figure[1,1],
    xlabel = "Feature",
    xticks = (collect(1:nrow(supptblA))),
    ylabel = "Cumulative relative importance",
    yticks = 0:10:100,
    yaxisposition = :right)

ylims!(ax1L, [0.0, 10.0])
ylims!(ax1R, [0.0, 100.0])
hidexdecorations!(ax1R)
# hidespines!(ax1R)
linkxaxes!(ax1L, ax1R)

barplot!(
    ax1L,
    collect(1:nrow(supptblA)), (100 .* supptblA.relativeWeightedImportance),
    color = barcolor, strokecolor = :black, strokewidth = 1)
scatter!(
    ax1R,
    collect(1:nrow(supptblA)), (100 .* supptblA.cumulativeWeightedImportance),
    color = curvecolor, markersize = 15)
lines!(
    ax1R,
    collect(1:nrow(supptblA)), (100 .* supptblA.cumulativeWeightedImportance),
    color = curvecolor, linewidth = 5)

ax2L = Axis(
    figure[2,1],
    xlabel = "Feature",
    xticks = (collect(1:nrow(supptblB)), replace.(supptblB.variable, "_"=>" ")),
    xticklabelsize=8,
    xticklabelrotation= pi/4,
    xticklabelfont="TeX Gyre Heros Makie Italic",
    title = "18 to 120 months",
    ylabel = "Individual relative importance")
ax2R = Axis(
    figure[2,1],
    xlabel = "Feature",
    xticks = (collect(1:nrow(supptblB))),
    ylabel = "Cumulative relative importance",
    yticks = 0:10:100,
    yaxisposition = :right)

ylims!(ax2L, [0.0, 10.0])
ylims!(ax2R, [0.0, 100.0])
hidexdecorations!(ax2R)
# hidespines!(ax1R)
linkxaxes!(ax2L, ax2R)

barplot!(
    ax2L,
    collect(1:nrow(supptblB)), (100 .* supptblB.relativeWeightedImportance),
    color = barcolor, strokecolor = :black, strokewidth = 1)
scatter!(
    ax2R,
    collect(1:nrow(supptblB)), (100 .* supptblB.cumulativeWeightedImportance),
    color = curvecolor, markersize = 15)
lines!(
    ax2R,
    collect(1:nrow(supptblB)), (100 .* supptblB.cumulativeWeightedImportance),
    color = curvecolor, linewidth = 5)

figure