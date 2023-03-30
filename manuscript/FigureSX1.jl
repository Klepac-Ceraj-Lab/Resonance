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

## 2. Building the tables
supptblA = singlemodel_importances_suppltable(regression_currentCogScores_00to06mo_onlytaxa)
supptblB = singlemodel_importances_suppltable(regression_currentCogScores_18to120mo_onlytaxa)

## 3. Building the Figure
figure = Figure(resolution = (1200, 1200))

this_barcolor = :lightblue
this_curvecolor = :orange

plot_importances_pareto!(figure[1,1], supptblA, "Pareto Plot - 0 to 6 months"; barcolor = this_barcolor, curvecolor = this_curvecolor)
plot_importances_pareto!(figure[2,1], supptblB[1:70, :], "Pareto Plot - 18 to 120 months"; barcolor = this_barcolor, curvecolor = this_curvecolor)

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

tightlimits(ax1L, Left(), Right())
tightlimits(ax1R, Left(), Right())

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

tightlimits(ax1L, Left(), Right())
tightlimits(ax1R, Left(), Right())

figure