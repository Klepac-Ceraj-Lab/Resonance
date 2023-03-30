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

Legend(
    figure[3, 1],
    [
        PolyElement(; color = this_barcolor, strokewidth = 1, strokecolor = :black),
        LineElement(; color = this_curvecolor, linewidth = 5)
    ] , [
        "Individual relative Importance",
        "Cumulative relative importance"

    ],
    tellheight = false,
    tellwidth = false,
    margin = (10, 10, 10, 10),
    halign = :right, valign = :bottom, orientation = :vertical
)

figure