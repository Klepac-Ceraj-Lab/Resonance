#####
# Notebook DX5 - Brain univariate result plots
#####

using Resonance
using CairoMakie
using Statistics
using MLJ
using DecisionTree
using JLD2
using StableRNGs
using GLM
using ProgressMeter

#####
# Functions
#####

function predict_bestsplit(regression_results::UnivariateRandomForestRegressor)

    selected_split = regression_results.selected_split[2]
    X, y = regression_results.inputs_outputs
    train, test = regression_results.dataset_partitions[selected_split]
    slope_correction = regression_results.slope_correction[selected_split]
    selected_machine = regression_results.models[selected_split]
    yhat = GLM.predict(slope_correction, DataFrame( :yhat => MLJ.predict(selected_machine, X)))

    return y, yhat, train, test
end

function singlemodel_merit_scatterplot!(
    figure::Figure,
    res::UnivariateRandomForestRegressor,
    pos::Tuple{Int64, Int64},
    plot_title::String)

    # Build the axis
    ax = Axis(
        figure[pos[1],pos[2]];
        xlabel = "Ground Truth",
        ylabel = "Prediction",
        title = plot_title
    )

    # Plot barplot
    y, yhat, train, test = predict_bestsplit(res)
    scatter!(ax, y[train], yhat[train]; color=:orange)
    scatter!(ax, y[test], yhat[test]; color=:purple)
    ablines!(ax, 0, 1; color=:grey)
    annotations!( ax, ["r = $(round(cor(y, yhat); digits = 2))"], [Point(1.1*minimum(y), 0.9*maximum(yhat))], textsize = 20)

    return figure

end # end function

function plot_all_results!(
    figure::Figure,
    result_set::Vector{T} where T <: ResonancePredictor,
    plot_positions::Vector{Tuple{Int64, Int64}}
    )

    for i in eachindex(result_set)
        singlemodel_merit_scatterplot!(figure, result_set[i], plot_positions[i], string(result_set[i].name))
    end

    return figure

end

#####
# Load Results
#####

RandomForestRegressor= MLJ.@load RandomForestRegressor pkg=DecisionTree
JLD2.@load "models/brain_univariate_regression.jld" brain_results

for m in brain_results
    println(m.name)
    @show report_merits(m)
end

using CairoMakie
fig = Figure(resolution = (4096, 4096))

cols = collect(1:10)
rows = collect(1:10)
plot_positions = vec(permutedims(collect(Base.product(rows, cols))))

plot_all_results!(fig, brain_results, plot_positions)

save("figures/Brain_univariate_scatterplots.png", fig)