using Resonance
using CairoMakie
using Statistics
using MLJ
using DecisionTree
using JLD2
using StableRNGs
using GLM

non_na_mean(a) = mean(a[.!(isnan.(a))])

RandomForestClassifier = MLJ.@load RandomForestClassifier pkg=DecisionTree
RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree
# concurrent cogScore classification from taxonomic profiles
JLD2.@load "/home/guilherme/Documents/models/classification_currentCogScores_00to06_fromtaxa_results.jld"
JLD2.@load "/home/guilherme/Documents/models/classification_currentCogScores_06to12_fromtaxa_results.jld"
JLD2.@load "/home/guilherme/Documents/models/classification_currentCogScores_12to18_fromtaxa_results.jld"
JLD2.@load "/home/guilherme/Documents/models/classification_currentCogScores_18to24_fromtaxa_results.jld"
# concurrent cogScore regression from taxonomic profiles
JLD2.@load "/home/guilherme/Documents/models/regression_currentCogScores_00to06_fromtaxa_results.jld"
JLD2.@load "/home/guilherme/Documents/models/regression_currentCogScores_06to12_fromtaxa_results.jld"
JLD2.@load "/home/guilherme/Documents/models/regression_currentCogScores_12to18_fromtaxa_results.jld"
JLD2.@load "/home/guilherme/Documents/models/regression_currentCogScores_18to24_fromtaxa_results.jld"
# future cogScore classification from taxonomic profiles
JLD2.@load "/home/guilherme/Documents/models/classification_futureCogScores_allselected_fromtaxa_results.jld"
# future cogScore regression from taxonomic profiles
JLD2.@load "/home/guilherme/Documents/models/regression_futureCogScores_allselected_fromtaxa_results.jld"

# Axis
function average_classification_importances(res::Dict)
    # map(non_na_mean, eachrow(reduce(hcat, [ impurity_importance(res[:models][i].fitresult[1]) for i in 1:5 ])))
    sort!(
        DataFrame(
            :Variable => names(res[:inputs_outputs][1]),
            :Importance => map(non_na_mean, eachrow(reduce(hcat, [ impurity_importance(res[:models][i].fitresult[1]) for i in 1:5 ])))
            ),
            :Importance, rev = true)
end

function average_regression_importances(res::Dict)
    # map(non_na_mean, eachrow(reduce(hcat, [ impurity_importance(res[:models][i].fitresult[1]) for i in 1:5 ])))
    sort!(
        DataFrame(
            :Variable => names(res[:inputs_outputs][1]),
            :Importance => map(non_na_mean, eachrow(reduce(hcat, [ impurity_importance(res[:models][i].fitresult) for i in 1:5 ])))
            ),
            :Importance, rev = true)
end

result_set = [
    classification_currentCogScores_00to06_fromtaxa_results,
    classification_currentCogScores_06to12_fromtaxa_results,
    classification_currentCogScores_12to18_fromtaxa_results,
    classification_currentCogScores_18to24_fromtaxa_results,
    classification_futureCogScores_allselected_fromtaxa_results,    
]


a = barplot_summary_importance(
    result_set,
    "Summary of variable importances throughout all classification models",
    "figures/importance_plots/Summarybarplot_classification_fromtaxa.png"
)

barplot_summary_importance = function(result_set::Vector{T} where T <: Dict, plot_title::String, outpath::String; n = 50)

    model_importances = map(x -> rename!(average_classification_importances(x), :Importance => :Clas_Importance), result_set)
    model_topns = map(x -> DataFrame(
        :Variable => x.Variable,
        :TopN => vcat(ones(n), zeros(nrow(x)-n))), model_importances
    )

    joined_importances = reduce((x, y) -> outerjoin(x, y, on = :Variable, makeunique = true), model_importances)
    joined_topns = reduce((x, y) -> outerjoin(x, y, on = :Variable, makeunique = true), model_topns)

    return(joined_importances)

    average_importances = DataFrame(
        :Variable => joined_importances.Variable,
        :AvgImportance => map(x -> mean(x[.!(ismissing.(x))]), eachrow(Matrix(joined_importances[:, 2:end])))
    )

    sum_topns = DataFrame(
        :Variable => joined_topns.Variable,
        :SumTopN => map(x -> floor(Int64, sum(x[.!(ismissing.(x))])), eachrow(Matrix(joined_topns[:, 2:end])))
    )

    plot_df = leftjoin(average_importances, sum_topns, on = :Variable);
    dropmissing!(plot_df);
    sort!(plot_df, :AvgImportance, rev = true)

    ######
    # Building Figure
    ######
    
    figure = Figure(resolution = (1000, 1000))
    bar_colors = ColorSchemes.viridis.colors[floor.(Int64, collect(range(256, 1, length = 1+maximum(plot_df.SumTopN))))]

    # Build the axis
    ax1_1 = Axis(
        figure[1,1];
        ylabel = "Most important Taxa",
        yticks = (reverse(collect(1:n)), plot_df.Variable[1:n]),
        xlabel = "Average MDI (Gini) Importance over all ($(length(result_set))) models",
        title = plot_title
    )
    
    # Plot barplot
    barplot!(ax1_1, reverse(collect(1:n)), plot_df.AvgImportance[1:n], color = bar_colors[plot_df.SumTopN[1:n] .+ 1], direction=:x)
    
    # Plot Legend
    labels = string.( collect(0:1:maximum(plot_df.SumTopN)) )
    elements = [ PolyElement(polycolor = bar_colors[i]) for i in 1:length(labels) ]
    Legend(figure[2,1], elements, labels, "Models for which this taxon is among the top $(n) important factors", orientation=:horizontal)

    # # Export
    save(outpath, figure)

end # end function

#####
# Outputting all plots
#####

## Classification, current, from taxa

scatterplot_cvsr_importance(
    classification_currentCogScores_00to06_fromtaxa_results,
    regression_currentCogScores_00to06_fromtaxa_results,
    "Importance comparison, concurrent cogScores, samples 0 to 6 months (n = 73)",
    "figures/importance_plots/Scatter_CvsR_currentCogScores_00to06_fromtaxa.png"
)

scatterplot_cvsr_importance(
    classification_currentCogScores_06to12_fromtaxa_results,
    regression_currentCogScores_06to12_fromtaxa_results,
    "Importance comparison, concurrent cogScores, samples 6 to 12 months (n = 61)",
    "figures/importance_plots/Scatter_CvsR_currentCogScores_06to12_fromtaxa.png"
)

plot_and_save_importance_barplot(
    classification_currentCogScores_12to18_fromtaxa_results,
    "Classification model, concurrent cogScores, samples 12 to 18 months (n = 39)",
    "figures/importance_plots/Importances_classification_currentCogScores_12to18_fromtaxa.png"
)

plot_and_save_importance_barplot(
    classification_currentCogScores_18to24_fromtaxa_results,
    "Classification model, concurrent cogScores, samples 18 to 24 months (n = 50)",
    "figures/importance_plots/Importances_regression_currentCogScores_18to24_fromtaxa.png"
)

## Regression, current, from taxa

plot_and_save_importance_barplot(
    regression_currentCogScores_00to06_fromtaxa_results,
    "Regression model, concurrent cogScores, samples 0 to 6 months (n = 73)",
    "figures/importance_plots/Importances_regression_currentCogScores_00to06_fromtaxa.png"
)

plot_and_save_importance_barplot(
    regression_currentCogScores_06to12_fromtaxa_results,
    "Regression model, concurrent cogScores, samples 6 to 12 months (n = 61)",
    "figures/importance_plots/Importances_regression_currentCogScores_06to12_fromtaxa.png"
)

plot_and_save_importance_barplot(
    regression_currentCogScores_12to18_fromtaxa_results,
    "Regression model, concurrent cogScores, samples 12 to 18 months (n = 39)",
    "figures/importance_plots/Importances_regression_currentCogScores_12to18_fromtaxa.png"
)

plot_and_save_importance_barplot(
    regression_currentCogScores_18to24_fromtaxa_results,
    "Regression model, concurrent cogScores, samples 18 to 24 months (n = 50)",
    "figures/importance_plots/Importances_regression_currentCogScores_18to24_fromtaxa.png"
)

## Classification, future, from taxa

plot_and_save_importance_barplot(
    classification_futureCogScores_allselected_fromtaxa_results,
    "Classification model, future cogScores, all qualifying samples\n(sample < 12mo, closest later assesment < 24mo)",
    "figures/importance_plots/Importances_classification_futureCogScores_allselected_fromtaxa.png"
)

## Regression, future, from taxa

plot_and_save_importance_barplot(
    regression_futureCogScores_allselected_fromtaxa_results,
    "Regression model, future cogScores, all qualifying samples\n(sample < 12mo, closest later assesment < 24mo)",
    "figures/importance_plots/Importances_regression_futureCogScores_allselected_fromtaxa.png"
)