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

# concurrent cogScore classification from taxonomic profiles
JLD2.@load "/home/guilherme/Documents/models/classification_currentCogScores_00to06_fromfunctions_results.jld"
JLD2.@load "/home/guilherme/Documents/models/classification_currentCogScores_06to12_fromfunctions_results.jld"
JLD2.@load "/home/guilherme/Documents/models/classification_currentCogScores_12to18_fromfunctions_results.jld"
JLD2.@load "/home/guilherme/Documents/models/classification_currentCogScores_18to24_fromfunctions_results.jld"
# concurrent cogScore regression from taxonomic profiles
JLD2.@load "/home/guilherme/Documents/models/regression_currentCogScores_00to06_fromfunctions_results.jld"
JLD2.@load "/home/guilherme/Documents/models/regression_currentCogScores_06to12_fromfunctions_results.jld"
JLD2.@load "/home/guilherme/Documents/models/regression_currentCogScores_12to18_fromfunctions_results.jld"
JLD2.@load "/home/guilherme/Documents/models/regression_currentCogScores_18to24_fromfunctions_results.jld"
# future cogScore classification from taxonomic profiles
#JLD2.@load "/home/guilherme/Documents/models/classification_futureCogScores_allselected_fromfunctions_results.jld"
# future cogScore regression from taxonomic profiles
#JLD2.@load "/home/guilherme/Documents/models/regression_futureCogScores_allselected_fromfunctions_results.jld"

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

scatterplot_cvsr_importance = function(clas_results::Dict, regr_results::Dict, plot_title::String, outpath::String; n = 50)

    clas_importances = rename!(average_classification_importances(clas_results), :Importance => :Clas_Importance)
    regr_importances = rename!(average_regression_importances(regr_results), :Importance => :Regr_Importance)

    clas_topvars = clas_importances.Variable[1:n]
    regr_topvars = regr_importances.Variable[1:n]
    # union_topvars = union(clas_topvars, regr_topvars)
    intersection_topvars = intersect(clas_topvars, regr_topvars)

    joined_importances = @chain clas_importances begin
        innerjoin(regr_importances, on = :Variable => :Variable)
        filter([:name, :grade_2020] => complex_filter)
        # subset(:Variable => x -> x .∈ Ref(union_topvars))
    end

    joined_importances.SetPertain = joined_importances.Variable .∈ Ref(intersection_topvars)

    # Create the figure
    figure = Figure(resolution = (800, 600))  

    # Build the axis
    ax1_1 = Axis(
        figure[1,1];
        xlabel = "Classification Importance",
        ylabel = "Regression Importance",
        title = plot_title
    )

    # Plot barplot
    sc1 = scatter!(
        ax1_1,
        joined_importances.Clas_Importance[joined_importances.SetPertain],
        joined_importances.Regr_Importance[joined_importances.SetPertain],
        color = :purple4
    )

    sc2 = scatter!(
        ax1_1,
        joined_importances.Clas_Importance[.!(joined_importances.SetPertain)],
        joined_importances.Regr_Importance[.!(joined_importances.SetPertain)],
        color = :goldenrod
    )

    # Plot Legend
    labels = ["Intersection", "Symmetric Difference"]
    # elements = [PolyElement(polycolor = el) for el in [:orange, :purple]]
    Legend(figure[2,1], [sc1, sc2], labels, orientation=:horizontal)

    # Export
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