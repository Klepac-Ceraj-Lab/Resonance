
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

# concurrent cogScore classification from taxonomic profiles
report_classification_merit(classification_currentCogScores_00to06_fromtaxa_results)
report_classification_merit(classification_currentCogScores_06to12_fromtaxa_results)
report_classification_merit(classification_currentCogScores_12to18_fromtaxa_results)
report_classification_merit(classification_currentCogScores_18to24_fromtaxa_results)
# concurrent cogScore regression from taxonomic profiles
report_regression_merit(regression_currentCogScores_00to06_fromtaxa_results)
report_regression_merit(regression_currentCogScores_06to12_fromtaxa_results)
report_regression_merit(regression_currentCogScores_12to18_fromtaxa_results)
report_regression_merit(regression_currentCogScores_18to24_fromtaxa_results)
# concurrent cogScore classification from functional profiles
report_classification_merit(classification_currentCogScores_00to06_fromfunctions_results)
report_classification_merit(classification_currentCogScores_06to12_fromfunctions_results)
report_classification_merit(classification_currentCogScores_12to18_fromfunctions_results)
report_classification_merit(classification_currentCogScores_18to24_fromfunctions_results)
# concurrent cogScore regression from functional profiles
report_regression_merit(regression_currentCogScores_00to06_fromfunctions_results)
report_regression_merit(regression_currentCogScores_06to12_fromfunctions_results)
report_regression_merit(regression_currentCogScores_12to18_fromfunctions_results)
report_regression_merit(regression_currentCogScores_18to24_fromfunctions_results)

#####
# Sanity check on Importances
#####

# classification_currentCogScores_12to18_fromtaxa_results[:importance_df] = sort!(
# DataFrame(
#     :Variable => names(classification_currentCogScores_12to18_fromtaxa_results[:inputs_outputs][1]),
#     :Importance => map(non_na_mean, eachrow(reduce(hcat, [ impurity_importance(classification_currentCogScores_12to18_fromtaxa_results[:models][i].fitresult[1]) for i in 1:5 ])))
#     ),
#     :Importance, rev = true)

# a = classification_currentCogScores_12to18_fromtaxa_results[:models]
# b = classification_currentCogScores_12to18_fromtaxa_results[:selected_trial][2]

# impurity_importance(a[5].fitresult[1])

# m = a[b]

### Plot group 1 - Concurrent cogScore classification

#### Plot 1.1 - Confusion Matrix / Accuracies for 00 to 06 months

# Axis

plot_and_save_importance_barplot = function(results::Dict, plot_title::String, outpath::String; n = 30)
    
    # Create the figure
    figure = Figure(resolution = (1000, 800))

    # Collect the importances
    res = results[:importance_df]

    # Build the axis
    ax1_1 = Axis(
        figure[1,1];
        ylabel = "Most important Taxa",
        yticks = (reverse(collect(1:n)), res.Variable[1:n]),
        #yticklabelrotation=-pi/4,
        xlabel = "Mean Decrease in Impurity (Gini) Importance",
        title = plot_title
    )

    # Plot barplot
    barplot!(ax1_1, reverse(collect(1:n)), res.Importance[1:n], color = :blue, direction=:x)

    # Export
    save(outpath, figure)

end # end function

#####
# Outputting all plots
#####

## Classification, current, from taxa

plot_and_save_importance_barplot(
    classification_currentCogScores_00to06_fromtaxa_results,
    "Classification model, concurrent cogScores, samples 0 to 6 months (n = 73)",
    "figures/importance_plots/Importances_classification_currentCogScores_00to06_fromtaxa.png"
)

plot_and_save_importance_barplot(
    classification_currentCogScores_06to12_fromtaxa_results,
    "Classification model, concurrent cogScores, samples 6 to 12 months (n = 61)",
    "figures/importance_plots/Importances_classification_currentCogScores_06to12_fromtaxa.png"
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