# Figure4 - RandomForest prediction of cogScores from Taxonomic Profiles

```julia
using Resonance
using CairoMakie
using Statistics
using MLJ
using DecisionTree
using JLD2
using StableRNGs
using GLM
```

## Loading the pretrained models

```julia
RandomForestClassifier = MLJ.@load RandomForestClassifier pkg=DecisionTree
# concurrent cogScore classification
JLD2.@load "/home/guilherme/Documents/models/classification_currentCogScores_00to06_fromtaxa_results.jld"
JLD2.@load "/home/guilherme/Documents/models/classification_currentCogScores_06to12_fromtaxa_results.jld"
JLD2.@load "/home/guilherme/Documents/models/classification_currentCogScores_12to18_fromtaxa_results.jld"
JLD2.@load "/home/guilherme/Documents/models/classification_currentCogScores_18to24_fromtaxa_results.jld"
# concurrent cogScore regression
JLD2.@load "/home/guilherme/Documents/models/regression_currentCogScores_00to06_fromtaxa_results.jld"
JLD2.@load "/home/guilherme/Documents/models/regression_currentCogScores_06to12_fromtaxa_results.jld"
JLD2.@load "/home/guilherme/Documents/models/regression_currentCogScores_12to18_fromtaxa_results.jld"
JLD2.@load "/home/guilherme/Documents/models/regression_currentCogScores_18to24_fromtaxa_results.jld"
# future cogScore classification
JLD2.@load "/home/guilherme/Documents/models/regression_currentCogScores_18to24_fromtaxa_results.jld"
# future cogScore regression
JLD2.@load "/home/guilherme/Documents/models/results_regression_futureCogScores_combined_onlytaxa.jld"
```

## Reporting the Figures of Merit for Classification

```julia
# concurrent cogScore classification
report_classification_merit(classification_currentCogScores_00to06_fromtaxa_results)
report_classification_merit(classification_currentCogScores_06to12_fromtaxa_results)
report_classification_merit(classification_currentCogScores_12to18_fromtaxa_results)
report_classification_merit(classification_currentCogScores_18to24_fromtaxa_results)
# concurrent cogScore regression
report_regression_merit(regression_currentCogScores_00to06_fromtaxa_results)
report_regression_merit(regression_currentCogScores_06to12_fromtaxa_results)
report_regression_merit(regression_currentCogScores_12to18_fromtaxa_results)
report_regression_merit(regression_currentCogScores_18to24_fromtaxa_results)
```

## Building Figure

```julia

figure = Figure(resolution = (1200, 800))
confplot_colors = [:blue3, :red3, :tomato2, :dodgerblue2]
scatterplot

```

### Plot group 1 - Concurrent cogScore classification

#### Plot 1.1 - Confusion Matrix / Accuracies for 00 to 06 months

```julia
# Axis
ax1_1 = Axis(
    figure[1,1];
    xlabel = "Classifications",
    xticks = (1:2, ["Correct", "Incorrect"]),
    ylabel = "Proportion",
    title = "0 to 6 months"
)
ylims!(ax1_1, [0.0, 1.0])

# Plot
tbl1_1 = confmatrix2barplot(classification_00to06_results)
barplot!(ax1_1, tbl1_1.x, tbl1_1.value,
        stack = tbl1_1.grp,
        color = confplot_colors[tbl1_1.color])
```

#### Plot 1.2 - Confusion Matrix / Accuracies for 00 to 06 months
```julia
# Axis
ax1_2 = Axis(
    figure[1,2];
    xlabel = "Classifications",
    xticks = (1:2, ["Correct", "Incorrect"]),
    ylabel = "Proportion",
    title = "6 to 12 months"
)
ylims!(ax1_2, [0.0, 1.0])

# Plot
tbl1_2 = confmatrix2barplot(classification_06to12_results)
barplot!(ax1_2, tbl1_2.x, tbl1_2.value,
        stack = tbl1_2.grp,
        color = confplot_colors[tbl1_2.color])
```

#### Plot 1.3 - Confusion Matrix / Accuracies for 00 to 06 months
```julia
# Axis
ax1_3 = Axis(
    figure[1,3];
    xlabel = "Classifications",
    xticks = (1:2, ["Correct", "Incorrect"]),
    ylabel = "Proportion",
    title = "12 to 18 months"
)
ylims!(ax1_3, [0.0, 1.0])

# Plot
tbl1_3 = confmatrix2barplot(classification_12to18_results)
barplot!(ax1_3, tbl1_3.x, tbl1_3.value,
        stack = tbl1_3.grp,
        color = confplot_colors[tbl1_3.color])
```

#### Plot 1.4 - Confusion Matrix / Accuracies for 00 to 06 months
```julia
# Axis
ax1_4 = Axis(
    figure[1,4];
    xlabel = "Classifications",
    xticks = (1:2, ["Correct", "Incorrect"]),
    ylabel = "Proportion",
    title = "18 to 24 months"
)
ylims!(ax1_4, [0.0, 1.0])

# Plot
tbl1_4 = confmatrix2barplot(classification_18to24_results)
barplot!(ax1_4, tbl1_4.x, tbl1_4.value,
        stack = tbl1_4.grp,
        color = confplot_colors[tbl1_4.color])
```

#### Plot group 1 - Legend and title
```julia
Label(fig[1, :, Top()], "Average figures of merit for binary classification of cogScore (above/below average for age bracket)\nat the time of stool collection, for the independent test set", valign = :bottom,
    # font = "TeX Gyre Heros Bold",
    padding = (0, 0, 5, 0))

labels = ["True Negatives", "False Positives", "False Negatives", "True Positives"]
elements = [PolyElement(polycolor = confplot_colors[i]) for i in 1:length(labels)]
Legend(figure[2,1:4], elements, labels, "Classification result", orientation=:horizontal)
```

### plot group 2 - Concurrent cogScore regression

#### Plot 2.1 - Predictions vs ground truth for 00 to 06 months

```julia
# Axis
ax3_1 = Axis(
    figure[3,1];
    xlabel = "Ground Truth cogScore",
    ylabel = "Predicted cogScore",
    title = "0 to 6 months"
)

# Plot
y, yhat, train, test = regression_bestprediction(regression_currentCogScores_00to06_fromtaxa_results)
sc_train = scatter!(ax, y[train], train_y_hat; color=:orange)
sc_test = scatter!(ax, y[test], test_y_hat; color=:purple)
ablines!(ax, 0, 1; color=:grey)
annotations!( ax, ["r = $(round(cor(y, yhat); digits = 2))"], [Point(60, 95)], textsize = 40)
```

#### Plot Group 2 - Legend and Title
```julia
Label(fig[1, :, Top()], "Predicted vs ground truth values for regression of cogScore\nat the time of stool collection for the best validated model", valign = :bottom,
    # font = "TeX Gyre Heros Bold",
    padding = (0, 0, 5, 0))

labels = ["True Negatives", "False Positives", "False Negatives", "True Positives"]
elements = [PolyElement(polycolor = confplot_colors[i]) for i in 1:length(labels)]
Legend(figure[4,1:4], elements, labels, "Classification result", orientation=:horizontal)
save("figure.png", figure)
```