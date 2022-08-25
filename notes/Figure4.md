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
JLD2.@load "/home/guilherme/Documents/models/classification_futureCogScores_allselected_fromtaxa_results.jld"
# future cogScore regression
JLD2.@load "/home/guilherme/Documents/models/regression_futureCogScores_allselected_fromtaxa_results.jld"
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

## Building Figure 4.1 - current cogscores from taxonomic profiles

```julia
figure = Figure(resolution = (1200, 900))
confplot_colors = [:blue3, :red3, :tomato2, :dodgerblue2]
```

### Plot group 1 - Concurrent cogScore classification

#### Plot 1.1 - Confusion Matrix / Accuracies for 00 to 06 months

```julia
# Axis
ax1_1 = Axis(
    figure[1,1];
    # xlabel = "Classifications",
    xticks = (1:2, ["Correct", "Incorrect"]),
    ylabel = "Proportion",
    title = "0 to 6 months (n = 73)"
)
ylims!(ax1_1, [0.0, 1.0])

# Plot
tbl1_1 = confmatrix2barplot(classification_currentCogScores_00to06_fromtaxa_results)
barplot!(ax1_1, tbl1_1.x, tbl1_1.value,
        stack = tbl1_1.grp,
        color = confplot_colors[tbl1_1.color])
```

#### Plot 1.2 - Confusion Matrix / Accuracies for 00 to 06 months
```julia
# Axis
ax1_2 = Axis(
    figure[1,2];
    # xlabel = "Classifications",
    xticks = (1:2, ["Correct", "Incorrect"]),
    ylabel = "Proportion",
    title = "6 to 12 months (n = 61)"
)
ylims!(ax1_2, [0.0, 1.0])

# Plot
tbl1_2 = confmatrix2barplot(classification_currentCogScores_06to12_fromtaxa_results)
barplot!(ax1_2, tbl1_2.x, tbl1_2.value,
        stack = tbl1_2.grp,
        color = confplot_colors[tbl1_2.color])
```

#### Plot 1.3 - Confusion Matrix / Accuracies for 00 to 06 months
```julia
# Axis
ax1_3 = Axis(
    figure[1,3];
    # xlabel = "Classifications",
    xticks = (1:2, ["Correct", "Incorrect"]),
    ylabel = "Proportion",
    title = "12 to 18 months (n = 39)"
)
ylims!(ax1_3, [0.0, 1.0])

# Plot
tbl1_3 = confmatrix2barplot(classification_currentCogScores_12to18_fromtaxa_results)
barplot!(ax1_3, tbl1_3.x, tbl1_3.value,
        stack = tbl1_3.grp,
        color = confplot_colors[tbl1_3.color])
```

#### Plot 1.4 - Confusion Matrix / Accuracies for 00 to 06 months
```julia
# Axis
ax1_4 = Axis(
    figure[1,4];
    # xlabel = "Classifications",
    xticks = (1:2, ["Correct", "Incorrect"]),
    ylabel = "Proportion",
    title = "18 to 24 months (n = 50)"
)
ylims!(ax1_4, [0.0, 1.0])

# Plot
tbl1_4 = confmatrix2barplot(classification_currentCogScores_18to24_fromtaxa_results)
barplot!(ax1_4, tbl1_4.x, tbl1_4.value,
        stack = tbl1_4.grp,
        color = confplot_colors[tbl1_4.color])
```

#### Plot group 1 - Legend and title
```julia
Label(figure[1, :, Top()], "Average figures of merit for binary classification of concurrent cogScore (above/below average for age bracket)\nat the time of stool collection, for the independent test set of each split", valign = :bottom,
    # font = "TeX Gyre Heros Bold",
    padding = (0, 0, 30, 0))

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
    title = "0 to 6 months (n = 73)"
)

# Plot
y, yhat, train, test = regression_bestprediction(regression_currentCogScores_00to06_fromtaxa_results)
sc_train = scatter!(ax3_1, y[train], yhat[train]; color=:orange)
sc_test = scatter!(ax3_1, y[test], yhat[test]; color=:purple)
ablines!(ax3_1, 0, 1; color=:grey)
annotations!( ax3_1, ["r = $(round(cor(y, yhat); digits = 2))"], [Point(60, 110)], textsize = 20)
```

#### Plot 2.2 - Predictions vs ground truth for 06 to 12 months

```julia
# Axis
ax3_2 = Axis(
    figure[3,2];
    xlabel = "Ground Truth cogScore",
    ylabel = "Predicted cogScore",
    title = "6 to 12 months (n = 61)"
)

# Plot
y, yhat, train, test = regression_bestprediction(regression_currentCogScores_06to12_fromtaxa_results)
sc_train = scatter!(ax3_2, y[train], yhat[train]; color=:orange)
sc_test = scatter!(ax3_2, y[test], yhat[test]; color=:purple)
ablines!(ax3_2, 0, 1; color=:grey)
annotations!( ax3_2, ["r = $(round(cor(y, yhat); digits = 2))"], [Point(60, 110)], textsize = 20)
```

#### Plot 2.3 - Predictions vs ground truth for 12 to 18 months

```julia
# Axis
ax3_3 = Axis(
    figure[3,3];
    xlabel = "Ground Truth cogScore",
    ylabel = "Predicted cogScore",
    title = "12 to 18 months (n = 39)"
)

# Plot
y, yhat, train, test = regression_bestprediction(regression_currentCogScores_12to18_fromtaxa_results)
sc_train = scatter!(ax3_3, y[train], yhat[train]; color=:orange)
sc_test = scatter!(ax3_3, y[test], yhat[test]; color=:purple)
ablines!(ax3_3, 0, 1; color=:grey)
annotations!( ax3_3, ["r = $(round(cor(y, yhat); digits = 2))"], [Point(60, 98)], textsize = 20)
```

#### Plot 2.4 - Predictions vs ground truth for 18 to 24 months

```julia
# Axis
ax3_4 = Axis(
    figure[3,4];
    xlabel = "Ground Truth cogScore",
    ylabel = "Predicted cogScore",
    title = "18 to 24 months (n = 50)"
)

# Plot
y, yhat, train, test = regression_bestprediction(regression_currentCogScores_18to24_fromtaxa_results)
sc_train = scatter!(ax3_4, y[train], yhat[train]; color=:orange)
sc_test = scatter!(ax3_4, y[test], yhat[test]; color=:purple)
ablines!(ax3_4, 0, 1; color=:grey)
annotations!( ax3_4, ["r = $(round(cor(y, yhat); digits = 2))"], [Point(50, 102)], textsize = 20)
```

#### Plot Group 2 - Legend and Title
```julia
Label(figure[3, :, Top()], "Predicted vs ground truth values for regression of concurrent cogScore at the time of stool collection for the best (split + validated model) combination", valign = :bottom,
    # font = "TeX Gyre Heros Bold",
    padding = (0, 0, 30, 0))

labels = ["Training samples", "Test samples"]
elements = [PolyElement(polycolor = el) for el in [:orange, :purple]]
Legend(figure[4,1:4], elements, labels, "Sample set", orientation=:horizontal)
save("figures/Figure4_current_taxa.png", figure); save("figures/Figure4_current_taxa.svg", figure)
```

## Building Figure 4.2 - future cogscores from taxonomic profiles

```julia
figure = Figure(resolution = (900, 600))
confplot_colors = [:blue3, :red3, :tomato2, :dodgerblue2]
```

#### Plot 2.1 - Confusion Matrix / Accuracies for all ages

```julia
# Axis
ax1_1 = Axis(
    figure[1,1];
    # xlabel = "Classifications",
    xticks = (1:2, ["Correct", "Incorrect"]),
    ylabel = "Proportion",
    title = "All qualifying samples (n = 112)"
)
ylims!(ax1_1, [0.0, 1.0])

# Plot
tbl1_1 = confmatrix2barplot(classification_futureCogScores_allselected_fromtaxa_results)
barplot!(ax1_1, tbl1_1.x, tbl1_1.value,
        stack = tbl1_1.grp,
        color = confplot_colors[tbl1_1.color])
```

#### Plot 2.2 - Predictions vs ground truth for 06 to 12 months

```julia
# Axis
ax3_2 = Axis(
    figure[1,2];
    xlabel = "Ground Truth cogScore",
    ylabel = "Predicted cogScore",
    title = "All qualifying samples (n = 112)"
)

# Plot
y, yhat, train, test = regression_bestprediction(regression_futureCogScores_allselected_fromtaxa_results)
sc_train = scatter!(ax3_2, y[train], yhat[train]; color=:orange)
sc_test = scatter!(ax3_2, y[test], yhat[test]; color=:purple)
ablines!(ax3_2, 0, 1; color=:grey)
annotations!( ax3_2, ["r = $(round(cor(y, yhat); digits = 2))"], [Point(60, 110)], textsize = 20)
```

#### Legends and title
```julia
labels = ["True Negatives", "False Positives", "False Negatives", "True Positives"]
elements = [PolyElement(polycolor = confplot_colors[i]) for i in 1:length(labels)]
Legend(figure[2,1], elements, labels, "Classification result", nbanks = 2, orientation = :horizontal, valign = :top)

labels = ["Training samples", "Test samples"]
elements = [PolyElement(polycolor = el) for el in [:orange, :purple]]
Legend(figure[2,2], elements, labels, "Sample set", nbanks = 1, orientation=:horizontal, valign = :top)

Label(figure[1, :, Top()], "Binary classification and regression results for prediction of future cognitive scores from present taxonomic profiles for all\nqualifying samples with stool collection under 12 months and at least one future cognitive assessment before 24 months", valign = :bottom,
    # font = "TeX Gyre Heros Bold",
    padding = (0, 50, 40, 0))

save("figures/Figure4_future_taxa.png", figure); save("figures/Figure4_future_taxa.svg", figure)
```