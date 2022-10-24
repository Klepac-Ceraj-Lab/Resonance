# Figure3 - RandomForest prediction of cogScores from SES, Taxonomic and Functional Profiles

#```julia
using Resonance
using CairoMakie
using Statistics
using MLJ
using DecisionTree
using JLD2
using StableRNGs
using GLM
#```

## Loading the pretrained models

#```julia
RandomForestClassifier = MLJ.@load RandomForestClassifier pkg=DecisionTree
RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree
# concurrent cogScore classification from taxonomic profiles
JLD2.@load "models/classification_currentCogScores_00to06mo_onlyses.jld"
JLD2.@load "models/classification_currentCogScores_00to06mo_onlytaxa.jld"
JLD2.@load "models/classification_currentCogScores_00to06mo_onlyecs.jld"
JLD2.@load "models/classification_currentCogScores_00to06mo_sesplusecs.jld"
JLD2.@load "models/classification_currentCogScores_00to06mo_sesplustaxa.jld"
#```

## Building Figure 4.1 - current cogscores from taxonomic profiles

#```julia
figure = Figure(resolution = (1200, 1800))
confplot_colors = [:blue3, :red3, :tomato2, :dodgerblue2]
#```

### Plot group 1 - Concurrent cogScore classification from taxonomic profiles

#```julia
singlemodel_merit_barplot!(figure, classification_currentCogScores_00to06_fromtaxa_results, (1,1), "0 to 6 months (n = 73)"; confplot_colors = confplot_colors)
singlemodel_merit_barplot!(figure, classification_currentCogScores_06to12_fromtaxa_results, (1,2), "6 to 12 months (n = 61)"; confplot_colors = confplot_colors)
singlemodel_merit_barplot!(figure, classification_currentCogScores_12to18_fromtaxa_results, (1,3), "12 to 18 months (n = 39)"; confplot_colors = confplot_colors)
singlemodel_merit_barplot!(figure, classification_currentCogScores_18to24_fromtaxa_results, (1,4), "18 to 24 months (n = 50)"; confplot_colors = confplot_colors)
#```

#### Plot group 1 - Legend and title
#```julia
Label(figure[1, :, Top()], "Average figures of merit for binary classification of concurrent cogScore (above/below average for age bracket)\nat the time of stool collection, for the independent test set of each split", valign = :bottom,
    # font = "TeX Gyre Heros Bold",
    padding = (0, 0, 30, 0))

labels = ["True Negatives", "False Positives", "False Negatives", "True Positives"]
elements = [PolyElement(polycolor = confplot_colors[i]) for i in 1:length(labels)]
Legend(figure[2,1:4], elements, labels, "Classification result", orientation=:horizontal)
#```

### plot group 2 - Concurrent cogScore regression from taxonomic profiles

#```julia
singlemodel_merit_scatterplot!(figure, regression_currentCogScores_00to06_fromtaxa_results, (3,1), "0 to 6 months (n = 73)")
singlemodel_merit_scatterplot!(figure, regression_currentCogScores_06to12_fromtaxa_results, (3,2), "6 to 12 months (n = 61)")
singlemodel_merit_scatterplot!(figure, regression_currentCogScores_12to18_fromtaxa_results, (3,3), "12 to 18 months (n = 39)")
singlemodel_merit_scatterplot!(figure, regression_currentCogScores_18to24_fromtaxa_results, (3,4), "18 to 24 months (n = 50)")
#```

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

### Plot group 3 - Concurrent cogScore classification from functional profiles

```julia
singlemodel_merit_barplot!(figure, classification_currentCogScores_00to06_fromfunctions_results, (5,1), "0 to 6 months (n = 73)"; confplot_colors = confplot_colors)
singlemodel_merit_barplot!(figure, classification_currentCogScores_06to12_fromfunctions_results, (5,2), "6 to 12 months (n = 61)"; confplot_colors = confplot_colors)
singlemodel_merit_barplot!(figure, classification_currentCogScores_12to18_fromfunctions_results, (5,3), "12 to 18 months (n = 39)"; confplot_colors = confplot_colors)
singlemodel_merit_barplot!(figure, classification_currentCogScores_18to24_fromfunctions_results, (5,4), "18 to 24 months (n = 50)"; confplot_colors = confplot_colors)
```

#### Plot group 3 - Legend and title

```julia
Label(figure[5, :, Top()], "Average figures of merit for binary classification of concurrent cogScore (above/below average for age bracket)\nat the time of stool collection, for the independent test set of each split", valign = :bottom,
    # font = "TeX Gyre Heros Bold",
    padding = (0, 0, 30, 0))

labels = ["True Negatives", "False Positives", "False Negatives", "True Positives"]
elements = [PolyElement(polycolor = confplot_colors[i]) for i in 1:length(labels)]
Legend(figure[6,1:4], elements, labels, "Classification result", orientation=:horizontal)
```

### plot group 4 - Concurrent cogScore regression

```julia
singlemodel_merit_scatterplot!(figure, regression_currentCogScores_00to06_fromfunctions_results, (7,1), "0 to 6 months (n = 73)")
singlemodel_merit_scatterplot!(figure, regression_currentCogScores_06to12_fromfunctions_results, (7,2), "6 to 12 months (n = 61)")
singlemodel_merit_scatterplot!(figure, regression_currentCogScores_12to18_fromfunctions_results, (7,3), "12 to 18 months (n = 39)")
singlemodel_merit_scatterplot!(figure, regression_currentCogScores_18to24_fromfunctions_results, (7,4), "18 to 24 months (n = 50)")
```

#### Plot Group 4 - Legend and Title
```julia
Label(figure[7, :, Top()], "Predicted vs ground truth values for regression of concurrent cogScore at the time of stool collection for the best (split + validated model) combination", valign = :bottom,
    # font = "TeX Gyre Heros Bold",
    padding = (0, 0, 30, 0))

labels = ["Training samples", "Test samples"]
elements = [PolyElement(polycolor = el) for el in [:orange, :purple]]
Legend(figure[8,1:4], elements, labels, "Sample set", orientation=:horizontal)
save("figures/Figure4_current.png", figure); save("figures/Figure4_current.svg", figure)
```

#####
# FUTURE
#####

## Building Figure 4.2 - future cogscores from taxonomic profiles

```julia
figure = Figure(resolution = (1200, 600))
confplot_colors = [:blue3, :red3, :tomato2, :dodgerblue2]
```

#### Plot 2.1 - Confusion Matrix / Accuracies for all ages

```julia
singlemodel_merit_barplot!(figure, classification_futureCogScores_allselected_fromtaxa_results, (1,1), "All qualifying samples (n = 112)"; confplot_colors = confplot_colors)
singlemodel_merit_barplot!(figure, classification_futureCogScores_allselected_fromtaxa_results, (1,2), "All qualifying samples (n = 112)"; confplot_colors = confplot_colors)
singlemodel_merit_scatterplot!(figure, regression_futureCogScores_allselected_fromfunctions_results, (1,3), "All qualifying samples (n = 112)")
singlemodel_merit_scatterplot!(figure, regression_futureCogScores_allselected_fromfunctions_results, (1,4), "All qualifying samples (n = 112)")
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

save("figures/Figure4_future.png", figure); save("figures/Figure4_future.svg", figure)
```