# Figure4 - RandomForest prediction of cogScores from Taxonomic Profiles

```julia
using Resonance
using CairoMakie
using Statistics
using MLJ
using DecisionTree
using JLD2
using StableRNGs
```

## Loading the pretrained models

```julia

RandomForestClassifier = MLJ.@load RandomForestClassifier pkg=DecisionTree

#JLD2.@load "/home/guilherme/Documents/models/results_classification_currentCogScores_00to06_onlytaxa.jld" # Wrong name! redo!
JLD2.@load "/home/guilherme/Documents/models/results_classification_currentCogScores_06to12_onlytaxa.jld"
JLD2.@load "/home/guilherme/Documents/models/results_classification_currentCogScores_12to18_onlytaxa.jld"
JLD2.@load "/home/guilherme/Documents/models/results_classification_currentCogScores_18to24_onlytaxa.jld"
JLD2.@load "/home/guilherme/Documents/models/results_regression_currentCogScores_00to06_onlytaxa.jld"
JLD2.@load "/home/guilherme/Documents/models/results_regression_currentCogScores_06to12_onlytaxa.jld"
JLD2.@load "/home/guilherme/Documents/models/results_regression_currentCogScores_12to18_onlytaxa.jld"
JLD2.@load "/home/guilherme/Documents/models/results_regression_currentCogScores_18to24_onlytaxa.jld"
JLD2.@load "/home/guilherme/Documents/models/results_classification_futureCogScores_combined_onlytaxa.jld"
JLD2.@load "/home/guilherme/Documents/models/results_regression_futureCogScores_combined_onlytaxa.jld"
```

## Reporting the Figures of Merit for Classification

```julia

function report_classification_merit(classification_results::Dict)

    classification_merits = DataFrame(
    :Trial => collect(1:classification_results[:n_trials]),
    :Train_ACC => classification_results[:train_accuracies],
    :Test_ACC => classification_results[:test_accuracies],
    )

    push!(classification_merits, Dict(
        :Trial => 0,
        :Train_ACC => mean(classification_merits.Train_ACC),
        :Test_ACC => mean(classification_merits.Test_ACC)
    ))

    return classification_merits
    
end

function report_regression_merit(regression_results::Dict)

    regression_merits = DataFrame(
        :Trial => collect(1:regression_results[:n_trials]),
        :Train_MAE => regression_results[:train_maes],
        :Test_MAE => regression_results[:test_maes],
        :Train_MAPE => regression_results[:train_mapes],
        :Test_MAPE => regression_results[:test_mapes],
        :Train_COR => regression_results[:train_correlations],
        :Test_COR => regression_results[:test_correlations]
        )

    push!(regression_merits, Dict(
        :Trial => 0,
        :Train_MAE => mean(regression_merits.Train_MAE),
        :Test_MAE => mean(regression_merits.Test_MAE),
        :Train_MAPE => mean(regression_merits.Train_MAPE),
        :Test_MAPE => mean(regression_merits.Test_MAPE),
        :Train_COR => mean(regression_merits.Train_COR),
        :Test_COR => mean(regression_merits.Test_COR)
    ))

    return regression_merits
    
end

#report_classification_merit(classification_currentCogScores_00to06_fromtaxa_results)
# report_classification_merit(classification_currentCogScores_06to12_fromtaxa_results)
# report_classification_merit(classification_currentCogScores_12to18_fromtaxa_results)
# report_classification_merit(classification_currentCogScores_18to24_fromtaxa_results)
# report_classification_merit(classification_futureCogScores_allages_fromtaxa_results)

# report_regression_merit(regression_currentCogScores_00to06_fromtaxa_results)
# report_regression_merit(regression_currentCogScores_06to12_fromtaxa_results)
# report_regression_merit(regression_currentCogScores_12to18_fromtaxa_results)
# report_regression_merit(regression_currentCogScores_18to24_fromtaxa_results)
# report_regression_merit(regression_futureCogScores_allages_fromtaxa_results)

function build_confusion_matrix(classification_results::Dict, trial::Int)

    y = coerce(classification_results[:input_data][classification_results[:dataset_partitions][trial][2], :cogScore] .>= mean(classification_results[:input_data].cogScore), OrderedFactor)
    yhat = MLJ.predict_mode(classification_results[:models][trial], classification_results[:input_data][classification_results[:dataset_partitions][trial][2],6:end])
    confmat = ConfusionMatrix()(yhat, y)

end

function average_confusion_matrix(classification_results::Dict)

    matrix_vector = [ build_confusion_matrix(classification_results, i).mat for i in 1:classification_results[:n_trials] ]
    concatenated_matrix = vcat([ vec(mat)' for mat in matrix_vector ]...) # Column order is TN, FP, FN, TP !
    averages = [ sum(concatenated_matrix[:,i])/sum(concatenated_matrix) for i in 1:4]

end

function confmatrix2barplot(classification_results::Dict)

    confmat_inputs = (
        x = [1, 2, 2, 1], # Column order is TN, FP, FN, TP !
        value = average_confusion_matrix(classification_results),
        grp = [1, 2, 2, 1],
        color = [1, 2, 3, 4]
    )

end
```

## Building Figure

```julia

figure = Figure(resolution = (1200, 400))
confplot_colors = [:blue3, :red3, :tomato2, :dodgerblue2]
scatterplot

```

### Row 1 - Current cogScore classification

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
tbl1_1 = confmatrix2barplot(classification_06to12_results)
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

#### Plot 1 - Legend
```julia
labels = ["True Negatives", "False Positives", "False Negatives", "True Positives"]
elements = [PolyElement(polycolor = confplot_colors[i]) for i in 1:length(labels)]
Legend(figure[2,1:4], elements, labels, "Classification result", orientation=:horizontal)
```

### Row 2 - Current cogScore regression

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
yhat = GLM.predict(slope_correction, DataFrame( :yhat => MLJ.predict(selected_machine, X)))
sc_train = scatter!(ax, y[train], train_y_hat; color=:orange)
sc_test = scatter!(ax, y[test], test_y_hat; color=:blue)
ablines!(ax, 0, 1; color=:grey)
annotations!( ax, ["r = $(round(cor(y, yhat); digits = 2))"], [Point(60, 95)], textsize = 40)
```

#### Plot 1 - Legend
```julia
labels = ["True Negatives", "False Positives", "False Negatives", "True Positives"]
elements = [PolyElement(polycolor = confplot_colors[i]) for i in 1:length(labels)]
Legend(figure[4,1:4], elements, labels, "Classification result", orientation=:horizontal)
save("figure.png", figure)
```