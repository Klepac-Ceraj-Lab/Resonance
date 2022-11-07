# Figure4 - RandomForest prediction of cogScores from DEM, Taxonomic and Functional Profiles

#```julia
using Resonance
using CairoMakie
using Colors
using Statistics
using MultivariateStats
using MLJ
using DecisionTree
using JLD2
using StableRNGs
using GLM
#```

## Defining relevant functions

### Calculate AUROC
function res_auroc(res::UnivariateRandomForestClassifier, split_index::Int; sample_set = "test")

    (split_index > length(res.models)) && error("Out of Bounds model/split index selected")
    (split_index == 0) && (split_index = res.selected_split[2])

    train, test = res.dataset_partitions[split_index]

    if sample_set == "test"
        X = res.inputs_outputs[1][test, :]
        y = res.inputs_outputs[2][test]
    elseif sample_set == "train"
        X = res.inputs_outputs[1][train, :]
        y = res.inputs_outputs[2][train]
    else
        X = res.inputs_outputs[1]
        y = res.inputs_outputs[2]
    end

    yhat = Resonance.predict_proba(res, X; split_index=split_index)
    auroc = AreaUnderCurve()(yhat, y)

    return auroc

end

### Calculate classification measures
function classification_measures(resvec::Vector{UnivariateRandomForestClassifier}; sample_set = "test")

    allmodels_measures = DataFrame(
        model = String[],
        accuracy = Float64[],
        AUROC = Float64[],
        TPR = Float64[],
        TNR = Float64[],
        FPR = Float64[],
        FNR = Float64[]
    )

    for (i, res) in enumerate(resvec)

        thismodel_measures = DataFrame(
            accuracy = Float64[],
            AUROC = Float64[],
            TPR = Float64[],
            TNR = Float64[],
            FPR = Float64[],
            FNR = Float64[]
        )

        for this_split in 1:res.n_splits

            auroc = res_auroc(res, this_split; sample_set = sample_set)
            
            confmat = build_confusion_matrix(res, this_split).mat
            tns = confmat[1,1]
            fps = confmat[2,1]
            fns = confmat[1,2]
            tps = confmat[2,2]

            push!(
                thismodel_measures, 
                (
                    accuracy = (tns + tps)/( tns + tps + fns + fps ),
                    AUROC = auroc,
                    TPR = ( tps )/( tps + fns ),
                    TNR = ( tns )/( tns + fps ),
                    FPR = ( fps )/( tns + fps ),
                    FNR = ( fns )/( tps + fns ),
                )
            )
        end

        thismodel_summary = map(mean, eachcol(thismodel_measures))

        push!(
            allmodels_measures, 
            (
                model = res.name,
                accuracy = thismodel_summary[1],
                AUROC = thismodel_summary[2],
                TPR = thismodel_summary[3],
                TNR = thismodel_summary[4],
                FPR = thismodel_summary[5],
                FNR = thismodel_summary[6],
            )
        )
    end

    return(allmodels_measures)

end

### Multimodel stacked bar plot
function confmatrices2barplots(resvec::Vector{UnivariateRandomForestClassifier})

    xs = Int64[]
    values = Float64[]
    grps = Int64[]
    colors = Int64[]

    for (i, res) in enumerate(resvec)

        append!(xs, repeat( [ i ], 4))
        append!(values, average_confusion_matrix(res) .- 0.25)
        append!(grps, [ 1, 2, 2, 1 ] )
        append!(colors, [ 1, 2, 3, 4] )
    end

    return xs, values, grps, colors

end

### Generate ROC curve
function res_roc_curve(res::UnivariateRandomForestClassifier, split_index::Int; sample_set = "test")

    (split_index > length(res.models)) && error("Out of Bounds model/split index selected")
    (split_index == 0) && (split_index = res.selected_split[2])

    train, test = res.dataset_partitions[split_index]

    if sample_set == "test"
        X = res.inputs_outputs[1][test, :]
        y = res.inputs_outputs[2][test]
    elseif sample_set == "train"
        X = res.inputs_outputs[1][train, :]
        y = res.inputs_outputs[2][train]
    else
        X = res.inputs_outputs[1]
        y = res.inputs_outputs[2]
    end

    yhat = Resonance.predict_proba(res, X; split_index=split_index)

    fprs, tprs, ts = roc_curve(yhat, y)

    return fprs, tprs, ts

end

## Loading the pretrained models

#```julia
RandomForestClassifier = MLJ.@load RandomForestClassifier pkg=DecisionTree
RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree

# future cogScore classification from taxonomic profiles
JLD2.@load "models/classification_futureCogScores_00to12mo_onlydemo.jld"
JLD2.@load "models/classification_futureCogScores_00to12mo_onlytaxa.jld"
JLD2.@load "models/classification_futureCogScores_00to12mo_demoplustaxa.jld"
JLD2.@load "models/classification_futureCogScores_00to12mo_onlyecs.jld"
JLD2.@load "models/classification_futureCogScores_00to12mo_demoplusecs.jld"
# future cogScore reression from taxonomic profiles
JLD2.@load "models/regression_futureCogScores_00to12mo_onlydemo.jld"
JLD2.@load "models/regression_futureCogScores_00to12mo_onlytaxa.jld"
JLD2.@load "models/regression_futureCogScores_00to12mo_demoplustaxa.jld"
JLD2.@load "models/regression_futureCogScores_00to12mo_onlyecs.jld"
JLD2.@load "models/regression_futureCogScores_00to12mo_demoplusecs.jld"

classification_results = [
    classification_futureCogScores_00to12mo_onlydemo,
    classification_futureCogScores_00to12mo_onlytaxa,
    classification_futureCogScores_00to12mo_demoplustaxa,
    classification_futureCogScores_00to12mo_onlyecs,
    classification_futureCogScores_00to12mo_demoplusecs
]

regression_results = [
    regression_futureCogScores_00to12mo_onlydemo,
    regression_futureCogScores_00to12mo_onlytaxa,
    regression_futureCogScores_00to12mo_demoplustaxa,
    regression_futureCogScores_00to12mo_onlyecs,
    regression_futureCogScores_00to12mo_demoplusecs
]

#```

## Initializing Figure 3

#```julia
figure = Figure(resolution = (1920, 1080))
left_subfig = GridLayout(figure[1,1])
right_subfig = GridLayout(figure[1,2])

A_subfig = GridLayout(left_subfig[1,1])
B_subfig = GridLayout(left_subfig[2,1])

CDEF_subfig = GridLayout(right_subfig[1,1])
Legends_subfig = GridLayout(right_subfig[2,1])

colsize!(figure.layout, 1, Relative(0.4))
colsize!(figure.layout, 2, Relative(0.6))

confplot_colors = [:blue3, :red3, :tomato2, :dodgerblue2]
#```

### Plot panel A - Concurrent cogScore classification

#```julia
xs, values, grps, colors = confmatrices2barplots(classification_results)

axA1 = Axis(
    A_subfig[1,1];
    xticks = (1:5, ["DEM", "TAXA", "TAXA+DEM", "ECS", "ECS+DEM"]),
    ylabel = "Confusio matrix enrichment raletive to a random classifier",
    title = "A - Classification - above/below 50th percentile"
)

ylims!(axA1, [-0.5, 0.5])

barplot!(
    axA1, xs, values,
    dodge = grps,
    stack = grps,
    color = confplot_colors[colors]
)

labels = ["True Negatives", "False Positives", "False Negatives", "True Positives"]
elements = [PolyElement(polycolor = confplot_colors[i]) for i in 1:length(labels)]
Legend(A_subfig[2,1], elements, labels, "Classification result", orientation=:horizontal)

#```

### Plot panel B - Concurrent cogScore regression
#```julia

axB1 = Axis(
    B_subfig[1,1];
    xticks = (1:5, ["DEM", "TAXA", "TAXA+DEM", "ECS", "ECS+DEM"]),
    ylabel = "Test set Mean Abssolut Error\nTest set Person's R coefficient",
    title = "B - Regression - numeric percentile"
)

# ## Median Version
plot_points = DataFrame(
    [
        report_merits(m)[end, :] for m in regression_results
    ]
)
selected_rmses = plot_points.Test_RMSE
selected_correlations = plot_points.Test_COR .^ 2
selected_indexes = Vector{Int64}()
for i in eachindex(regression_results)
    bottomdist_models = findall(
        report_merits(regression_results[i])[1:end-2, :Test_RMSE] .- selected_rmses[i] .< 0.0
    )
    plot_model = findmin(
        abs.(report_merits(regression_results[i])[bottomdist_models, :Test_RMSE] .- selected_rmses[i])
    )
    push!(selected_indexes, plot_model[2])
end

# ## Mean Version
# plot_points = DataFrame(
#     [
#         report_merits(m)[end-1, :] for m in regression_results
#     ]
# )
# selected_maes = plot_points.Test_RMSE
# selected_correlations = plot_points.Test_COR .^ 2
# minmse_splits = [ findmin(report_merits(m)[1:end-1, :Test_RMSE]) for m in regression_results ]
# selected_indexes = [ el[2] for el in minmse_splits]

lines!(
    axB1,
    1:5,
    selected_rmses,
    color = :pink
)
scatter!(
    axB1,
    1:5,
    selected_rmses,
    color = :pink,
    size = 5
)

lines!(
    axB1,
    1:5,
    selected_correlations,
    color = :green
)
scatter!(
    axB1,
    1:5,
    selected_correlations,
    color = :green,
    size = 5
)

labels = ["Determination coefficient (R2)", "Mean Absolute Error"]
elements = [PolyElement(polycolor = c) for c in [:green, :pink]]
Legend(B_subfig[2,1], elements, labels, "Classification result", orientation=:horizontal)

#```


### Plot panel C - Taxonomic profile ROC curves, 0-6 mo, color by Microbiome vs pure DEM
#```julia

axC = Axis(
    CDEF_subfig[1,1];
    xlabel = "False-positive rate",
    ylabel = "True-positive rate",
    title = "C - ROC Curves, TAXA"
)

for i in 1:classification_futureCogScores_00to12mo_onlydemo.n_splits
    fprs, tprs, ts = res_roc_curve(classification_futureCogScores_00to12mo_onlydemo, i; sample_set = "all")
    lines!(axC, fprs, tprs; color = RGBA(1.0, 0.0, 0.0, 0.4), transparency = true)
end

for i in 1:classification_futureCogScores_00to12mo_onlytaxa.n_splits
    fprs, tprs, ts = res_roc_curve(classification_futureCogScores_00to12mo_onlytaxa, i; sample_set = "all")
    lines!(axC, fprs, tprs; color = RGBA(0.0, 0.0, 1.0, 0.4), transparency = true)
end

#```

### Plot panel D - Taxonomic profile regression scatterplots, 0-6 mo, color by train/test
#```julia

axD = Axis(
    CDEF_subfig[1,2];
    xlabel = "Ground Truth Percentile",
    ylabel = "Predicted percentile",
    title = "D - Scatterplots, TAXA"
)

xlims!(axD, [0.0, 1.0])
ylims!(axD, [0.0, 1.0])

singlemodel_merit_scatterplot!(
    axD, regression_futureCogScores_00to12mo_onlytaxa;
    split_index = selected_indexes[2],
    traincorr_pos = Point(0.05, 0.85),
    testcorr_pos = Point(0.05, 0.7)
)

#```

### Plot panel E - Functional profile ROC curves, 0-6 mo, color by Microbiome vs pure DEM
#```julia

axE = Axis(
    CDEF_subfig[2,1];
    xlabel = "False-positive rate",
    ylabel = "True-positive rate",
    title = "G - ROC Curves, ECS"
)

for i in 1:classification_futureCogScores_00to12mo_onlydemo.n_splits
    fprs, tprs, ts = res_roc_curve(classification_futureCogScores_00to12mo_onlydemo, i; sample_set = "all")
    lines!(axE, fprs, tprs; color = RGBA(1.0, 0.0, 0.0, 0.4), transparency = true)
end

for i in 1:classification_futureCogScores_00to12mo_onlyecs.n_splits
    fprs, tprs, ts = res_roc_curve(classification_futureCogScores_00to12mo_onlyecs, i; sample_set = "all")
    lines!(axE, fprs, tprs; color = RGBA(0.0, 0.0, 1.0, 0.4), transparency = true)
end

figure
#```

### Plot panel F - Regression scatterplots, 0-6 mo, functional profiles
#```julia

axF = Axis(
    CDEF_subfig[2,2];
    xlabel = "Ground Truth Percentile",
    ylabel = "Predicted percentile",
    title = "H - Scatterplots, 0-6, ECS"
)

xlims!(axF, [0.0, 1.0])
ylims!(axF, [0.0, 1.0])

singlemodel_merit_scatterplot!(
    axF, regression_futureCogScores_00to12mo_onlyecs;
    split_index = selected_indexes[4],
    traincorr_pos = Point(0.05, 0.85),
    testcorr_pos = Point(0.05, 0.7)
)
#```

### Legends

labels = [ "Microbiome", "Demographics" ]
elements = [PolyElement(polycolor = el) for el in [ :blue, :red ] ]
Legend(Legends_subfig[1, 1], elements, labels, "Model ROC", orientation = :horizontal)

labels = [ "Train set", "Test set" ]
elements = [PolyElement(polycolor = el) for el in [ :orange, :purple ] ]
Legend(Legends_subfig[1, 2], elements, labels, "Sample set", orientation = :horizontal)

figure

save("figures/Figure4.png", figure)

# Accessory figures

acc_fig3 = Figure(resolution = (1000, 400))
axO = Axis(
    acc_fig3[1,1];
    xlabel = "Ground Truth Percentile",
    ylabel = "Predicted percentile",
    title = "O - Scatterplots, 00-06, Demographics"
)

xlims!(axO, [0.0, 1.0])
ylims!(axO, [0.0, 1.0])

singlemodel_merit_scatterplot!(
    axO, regression_currentCogScores_00to06mo_onlydemo;
    split_index = selected_indexes[1],
    traincorr_pos = Point(0.05, 0.85),
    testcorr_pos = Point(0.05, 0.7)
)

axP = Axis(
    acc_fig3[1,2];
    xlabel = "Ground Truth Percentile",
    ylabel = "Predicted percentile",
    title = "P - Scatterplots, 18-120, Demographics"
)

xlims!(axP, [0.0, 1.0])
ylims!(axP, [0.0, 1.0])

singlemodel_merit_scatterplot!(
    axP, regression_currentCogScores_18to120mo_onlydemo;
    split_index = selected_indexes[6],
    traincorr_pos = Point(0.05, 0.85),
    testcorr_pos = Point(0.05, 0.7)
)