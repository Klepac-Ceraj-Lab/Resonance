# Figure3 - RandomForest prediction of cogScores from DEM, Taxonomic and Functional Profiles

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
        FNR = Float64[],
        PPV = Float64[],
        NPV = Float64[],
        F1 = Float64[]
    )

    for (i, res) in enumerate(resvec)

        thismodel_measures = DataFrame(
            accuracy = Float64[],
            AUROC = Float64[],
            TPR = Float64[],
            TNR = Float64[],
            FPR = Float64[],
            FNR = Float64[],
            PPV = Float64[],
            NPV = Float64[],
            F1 = Float64[]
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
                    PPV = ( tps )/( tps + fps ),
                    NPV = ( tns )/( tns + fns ),
                    F1 = ( tps ) / ( tps + 0.5 * (fps + fns) )
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
                PPV = thismodel_summary[7],
                NPV = thismodel_summary[8],
                F1 = thismodel_summary[9],
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
# concurrent cogScore classification from taxonomic profiles
JLD2.@load "models/classification_currentCogScores_00to06mo_onlydemo.jld"
JLD2.@load "models/classification_currentCogScores_00to06mo_onlytaxa.jld"
JLD2.@load "models/classification_currentCogScores_00to06mo_demoplustaxa.jld"
JLD2.@load "models/classification_currentCogScores_00to06mo_onlyecs.jld"
JLD2.@load "models/classification_currentCogScores_00to06mo_demoplusecs.jld"
JLD2.@load "models/classification_currentCogScores_06to12mo_onlydemo.jld"
JLD2.@load "models/classification_currentCogScores_06to12mo_onlytaxa.jld"
JLD2.@load "models/classification_currentCogScores_06to12mo_demoplustaxa.jld"
JLD2.@load "models/classification_currentCogScores_06to12mo_onlyecs.jld"
JLD2.@load "models/classification_currentCogScores_06to12mo_demoplusecs.jld"
JLD2.@load "models/classification_currentCogScores_12to18mo_onlydemo.jld"
JLD2.@load "models/classification_currentCogScores_12to18mo_onlytaxa.jld"
JLD2.@load "models/classification_currentCogScores_12to18mo_demoplustaxa.jld"
JLD2.@load "models/classification_currentCogScores_12to18mo_onlyecs.jld"
JLD2.@load "models/classification_currentCogScores_12to18mo_demoplusecs.jld"
JLD2.@load "models/classification_currentCogScores_18to120mo_onlydemo.jld"
JLD2.@load "models/classification_currentCogScores_18to120mo_onlytaxa.jld"
JLD2.@load "models/classification_currentCogScores_18to120mo_demoplustaxa.jld"
JLD2.@load "models/classification_currentCogScores_18to120mo_onlyecs.jld"
JLD2.@load "models/classification_currentCogScores_18to120mo_demoplusecs.jld"
# concurrent cogScore reression from taxonomic profiles
JLD2.@load "models/regression_currentCogScores_00to06mo_onlydemo.jld"
JLD2.@load "models/regression_currentCogScores_00to06mo_onlytaxa.jld"
JLD2.@load "models/regression_currentCogScores_00to06mo_demoplustaxa.jld"
JLD2.@load "models/regression_currentCogScores_00to06mo_onlyecs.jld"
JLD2.@load "models/regression_currentCogScores_00to06mo_demoplusecs.jld"
JLD2.@load "models/regression_currentCogScores_06to12mo_onlydemo.jld"
JLD2.@load "models/regression_currentCogScores_06to12mo_onlytaxa.jld"
JLD2.@load "models/regression_currentCogScores_06to12mo_demoplustaxa.jld"
JLD2.@load "models/regression_currentCogScores_06to12mo_onlyecs.jld"
JLD2.@load "models/regression_currentCogScores_06to12mo_demoplusecs.jld"
JLD2.@load "models/regression_currentCogScores_12to18mo_onlydemo.jld"
JLD2.@load "models/regression_currentCogScores_12to18mo_onlytaxa.jld"
JLD2.@load "models/regression_currentCogScores_12to18mo_demoplustaxa.jld"
JLD2.@load "models/regression_currentCogScores_12to18mo_onlyecs.jld"
JLD2.@load "models/regression_currentCogScores_12to18mo_demoplusecs.jld"
JLD2.@load "models/regression_currentCogScores_18to120mo_onlydemo.jld"
JLD2.@load "models/regression_currentCogScores_18to120mo_onlytaxa.jld"
JLD2.@load "models/regression_currentCogScores_18to120mo_demoplustaxa.jld"
JLD2.@load "models/regression_currentCogScores_18to120mo_onlyecs.jld"
JLD2.@load "models/regression_currentCogScores_18to120mo_demoplusecs.jld"
#```

## Initializing Figure 3

#```julia
figure = Figure(resolution = (1920, 1080))
left_subfig = GridLayout(figure[1,1])
right_subfig = GridLayout(figure[1,2])

A_subfig = GridLayout(left_subfig[1,1])
B_subfig = GridLayout(left_subfig[2,1])

CDEFGH_subfig = GridLayout(right_subfig[1,1])
IJKLMN_subfig = GridLayout(right_subfig[2,1])
Legends_subfig = GridLayout(right_subfig[3,1])

colsize!(figure.layout, 1, Relative(0.4))
colsize!(figure.layout, 2, Relative(0.6))

confplot_colors = [:blue3, :red3, :tomato2, :dodgerblue2]
#```

### Plot panel A - Concurrent cogScore classification

#```julia
classification_results = [
    classification_currentCogScores_00to06mo_onlydemo,
    classification_currentCogScores_00to06mo_onlytaxa,
    classification_currentCogScores_00to06mo_demoplustaxa,
    classification_currentCogScores_00to06mo_onlyecs,
    classification_currentCogScores_00to06mo_demoplusecs,
    classification_currentCogScores_06to12mo_onlydemo,
    classification_currentCogScores_06to12mo_onlytaxa,
    classification_currentCogScores_06to12mo_demoplustaxa,
    classification_currentCogScores_06to12mo_onlyecs,
    classification_currentCogScores_06to12mo_demoplusecs,
    classification_currentCogScores_12to18mo_onlydemo,
    classification_currentCogScores_12to18mo_onlytaxa,
    classification_currentCogScores_12to18mo_demoplustaxa,
    classification_currentCogScores_12to18mo_onlyecs,
    classification_currentCogScores_12to18mo_demoplusecs,
    classification_currentCogScores_18to120mo_onlydemo,
    classification_currentCogScores_18to120mo_onlytaxa,
    classification_currentCogScores_18to120mo_demoplustaxa,
    classification_currentCogScores_18to120mo_onlyecs,
    classification_currentCogScores_18to120mo_demoplusecs
]

xs, values, grps, colors = confmatrices2barplots(classification_results)

axA1 = Axis(
    A_subfig[1,1];
    xticks = (1:20, ["DEM", "TAXA", "TAXA+DEM", "ECS", "ECS+DEM", "DEM", "TAXA", "TAXA+DEM", "ECS", "ECS+DEM", "DEM", "TAXA", "TAXA+DEM", "ECS", "ECS+DEM", "DEM", "TAXA", "TAXA+DEM", "ECS", "ECS+DEM"]),
    ylabel = "Test set average Proportion",
    title = "A - Classification - above/below 50th percentile"
)

ylims!(axA1, [-0.4, 0.4])

axA2 = Axis(
    A_subfig[1,1];
    xticks = (1:20, [
        "", "", "Subjects from birth to 6 months", "", "",
        "", "", "Subjects from 6 to 12 months", "", "",
        "", "", "Subjects from 12 to 18 months", "", "",
        "", "", "Subjects from 18 months to 120 months", "", ""
        ]),
    xticklabelpad = 30.0
)

linkxaxes!(axA1, axA2)

hideydecorations!(axA2)
hidexdecorations!(axA2, ticklabels = false)

barplot!(
    axA1, xs, values,
    dodge = grps,
    stack = grps,
    color = confplot_colors[colors],
    width = 0.7
)

labels = ["True Negatives", "False Positives", "False Negatives", "True Positives"]
elements = [PolyElement(polycolor = confplot_colors[i]) for i in 1:length(labels)]
Legend(A_subfig[2,1], elements, labels, "Classification result", orientation=:horizontal)

#```

### Plot panel B - Concurrent cogScore regression
#```julia
regression_results = [
    regression_currentCogScores_00to06mo_onlydemo,
    regression_currentCogScores_00to06mo_onlytaxa,
    regression_currentCogScores_00to06mo_demoplustaxa,
    regression_currentCogScores_00to06mo_onlyecs,
    regression_currentCogScores_00to06mo_demoplusecs,
    regression_currentCogScores_06to12mo_onlydemo,
    regression_currentCogScores_06to12mo_onlytaxa,
    regression_currentCogScores_06to12mo_demoplustaxa,
    regression_currentCogScores_06to12mo_onlyecs,
    regression_currentCogScores_06to12mo_demoplusecs,
    regression_currentCogScores_12to18mo_onlydemo,
    regression_currentCogScores_12to18mo_onlytaxa,
    regression_currentCogScores_12to18mo_demoplustaxa,
    regression_currentCogScores_12to18mo_onlyecs,
    regression_currentCogScores_12to18mo_demoplusecs,
    regression_currentCogScores_18to120mo_onlydemo,
    regression_currentCogScores_18to120mo_onlytaxa,
    regression_currentCogScores_18to120mo_demoplustaxa,
    regression_currentCogScores_18to120mo_onlyecs,
    regression_currentCogScores_18to120mo_demoplusecs
]

axB1 = Axis(
    B_subfig[1,1];
        xticks = (1:20, ["DEM", "TAXA", "TAXA+DEM", "ECS", "ECS+DEM", "DEM", "TAXA", "TAXA+DEM", "ECS", "ECS+DEM", "DEM", "TAXA", "TAXA+DEM", "ECS", "ECS+DEM", "DEM", "TAXA", "TAXA+DEM", "ECS", "ECS+DEM"]),
    ylabel = "Test set Mean Abssolut Error\nTest set Person's R coefficient",
    title = "B - Regression - numeric percentile"
)

#ylims!(axB1, [0.0, 1.0])

axB2 = Axis(
    B_subfig[1,1];
        xticks = (1:20, [
        "", "", "Subjects from birth to 6 months", "", "",
        "", "", "Subjects from 6 to 12 months", "", "",
        "", "", "Subjects from 12 to 18 months", "", "",
        "", "", "Subjects from 18 months to 120 months", "", ""
        ]),
    xticklabelpad = 30.0
)

linkxaxes!(axB1, axB2)

hideydecorations!(axB2)
hidexdecorations!(axB2, ticklabels = false)

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
    selected_rmses[1:5],
    color = :pink
)
lines!(
    axB1,
    6:10,
    selected_rmses[6:10],
    color = :pink
)
lines!(
    axB1,
    11:15,
    selected_rmses[11:15],
    color = :pink
)
lines!(
    axB1,
    16:20,
    selected_rmses[16:20],
    color = :pink
)
scatter!(
    axB1,
    1:20,
    selected_rmses,
    color = :pink,
    size = 5
)

lines!(
    axB1,
    1:5,
    selected_correlations[1:5],
    color = :green
)
lines!(
    axB1,
    6:10,
    selected_correlations[6:10],
    color = :green
)
lines!(
    axB1,
    11:15,
    selected_correlations[11:15],
    color = :green
)
lines!(
    axB1,
    16:20,
    selected_correlations[16:20],
    color = :green
)
scatter!(
    axB1,
    1:20,
    selected_correlations,
    color = :green,
    size = 5
)

labels = ["Correlation coefficient (r)", "Mean Absolute Error"]
elements = [PolyElement(polycolor = c) for c in [:green, :pink]]
Legend(B_subfig[2,1], elements, labels, "Classification result", orientation=:horizontal)

#```

### Plot panel C - Taxonomic profile PCA, 0-6 mo, color by percentile
#```julia

axC = Axis(
    CDEFGH_subfig[1,1];
    xlabel = "PCo1",
    ylabel = "PCo2",
    title = "C - TAXA Percentile PCA - 0 to 6 months"
)

hidedecorations!(axC, grid = false, label = false)

C_pca_model = fit(
    PCA,
    transpose(Matrix(classification_currentCogScores_00to06mo_onlytaxa.inputs_outputs[1])),
    maxoutdim=2
)

C_prcomps = permutedims(MultivariateStats.predict(
    C_pca_model,
    transpose(Matrix(classification_currentCogScores_00to06mo_onlytaxa.inputs_outputs[1]))
))

scatter!(
    axC,
    C_prcomps[:,1], C_prcomps[:,2];
    color = classification_currentCogScores_00to06mo_onlytaxa.original_data.cogScorePercentile
)

#```

### Plot panel D - Taxonomic profile ROC curves, 0-6 mo, color by Microbiome vs pure DEM
#```julia

axD = Axis(
    CDEFGH_subfig[1,2];
    xlabel = "False-positive rate",
    ylabel = "True-positive rate",
    title = "D - ROC Curves, 0-6, TAXA"
)

for i in 1:classification_currentCogScores_00to06mo_onlydemo.n_splits
    fprs, tprs, ts = res_roc_curve(classification_currentCogScores_00to06mo_onlydemo, i; sample_set = "all")
    lines!(axD, fprs, tprs; color = RGBA(1.0, 0.0, 0.0, 0.3), transparency = true)
end

for i in 1:classification_currentCogScores_00to06mo_onlytaxa.n_splits
    fprs, tprs, ts = res_roc_curve(classification_currentCogScores_00to06mo_onlytaxa, i; sample_set = "all")
    lines!(axD, fprs, tprs; color = RGBA(0.0, 0.0, 1.0, 0.3), transparency = true)
end

#```

### Plot panel E - Taxonomic profile regression scatterplots, 0-6 mo, color by train/test
#```julia

axE = Axis(
    CDEFGH_subfig[1,3];
    xlabel = "Ground Truth Percentile",
    ylabel = "Predicted percentile",
    title = "E - Scatterplots, 0-6, TAXA"
)

xlims!(axE, [0.0, 1.0])
ylims!(axE, [0.0, 1.0])

singlemodel_merit_scatterplot!(
    axE, regression_currentCogScores_00to06mo_onlytaxa;
    split_index = selected_indexes[2],
    traincorr_pos = Point(0.05, 0.85),
    testcorr_pos = Point(0.05, 0.7)
)

#```

### Plot panel F - Functional profile PCA, 0-6 mo, color by percentile
#```julia

axF = Axis(
    CDEFGH_subfig[2,1];
    xlabel = "PCo1",
    ylabel = "PCo2",
    title = "F - ECS Percentile PCA - 0 to 6 months"
)

hidedecorations!(axF, grid = false, label = false)

F_pca_model = fit(
    PCA,
    transpose(Matrix(classification_currentCogScores_00to06mo_onlyecs.inputs_outputs[1]));
    maxoutdim=2,
    pratio = 1.0
)

F_prcomps = permutedims(MultivariateStats.predict(
    F_pca_model,
    transpose(Matrix(classification_currentCogScores_00to06mo_onlyecs.inputs_outputs[1]))
))

scatter!(
    axF,
    F_prcomps[:,1], F_prcomps[:,2];
    color = classification_currentCogScores_00to06mo_onlyecs.original_data.cogScorePercentile
)

#```

### Plot panel G - Functional profile ROC curves, 0-6 mo, color by Microbiome vs pure DEM
#```julia

axG = Axis(
    CDEFGH_subfig[2,2];
    xlabel = "False-positive rate",
    ylabel = "True-positive rate",
    title = "G - ROC Curves, 0-6, ECS"
)

for i in 1:classification_currentCogScores_00to06mo_onlydemo.n_splits
    fprs, tprs, ts = res_roc_curve(classification_currentCogScores_00to06mo_onlydemo, i; sample_set = "all")
    lines!(axG, fprs, tprs; color = RGBA(1.0, 0.0, 0.0, 0.3), transparency = true)
end

for i in 1:classification_currentCogScores_00to06mo_onlyecs.n_splits
    fprs, tprs, ts = res_roc_curve(classification_currentCogScores_00to06mo_onlyecs, i; sample_set = "all")
    lines!(axG, fprs, tprs; color = RGBA(0.0, 0.0, 1.0, 0.3), transparency = true)
end

figure
#```

### Plot panel H - Regression scatterplots, 0-6 mo, functional profiles
#```julia

axH = Axis(
    CDEFGH_subfig[2,3];
    xlabel = "Ground Truth Percentile",
    ylabel = "Predicted percentile",
    title = "H - Scatterplots, 0-6, ECS"
)

xlims!(axH, [0.0, 1.0])
ylims!(axH, [0.0, 1.0])

singlemodel_merit_scatterplot!(
    axH, regression_currentCogScores_00to06mo_onlyecs;
    split_index = selected_indexes[4],
    traincorr_pos = Point(0.05, 0.85),
    testcorr_pos = Point(0.05, 0.7)
)

figure
#```


### Plot panel I - Taxonomic profile PCA, 18-120mo, color by percentile
#```julia

axI = Axis(
    IJKLMN_subfig[1,1];
    xlabel = "PCo1",
    ylabel = "PCo2",
    title = "I - TAXA Percentile PCA - 18-120 months"
)

hidedecorations!(axI, grid = false, label = false)

I_pca_model = fit(
    PCA,
    transpose(Matrix(classification_currentCogScores_18to120mo_onlytaxa.inputs_outputs[1]));
    maxoutdim=2,
    pratio = 1.0
)

I_prcomps = permutedims(MultivariateStats.predict(
    I_pca_model,
    transpose(Matrix(classification_currentCogScores_18to120mo_onlytaxa.inputs_outputs[1]))
))

scatter!(
    axI,
    I_prcomps[:,1], I_prcomps[:,2];
    color = classification_currentCogScores_18to120mo_onlytaxa.original_data.cogScorePercentile
)

#```

### Plot panel J - ROC curves, 18-120 mo, taxonomic profiles
#```julia

axJ = Axis(
    IJKLMN_subfig[1,2];
    xlabel = "False-positive rate",
    ylabel = "True-positive rate",
    title = "J - ROC Curves, 18-120, TAXA"
)

for i in 1:classification_currentCogScores_18to120mo_onlydemo.n_splits
    fprs, tprs, ts = res_roc_curve(classification_currentCogScores_18to120mo_onlydemo, i; sample_set = "all")
    lines!(axJ, fprs, tprs; color = RGBA(1.0, 0.0, 0.0, 0.3), transparency = true)
end

for i in 1:classification_currentCogScores_18to120mo_onlytaxa.n_splits
    fprs, tprs, ts = res_roc_curve(classification_currentCogScores_18to120mo_onlytaxa, i; sample_set = "all")
    lines!(axJ, fprs, tprs; color = RGBA(0.0, 0.0, 1.0, 0.3), transparency = true)
end
#```

### Plot panel K - Regression scatterplots, 18-120 mo, taxonomic profiles
#```julia

axK = Axis(
    IJKLMN_subfig[1,3];
    xlabel = "Ground Truth Percentile",
    ylabel = "Predicted percentile",
    title = "K - Scatterplots, 18-120, TAXA"
)

xlims!(axK, [0.0, 1.0])
ylims!(axK, [0.0, 1.0])

singlemodel_merit_scatterplot!(
    axK, regression_currentCogScores_18to120mo_onlytaxa;
    split_index = selected_indexes[7],
    traincorr_pos = Point(0.05, 0.85),
    testcorr_pos = Point(0.05, 0.7)
)

#```

### Plot panel L - Functional profile PCA, 18-120mo, color by percentile
#```julia

axL = Axis(
    IJKLMN_subfig[2,1];
    xlabel = "PCo1",
    ylabel = "PCo2",
    title = "L - ECS Percentile PCA - 18-120 months"
)

hidedecorations!(axL, grid = false, label = false)

L_pca_model = fit(
    PCA,
    transpose(Matrix(classification_currentCogScores_18to120mo_onlyecs.inputs_outputs[1]));
    maxoutdim=2,
    pratio = 1.0
)

L_prcomps = permutedims(MultivariateStats.predict(
    L_pca_model,
    transpose(Matrix(classification_currentCogScores_18to120mo_onlyecs.inputs_outputs[1]))
))

scatter!(
    axL,
    L_prcomps[:,1], L_prcomps[:,2];
    color = classification_currentCogScores_18to120mo_onlytaxa.original_data.cogScorePercentile
)

#```

### Plot panel M - ROC curves, 18-120 mo, functional profiles
#```julia

axM = Axis(
   IJKLMN_subfig[2,2];
    xlabel = "False-positive rate",
    ylabel = "True-positive rate",
    title = "M - ROC Curves, 18-120, ECS"
)

for i in 1:classification_currentCogScores_18to120mo_onlydemo.n_splits
    fprs, tprs, ts = res_roc_curve(classification_currentCogScores_18to120mo_onlydemo, i; sample_set = "all")
    lines!(axM, fprs, tprs; color = RGBA(1.0, 0.0, 0.0, 0.3), transparency = true)
end

for i in 1:classification_currentCogScores_18to120mo_onlyecs.n_splits
    fprs, tprs, ts = res_roc_curve(classification_currentCogScores_18to120mo_onlyecs, i; sample_set = "all")
    lines!(axM, fprs, tprs; color = RGBA(0.0, 0.0, 1.0, 0.3), transparency = true)
end

#```

### Plot panel N - Regression scatterplots, 18-120 mo, functional profiles
#```julia

axN = Axis(
    IJKLMN_subfig[2,3];
    xlabel = "Ground Truth Percentile",
    ylabel = "Predicted percentile",
    title = "N - Scatterplots, 18-120, ECS"
)

xlims!(axN, [0.0, 1.0])
ylims!(axN, [0.0, 1.0])

singlemodel_merit_scatterplot!(
    axN, regression_currentCogScores_18to120mo_onlyecs;
    split_index = selected_indexes[7],
    traincorr_pos = Point(0.05, 0.85),
    testcorr_pos = Point(0.05, 0.7)
)

#```

### Legends

Colorbar(Legends_subfig[1, 1], colorrange = (0.0, 1.0), colormap = :viridis, label="Cogscore Percentile", vertical = false)

labels = [ "Microbiome", "Demographics" ]
elements = [PolyElement(polycolor = el) for el in [ :blue, :red ] ]
Legend(Legends_subfig[1, 2], elements, labels, "Model ROC", orientation = :horizontal)

labels = [ "Train set", "Test set" ]
elements = [PolyElement(polycolor = el) for el in [ :orange, :purple ] ]
Legend(Legends_subfig[1, 3], elements, labels, "Sample set", orientation = :horizontal)

colsize!(Legends_subfig, 1, Relative(0.3))
colsize!(Legends_subfig, 2, Relative(0.4))
colsize!(Legends_subfig, 3, Relative(0.3))

figure

save("figures/Figure3.png", figure)

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