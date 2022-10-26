# Figure3 - RandomForest prediction of cogScores from SES, Taxonomic and Functional Profiles

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
        append!(values, average_confusion_matrix(res))
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
JLD2.@load "models/classification_currentCogScores_00to06mo_onlyses.jld"
JLD2.@load "models/classification_currentCogScores_00to06mo_onlytaxa.jld"
JLD2.@load "models/classification_currentCogScores_00to06mo_sesplustaxa.jld"
JLD2.@load "models/classification_currentCogScores_00to06mo_onlyecs.jld"
JLD2.@load "models/classification_currentCogScores_00to06mo_sesplusecs.jld"
JLD2.@load "models/classification_currentCogScores_18to120mo_onlyses.jld"
JLD2.@load "models/classification_currentCogScores_18to120mo_onlytaxa.jld"
JLD2.@load "models/classification_currentCogScores_18to120mo_sesplustaxa.jld"
JLD2.@load "models/classification_currentCogScores_18to120mo_onlyecs.jld"
JLD2.@load "models/classification_currentCogScores_18to120mo_sesplusecs.jld"
# concurrent cogScore reression from taxonomic profiles
JLD2.@load "models/regression_currentCogScores_00to06mo_onlyses.jld"
JLD2.@load "models/regression_currentCogScores_00to06mo_onlytaxa.jld"
JLD2.@load "models/regression_currentCogScores_00to06mo_sesplustaxa.jld"
JLD2.@load "models/regression_currentCogScores_00to06mo_onlyecs.jld"
JLD2.@load "models/regression_currentCogScores_00to06mo_sesplusecs.jld"
JLD2.@load "models/regression_currentCogScores_18to120mo_onlyses.jld"
JLD2.@load "models/regression_currentCogScores_18to120mo_onlytaxa.jld"
JLD2.@load "models/regression_currentCogScores_18to120mo_sesplustaxa.jld"
JLD2.@load "models/regression_currentCogScores_18to120mo_onlyecs.jld"
#```

## Initializing Figure 3

#```julia
figure = Figure(resolution = (1920, 1080))

A_subfig = GridLayout(figure[1,1])
B_subfig = GridLayout(figure[2,1])
CDEFGH_subfig = GridLayout(figure[1,2])
IJKLMN_subfig = GridLayout(figure[2,2])

colsize!(figure.layout, 1, Relative(0.4))
colsize!(figure.layout, 2, Relative(0.6))

confplot_colors = [:blue3, :red3, :tomato2, :dodgerblue2]
#```

### Plot panel A - Concurrent cogScore classification

#```julia
classification_results = [
    classification_currentCogScores_00to06mo_onlyses,
    classification_currentCogScores_00to06mo_onlytaxa,
    classification_currentCogScores_00to06mo_sesplustaxa,
    classification_currentCogScores_00to06mo_onlyecs,
    classification_currentCogScores_00to06mo_sesplusecs,
    classification_currentCogScores_18to120mo_onlyses,
    classification_currentCogScores_18to120mo_onlytaxa,
    classification_currentCogScores_18to120mo_sesplustaxa,
    classification_currentCogScores_18to120mo_onlyecs,
    classification_currentCogScores_18to120mo_sesplusecs
]

xs, values, grps, colors = confmatrices2barplots(classification_results)

axA1 = Axis(
    A_subfig[1,1];
    xticks = (1:10, ["SES", "TAXA", "TAXA+SES", "ECS", "ECS+SES", "SES", "TAXA", "TAXA+SES", "ECS", "ECS+SES"]),
    ylabel = "Proportion",
    title = "Classification - above/below 50th percentile"
)

ylims!(axA1, [0.0, 1.0])

axA2 = Axis(
    A_subfig[1,1];
    xticks = (1:10, ["", "", "Subjects from birth to 6 months", "", "", "", "", "Subjects from 18 months to 120 months", "", ""]),
    xticklabelpad = 30.0
)

linkxaxes!(axA1, axA2)

hideydecorations!(axA2)
hidexdecorations!(axA2, ticklabels = false)

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
regression_results = [
    regression_currentCogScores_00to06mo_onlyses,
    regression_currentCogScores_00to06mo_onlytaxa,
    regression_currentCogScores_00to06mo_sesplustaxa,
    regression_currentCogScores_00to06mo_onlyecs,
    regression_currentCogScores_00to06mo_sesplusecs,
    regression_currentCogScores_18to120mo_onlyses,
    regression_currentCogScores_18to120mo_onlytaxa,
    regression_currentCogScores_18to120mo_sesplustaxa,
    regression_currentCogScores_18to120mo_onlyecs,
    regression_currentCogScores_00to06mo_sesplusecs
]

merits = DataFrame( [ report_merits(regression_results[i])[end, :] for i in eachindex(regression_results) ] )

axB1 = Axis(
    B_subfig[1,1];
    xticks = (1:10, ["SES", "TAXA", "TAXA+SES", "ECS", "ECS+SES", "SES", "TAXA", "TAXA+SES", "ECS", "ECS+SES"]),
    ylabel = "Proportion",
    title = "Classification - above/below 50th percentile"
)

ylims!(axB1, [0.0, 1.0])

axB2 = Axis(
    B_subfig[1,1];
    xticks = (1:10, ["", "", "Subjects from birth to 6 months", "", "", "", "", "Subjects from 18 months to 120 months", "", ""]),
    xticklabelpad = 30.0
)

linkxaxes!(axB1, axB2)

hideydecorations!(axB2)
hidexdecorations!(axB2, ticklabels = false)

lines!(
    1:10,
    [ 0.4, 0.2, 0.2, 0.2, 0.2, 0.3, 0.2, 0.2, 0.2, 0.2 ];
    color = :gray
)
scatter!(
    1:10,
    [ 0.4, 0.2, 0.2, 0.2, 0.2, 0.3, 0.2, 0.2, 0.2, 0.2 ];
    color = :gray
)

lines!(
    1:10,
    [ 0.4, 0.2, 0.2, 0.2, 0.2, 0.3, 0.2, 0.2, 0.2, 0.2 ];
    color = :gray
)
scatter!(
    1:10,
    [ 0.4, 0.2, 0.2, 0.2, 0.2, 0.3, 0.2, 0.2, 0.2, 0.2 ];
    color = :gray
)

labels = ["Correlation coefficient (r)", "Mean Absolute Error"]
elements = [PolyElement(polycolor = c) for c in [:green, :yellow]]
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

hidedecorations!(axC, grid = false)

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

### Plot panel D - Taxonomic profile ROC curves, 0-6 mo, color by Microbiome vs pure SES
#```julia

axD = Axis(
    CDEFGH_subfig[1,2];
    xlabel = "False-positive rate",
    ylabel = "True-positive rate",
    title = "D - ROC Curves, 0-6, TAXA"
)

for i in 1:classification_currentCogScores_00to06mo_onlyses.n_splits
    fprs, tprs, ts = res_roc_curve(classification_currentCogScores_00to06mo_onlyses, i; sample_set = "all")
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

singlemodel_merit_scatterplot!(axE, regression_currentCogScores_00to06mo_onlytaxa)

#```

### Plot panel F - Functional profile PCA, 0-6 mo, color by percentile
#```julia

axF = Axis(
    CDEFGH_subfig[2,1];
    xlabel = "PCo1",
    ylabel = "PCo2",
    title = "F - ECS Percentile PCA - 0 to 6 months"
)

hidedecorations!(axF, grid = false)

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

### Plot panel G - Functional profile ROC curves, 0-6 mo, color by Microbiome vs pure SES
#```julia

axG = Axis(
    CDEFGH_subfig[2,2];
    xlabel = "False-positive rate",
    ylabel = "True-positive rate",
    title = "G - ROC Curves, 0-6, ECS"
)

for i in 1:classification_currentCogScores_00to06mo_onlyses.n_splits
    fprs, tprs, ts = res_roc_curve(classification_currentCogScores_00to06mo_onlyses, i; sample_set = "all")
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

singlemodel_merit_scatterplot!(axH, regression_currentCogScores_00to06mo_onlyecs; split_index = 5)

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

hidedecorations!(axI, grid = false)

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

for i in 1:classification_currentCogScores_18to120mo_onlyses.n_splits
    fprs, tprs, ts = res_roc_curve(classification_currentCogScores_18to120mo_onlyses, i; sample_set = "all")
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

singlemodel_merit_scatterplot!(axK, regression_currentCogScores_00to06mo_onlytaxa)

#```

### Plot panel L - Functional profile PCA, 18-120mo, color by percentile
#```julia

axL = Axis(
    IJKLMN_subfig[2,1];
    xlabel = "PCo1",
    ylabel = "PCo2",
    title = "L - ECS Percentile PCA - 18-120 months"
)

hidedecorations!(axL, grid = false)

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

for i in 1:classification_currentCogScores_18to120mo_onlyses.n_splits
    fprs, tprs, ts = res_roc_curve(classification_currentCogScores_18to120mo_onlyses, i; sample_set = "all")
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

singlemodel_merit_scatterplot!(axN, regression_currentCogScores_00to06mo_onlyecs; split_index = 5)

figure
#```

save("figures/Figure3.png", figure)