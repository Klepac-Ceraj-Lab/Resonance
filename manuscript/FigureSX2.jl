#####
# Supplementary Figure: Heatmap from Figure 4B, but for the entire input set.
#####
using Chain
using CSV
using Statistics
using Random
using DataFrames
using StableRNGs
using MLJ
using CairoMakie
using DecisionTree
using JLD2
using Resonance
using Distributions
using Clustering, Distances
using ColorSchemes
ml_rng = StableRNG(0)

## Loading the pretrained models
RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree
brain_models = JLD2.load(modelfiles("manuscript","brain_models.jld"))["brain_models"]
include(joinpath("manuscript", "Figure4-definitions.jl"))

## Compute the mean figures of merit for the Brain segment regressions
mean_brain_merits = reduce(
    vcat,
    [ DataFrame(:variable => ordered_brain_segments_list[i], :Train_Cor => mean(brain_models[ordered_brain_segments_list[i]].merits.Train_Cor), :Test_Cor => mean(brain_models[ordered_brain_segments_list[i]].merits.Test_Cor)) for i in eachindex(ordered_brain_segments_list) ]
)

## Compute the mean importances loading of taxa on Brain segment regressions
weighted_brain_importances = reduce(
    (x, y) -> outerjoin(x, y, on = :variable, makeunique=true),
    [ rename!(weighted_hpimportances(j; normalize_importances=false), :weightedImportance => Symbol(split(j.name, '_')[2]))
      for (i, j) in brain_models ]
)
relative_brain_importances = reduce(
    (x, y) -> outerjoin(x, y, on = :variable, makeunique=true),
    [ rename!(weighted_hpimportances(j; normalize_importances=true), :weightedImportance => Symbol(split(j.name, '_')[2]))
      for (i, j) in brain_models ]
)

## Filter the merits and importances tables according to the list of interesting segments
interesting_segments_idxes = mean_brain_merits.variable .∈ Ref(interesting_segments)
mean_brain_merits = mean_brain_merits[interesting_segments_idxes, :]
relative_brain_importances = dropmissing(relative_brain_importances[:, vcat(["variable"], interesting_segments)]) # no row filter because the purpose of this is full-dataset

## Perform hierarchical clustering
transposedImportances = permutedims(relative_brain_importances, :variable)
insertcols!(transposedImportances, 1, :symmetricSegment => symmetric_segment_strings)

combinedTransposedImportances = DataFrames.combine(
    groupby(transposedImportances, :symmetricSegment),
    Symbol.(names(transposedImportances))[3:end] .=> mean .=> Symbol.(names(transposedImportances))[3:end]
)

dist_taxa = pairwise(Euclidean(), Matrix(combinedTransposedImportances[:, 2:end]); dims=2)
dist_segments = pairwise(Euclidean(), Matrix(combinedTransposedImportances[:, 2:end]); dims=1)
hcl_taxa = hclust(dist_taxa; linkage=:average, branchorder=:optimal)
hcl_segments = hclust(dist_segments; linkage=:average, branchorder=:optimal)
hclust_taxa_order = hcl_taxa.order
hclust_symmetric_segment_order = hcl_segments.order

## Reordering the tables according to the HClust analysis

reorder_segments_df = innerjoin(
    DataFrame(
        :original_segment => interesting_segments,
        :left_or_unique => vcat(repeat([1, 0], 24), repeat([ 1 ], 3)),
        :symmetric_segment => symmetric_segment_strings),
    DataFrame(
        :symmetric_segment => combinedTransposedImportances.symmetricSegment[hclust_symmetric_segment_order],
        :symmetric_order => hclust_symmetric_segment_order,
        :original_order => collect(1:length(combinedTransposedImportances.symmetricSegment[hclust_symmetric_segment_order]))),
    on = :symmetric_segment
)
insertcols!(reorder_segments_df, 1, :plot_order => collect(1:nrow(reorder_segments_df)))
sort!(reorder_segments_df, [:original_order, :left_or_unique])

plot_segments_order = reorder_segments_df.plot_order

#####
# Plotting
#####

## Initialize figure
figure = Figure(resolution = (1920, 1080))

mean_brain_merits = mean_brain_merits[plot_segments_order, :]
relative_brain_importances = relative_brain_importances[hclust_taxa_order, vcat( [ 1 ], (plot_segments_order .+ 1))]

## Importance Heatmap
ax = Axis(
    figure[1, 1];
    ylabel = "Target brain segment",
    xticks = (collect(1:nrow(relative_brain_importances)), relative_brain_importances.variable),
    yticks = (collect(1:length(interesting_segments)), names(relative_brain_importances)[2:end]),
    yticklabelsize=14,
    yreversed=true,
    title = "Brain segmentation data variable importances"
)
hidexdecorations!(ax; ticks=false)
ax_ticks = let
    bugs = relative_brain_importances.variable
    Axis(
        figure[1, 1];
        xlabel = "Predictor",
        xticks = (collect(1:nrow(relative_brain_importances)),
                replace.(bugs, "_"=>" ")),
        xticklabelsize=10,
        xticklabelrotation= pi/4,
        xticklabelfont="TeX Gyre Heros Makie Italic"
    )
end

tightlimits!.([ax, ax_ticks])
linkxaxes!(ax, ax_ticks)
hidexdecorations!(ax_ticks; ticklabels=false, label=false)
hideydecorations!(ax_ticks)
hidespines!(ax_ticks)

hm = CairoMakie.heatmap!(ax, Matrix(relative_brain_importances[:, 2:end]), yflip=true)
Colorbar(figure[1,2], hm; label= "Relative feature importance")

save(figurefiles("FigureSX2.svg"), figure)
save("manuscript/assets/FigureSX2.png", figure)