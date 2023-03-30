# Figure4 - RandomForest prediction of Brain Segments from DEMO + Taxonomic Profiles
```julia
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
```

## Loading the pretrained models

```julia
RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree
# concurrent brain regression from taxonomic profiles
brain_models = JLD2.load(modelfiles("brain_models.jld"))["brain_models"]
include("Figure4-definitions.jl")
```

### Plot figure 4 - FULL version without filtering segments or taxa, taxa ordered by mean importance, segments by hclust

#### Compute the mean figures of merit for the Brain segment regressions

```julia
mean_brain_merits = reduce(
    vcat,
    [ DataFrame(:variable => ordered_brain_segments_list[i], :Train_Cor => mean(brain_models[ordered_brain_segments_list[i]].merits.Train_Cor), :Test_Cor => mean(brain_models[ordered_brain_segments_list[i]].merits.Test_Cor)) for i in eachindex(ordered_brain_segments_list) ]
)
```

#### Compute the mean importances loading of taxa on Brain segment regressions

```julia
weighted_brain_importances = reduce(
    (x, y) -> outerjoin(x, y, on = :variable, makeunique=true),
    [ rename!(weighted_hpimportances(j; normalize_importances=false), :weightedImportance => Symbol(split(j.name, '_')[2]))
      for (i, j) in brain_models ] )

relative_brain_importances = reduce(
    (x, y) -> outerjoin(x, y, on = :variable, makeunique=true),
    [ rename!(weighted_hpimportances(j; normalize_importances=true), :weightedImportance => Symbol(split(j.name, '_')[2]))
      for (i, j) in brain_models ] )
```

#### Filter the merits and importances tables according to the list of interesting segments and taxa

```julia
interesting_segments_idxes = mean_brain_merits.variable .∈ Ref(interesting_segments)
mean_brain_merits = mean_brain_merits[interesting_segments_idxes, :]

interesting_taxa_idxes = relative_brain_importances.variable .∈ Ref(interesting_taxa)
relative_brain_importances = dropmissing(relative_brain_importances[interesting_taxa_idxes, vcat(["variable"], interesting_segments) ])
```

#### Perform hierarchical clustering

```julia
## Transpose the matrix to add information about the symmetric segments
transposedImportances = permutedims(relative_brain_importances, :variable)
insertcols!(transposedImportances, 1, :symmetricSegment => symmetric_segment_strings)
## Combine left and right of symmetric segments
combinedTransposedImportances = DataFrames.combine(
    groupby(transposedImportances, :symmetricSegment),
    Symbol.(names(transposedImportances))[3:end] .=> mean .=> Symbol.(names(transposedImportances))[3:end]
)
## Perform the actual HCA and store the order
dist_taxa = pairwise(Euclidean(), Matrix(combinedTransposedImportances[:, 2:end]); dims=2)
dist_segments = pairwise(Euclidean(), Matrix(combinedTransposedImportances[:, 2:end]); dims=1)
hcl_taxa = hclust(dist_taxa; linkage=:average, branchorder=:optimal)
hcl_segments = hclust(dist_segments; linkage=:average, branchorder=:optimal)
hclust_taxa_order = hcl_taxa.order
hclust_symmetric_segment_order = hcl_segments.order
```
#### Reordering the tables according to the HClust analysis

```julia
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
```

#### using the hclust to reorder all we need

```julia
mean_brain_merits = mean_brain_merits[plot_segments_order, :]
relative_brain_importances = relative_brain_importances[hclust_taxa_order, vcat( [ 1 ], (plot_segments_order .+ 1))]
#relative_brain_importances = relative_brain_importances[sortperm(map(mean, eachrow(relative_brain_importances[:, 2:end])); rev=true), vcat( [ 1 ], (plot_segments_order .+ 1))]
```

#### Some other subsets

```julia
weighted_brain_importances = reduce(
    (x, y) -> outerjoin(x, y, on = :variable, makeunique=true),
        [ rename!(weighted_hpimportances(j; normalize_importances=true), 
                :weightedImportance => Symbol(split(j.name, '_')[2])) for (i, j) in brain_models ]
)

weighted_brain_importances = dropmissing(weighted_brain_importances[:, vcat(["variable"], interesting_segments) ])
weighted_noage_importances = weighted_brain_importances[2:end, :]
```

### Plotting

#### Initialize figure

```julia
figure = Figure(resolution = (1920, 1536))


AB_Subfig = GridLayout(figure[1,1], alignmode=Outside()) 
# A_subfig = GridLayout(AB_Subfig[1,1])
# B_subfig = GridLayout(AB_Subfig[1,2])
C_subfig = GridLayout(figure[2,1], alignmode=Outside())
```

```julia
highlight_bugs = [
    "Anaerostipes_hadrus",
    "Bacteroides_vulgatus",
    "Blautia_wexlerae",
    "Eubacterium_eligens",
    "Coprococcus_comes",
    "Ruminococcus_torques",
    "Adlercreutzia_equolifaciens",
    "Asaccharobacter_celatus",
    "Bacteroides_ovatus",
    "Coprococcus_eutactus"
]

hlbugs_color = Dict(b=> c for (c,b) in zip(ColorSchemes.Set3_12, highlight_bugs))
```
#### Calling the plot functions with the ordered data

```julia
axA = Axis(
    AB_Subfig[1, 1];
    xlabel = "Correlation",
    yticks = (reverse(collect(1:length(interesting_segments))), mean_brain_merits.variable),
    ylabel = "Target Variable",
    xticks = [-0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
    title = "Mean Random Forest regression correlations",
    yticklabelsize=16,
    alignmode=Inside(),
    #yticklabelrotation= -pi/2
)

tightlimits!(axA, Top())
tightlimits!(axA, Bottom())

barplot!(
    axA,
    reverse(collect(1:length(interesting_segments))),
    mean_brain_merits.Test_Cor,
    color = vcat( repeat( [ "blue", "red" ], 24),repeat( [ "purple" ], 3) )[plot_segments_order],
    direction=:x
)

Legend(AB_Subfig[2,1],
    [MarkerElement(; marker=:rect, color=c) for c in ("blue", "red", "purple")],
    ["Left hemisphere", "Right hemisphere", "non-lateral"];
)
```

#### Importance Heatmaps

```julia
nbugs_toplot = length(interesting_taxa) # the top 1/3 bugs (out of 129) From intersecting the top mean importances and top max importances

axB = Axis(
    AB_Subfig[1, 2];
    ylabel = "Target brain segment",
    xticks = (collect(1:nbugs_toplot), relative_brain_importances.variable[1:nbugs_toplot]),
    yticks = (collect(1:length(interesting_segments)), names(relative_brain_importances)[2:end]),
    yticklabelsize=16,
    yreversed=true,
    title = "Brain segmentation data variable importances"
)
hidexdecorations!(axB; ticks=false)
axBticks = let
    bugs = relative_brain_importances.variable[1:nbugs_toplot]
    Axis(
        AB_Subfig[2,2];
        xlabel = "Predictor",
        xticks = (collect(1:nbugs_toplot),
                replace.(bugs, "_"=>" ")),
        xticklabelsize=16,
        xticklabelrotation= pi/4,
        xticklabelfont="TeX Gyre Heros Makie Italic",
        alignmode=Outside()
    )
end
tightlimits!.([axB, axBticks])

linkxaxes!(axB, axBticks)
hidexdecorations!(axBticks; ticklabels=false, label=false)
hideydecorations!(axBticks)
hidespines!(axBticks)
hm = CairoMakie.heatmap!(axB, Matrix(relative_brain_importances[1:nbugs_toplot, 2:end]), yflip=true)
```

```julia
let
    bugs = relative_brain_importances.variable[1:nbugs_toplot]

    for (i, bug) in enumerate(bugs)
        if bug in highlight_bugs
            poly!(axBticks, Point2f[(i-0.5, 0), (i-0.5, 1), (i+0.5, 1), (i+0.5, 0)];
                color=hlbugs_color[bug])
        end
    end
end
Colorbar(AB_Subfig[1,3], hm; label= "Relative feature importance")
```


```julia
idx = sortperm([median(row[2:end]) for row in eachrow(weighted_noage_importances)])
ys = reduce(vcat, [values(weighted_noage_importances[i, 2:end])...] for i in idx)
xs = repeat(1:length(idx); inner=ncol(weighted_noage_importances)-1)

axC = Axis(
    C_subfig[1,1];
    xlabel="Bugs (rank median)",
    ylabel = "importances",
    title = "Importances on all brain segments, by taxa")

CairoMakie.xlims!(axC, [0, nrow(weighted_noage_importances)+1])

maximp = maximum(Matrix(weighted_noage_importances[!, 2:end]))
   
for bug in highlight_bugs
    i = findfirst(==(bug), weighted_noage_importances.variable[idx])
    poly!(axC, Point2f[(i-0.5, 0), (i+0.5, 0), (i+0.5, maximp), (i-0.5, maximp)]; color=hlbugs_color[bug])
end
CairoMakie.scatter!(axC, xs .+ rand(Normal(0, 0.1), length(xs)), ys;)
Legend(C_subfig[1,2], [MarkerElement(; marker=:rect, color = hlbugs_color[bug]) for bug in highlight_bugs],
        replace.(highlight_bugs, "_"=> " ");
        labelfont="TeX Gyre Heros Makie Italic")

Label(AB_Subfig[1, 1, TopLeft()], "A", fontsize = 26,font = :bold, padding = (0, 240, 5, 0), halign = :right)
Label(AB_Subfig[1, 2, TopLeft()], "B", fontsize = 26,font = :bold, padding = (0, 200, 5, 0), halign = :right)
Label(C_subfig[1, 1, TopLeft()], "C", fontsize = 26,font = :bold, padding = (0, 40, 5, 0), halign = :right)
```

Then resize

```julia
colsize!(AB_Subfig, 1, Relative(0.2))
colsize!(AB_Subfig, 2, Relative(0.8))
rowsize!(AB_Subfig, 2, Relative(0.22))
axA.alignmode=Mixed(; bottom=-60)
axB.alignmode=Mixed(; bottom=-20)
rowsize!(figure.layout, 1, Relative(0.75))
rowsize!(figure.layout, 2, Relative(0.25))
figure
```

```julia
save(figurefiles("Figure4.svg"), figure)
save("manuscript/assets/Figure4.png", figure)
```

## Supplementary Figures

### Supplementary Figure 5

```julia
relative_brain_importances = reduce(
    (x, y) -> outerjoin(x, y, on = :variable, makeunique=true),
    [ rename!(weighted_hpimportances(j; normalize_importances=true), :weightedImportance => Symbol(split(j.name, '_')[2]))
      for (i, j) in brain_models ]
)

## Filter the merits and importances tables according to the list of interesting segments
interesting_segments_idxes = mean_brain_merits.variable .∈ Ref(interesting_segments)
relative_brain_importances = dropmissing(relative_brain_importances[:, vcat(["variable"], interesting_segments)]) # no row filter because the purpose of this is full-dataset

## Perform hierarchical clustering
transposedImportances = permutedims(relative_brain_importances, :variable)
dist_taxa = pairwise(Euclidean(), Matrix(transposedImportances[:, 2:end]); dims=2)
dist_segments = pairwise(Euclidean(), Matrix(transposedImportances[:, 2:end]); dims=1)
hcl_taxa = hclust(dist_taxa; linkage=:average, branchorder=:optimal)
hcl_segments = hclust(dist_segments; linkage=:average, branchorder=:optimal)
hclust_taxa_order = hcl_taxa.order
hclust_segments_order = hcl_segments.order

## Initialize figure
figure = Figure(resolution = (1920, 1080))
relative_brain_importances = relative_brain_importances[hclust_taxa_order, vcat( [ 1 ], (hclust_segments_order .+ 1))]

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

hm = CairoMakie.heatmap!(ax, Matrix(relative_brain_importances[:, 2:end]), yflip=true, colorrange = (0.0, 0.04), highclip = :yellow)
Colorbar(figure[1,2], hm; label= "Relative feature importance")

## Saving
save(figurefiles("Supp_Figure5.svg"), figure)
save("manuscript/assets/Supp_Figure5.png", figure)
```

### Supp Figure 6

```julia
let fig = Figure(; resolution=(2000, 500))
    importances = reduce(
        (x, y) -> outerjoin(x, y, on = :variable, makeunique=true),
            [ rename!(weighted_hpimportances(j; normalize_importances=true), 
                    :weightedImportance => Symbol(split(j.name, '_')[2])) for (i, j) in brain_models ]
    )
    subset!(importances, "variable" => ByRow(v-> v != "ageMonths"))

    idx = sortperm([median(row[2:end]) for row in eachrow(importances)])
    ys = reduce(vcat, [values(importances[i, 2:end])...] for i in idx)
    xs = repeat(1:length(idx); inner=ncol(importances)-1)
    @info length(ys) length(xs)
    
    ax = Axis(fig[1,1]; 
        xlabel="Bugs (rank median)",
        ylabel = "importances",
        xticks = (1:nrow(importances), replace.(importances[idx, "variable"], "_"=>"")),
        xticklabelrotation = π/4,
        xticklabelfont="TeX Gyre Heros Makie Italic",
        xticklabelsize=10

    )
    scatter!(ax, xs .+ rand(Normal(0, 0.1), length(xs)), ys;)
    xlims!(ax, 0, nrow(importances)+1)
    save(figurefiles("Supp_Figure6.svg"), fig)
    save("manuscript/assets/Supp_Figure6.png", fig)
    fig
end

```



### Supplementary Figure 7

```julia
relative_brain_importances = reduce(
    (x, y) -> outerjoin(x, y, on = :variable, makeunique=true),
    [ rename!(weighted_hpimportances(j; normalize_importances=true), :weightedImportance => Symbol(split(j.name, '_')[2]))
      for (i, j) in brain_models ]
)

interesting_segments_idxes = mean_brain_merits.variable .∈ Ref(interesting_segments)
relative_brain_importances = dropmissing(relative_brain_importances[:, vcat(["variable"], interesting_segments)]) # no row filter 

supptbl_like_braincombine = DataFrame(
    :variable => relative_brain_importances.variable,
    :relativeWeightedImportance => map(mean, eachrow(Matrix(relative_brain_importances[:, 2:end])))
)
sort!(supptbl_like_braincombine, :relativeWeightedImportance; rev = true)
supptbl_like_braincombine.cumulativeWeightedImportance = cumsum(supptbl_like_braincombine.relativeWeightedImportance)

## Initialize figure and plot Pareto
figure = Figure(resolution = (1920, 1080))
this_barcolor = :lightblue
this_curvecolor = :orange
plot_importances_pareto!(figure[1,1], supptbl_like_braincombine, "Pareto Plot - average over all regions"; barcolor = this_barcolor, curvecolor = this_curvecolor)

Legend(
    figure[2, 1],
    [ PolyElement(; color = this_barcolor, strokewidth = 1, strokecolor = :black), LineElement(; color = this_curvecolor, linewidth = 5)],
    [ "Individual relative Importance", "Cumulative relative importance" ],
    tellheight = true, tellwidth = false,
    margin = (10, 10, 10, 10),
    orientation = :horizontal
)

## Saving
save(figurefiles("Supp_Figure7.svg"), figure)
save("manuscript/assets/Supp_Figure7.png", figure)
```