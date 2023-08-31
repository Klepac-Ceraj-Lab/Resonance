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

RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree

include("manuscript/Figure4-definitions.jl")
```


## Plot figure 4 panels A and B

#### Loading the pretrained models

```julia
# Composite cognitive scores and Mullen subscales
regression_currentCogScores_18to120mo_onlydemo = JLD2.load(
    modelfiles("regression_currentCogScores_18to120mo_onlydemo.jld")
    )["regression_currentCogScores_18to120mo_onlydemo"];
regression_currentCogScores_18to120mo_demoplustaxa = JLD2.load(
    modelfiles("regression_currentCogScores_18to120mo_demoplustaxa.jld")
    )["regression_currentCogScores_18to120mo_demoplustaxa"];
regression_currentCogScores_18to120mo_demoplusecs = JLD2.load(
    modelfiles("regression_currentCogScores_18to120mo_demoplusecs.jld")
    )["regression_currentCogScores_18to120mo_demoplusecs"];
regression_currentExpressiveLanguages_18to120mo_onlydemo = JLD2.load(
    modelfiles("regression_currentExpressiveLanguages_18to120mo_onlydemo.jld")
    )["regression_currentExpressiveLanguages_18to120mo_onlydemo"];
regression_currentExpressiveLanguages_18to120mo_demoplustaxa = JLD2.load(
    modelfiles("regression_currentExpressiveLanguages_18to120mo_demoplustaxa.jld")
    )["regression_currentExpressiveLanguages_18to120mo_demoplustaxa"];
regression_currentExpressiveLanguages_18to120mo_demoplusecs = JLD2.load(
    modelfiles("regression_currentExpressiveLanguages_18to120mo_demoplusecs.jld")
    )["regression_currentExpressiveLanguages_18to120mo_demoplusecs"];
regression_currentGrossMotors_18to120mo_onlydemo = JLD2.load(
    modelfiles("regression_currentGrossMotors_18to120mo_onlydemo.jld")
    )["regression_currentGrossMotors_18to120mo_onlydemo"];
regression_currentGrossMotors_18to120mo_demoplustaxa = JLD2.load(
    modelfiles("regression_currentGrossMotors_18to120mo_demoplustaxa.jld")
    )["regression_currentGrossMotors_18to120mo_demoplustaxa"];
regression_currentGrossMotors_18to120mo_demoplusecs = JLD2.load(
    modelfiles("regression_currentGrossMotors_18to120mo_demoplusecs.jld")
    )["regression_currentGrossMotors_18to120mo_demoplusecs"];
regression_currentVisualReceptions_18to120mo_onlydemo = JLD2.load(
    modelfiles("regression_currentVisualReceptions_18to120mo_onlydemo.jld")
    )["regression_currentVisualReceptions_18to120mo_onlydemo"];
regression_currentVisualReceptions_18to120mo_demoplustaxa = JLD2.load(
    modelfiles("regression_currentVisualReceptions_18to120mo_demoplustaxa.jld")
    )["regression_currentVisualReceptions_18to120mo_demoplustaxa"];
regression_currentVisualReceptions_18to120mo_demoplusecs = JLD2.load(
    modelfiles("regression_currentVisualReceptions_18to120mo_demoplusecs.jld")
    )["regression_currentVisualReceptions_18to120mo_demoplusecs"];
# concurrent brain regression from taxonomic profiles
brain_models = JLD2.load(modelfiles("brain_models.jld"))["brain_models"]
```

#### Compute the mean figures of merit for the Brain segment regressions

```julia
mean_brain_merits = reduce(
    vcat,
    [ DataFrame(:variable => ordered_brain_segments_list[i], :Train_Cor => mean(brain_models[ordered_brain_segments_list[i]].merits.Train_Cor), :Test_Cor => mean(brain_models[ordered_brain_segments_list[i]].merits.Test_Cor)) for i in eachindex(ordered_brain_segments_list) ]
)
```

#### Compute the mean importances loading of taxa on Brain segment regressions

```julia
weighted_brain_importances = dropmissing(reduce(
    (x, y) -> outerjoin(x, y, on = :variable, makeunique=true),
    [ rename!(weighted_hpimportances(j; normalize_importances=false), :weightedImportance => Symbol(split(j.name, '_')[2]))
      for (i, j) in brain_models ] ) )

relative_brain_importances = dropmissing(reduce(
    (x, y) -> outerjoin(x, y, on = :variable, makeunique=true),
    [ rename!(weighted_hpimportances(j; normalize_importances=true), :weightedImportance => Symbol(split(j.name, '_')[2]))
      for (i, j) in brain_models ] ))
```


#### Getting a table with figures of merit to plot on Panel A

```julia
function generate_cols(model_df, xs, grp, color)
    insertcols(DataFrames.combine(
        groupby(model_df.merits, :Hyperpar_Idx),
            "Test_Cor" => (x -> round(mean(x); digits=2))=> "Test_Cor_mean"), 1,
            "xs"=> xs, "grp"=> grp, "color"=> color
    )

end 
cm = [ColorSchemes.Spectral_10[i] for i in (10,2,4,8)]
panelA_plot_df = let
    vcat(
        generate_cols(regression_currentCogScores_18to120mo_onlydemo,               1,  1, cm[1] ),
        generate_cols(regression_currentCogScores_18to120mo_demoplustaxa,           2,  1, cm[1] ),
        generate_cols(regression_currentCogScores_18to120mo_demoplusecs,            3,  1, cm[1] ),
        generate_cols(regression_currentExpressiveLanguages_18to120mo_onlydemo,     1,  2, cm[2] ),
        generate_cols(regression_currentExpressiveLanguages_18to120mo_demoplustaxa, 2,  2, cm[2] ),
        generate_cols(regression_currentExpressiveLanguages_18to120mo_demoplusecs,  3,  2, cm[2] ),
        generate_cols(regression_currentGrossMotors_18to120mo_onlydemo,             1,  3, cm[3] ),
        generate_cols(regression_currentGrossMotors_18to120mo_demoplustaxa,         2,  3, cm[3] ),
        generate_cols(regression_currentGrossMotors_18to120mo_demoplusecs,          3,  3, cm[3] ),
        generate_cols(regression_currentVisualReceptions_18to120mo_onlydemo,        1,  4, cm[4] ),
        generate_cols(regression_currentVisualReceptions_18to120mo_demoplustaxa,    2,  4, cm[4] ),
        generate_cols(regression_currentVisualReceptions_18to120mo_demoplusecs,     3,  4, cm[4] ),
    )
end
```

#### Getting the individual and combined importances, finding the bugs to plot on Panel B
```julia
composite_importances = rename(
    weighted_hpimportances(regression_currentCogScores_18to120mo_demoplustaxa;
        normalize_importances=true),
        :weightedImportance => :Composite
)
el_importances = rename(
    weighted_hpimportances(regression_currentExpressiveLanguages_18to120mo_demoplustaxa;
        normalize_importances=true),
        :weightedImportance => :ExpressiveLanguage
)
gm_importances = rename(
    weighted_hpimportances(regression_currentGrossMotors_18to120mo_demoplustaxa;
        normalize_importances=true),
        :weightedImportance => :GrossMotor
)
vr_importances = rename(
    weighted_hpimportances(regression_currentVisualReceptions_18to120mo_demoplustaxa;
        normalize_importances=true),
        :weightedImportance => :VisualReception
)

sidebyside_importances = reduce( (x, y) -> outerjoin(x, y, on = :variable; makeunique = true), [ composite_importances, el_importances, gm_importances, vr_importances ])
sidebyside_importances.avg = [ mean(skipmissing(collect(rr))) for rr in eachrow(sidebyside_importances[:, 2:end]) ]
important_bugs_composite = sort(dropmissing(sidebyside_importances, :Composite), :Composite; rev = true).variable[1:12]
important_bugs_ExpressiveLanguage = sort(dropmissing(sidebyside_importances, :ExpressiveLanguage), :ExpressiveLanguage; rev = true).variable[1:12]
important_bugs_GrossMotor = sort(dropmissing(sidebyside_importances, :GrossMotor), :GrossMotor; rev = true).variable[1:12]
important_bugs_VisualReception = sort(dropmissing(sidebyside_importances, :VisualReception), :VisualReception; rev = true).variable[1:12]

@show union(important_bugs_composite, important_bugs_ExpressiveLanguage, important_bugs_GrossMotor, important_bugs_VisualReception)

panelB_taxa = [
    "Faecalibacterium_prausnitzii",
    "Blautia_wexlerae",
    "Eubacterium_eligens",
    "Bifidobacterium_pseudocatenulatum",
    "Bifidobacterium_longum",
    "Parasutterella_excrementihominis",
    "Veillonella_dispar",
    "Fusicatenibacter_saccharivorans",
    "Ruminococcus_gnavus",
    "Roseburia_inulinivorans",
    "Flavonifractor_plautii",
    "Roseburia_faecis",
    "Streptococcus_salivarius",
    "Clostridium_symbiosum",
    "Bacteroides_uniformis",
    "Clostridium_innocuum",
    "Bacteroides_vulgatus",
    "Parabacteroides_merdae",
    "Intestinibacter_bartlettii",
]

panelB_taxa_idxer = Dict( [ panelB_taxa[i] => i for i in eachindex(panelB_taxa) ] )
```

#### Preparing data to plot on panel B
```julia
panelB_plot_df = @chain vcat(
    insertcols(
        weighted_hpimportances(
            regression_currentCogScores_18to120mo_demoplustaxa;
            normalize_importances=true
        ), 1, :grp => 1, :color => cm[1]),
    insertcols(
        weighted_hpimportances(
            regression_currentExpressiveLanguages_18to120mo_demoplustaxa;
            normalize_importances=true
        ), 1, :grp => 2, :color => cm[2]),
    insertcols(
        weighted_hpimportances(
            regression_currentGrossMotors_18to120mo_demoplustaxa;
            normalize_importances=true
        ), 1, :grp => 3, :color => cm[3]),
    insertcols(
        weighted_hpimportances(
            regression_currentVisualReceptions_18to120mo_demoplustaxa;
            normalize_importances=true
        ), 1, :grp => 4, :color => cm[4])
) begin
    subset(:variable => ( x -> x .∈ Ref(panelB_taxa) ))
    transform!(:variable => (x -> [ panelB_taxa_idxer[el] for el in x ]) => :xs; renamecols = false)
end
```
#### Filter the merits tables according to the list of interesting segments

```julia
interesting_segments_idxes = mean_brain_merits.variable .∈ Ref(interesting_segments)
mean_brain_merits = mean_brain_merits[interesting_segments_idxes, :]
```

#### Using the filtered merits to explore the importances
```julia
explore_importances_df = dropmissing(relative_brain_importances[:, vcat(["variable"], mean_brain_merits.variable) ])
rename!(explore_importances_df, :variable => :predictor)
explore_importances_df = stack(explore_importances_df, 2:ncol(explore_importances_df))
subset!(explore_importances_df, :predictor => (x -> x .!= "ageMonths"))
sort!(explore_importances_df, :value; rev = true)
@show unique(explore_importances_df.predictor)[1:20]
```

#### Filter the importances tables according to the list of interesting segments and taxa
```julia
interesting_taxa_idxes = relative_brain_importances.variable .∈ Ref(interesting_taxa)
relative_brain_importances = dropmissing(relative_brain_importances[interesting_taxa_idxes, vcat(["variable"], mean_brain_merits.variable) ])
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
        :left_or_unique => repeat([1, 0], 16),
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

### Declaring the Figure and its subdivisions

```julia
figure = Figure(resolution = (1920, 1536));

AB_Subfig = GridLayout(figure[1,1], alignmode=Outside()) 
A_Subfig = GridLayout(AB_Subfig[1,1])
B_Subfig = GridLayout(AB_Subfig[1,2])
CD_Subfig = GridLayout(figure[2,1], alignmode=Outside())
E_Subfig  = GridLayout(figure[1:2,2], alignmode=Outside())

```

#### Plotting Panel A
```julia
axA = Axis(
    A_Subfig[1, 1];
    xlabel = "Correlation",
    yticks = (collect(1:3), ["Only demographics", "demographics and\ntaxonomic profiles", "demographics and\nfunctional profiles (ECs)"]),
    ylabel = "Model input composition",
    title = "Average correlation of composite and subscale predictor model",
    yticklabelsize=16,
    alignmode=Outside(),
)

barplot!(axA,
    panelA_plot_df.xs,
    panelA_plot_df.Test_Cor_mean,
    dodge = panelA_plot_df.grp,
    color = panelA_plot_df.color,
    direction = :x
)
```
```julia
highlight_bugs = [
    "Anaerostipes_hadrus",
    "Bacteroides_vulgatus",
    "Blautia_wexlerae",
    "Eubacterium_eligens",
    "Coprococcus_comes",
    "Ruminococcus_torques",
    "Bifidobacterium_pseudocatenulatum",
    "Asaccharobacter_celatus",
    "Bacteroides_ovatus",
    "Coprococcus_eutactus"
]

hlbugs_color = Dict(b=> c for (c,b) in zip(ColorSchemes.Set3_12, highlight_bugs))
@assert length(interesting_segments) % 2 == 0
@assert all(1:2:length(interesting_segments)) do i
    m1 = match(r"^(left)-(.+)$", interesting_segments[i])
    m2 = match(r"^(right)-(.+)$", interesting_segments[i+1])
    any(isnothing, (m1,m2)) && return false
    m1[2] == m2[2]
end
```


#### Plotting Panel B
```julia
axB = Axis(
    B_Subfig[1, 1];
    ylabel = "Relative weighted Importance",
    xticks = (collect(eachindex(panelB_taxa)), replace.(panelB_taxa, "_"=>" ")),
    xticklabelfont="TeX Gyre Heros Makie Italic",
    xticklabelsize=16,
    xticklabelrotation= pi/4,
    yreversed=false,
    alignmode = Outside(),
    title = "Importance of selected taxa on composite and subscores"
)

hideydecorations!(axB)
hidexdecorations!(axB; ticks=false)
axBticks = let
    bugs = panelB_taxa
    Axis(
        B_Subfig[2,1];
        xlabel = "Predictor",
        xticks = (collect(1:length(bugs)),
                format_species_labels(bugs)),
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

barplot!(axB,
    panelB_plot_df.xs,
    panelB_plot_df.weightedImportance,
    dodge = panelB_plot_df.grp,
    color = panelB_plot_df.color
    )

Legend(A_Subfig[2,1],
    [
        MarkerElement(; marker=:rect, color=cm[1]),
        MarkerElement(; marker=:rect, color=cm[2]),
        MarkerElement(; marker=:rect, color=cm[3]),
        MarkerElement(; marker=:rect, color=cm[4])
    ],
    [ "Composite Score", "Expressive Language", "Gross Motor", "Visual Reception"];
    orientation=:horizontal,
    nbanks=2
)
let
    bugs = panelB_taxa

    for (i, bug) in enumerate(bugs)
        if bug in highlight_bugs
            poly!(axBticks, Point2f[(i-0.5, 0), (i-0.5, 1), (i+0.5, 1), (i+0.5, 0)];
                color=hlbugs_color[bug])
        end
    end
end
```

### Plot figure 4 panels C, d and E

### Plotting

#### Calling the plot functions with the ordered data

```julia
axC = Axis(
    CD_Subfig[1, 1];
    xlabel = "Correlation",
    yticks = (reverse(collect(1.5:2:length(interesting_segments)+0.5)),
              replace.(mean_brain_merits.variable[1:2:end], r"(right|left)-"=> "")),
    ylabel = "Target Variable",
    xticks = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
    title = "Mean RF correlations",
    yticklabelsize=16,
    alignmode=Inside(),
    #yticklabelrotation= -pi/2
)

tightlimits!(axC, Top())
tightlimits!(axC, Bottom())

barplot!(
    axC,
    reverse(collect(1:length(interesting_segments))),
    mean_brain_merits.Test_Cor,
    color = repeat( [ "blue", "red" ], 16)[plot_segments_order],
    direction=:x
)

Legend(CD_Subfig[2,1],
    [MarkerElement(; marker=:rect, color=c) for c in ("blue", "red")],
    ["Left hemisphere", "Right hemisphere"];
)
```

#### Importance Heatmaps

```julia
nbugs_toplot = length(interesting_taxa) # the top 1/3 bugs (out of 129) From intersecting the top mean importances and top max importances

axD = Axis(
    CD_Subfig[1, 2];
    xticks = (collect(1:nbugs_toplot), relative_brain_importances.variable[1:nbugs_toplot]),
    yreversed=true,
    title = "Brain segmentation data variable importances"
)
hideydecorations!(axD)
hidexdecorations!(axD; ticks=false)
axDticks = let
    bugs = relative_brain_importances.variable[1:nbugs_toplot]
    Axis(
        CD_Subfig[2,2];
        xlabel = "Predictor",
        xticks = (collect(1:nbugs_toplot),
                replace.(bugs, "_"=>" ")),
        xticklabelsize=16,
        xticklabelrotation= pi/4,
        xticklabelfont="TeX Gyre Heros Makie Italic",
        alignmode=Outside()
    )
end
tightlimits!.([axD, axDticks])

linkxaxes!(axD, axDticks)
hidexdecorations!(axDticks; ticklabels=false, label=false)
hideydecorations!(axDticks)
hidespines!(axDticks)
hm = CairoMakie.heatmap!(axD, Matrix(relative_brain_importances[1:nbugs_toplot, 2:end]), yflip=true)
```

```julia
let
    bugs = relative_brain_importances.variable[1:nbugs_toplot]

    for (i, bug) in enumerate(bugs)
        if bug in highlight_bugs
            poly!(axDticks, Point2f[(i-0.5, 0), (i-0.5, 1), (i+0.5, 1), (i+0.5, 0)];
                color=hlbugs_color[bug])
        end
    end
end
Colorbar(CD_Subfig[1,3], hm; label= "Relative feature importance", ticks=0:0.01:0.04, minorticksvisible=true)
```


```julia
idx = sortperm([median(row[2:end]) for row in eachrow(weighted_noage_importances)])
xs = reduce(vcat, [values(weighted_noage_importances[i, 2:end])...] for i in idx)
ys = repeat(1:length(idx); inner=ncol(weighted_noage_importances)-1)

axE = Axis(
    E_Subfig[1,1];
    ylabel="Bugs (rank median)",
    xlabel = "importances",
)
minidx = 75
CairoMakie.ylims!(axE, [minidx, nrow(weighted_noage_importances)+1])

maximp = maximum(Matrix(weighted_noage_importances[!, 2:end]))
   
for bug in highlight_bugs
    i = findfirst(==(bug), weighted_noage_importances.variable[idx])
    poly!(axE, Point2f[(0, i-0.5), (0, i+0.5), (maximp, i+0.5), (maximp, i-0.5)]; color=hlbugs_color[bug])
end
CairoMakie.scatter!(axE, xs, ys .+ rand(Normal(0, 0.1), length(xs));)
Legend(E_Subfig[1,2], [MarkerElement(; marker=:rect, color = hlbugs_color[bug]) for bug in highlight_bugs],
        replace.(highlight_bugs, "_"=> " ");
        labelfont="TeX Gyre Heros Makie Italic")

```

### Final labeling and layout resizing

```julia
Label(A_Subfig[1, 1, TopLeft()], "A", fontsize = 26,font = :bold, halign = :right)
Label(B_Subfig[1, 1, TopLeft()], "B", fontsize = 26,font = :bold, padding = (0, 10, 5, 0), halign = :right)
Label(CD_Subfig[1, 1, TopLeft()], "C", fontsize = 26,font = :bold, padding = (0, 10, 5, 0), halign = :right)
Label(CD_Subfig[1, 2, TopLeft()], "D", fontsize = 26,font = :bold, padding = (0, 10, 5, 0), halign = :right)
Label(E_Subfig[1, 1, TopLeft()], "E", fontsize = 26,font = :bold, padding = (0, 40, 5, 0), halign = :right)

#colsize!(AB_Subfig, 1, Relative(0.4))
#colsize!(AB_Subfig, 2, Relative(0.6))
#colsize!(CD_Subfig, 1, Relative(0.2))
#colsize!(CD_Subfig, 2, Relative(0.8))
#rowsize!(CD_Subfig, 2, Relative(0.40))
#axC.alignmode=Mixed(; bottom=-60)
#axD.alignmode=Mixed(; bottom=-20)
#rowsize!(figure.layout, 1, Relative(0.35))
#rowsize!(figure.layout, 2, Relative(0.45))
#rowsize!(figure.layout, 3, Relative(0.2))
#colgap!(CD_Subfig, Fixed(5))
figure
```

```julia
# save(figurefiles("Figure4.svg"), figure)
save(figurefiles("Figure4.png"), figure)
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

hm = CairoMakie.heatmap!(ax, Matrix(relative_brain_importances[:, 2:end]), yflip=true, colorrange = (0.0, 0.04), highclip = :white)
Colorbar(figure[1,2], hm; label= "Relative feature importance", ticks=0:0.01:0.04, minorticksvisible=true)

## Saving
save(figurefiles("Supp_Figure5.svg"), figure)
save("manuscript/assets/Supp_Figure5.png", figure)
figure
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