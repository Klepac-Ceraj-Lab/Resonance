# Figure3 - RandomForest prediction of cogScores from DEM, Taxonomic and Functional Profiles

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
ml_rng = StableRNG(0)
```

## Defining relevant functions

### Calculating fitness-weighted importances
```julia
function calculate_fitness(train_cor::Vector{Float64}, test_cor::Vector{Float64})

    positive_train_cor = map( x -> maximum([x, 0.0]), train_cor)
    positive_test_cor = map( x -> maximum([x, 0.0]), test_cor)

    return positive_train_cor .* positive_test_cor

end

function weightedimportances(m; change_hashnames=false, hashnamestable::Union{Nothing, DataFrame} = nothing)

    merits = m.merits
    importances = m.importances
    fitnesses = calculate_fitness(merits.Train_Cor, merits.Test_Cor)
    weighted_importances = map(x -> sum(x .* fitnesses)/sum(fitnesses .> 0.0), collect(eachrow(Matrix(importances[:, 2:end]))))
    weighted_importances_df = sort( DataFrame(:variable => importances[:,1], :weightedImportance => weighted_importances), :weightedImportance; rev=true)
    sum(fitnesses .> 0.0) / nrow(merits)
    contributing_rows = merits[fitnesses .> 0.0, :]
    sum(contributing_rows.Train_Cor .> contributing_rows.Test_Cor)
    sum(contributing_rows.Train_Cor .< contributing_rows.Test_Cor)

    if change_hashnames
        newnames = leftjoin(weighted_importances_df, hashnamestable, on = :variable => :hashname).longname
        weighted_importances_df.variable = newnames
    end

    return weighted_importances_df
end
```

### Comparative Importance plots
```julia
function compute_joined_importances(modelwithoutdemo, modelwithdemo; imp_fun = hpimportances)

    importanceswithoutdemo = imp_fun(modelwithoutdemo)
    importanceswithdemo = imp_fun(modelwithdemo)
    rename!(importanceswithoutdemo, "weightedImportance" => "ImportanceWithoutDemo")
    rename!(importanceswithdemo, "weightedImportance" => "ImportanceWithDemo")

    joined_importances = dropmissing(outerjoin(importanceswithdemo, importanceswithoutdemo, on = [ :variable => :variable ]))
    sort!(joined_importances, :ImportanceWithoutDemo, rev = true)
    insertcols!(joined_importances, 1, :rank => collect(1:nrow(joined_importances)))

    return joined_importances

end

function plot_comparativedemo_importance_barplots!(plot_axis::Axis, joined_importances::DataFrame; n_rows = 40)

    # Plot barplot
    barplot!(plot_axis, reverse(collect(1:n_rows)), joined_importances.ImportanceWithoutDemo[1:n_rows], color = "medium turquoise", direction = :x)
    scatter!(plot_axis, joined_importances.ImportanceWithDemo[1:n_rows], reverse(collect(1:n_rows)), marker = :star8, color = "purple", direction = :x)
    
    return plot_axis

end
```

### Comparative Linear/RF importances scatterplots
```julia
function plot_comparative_lmvsrf_scatterplots!(plot_axis::Axis, rf_model::ProbeData, lm_path::String;)

    rf_model_importances = weightedimportances(rf_model; change_hashnames=false)
    lm_coefs = CSV.read(lm_path, DataFrame)
    lm_coefs = select(lm_coefs, [:feature, :coef, :pvalue])
    lm_coefs.coef = sqrt.(lm_coefs.coef .^ 2)
    plot_comaprative_df = outerjoin(rf_model_importances, lm_coefs, on = [ :variable => :feature ])

    # Plot scatterplot
    scatter!(plot_axis, log.(plot_comaprative_df.pvalue).*(-1), plot_comaprative_df.weightedImportance, color = "blue")

    return plot_axis
end
```

### Taxon deepdive quantile plots
```julia
function plot_taxon_deepdive!(figure_layout::GridLayout, figure_col::Int, taxonomic_profile::CommunityProfile, filter_row::Symbol, taxon_to_dive::String;)

    filtered_spec = taxonomic_profile[:, get(spec, filter_row)]
    (lq, uq) = quantile(get(filtered_spec, :cogScore), [0.25, 0.75])
    lqidx = findall(x-> x <= lq, get(filtered_spec, :cogScore))
    mqidx = findall(x-> lq < x <= uq, get(filtered_spec, :cogScore))
    uqidx = findall(x-> x > uq, get(filtered_spec, :cogScore))

    taxab = vec(abundances(filtered_spec[Regex(taxon_to_dive), :]))
    g1 = filter(!=(0), taxab[lqidx])
    g2 = filter(!=(0), taxab[mqidx])
    g3 = filter(!=(0), taxab[uqidx])

    g1p = count(!=(0), taxab[lqidx]) / length(lqidx)
    g2p = count(!=(0), taxab[mqidx]) / length(mqidx)
    g3p = count(!=(0), taxab[uqidx]) / length(uqidx)

    ax1 = Axis(figure_layout[1,figure_col]; ylabel="relative abundance\n$taxon_to_dive")
    ax2 = Axis(figure_layout[2,figure_col]; xticks = (1:3, ["lower", "mid", "upper"]), xlabel="quartile", ylabel = "Intra-q\nPrevalence")
    
    Random.seed!(0)
    scatter!(ax1, 1 .+ rand(Normal(0, 0.1), length(g1)), g1)
    Random.seed!(0)
    scatter!(ax1, 2 .+ rand(Normal(0, 0.1), length(g2)), g2)
    Random.seed!(0)
    scatter!(ax1, 3 .+ rand(Normal(0, 0.1), length(g3)), g3)
    barplot!(ax2, [1,2,3], [g1p, g2p, g3p]; color=Makie.wong_colors()[1:3])
    
    hidexdecorations!(ax1; grid=false)
    rowsize!(figure_layout, 2, Relative(1/3))

    return figure_layout
end
```

## Loading the pretrained models

```julia
RandomForestClassifier = MLJ.@load RandomForestClassifier pkg=DecisionTree
RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree
# concurrent cogScore reression from taxonomic profiles
JLD2.@load "models/2023-02-11/brain_models.jld"

ordered_brain_segments_list = [
    "left-temporal", "right-temporal",
    "left-orbitofrontal", "right-orbitofrontal",
    "left-parietal", "right-parietal",
    "left-middle-frontal", "right-middle-frontal",
    "left-anterior-cingulate", "right-anterior-cingulate",
    "left-lateral-occipital", "right-lateral-occipital",
    "left-cerebellum-white-matter", "right-cerebellum-white-matter",
    "left-thalamus-proper", "right-thalamus-proper",
    "left-caudate", "right-caudate",
    "left-putamen", "right-putamen",
    "left-pallidum", "right-pallidum",
    "left-hippocampus", "right-hippocampus",
    "left-amygdala", "right-amygdala",
    "left-accumbens-area", "right-accumbens-area",
    "left-basal-forebrain", "right-basal-forebrain", 
    "left-cuneus", "right-cuneus", 
    "left-entorhinal", "right-entorhinal",
    "left-fusiform", "right-fusiform",
    "left-isthmus-cingulate", "right-isthmus-cingulate",
    "left-lingual", "right-lingual",
    "left-parahippocampal", "right-parahippocampal",
    "left-paracentral", "right-paracentral",
    "left-pars-opercularis", "right-pars-opercularis",
    "left-pars-orbitalis", "right-pars-orbitalis",
    "left-pars-triangularis", "right-pars-triangularis",
    "left-pericalcarine", "right-pericalcarine",
    "left-postcentral", "right-postcentral",
    "left-posterior-cingulate", "right-posterior-cingulate",
    "left-precentral", "right-precentral",
    "left-precuneus", "right-precuneus",
    "left-superior-frontal", "right-superior-frontal",
    "left-supramarginal", "right-supramarginal",
    "left-insula", "right-insula",
    "cerebellar-vermal-lobules-I-V", "cerebellar-vermal-lobules-VI-VII", "cerebellar-vermal-lobules-VIII-X",
    "Brain-stem", "CSF"
]

symmetric_segment_strings = [
    "temporal", "temporal",
    "orbitofrontal", "orbitofrontal",
    "parietal", "parietal",
    "middle-frontal", "middle-frontal",
    "anterior-cingulate", "anterior-cingulate",
    "lateral-occipital", "lateral-occipital",
    "cerebellum-white-matter", "cerebellum-white-matter",
    "thalamus-proper", "thalamus-proper",
    "caudate", "caudate",
    "putamen", "putamen",
    "pallidum", "pallidum",
    "hippocampus", "hippocampus",
    "amygdala", "amygdala",
    "accumbens-area", "accumbens-area",
    "basal-forebrain", "basal-forebrain", 
    "cuneus", "cuneus", 
    "entorhinal", "entorhinal",
    "fusiform", "fusiform",
    "isthmus-cingulate", "isthmus-cingulate",
    "lingual", "lingual",
    "parahippocampal", "parahippocampal",
    "paracentral", "paracentral",
    "pars-opercularis", "pars-opercularis",
    "pars-orbitalis", "pars-orbitalis",
    "pars-triangularis", "pars-triangularis",
    "pericalcarine", "pericalcarine",
    "postcentral", "postcentral",
    "posterior-cingulate", "posterior-cingulate",
    "precentral", "precentral",
    "precuneus", "precuneus",
    "superior-frontal", "superior-frontal",
    "supramarginal", "supramarginal",
    "insula", "insula",
    "cerebellar-vermal-lobules-I-V", "cerebellar-vermal-lobules-VI-VII", "cerebellar-vermal-lobules-VIII-X",
    "Brain-stem", "CSF"
]
```

## Initializing Figure 4
```julia
figure = Figure(resolution = (1920, 1080))

A_subfig = GridLayout(figure[1,1])
B_subfig = GridLayout(figure[1,2])

colsize!(figure.layout, 1, Relative(0.2))
colsize!(figure.layout, 2, Relative(0.8))
```

### Plot figure 4 - FULL version without filtering segments or taxa, taxa ordered by mean importance, segments by hclust
```julia
mean_brain_merits = reduce(
    vcat,
    [ DataFrame(:variable => ordered_brain_segments_list[i], :Train_Cor => mean(brain_models[ordered_brain_segments_list[i]].merits.Train_Cor), :Test_Cor => mean(brain_models[ordered_brain_segments_list[i]].merits.Test_Cor)) for i in eachindex(ordered_brain_segments_list) ]
)

mean_brain_importances = reduce(
    (x, y) -> outerjoin(x, y, on = :variable, makeunique=true),
    [ 
        DataFrame(:variable => j.importances.variable, Symbol(split(j.name, '_')[2]) => map(mean, eachrow(j.importances[:, 2:end]))) 
        for (i, j) in brain_models
    ]
)

mean_brain_importances = mean_brain_importances[:, vcat([ "variable" ], ordered_brain_segments_list)]
mean_brain_importances.averageImportance = [ sum(row[.!(isnan.(collect(row)))])/length(collect(row)) for row in eachrow(mean_brain_importances[:, 2:end]) ]
sort!(mean_brain_importances, :averageImportance, rev=true)

transposedImportances = DataFrame(permutedims(Matrix(mean_brain_importances)), :auto)
transposedImportances = rename!(transposedImportances, Symbol.(Vector(transposedImportances[1,:])))[2:end-1,:]
insertcols!(transposedImportances, 1, :symmetricSegment => symmetric_segment_strings)
insertcols!(transposedImportances, 1, :longSegment => ordered_brain_segments_list)

combinedTransposedImportances = combine(
    groupby(transposedImportances, :symmetricSegment),
    Symbol.(names(transposedImportances))[3:end] .=> mean .=> Symbol.(names(transposedImportances))[3:end]
)

dist_taxa = pairwise(Euclidean(), Matrix(combinedTransposedImportances[:, 2:end]); dims=2)
dist_segments = pairwise(Euclidean(), Matrix(combinedTransposedImportances[:, 2:end]); dims=1)

hcl_taxa = hclust(dist_taxa; linkage=:average, branchorder=:optimal)
hcl_segments = hclust(dist_segments; linkage=:average, branchorder=:optimal)

hclust_taxa_order = hcl_taxa.order
hclust_symmetric_segment_order = hcl_segments.order

reorder_segments_df = innerjoin(
    DataFrame(
        :original_segment => ordered_brain_segments_list,
        :left_or_unique => vcat(repeat([1, 0], 33), repeat([ 1 ], 5)),
        :symmetric_segment => symmetric_segment_strings),
    DataFrame(
        :symmetric_segment => combinedTransposedImportances.symmetricSegment,
        :hclust_order => hclust_symmetric_segment_order),
    on = :symmetric_segment
)
reorder_segments_df.plot_order = collect(1:nrow(reorder_segments_df))
sort!(reorder_segments_df, [:hclust_order, :left_or_unique])

plot_segments_order = reorder_segments_df.plot_order
## using the hclust to order it
mean_brain_merits = mean_brain_merits[plot_segments_order, :]

axA = Axis(
    A_subfig[1, 1];
    xlabel = "Correlation",
    yticks = (reverse(collect(1:nrow(mean_brain_merits))), mean_brain_merits.variable),
    ylabel = "Target Variable",
    xticks = [-0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
    title = "Mean Random Forest regression correlations",
    yticklabelsize=16
    #yticklabelrotation= -pi/2
)

tightlimits!(axA, Top())
tightlimits!(axA, Bottom())
ylims!(axA, 0, length(ordered_brain_segments_list)+1)

# Plot barplot
barplot!(
    axA,
    reverse(collect(1:nrow(mean_brain_merits))),
    mean_brain_merits.Test_Cor,
    color = vcat( repeat( [ "blue", "red" ], 33),repeat( [ "purple" ], 5) )[plot_segments_order],
    direction=:x
)

######
# Importance Heatmaps
######

mean_brain_importances = mean_brain_importances[:, vcat([1], plot_segments_order .+ 1)]

axB = Axis(
    B_subfig[1, 1];
    xlabel = "Predictor",
    xticks = (collect(1:40), mean_brain_importances.variable[2:41]),
    xticklabelsize=16,
    xticklabelrotation= -pi/4,
    ylabel = "Target brain segment",
    yticks = (collect(1:ncol(mean_brain_importances)-2), names(mean_brain_importances)[2:end-1]),
    yticklabelsize=16,
    yreversed=true,
    title = "Brain segmentation data variable importances"
)

heatmap!(axB, Matrix(mean_brain_importances[2:41, 2:end-1]), yflip=true)
figure
```
### Plot figure 4 again, now with filtering segments and taxa based on predetermined lists, taxa ordered by mean importance

```julia
interesting_segments = [
    "left-temporal", "right-temporal",
    "left-orbitofrontal", "right-orbitofrontal",
    "left-parietal", "right-parietal",
    "left-middle-frontal", "right-middle-frontal",
    "left-anterior-cingulate", "right-anterior-cingulate",
    "left-lateral-occipital", "right-lateral-occipital",
    "left-thalamus-proper", "right-thalamus-proper",
    "left-caudate", "right-caudate",
    "left-putamen", "right-putamen",
    "left-pallidum", "right-pallidum",
    "left-hippocampus", "right-hippocampus",
    "left-amygdala", "right-amygdala",
    "left-accumbens-area", "right-accumbens-area",
    "left-cuneus", "right-cuneus", 
    "left-entorhinal", "right-entorhinal",
    "left-fusiform", "right-fusiform",
    "left-isthmus-cingulate", "right-isthmus-cingulate",
    "left-lingual", "right-lingual",
    "left-parahippocampal", "right-parahippocampal",
    "left-paracentral", "right-paracentral",
    "left-pars-opercularis", "right-pars-opercularis",
    "left-pars-orbitalis", "right-pars-orbitalis",
    "left-pericalcarine", "right-pericalcarine",
    "left-posterior-cingulate", "right-posterior-cingulate",
    "left-precentral", "right-precentral",
    "left-precuneus", "right-precuneus",
    "left-superior-frontal", "right-superior-frontal",
    "cerebellar-vermal-lobules-I-V", "cerebellar-vermal-lobules-VI-VII", "cerebellar-vermal-lobules-VIII-X",
]

interesting_segments_idxes = mean_brain_merits.variable .∈ Ref(interesting_segments)

interesting_taxa = [
    "Anaerostipes_hadrus",
    "Fusicatenibacter_saccharivorans",
    "Eubacterium_rectale",
    "Eubacterium_hallii",
    "Blautia_wexlerae",
    "Bacteroides_vulgatus",
    "Faecalibacterium_prausnitzii",
    "Ruminococcus_torques",
    "Bifidobacterium_longum",
    "Bacteroides_ovatus",
    "Bacteroides_uniformis",
    "Coprococcus_comes",
    "Alistipes_finegoldii",
    "Flavonifractor_plautii",
    "Roseburia_intestinalis",
    "Streptococcus_salivarius",
    "Eubacterium_eligens",
    "Ruminococcus_bicirculans",
    "Collinsella_aerofaciens",
    "Ruminococcus_gnavus",
    "Dorea_formicigenerans",
    "Eggerthella_lenta",
    "Asaccharobacter_celatus",
    "Gordonibacter_pamelaeae"
]

figure = Figure(resolution = (1920, 1080))

A_subfig = GridLayout(figure[1,1])
B_subfig = GridLayout(figure[1,2])

colsize!(figure.layout, 1, Relative(0.2))
colsize!(figure.layout, 2, Relative(0.8))

mean_brain_merits = reduce(
    vcat,
    [ DataFrame(:variable => ordered_brain_segments_list[i], :Train_Cor => mean(brain_models[ordered_brain_segments_list[i]].merits.Train_Cor), :Test_Cor => mean(brain_models[ordered_brain_segments_list[i]].merits.Test_Cor)) for i in eachindex(ordered_brain_segments_list) ]
)

axA = Axis(
    A_subfig[1, 1];
    xlabel = "Correlation",
    yticks = (reverse(collect(1:length(interesting_segments))), mean_brain_merits.variable[interesting_segments_idxes]),
    ylabel = "Target Variable",
    xticks = [-0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
    title = "Mean Random Forest regression correlations",
    yticklabelsize=16
    #yticklabelrotation= -pi/2
)

ylims!(axA, 0, length(interesting_segments)+1)

# Plot barplot
barplot!(
    axA,
    reverse(collect(1:length(interesting_segments))),
    mean_brain_merits.Test_Cor[interesting_segments_idxes],
    color = vcat( repeat( [ "blue", "red" ], 27),repeat( [ "purple" ], 3) ),
    direction=:x
)

######
# Importance Heatmaps
######

mean_brain_importances = reduce(
    (x, y) -> outerjoin(x, y, on = :variable, makeunique=true),
    [ 
        DataFrame(:variable => j.importances.variable, Symbol(split(j.name, '_')[2]) => map(mean, eachrow(j.importances[:, 2:end]))) 
        for (i, j) in brain_models
    ]
)

mean_brain_importances = mean_brain_importances[:, vcat([ "variable" ], ordered_brain_segments_list)]

mean_brain_importances.averageImportance = [ sum(row[.!(isnan.(collect(row)))])/length(collect(row)) for row in eachrow(mean_brain_importances[:, 2:end]) ]
sort!(mean_brain_importances, :averageImportance, rev=true)

interesting_taxa_idxes = mean_brain_importances.variable .∈ Ref(interesting_taxa)

axB = Axis(
    B_subfig[1, 1];
    xlabel = "Predictor",
    xticks = (collect(1:length(interesting_taxa)), interesting_taxa),
    xticklabelsize=16,
    xticklabelrotation= -pi/4,
    ylabel = "Target brain segment",
    yticks = (collect(1:length(interesting_segments)), names(mean_brain_importances)[2:end-1][interesting_segments_idxes]),
    yticklabelsize=16,
    yreversed=true,
    title = "Brain segmentation data variable importances"
)

heatmap!(axB, Matrix(mean_brain_importances[interesting_taxa_idxes, 2:end-1][!, interesting_segments_idxes]), yflip=true)
figure
```