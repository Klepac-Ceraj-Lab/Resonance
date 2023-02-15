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

function weighted_hpimportances(m, hp = 1; change_hashnames=false, hashnamestable::Union{Nothing, DataFrame} = nothing)
    merits = m.merits
    importances = m.importances
    fitnesses = calculate_fitness(merits.Train_Cor, merits.Test_Cor)
    fitness_subset_idx = findall(merits.Hyperpar_Idx .== hp)

    mean_importances = map(x -> sum(x .* fitnesses[fitness_subset_idx])/sum(fitnesses[fitness_subset_idx] .> 0.0), collect(eachrow(Matrix(importances[:, fitness_subset_idx .+ 1]))))
    mean_importances_df = sort( DataFrame(:variable => importances[:,1], :weightedImportance => mean_importances), :weightedImportance; rev=true)

    if change_hashnames
        newnames = leftjoin(mean_importances_df, hashnamestable, on = :variable => :hashname).longname
        mean_importances_df.variable = newnames
    end

    return mean_importances_df
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

interesting_segments = [
    "left-temporal", "right-temporal",
    "left-orbitofrontal", "right-orbitofrontal",
    "left-parietal", "right-parietal",
    "left-middle-frontal", "right-middle-frontal",
    "left-anterior-cingulate", "right-anterior-cingulate",
    "left-lateral-occipital", "right-lateral-occipital",
    "left-thalamus-proper", "right-thalamus-proper",
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

symmetric_segment_strings = [
    "temporal", "temporal",
    "orbitofrontal", "orbitofrontal",
    "parietal", "parietal",
    "middle-frontal", "middle-frontal",
    "anterior-cingulate", "anterior-cingulate",
    "lateral-occipital", "lateral-occipital",
    "thalamus-proper", "thalamus-proper",
    "hippocampus", "hippocampus",
    "amygdala", "amygdala",
    "accumbens-area", "accumbens-area",
    "cuneus", "cuneus", 
    "entorhinal", "entorhinal",
    "fusiform", "fusiform",
    "isthmus-cingulate", "isthmus-cingulate",
    "lingual", "lingual",
    "parahippocampal", "parahippocampal",
    "paracentral", "paracentral",
    "pars-opercularis", "pars-opercularis",
    "pars-orbitalis", "pars-orbitalis",
    "pericalcarine", "pericalcarine",
    "posterior-cingulate", "posterior-cingulate",
    "precentral", "precentral",
    "precuneus", "precuneus",
    "superior-frontal", "superior-frontal",
    "cerebellar-vermal-lobules-I-V", "cerebellar-vermal-lobules-VI-VII", "cerebellar-vermal-lobules-VIII-X",
]

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

#### Compute the mean figures of merit for the Brain sement regressions
```julia
mean_brain_merits = reduce(
    vcat,
    [ DataFrame(:variable => ordered_brain_segments_list[i], :Train_Cor => mean(brain_models[ordered_brain_segments_list[i]].merits.Train_Cor), :Test_Cor => mean(brain_models[ordered_brain_segments_list[i]].merits.Test_Cor)) for i in eachindex(ordered_brain_segments_list) ]
)
```

#### Compute the mean importances loading of taxa on Brain sement regressions
```julia
mean_brain_importances = reduce(
    (x, y) -> outerjoin(x, y, on = :variable, makeunique=true),
    [ 
        DataFrame(:variable => j.importances.variable, Symbol(split(j.name, '_')[2]) => map(mean, eachrow(j.importances[:, 2:end]))) 
        for (i, j) in brain_models
    ]
)
```

#### Filter the tables accordign to the lsit of interesting taxa and segments
```julia
interesting_segments_idxes = mean_brain_merits.variable .∈ Ref(interesting_segments)
interesting_taxa_idxes = mean_brain_importances.variable .∈ Ref(interesting_taxa)

mean_brain_merits = mean_brain_merits[interesting_segments_idxes, :]
mean_brain_importances = mean_brain_importances[interesting_taxa_idxes, vcat(["variable"], interesting_segments) ]
```

#### Perform hierarchical clustering
```julia
transposedImportances = DataFrame(permutedims(Matrix(mean_brain_importances)), :auto)
transposedImportances = rename!(transposedImportances, Symbol.(Vector(transposedImportances[1,:])))[2:end,:]
insertcols!(transposedImportances, 1, :symmetricSegment => symmetric_segment_strings)
insertcols!(transposedImportances, 1, :longSegment => interesting_segments)

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

### using the hclust to reorder all we need
```julia
mean_brain_merits = mean_brain_merits[plot_segments_order, :]
mean_brain_importances = mean_brain_importances[hclust_taxa_order, vcat( [ 1 ], (plot_segments_order .+ 1))]
```

### Calling the plot functions with the ordered data
```julia
axA = Axis(
    A_subfig[1, 1];
    xlabel = "Correlation",
    yticks = (reverse(collect(1:length(interesting_segments))), mean_brain_merits.variable),
    ylabel = "Target Variable",
    xticks = [-0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
    title = "Mean Random Forest regression correlations",
    yticklabelsize=16
    #yticklabelrotation= -pi/2
)

tightlimits!(axA, Top())
tightlimits!(axA, Bottom())

# Plot barplot
barplot!(
    axA,
    reverse(collect(1:length(interesting_segments))),
    mean_brain_merits.Test_Cor,
    color = vcat( repeat( [ "blue", "red" ], 24),repeat( [ "purple" ], 3) )[plot_segments_order],
    direction=:x
)

######
# Importance Heatmaps
######
axB = Axis(
    B_subfig[1, 1];
    xlabel = "Predictor",
    xticks = (collect(1:length(interesting_taxa)), mean_brain_importances.variable),
    xticklabelsize=16,
    xticklabelrotation= pi/4,
    ylabel = "Target brain segment",
    yticks = (collect(1:length(interesting_segments)), names(mean_brain_importances)[2:end]),
    yticklabelsize=16,
    yreversed=true,
    title = "Brain segmentation data variable importances"
)

heatmap!(axB, Matrix(mean_brain_importances[:, 2:end]), yflip=true)

save("manuscript/assets/Figure4.png", figure)
```