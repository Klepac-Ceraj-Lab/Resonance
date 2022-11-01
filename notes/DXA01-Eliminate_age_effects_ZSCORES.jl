#####
# Notebook D21 - Independent regression of current Brain Segment from taxonomic data
#####

using Resonance
using Chain
using CSV
using Statistics
using Random
using DataFrames
using StableRNGs
using MLJ
using DecisionTree
using GLM
using JLD2
using CubicSplines
using CairoMakie
ml_rng = StableRNG(0)

#####
# Loading Data
#####

mdata = @chain Resonance.load(Metadata()) begin
    select( [ :subject, :timepoint, :ageMonths, :cogScore, :sex, :omni ] )
end # no dropmissing(:omni) because subjects/timepoints without sample but with assessment must be considered on the prediction routine

brain = Resonance.load(Resonance.Neuroimaging(); timepoint_metadata = mdata)
brain_df = DataFrame([eachcol(collect(brain.abundances'))...], map( x -> string(x), brain.features), copycols=true)

mdata_brain_df =  @chain DataFrame(metadata(brain)) begin 
    hcat( _ , brain_df)
    select( _ , Not([:sample, :hassample, :cogScore]))
    subset(:ageMonths => x -> x .<= 120.0)
    dropmissing()
end ## Columns 1:3 are metadata, columns 4:552 are taxonomic profile, columns 553:651 are brain data

#####
# Computing the curves for cogScore
#####

cogscore_df = @chain mdata begin
    dropmissing([:cogScore])
    subset(:ageMonths => x -> x .<= 120)
end

# cogscore_intervals = [ 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 15.0, 18.0, 21.0, collect(24.0:6.0:42.0)..., collect(48.0:12.0:120.0)... ] # OLD version
cogscore_intervals = [ 0.0, 6.0, 12.0, 18.0, 24.0, 30.0, 36.0, 42.0, collect(48.0:18.0:120.0)... ]
hist(cogscore_df.ageMonths, bins = cogscore_intervals; color = :aqua, strokewidth = 1, strokecolor = :black)

cogScoreGrowthCurve = compute_growth_curve(cogscore_df, "cogScore", cogscore_intervals, nothing)

fig = Figure(resolution = (1200, 600))
plot_multiple_growthcurves!(fig, (1,1), cogScoreGrowthCurve, "")

Label(fig[1, :, Top()], "cogScores growth curves", padding = (0, 50, 40, 0))

Legend( fig[end+1,:],[
    LineElement(color=:red, linewidth = 2.0),
    LineElement(color=:red3, linewidth = 3.0),
    LineElement(color=:black, linewidth = 6.0),
    LineElement(color=:blue3, linewidth = 3.0),
    LineElement(color=:blue, linewidth = 2.0) ],
    ["10th pctl", "25th pctl", "50th pctl", "75th pctl", "90th pctl"], orientation = :horizontal
)

fig
save("figures/cogscore_growthcurve.png", fig)

#####
# Computing the curves for Brain
#####

# brain_intervals = [ 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 21.0, collect(24.0:6.0:120.0)... ]
brain_intervals = [ 0.0, 6.0, 12.0, 18.0, 24.0, 36.0, collect(48.0:18.0:120.0)... ]
hist(mdata_brain_df.ageMonths, bins = brain_intervals)

original_prediction_segments = [
    "left-cerebellum-exterior", "right-cerebellum-exterior",
    "left-cerebellum-white-matter", "right-cerebellum-white-matter",
    "left-thalamus-proper", "right-thalamus-proper",
    "left-caudate", "right-caudate",
    "left-putamen", "right-putamen",
    "left-pallidum", "right-pallidum",
    "left-hippocampus", "right-hippocampus",
    "left-amygdala", "right-amygdala",
    "left-accumbens-area", "right-accumbens-area",
    "left-ventral-DC", "right-ventral-DC",
    "left-basal-forebrain", "right-basal-forebrain",
    "cerebellar-vermal-lobules-I-V", "cerebellar-vermal-lobules-VI-VII", "cerebellar-vermal-lobules-VIII-X",
    "left-caudal-anterior-cingulate", "right-caudal-anterior-cingulate",
    "left-caudal-middle-frontal", "right-caudal-middle-frontal",
    "left-cuneus", "right-cuneus", 
    "left-entorhinal", "right-entorhinal",
    "left-fusiform", "right-fusiform",
    "left-inferior-parietal", "right-inferior-parietal",
    "left-inferior-temporal", "right-inferior-temporal",
    "left-isthmus-cingulate", "right-isthmus-cingulate",
    "left-lateral-occipital", "right-lateral-occipital",
    "left-lateral-orbitofrontal", "right-lateral-orbitofrontal",
    "left-lingual", "right-lingual", 
    "left-medial-orbitofrontal", "right-medial-orbitofrontal",
    "left-middle-temporal", "right-middle-temporal",
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
    "left-rostral-anterior-cingulate", "right-rostral-anterior-cingulate",
    "left-rostral-middle-frontal", "right-rostral-middle-frontal",
    "left-superior-frontal", "right-superior-frontal",
    "left-superior-parietal", "right-superior-parietal",
    "left-superior-temporal", "right-superior-temporal",
    "left-supramarginal", "right-supramarginal",
    "left-transverse-temporal", "right-transverse-temporal",
    "left-insula", "right-insula",
    "Brain-stem"
]

maleBrainGrowthCurves = Dict{String, GrowthCurve}()
femaleBrainGrowthCurves = Dict{String, GrowthCurve}()

for segment in original_prediction_segments
    println(segment)
    gc = compute_growth_curve(mdata_brain_df, segment, brain_intervals, "Male")
    push!(maleBrainGrowthCurves, segment => gc)
    gc = compute_growth_curve(mdata_brain_df, segment, brain_intervals, "Female")
    push!(femaleBrainGrowthCurves, segment => gc)
end

#####
# Methods to plot growth curves
#####

cols = collect(1:12)
rows = collect(1:8)
plot_positions = vec(permutedims(collect(Base.product(rows, cols))))

fig = Figure(resolution = (8192, 4096))

plot_all_results!(
    fig,
    original_prediction_segments,
    maleBrainGrowthCurves,
    plot_positions,
    "Growth curves for % volume of different brain segments - MALE subjects"
    )

save("figures/male_braingrowthcurve.png", fig)

fig = Figure(resolution = (8000, 4096))

plot_all_results!(
    fig,
    original_prediction_segments,
    femaleBrainGrowthCurves,
    plot_positions,
    "Growth curves for % volume of different brain segments - FEMALE subjects"
    )

save("figures/female_braingrowthcurve.png", fig)

#####
# Methods to find out the percentiles of each subject
#####

cogscore_df.cogScorePercentile = [
    get_cogscore_percentile(cogScoreGrowthCurve, cogscore_intervals, a, v) for (a,v) in zip(cogscore_df.ageMonths, cogscore_df.cogScore)
]

combined_sex_growthcurves = Dict(
    "Male" => maleBrainGrowthCurves,
    "Female" => femaleBrainGrowthCurves
)

for segment in original_prediction_segments

    insertcols!(
        mdata_brain_df,
        segment*"-percentile" => [
            get_brain_percentile(combined_sex_growthcurves, segment, brain_intervals, s, a, v) for (s,a,v) in zip(mdata_brain_df.sex, mdata_brain_df.ageMonths, mdata_brain_df[:, segment])
        ]
    )
end

CSV.write("cogscore_percentiles.csv", cogscore_df)