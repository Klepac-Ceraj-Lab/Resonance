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
# Data structures and Functions
#####

struct SexAgeBracketPercentiles
    sex::String
    lb::Float64
    ub::Float64
    samples::Vector{Integer}
    quantiles::Dict{Float64, Float64}
    quant_fun::CubicSpline{Float64}
end

struct AgeBracketPercentiles
    idx::Integer
    sex_age_brckets::Dict{String, SexAgeBracketPercentiles}
end

struct GrowthCurve
    variable::String
    intervals::Vector{Float64}
    bracket_percentiles::Vector{AgeBracketPercentiles}
    percentile_growth_curves::Dict{String, Dict{Float64, CubicSpline{Float64}}}
end

# struct GrowthCurveEnsemble
#     growth_curves::Dict{String, AgeBracket}
# end

function compute_sex_age_bracket(df, variable, lb, ub, sex)
    
    bracket_samples = findall( ( lb .< df.ageMonths .<= ub ) .& (df.sex .== sex) )
    bracket_values = df[bracket_samples , variable]

    # CDC/WHO Growh curves compute the following percentiles:
    # (3rd, 5th, 10th, 25th, 50th, 75th, 90th, 95th, and 97th).
    # We do not have enough points to compute all of those without overlapping,
    # so maybe we can start with
    # [0.0, 0.10, 0.25, 0.50, 0.75, 0.90, 1.0]

    quantiles = Dict(
        0.00 => quantile(bracket_values, 0.00),
        0.10 => quantile(bracket_values, 0.10),
        0.25 => quantile(bracket_values, 0.25),
        0.50 => quantile(bracket_values, 0.50),
        0.75 => quantile(bracket_values, 0.75),
        0.90 => quantile(bracket_values, 0.90),
        1.00 => quantile(bracket_values, 1.00)
    )

    spline = CubicSpline(
        sort(collect(values(quantiles))),
        sort(collect(keys(quantiles))),
        extrapl=[1,], extrapr=[1,]
    )

    return SexAgeBracketPercentiles(
        sex, lb, ub, bracket_samples,
        quantiles, spline
    )

end

function compute_age_bracket(idx, df, variable, lb, ub)

    AgeBracketPercentiles(
        idx,
        Dict(
            "Male" => compute_sex_age_bracket(df, variable, lb, ub, "Male"),
            "Female" => compute_sex_age_bracket(df, variable, lb, ub, "Female")
        )
    )

end

function compute_growth_curve(df, variable, intervals)

    bracket_percentiles =  [
        compute_age_bracket(idx, df, variable, lb, ub)
        for (idx, lb, ub) in zip(1:(length(intervals)-1), intervals[1:end-1], intervals[2:end])
    ]

    xs = [ (intervals[i] + intervals[i+1])/2 for i in 1:(length(intervals)-1)]
    percentile_gcs = Dict{String, Dict{Float64, CubicSpline{Float64}}}()
    for sex in [ "Male", "Female" ]
        gc = Dict{Float64, CubicSpline{Float64}}()
        for qntl in [0.0, 0.10, 0.25, 0.50, 0.75, 0.90, 1.0] 
            ys = [ bracket_percentiles[i].sex_age_brckets[sex].quantiles[qntl] for i in eachindex(bracket_percentiles) ]
            curve = CubicSpline(
                xs,
                ys
            )

            println(xs, ys)

            push!(gc, qntl => curve)
        end # end for quantile
        push!(percentile_gcs, sex => gc)
    end # end for sex

    return GrowthCurve(
        variable,
        intervals,
        bracket_percentiles,
        percentile_gcs
    )

end # end function

function plot_multiple_growthcurves!(
    figure::Figure,
    pos::Tuple{Int64, Int64},
    gc::GrowthCurve,
    # variable::String,
    sex::String,
    plot_title::String;
    xs = 3.0:0.5:100.0)

    ax = Axis(fig[pos[1],pos[2]], title = plot_title)
    
    lins = [
        lines!(xs, gc.percentile_growth_curves[sex][0.0][collect(xs)];color=:firebrick1, linewidth = 1.0)
        lines!(xs, gc.percentile_growth_curves[sex][0.10][collect(xs)];color=:red, linewidth = 2.0)
        lines!(xs, gc.percentile_growth_curves[sex][0.25][collect(xs)];color=:red3, linewidth = 3.0)
        lines!(xs, gc.percentile_growth_curves[sex][0.50][collect(xs)];color=:black, linewidth = 6.0)
        lines!(xs, gc.percentile_growth_curves[sex][0.75][collect(xs)];color=:blue3, linewidth = 3.0)
        lines!(xs, gc.percentile_growth_curves[sex][0.90][collect(xs)];color=:blue, linewidth = 2.0)
        lines!(xs, gc.percentile_growth_curves[sex][1.00][collect(xs)];color=:dodgerblue2, linewidth = 1.0)
    ]
 
    return figure
end

function plot_all_results!(
    figure::Figure,
    sex::String,
    variables::Vector{String},
    gcs::Dict{String, GrowthCurve},
    plot_positions::Vector{Tuple{Int64, Int64}},
    general_title::String
    )

    for i in eachindex(variables)
        plot_multiple_growthcurves!(
            figure,
            plot_positions[i],
            gcs[variables[i]],
            sex,
            variables[i])
    end

    Label(figure[1, :, Top()], general_title,
    padding = (0, 50, 40, 0))

    Legend(
    fig[end+1,:],
    [
        LineElement(color=:firebrick1, linewidth = 1.0),
        LineElement(color=:red, linewidth = 2.0),
        LineElement(color=:red3, linewidth = 3.0),
        LineElement(color=:black, linewidth = 6.0),
        LineElement(color=:blue3, linewidth = 3.0),
        LineElement(color=:blue, linewidth = 2.0),
        LineElement(color=:dodgerblue2, linewidth = 1.0),
    ], [
        "min",
        "10th pctl",
        "25th pctl",
        "50th pctl",
        "75th pctl",
        "90th pctl",
        "max"
    ],
    orientation = :horizontal
)

    return figure

end

#####
# Loading Data
#####

mdata = @chain Resonance.load(Metadata()) begin
    select( [ :subject, :timepoint, :ageMonths, :cogScore, :omni ] )
end # no dropmissing(:omni) because subjects/timepoints without sample but with assessment must be considered on the prediction routine

brain = Resonance.load(Resonance.Neuroimaging(); timepoint_metadata = mdata)
brain_df = DataFrame([eachcol(collect(brain.abundances'))...], map( x -> string(x), brain.features), copycols=true)

mdata_brain_df =  @chain DataFrame(metadata(brain)) begin 
    hcat( _ , brain_df)
    select( _ , Not([:sample, :hassample, :cogScore]))
    dropmissing()
end ## Columns 1:3 are metadata, columns 4:552 are taxonomic profile, columns 553:651 are brain data

#####
# Computing
#####

tentative_intervals = [ 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 21.0, collect(24.0:6.0:120.0)... ]

brainGrowthCurves = Dict{String, GrowthCurve}()

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

for segment in original_prediction_segments
    println(segment)
    gc = compute_growth_curve(mdata_brain_df, segment, tentative_intervals)
    push!(brainGrowthCurves, segment => gc)
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
    "Male",
    original_prediction_segments,
    brainGrowthCurves,
    plot_positions,
    "Growth curves for % volume of different brain segments - MALE subjects"
    )

save("figures/male_braingrowthcurve.png", fig)

fig = Figure(resolution = (8000, 4096))

plot_all_results!(
    fig,
    "Female",
    original_prediction_segments,
    brainGrowthCurves,
    plot_positions,
    "Growth curves for % volume of different brain segments - FEMALE subjects"
    )

save("figures/female_braingrowthcurve.png", fig)