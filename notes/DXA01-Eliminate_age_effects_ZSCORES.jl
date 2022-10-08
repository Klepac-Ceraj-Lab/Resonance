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
using CairoMakie
using CubicSplines
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
    percentile_growth_curves::Dict{Float64, CubicSpline{Float64}}
end

# struct GrowthCurveEnsemble
#     growth_curves::Dict{String, AgeBracket}
# end

function compute_sex_age_bracket(df, variable, lb, ub, sex)
    
    bracket_samples = findall( ( lb .< df.ageMonths .<= ub ) .& (df.sex .== sex) )
    bracket_values = df[bracket_samples , variable]

    # CDC/WHO Growh curves compute the following percentiles:
    # (3rd, 5th, 10th, 25th, 50th, 75th, 90th, 95th, and 97th).
    # We do not have enough points to compute all of those.
    # So maybe we can start with
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

    return SexAgeBracket(
        sex, lb, ub, bracket_samples,
        quantiles, spline
    )

end

function compute_age_bracket(idx, df, variable, lb, ub)

    AgeBracket(
        idx,
        Dict(
            "Male" => compute_sex_age_bracket(df, variable, lb, ub, "Male"),
            "Female" => compute_sex_age_bracket(df, variable, lb, ub, "Female")
        )
    )

end

function compute_growth_curve(df, variable, intervals)

    GrowthCurve(
        variable,
        intervals,
        [
            compute_age_bracket(idx, df, variable, lb, ub)
            for (idx, lb, ub) in zip(1:(length(intervals)-1), intervals[1:end-1], intervals[2:end])
        ]
    )     

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

gc = compute_growth_curve(mdata_brain_df, "left-thalamus-proper", tentative_intervals)

# xs = range(0.005, 0.009, length=1000)
# ys = bracket_spline[xs]
# plot(xs, ys)
