#####
# Functions and structures for percentile computation
#####

struct AgeBracketPercentiles
    idx::Int64
    sex::String
    lb::Float64
    ub::Float64
    samples::Vector{Integer}
    quantiles::Dict{Float64, Float64}
    quant_fun::CubicSpline{Float64}
end

struct GrowthCurve
    variable::String
    intervals::Vector{Float64}
    bracket_percentiles::Vector{AgeBracketPercentiles}
    percentile_growth_curves::Dict{Float64, CubicSpline{Float64}}
end

function compute_age_bracket(idx, df, variable, lb, ub; sex::Union{Nothing, String}=nothing)

    if sex isa Nothing
        bracket_samples = findall( lb .< df.ageMonths .<= ub )
    else   
        bracket_samples = findall( ( lb .< df.ageMonths .<= ub ) .& (df.sex .== sex) )
    end
    
    bracket_values = df[bracket_samples , variable]

    # CDC/WHO Growh curves compute the following percentiles:
    # (3rd, 5th, 10th, 25th, 50th, 75th, 90th, 95th, and 97th).
    # We do not have enough points to compute all of those without overlapping,
    # so maybe we can start with
    # [0.0, 0.10, 0.25, 0.50, 0.75, 0.90, 1.0] # annotation: changed.

    quantiles = Dict(
        0.00 => quantile(bracket_values, 0.00),
        0.05 => quantile(bracket_values, 0.05),
        0.10 => quantile(bracket_values, 0.10),
        0.25 => quantile(bracket_values, 0.25),
        0.50 => quantile(bracket_values, 0.50),
        0.75 => quantile(bracket_values, 0.75),
        0.90 => quantile(bracket_values, 0.90),
        0.95 => quantile(bracket_values, 0.95),
        1.00 => quantile(bracket_values, 1.00)
    )

    spline = CubicSpline(
        sort(collect(values(quantiles))),
        sort(collect(keys(quantiles))),
        extrapl=[0,], extrapr=[0,]
    )

    if sex isa Nothing
        return AgeBracketPercentiles(
            idx, "All", lb, ub, bracket_samples,
            quantiles, spline
        )    
    else   
        return AgeBracketPercentiles(
            idx, sex, lb, ub, bracket_samples,
            quantiles, spline
        )    
    end


end

function compute_growth_curve(df, variable, intervals, sex::Union{Nothing, String}=nothing)

    bracket_percentiles =  [
        compute_age_bracket(idx, df, variable, lb, ub)
        for (idx, lb, ub) in zip(1:(length(intervals)-1), intervals[1:end-1], intervals[2:end])
    ]

    xs = [ (intervals[i] + intervals[i+1])/2 for i in 1:(length(intervals)-1)]
    percentile_gcs = Dict{Float64, CubicSpline{Float64}}()

    for qntl in [0.00, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 1.00]
        ys = [ bracket_percentiles[i].quantiles[qntl] for i in eachindex(bracket_percentiles) ]
        curve = CubicSpline(
            xs,
            ys
        )
        println(xs, ys) # for debugging

        push!(percentile_gcs, qntl => curve)

    end # end for qntl

    return GrowthCurve(
        variable,
        intervals,
        bracket_percentiles,
        percentile_gcs
    )

end # end function

function get_cogscore_percentile(gc::GrowthCurve, intervals::Vector{Float64}, ageMonths::Float64, value::Float64)
    subject_bracket_idx = findfirst(ageMonths .<= intervals[2:end])
    return gc.bracket_percentiles[subject_bracket_idx].quant_fun[value]
end

function get_brain_percentile(gc::Dict{String, Dict{String, GrowthCurve}}, variable, intervals::Vector{Float64}, sex::String, ageMonths::Float64, value::Float64)
    subject_bracket_idx = findfirst(ageMonths .<= intervals[2:end])
    return gc[sex][variable].bracket_percentiles[subject_bracket_idx].quant_fun[value]
end

function plot_multiple_growthcurves!(
    figure::Figure,
    pos::Tuple{Int64, Int64},
    gc::GrowthCurve,
    # variable::String,
    plot_title::String;
    xs = 3.0:0.5:100.0)

    ax = Axis(fig[pos[1],pos[2]], title = plot_title)
    
    lins = [
        # lines!(xs, gc.percentile_growth_curves[0.05][collect(xs)];color=:firebrick1, linewidth = 1.0)
        lines!(xs, gc.percentile_growth_curves[0.10][collect(xs)];color=:red, linewidth = 2.0)
        lines!(xs, gc.percentile_growth_curves[0.25][collect(xs)];color=:red3, linewidth = 3.0)
        lines!(xs, gc.percentile_growth_curves[0.50][collect(xs)];color=:black, linewidth = 6.0)
        lines!(xs, gc.percentile_growth_curves[0.75][collect(xs)];color=:blue3, linewidth = 3.0)
        lines!(xs, gc.percentile_growth_curves[0.90][collect(xs)];color=:blue, linewidth = 2.0)
        # lines!(xs, gc.percentile_growth_curves[0.95][collect(xs)];color=:dodgerblue2, linewidth = 1.0)
    ]
 
    return figure
end

function plot_all_results!(
    figure::Figure,
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
            variables[i])
    end

    Label(figure[1, :, Top()], general_title,
    padding = (0, 50, 40, 0))

    Legend(
    figure[end+1,:],
    [
        # LineElement(color=:firebrick1, linewidth = 1.0),
        LineElement(color=:red, linewidth = 2.0),
        LineElement(color=:red3, linewidth = 3.0),
        LineElement(color=:black, linewidth = 6.0),
        LineElement(color=:blue3, linewidth = 3.0),
        LineElement(color=:blue, linewidth = 2.0),
        # LineElement(color=:dodgerblue2, linewidth = 1.0),
    ], [
        # "min",
        "10th pctl",
        "25th pctl",
        "50th pctl",
        "75th pctl",
        "90th pctl",
        # "max"
    ],
    orientation = :horizontal
    )

    return figure

end