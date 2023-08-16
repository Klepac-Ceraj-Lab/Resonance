using Resonance
using VKCComputing
using CategoricalArrays
using MiniLoggers
using LoggingExtras

global_logger(TeeLogger(MiniLogger(; message_mode=:markdown), MiniLogger(; io="cohort_generation.log")))

newmeta = CSV.read("input/fmp_alltp.csv", DataFrame)

@info "Input DF (`fmp_alltp.csv`): size = $(size(newmeta))"
@info "Input DF (`fmp_alltp.csv`): n subjects = $(length(unique(newmeta.subject)))"

kids00to120 = @chain newmeta begin
    sort(["subject", "timepoint"])
    subset("ageMonths" => ByRow(a-> !ismissing(a) && a <= 120))
    groupby("subject")
    subset("omni" => (omni-> any(!ismissing, omni)))
    groupby("subject")
    
    transform(AsTable(["omni", "cogScore"]) => (nt-> any(.!ismissing.(nt.omni) .& any(.!ismissing.(nt.cogScore)))) => "has_concurrent_stool_cog"; ungroup=false)

    transform(AsTable(["ageMonths", "omni", "cogScore"]) => (nt-> begin
        isfirst = fill(false, length(nt.ageMonths))
        i = findfirst(i-> !any(ismissing, (nt.ageMonths[i], nt.omni[i], nt.cogScore[i])), eachindex(nt.ageMonths))
        !isnothing(i) && (isfirst[i] = true)
        isfirst
    end) => "filter_00to120";
    ungroup = false
    )
    transform(AsTable(["ageMonths", "omni", "cogScore"]) => (nt-> begin
        isfirst = fill(false, length(nt.ageMonths))
        i = findfirst(i-> !any(ismissing, (nt.ageMonths[i], nt.omni[i], nt.cogScore[i])) && (nt.ageMonths[i] < 6), eachindex(nt.ageMonths))
        !isnothing(i) && (isfirst[i] = true)
        isfirst
    end) => "filter_00to06";
    ungroup = false
    )
    transform(AsTable(["ageMonths", "omni", "cogScore"]) => (nt-> begin
        isfirst = fill(false, length(nt.ageMonths))
        i = findfirst(i-> !any(ismissing, (nt.ageMonths[i], nt.omni[i], nt.cogScore[i])) && nt.ageMonths[i] <= 120, eachindex(nt.ageMonths))
        !isnothing(i) && (isfirst[i] = true)
        isfirst
    end) => "filter_18to120"
    )
end

@info "Kids < 120 months: size = $(size(kids00to120))"
@info "Kids < 120 months: n subjects = $(length(unique(kids00to120.subject)))"
@info "N in 00to120 months: n samples = $(count(kids00to120.filter_00to120))"

@info "**All visits < 120 months**"

@chain kids00to120 begin
    groupby("subject")
    DataFrames.combine("omni"=> length => "nsamples")
    @aside @info "Subjects with 1 visit: $(count(==(1), _.nsamples))"
    @aside @info "Subjects with 2 visits: $(count(==(2), _.nsamples))"
    @aside @info "Subjects with 3 or more visits: $(count(>=(3), _.nsamples))"
end 

@info "**Subjects with at least 1 concurrent stool / cogScore**"

@chain kids00to120 begin
    subset("has_concurrent_stool_cog"=> identity)
    groupby("subject")
    DataFrames.combine("omni"=> length => "nsamples")
    @aside @info "Subjects with 1 visit: $(count(==(1), _.nsamples))"
    @aside @info "Subjects with 2 visits: $(count(==(2), _.nsamples))"
    @aside @info "Subjects with 3 or more visits: $(count(>=(3), _.nsamples))"
end 

@chain kids00to120 begin
    subset("has_concurrent_stool_cog"=> identity)
    subset("omni"=> ByRow(!ismissing))
    groupby("subject")
    DataFrames.combine("omni"=> length => "nsamples")
    @aside @info "Subjects with 1 sample: $(count(==(1), _.nsamples))"
    @aside @info "Subjects with 2 samples: $(count(==(2), _.nsamples))"
    @aside @info "Subjects with 3 or more samples: $(count(>=(3), _.nsamples))"
end 



transform!(AsTable(["ageMonths", "omni", "cogScore"]) => (nt-> begin
    isfirst = fill(false, length(nt.ageMonths))
    i = findfirst(i-> !any(ismissing, nt.ageMonths[i], nt.omni[i], nt.cogScore[i]), eachindex(nt.ageMonths))
    !isnothing(i) && (isfirst[i] = true)
    isfirst
end) => "filter_00to120";
ungroup = false
)
transform!(AsTable(["ageMonths", "omni", "cogScore"]) => (nt-> begin
    isfirst = fill(false, length(nt.ageMonths))
    i = findfirst(i-> !any(ismissing, nt.ageMonths[i], nt.omni[i], nt.cogScore[i]) && nt.ageMonths[i] < 6, eachindex(nt.ageMonths))
    !isnothing(i) && (isfirst[i] = true)
    isfirst
end) => "filter_00to06";
ungroup = false
)
transform!(AsTable(["ageMonths", "omni", "cogScore"]) => (nt-> begin
    isfirst = fill(false, length(nt.ageMonths))
    i = findfirst(i-> !any(ismissing, nt.ageMonths[i], nt.omni[i], nt.cogScore[i]) && nt.ageMonths[i] <= 120, eachindex(nt.ageMonths))
    !isnothing(i) && (isfirst[i] = true)
    isfirst
end) => "filter_18to120"
)
    
