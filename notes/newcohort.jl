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
    
    transform(AsTable(["omni", "cogScore"]) => (nt-> any(.!ismissing.(nt.omni) .& .!ismissing.(nt.cogScore))) => "has_concurrent_stool_cog"; ungroup=false)

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
        i = findfirst(i-> !any(ismissing, (nt.ageMonths[i], nt.omni[i], nt.cogScore[i])) && 18 < nt.ageMonths[i] <= 120, eachindex(nt.ageMonths))
        !isnothing(i) && (isfirst[i] = true)
        isfirst
    end) => "filter_18to120";
    ungroup = false
    )

    transform(AsTable(["ageMonths", "omni", "cogScore"]) => (nt-> begin
        future_omni = fill(false, length(nt.ageMonths))
        future_cog = fill(false, length(nt.ageMonths))
        nm_omni = findall(!ismissing, nt.omni)
        nm_cog = findall(!ismissing, nt.cogScore)

        for i in reverse(nm_omni), j in reverse(nm_cog)
            nt.ageMonths[i] <= 12 || continue
            if i < j && ((nt.ageMonths[j] - nt.ageMonths[i]) > 3)
                future_omni[i] = true
                future_cog[j] = true
                break
            end
        end
        return (; future_omni, future_cog)
    end) => ["filter_future_omni", "filter_future_cog"])
end

CSV.write("input/complete_filtered_dataset.csv", kids00to120)

@info "Kids < 120 months: size = $(size(kids00to120))"
@info "Kids < 120 months: n subjects = $(length(unique(kids00to120.subject)))"
@info "With concurrent stool / cog: n subjects = $(subset(kids00to120, "has_concurrent_stool_cog" => identity).subject |> unique |> length)"
@info "With concurrent stool / cog: n samples = $(size(subset(kids00to120, AsTable(["omni", "cogScore"]) => ByRow(row-> !any(ismissing, row))),1))"

 
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
    @aside let nsubjects = length(_.subject |> unique)
        @info "Subjects with 1 sample: $(count(==(1), _.nsamples)) ($(round(count(==(1), _.nsamples) / nsubjects * 100; digits=2))%)"
        @info "Subjects with 2 samples: $(count(==(2), _.nsamples)) ($(round(count(==(2), _.nsamples) / nsubjects * 100; digits=2))%)"
        @info "Subjects with 3 or more samples: $(count(>=(3), _.nsamples)) ($(round(count(>=(3), _.nsamples) / nsubjects * 100; digits=2))%)"
    end
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
    
