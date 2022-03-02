include("startup.jl")

for col in names(tps)
    tps[!, col] = collect(tps[!, col])
end

tps_specific = stack(tps, ["ageMonths", "cogScore", "ECHOTPCoded", brainmeta...], [:subject, :timepoint])
subj_specific = stack(tps, [:mother_HHS_Education, :simple_race], [:subject])
grp = groupby(subj_specific, [:subject, :variable])
subj_specific = DataFrames.combine(grp, :value=> (v-> coalesce(v...))=> :value)

##

tps_specific = vcat(tps_specific, DataFrame(subject   = get(species, :subject), 
                                timepoint = get(species, :timepoint),
                                variable  = fill("sample", nsamples(species)),
                                value     = samplenames(species))
)

for s in samples(species)
    subject   = s.subject
    timepoint = s.timepoint
    for f in features(species)
        variable = name(f)
        value = species[f, s]
        push!(tps_specific, (; subject, timepoint, variable, value))
    end
end

CSV.write("data/wrangled/tidy_timepoints_with_brain.csv", tps_specific)
CSV.write("data/wrangled/tidy_subjects.csv", subj_specific)


##

