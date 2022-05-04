# Vanja + Aditi CBCL

using Resonance
omni, etoh, tps, complete_brain, metabolites, species = startup()
genes = Resonance.read_gfs_arrow()

aditi_sids = readlines("data/aditi_samples.txt")

CSV.write(joinpath(ENV["SCRATCH_SPACE"], "aditi_taxa.csv"), species[:, aditi_sids])
CSV.write(joinpath(ENV["SCRATCH_SPACE"], "aditi_genefamilies.csv"), genes[:, aditi_sids])

 

## Korean group, blood pressure

using Resonance

fecalsamples = CSV.read("data/wrangled/samples.csv", DataFrame)
@rsubset! fecalsamples :Fecal_EtOH == "F"

timepoints = CSV.read("data/wrangled/timepoints.csv", DataFrame)

function mapbp(label)
    bp = Dict(
        "E2a" => "Systolic",
        "E2b" => "Diastolic",
        "E2c" => "Heart rate"
    )
    for pair in pairs(bp)
        bp = replace(label, pair)
    end

    # bp = replace(bp)
end


@assert size(subset(timepoints, "ChildAnthropometryOverTwo::E2a1"=> ByRow(!ismissing), "ChildAnthropometryUnderTwo::E2a1" => ByRow(!ismissing)), 1) == 0

timepoints.Systolic_1 = map(eachrow(timepoints)) do row
    coalesce(row."ChildAnthropometryOverTwo::E2a1", row."ChildAnthropometryUnderTwo::E2a1")
end
timepoints.Systolic_2 = map(eachrow(timepoints)) do row
    coalesce(row."ChildAnthropometryOverTwo::E2a2", row."ChildAnthropometryUnderTwo::E2a2")
end
timepoints.Systolic_3 = map(eachrow(timepoints)) do row
    coalesce(row."ChildAnthropometryOverTwo::E2a3", row."ChildAnthropometryUnderTwo::E2a3")
end
timepoints.Diastolic_1 = map(eachrow(timepoints)) do row
    coalesce(row."ChildAnthropometryOverTwo::E2b1", row."ChildAnthropometryUnderTwo::E2b1")
end
timepoints.Diastolic_2 = map(eachrow(timepoints)) do row
    coalesce(row."ChildAnthropometryOverTwo::E2b2", row."ChildAnthropometryUnderTwo::E2b2")
end
timepoints.Diastolic_3 = map(eachrow(timepoints)) do row
    coalesce(row."ChildAnthropometryOverTwo::E2b3", row."ChildAnthropometryUnderTwo::E2b3")
end
timepoints.HeartRate_1 = map(eachrow(timepoints)) do row
    coalesce(row."ChildAnthropometryOverTwo::E2c1", row."ChildAnthropometryUnderTwo::E2c1")
end
timepoints.HeartRate_2 = map(eachrow(timepoints)) do row
    coalesce(row."ChildAnthropometryOverTwo::E2c2", row."ChildAnthropometryUnderTwo::E2c2")
end
timepoints.HeartRate_3 = map(eachrow(timepoints)) do row
    coalesce(row."ChildAnthropometryOverTwo::E2c3", row."ChildAnthropometryUnderTwo::E2c3")
end

df = select(timepoints, ["subject", "timepoint",  "childGender", "ageMonths", "Systolic_1", "Systolic_2", "Systolic_3", "Diastolic_1", "Diastolic_2", "Diastolic_3", "HeartRate_1", "HeartRate_2", "HeartRate_3"])
df = leftjoin(df, select(fecalsamples, 
                    ["subject", "timepoint", "sample", "sid_old", "Mgx_batch"]
              ),
              on=["subject", "timepoint"])

count(groupby(df, :subject)) do grp
    !any(row-> (!ismissing(row.sample) && !ismissing(row.ageMonths) && (6 < row.ageMonths < 18)), eachrow(grp)) && return false
    any(row-> (!ismissing(row.Systolic_1) && !ismissing(row.ageMonths) && (36 < row.ageMonths < 84)), eachrow(grp))
end

count(groupby(df, :subject)) do grp
    !any(row-> (!ismissing(row.sample) && !ismissing(row.ageMonths) && (6 < row.ageMonths < 18)), eachrow(grp)) && return false
    any(row-> !ismissing(row.Systolic_1), eachrow(grp))
end

count(groupby(df, :subject)) do grp
    !any(row-> !ismissing(row.sample), eachrow(grp)) && return false
    any(row-> !ismissing(row.Systolic_1), eachrow(grp))
end

count(groupby(df, :subject)) do grp
    !any(row-> (!ismissing(row.sample) && !ismissing(row.ageMonths) && row.ageMonths < 24), eachrow(grp)) && return false
    any(row-> !ismissing(row.Systolic_1), eachrow(grp))
end

CSV.write("/home/kevin/Desktop/RESONANCE_bp.csv", sort(df, ["subject", "timepoint"]))
