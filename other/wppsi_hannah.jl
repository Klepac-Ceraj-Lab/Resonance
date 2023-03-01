using Resonance

mdata = Resonance.load_raw_metadata()
taxa = Resonance.load_raw_metaphlan()

#-

wide = comm2wide(filter(t-> taxrank(t) == :species, taxa))

#-

smeta = DataFrame(get(taxa))
rename!(smeta, "sample_base"=> "omni")
smeta = leftjoin(smeta, subset(mdata, "omni"=> ByRow(!ismissing)); on="omni")

smeta.MARRIED = map(eachrow(smeta)) do row
    old = row."Marital_Status_old"
    mom = row."Marital_Status_Mother"
    ismissing(old) && ismissing(mom) && return missing
    !ismissing(old) && old == "Married" && return true
    !ismissing(mom) && occursin(r"^Married", mom) && return true
    return false
end

smeta.PAROUS = map(smeta."BasicFamilyAndChild::siblingInStudy") do sib
    ismissing(sib) && return false
    return occursin("yes", lowercase(sib))
end

smeta.kidbday = smeta."BasicFamilyAndChild::childBirthday"
# smeta.mombday = [ismissing(bday) ? missing : Date(bday, "m/d/y") for bday in smeta."BasicFamilyAndChild::motherBirthday"]
smeta.dadbday = smeta."BasicFamilyAndChild::fatherBirthday" 

smeta.mombday = smeta."BasicFamilyAndChild::motherBirthday"

smeta.MAGE = map(eachrow(smeta)) do row
    mage = row.kidbday - row.mombday
    ismissing(mage) && return missing
    return mage.value / 365
end

smeta.FAGE = map(eachrow(smeta)) do row
    fage = row.kidbday - row.dadbday
    ismissing(fage) && return missing
    return fage.value / 365
end

smeta.SMOKE = [ismissing(smoke) ? missing : smoke=="Yes" for smoke in smeta."Prenatal::smoke"]

smeta.GESTAGE = smeta."NewbornInfo::childGestationalPeriodWeeks"

# select!(smeta, r"^(subject|education|BreastfeedingDone|MARRIED|PAROUS|MAGE|FAGE|SMOKE|GESTAGE)")

smeta.AGE = map(m-> m / 12 * 365, smeta.ageMonths)

smeta.SCORE = map(s-> ismissing(s) ? missing : parse(Float64, s), smeta.cogScore)

smeta.MED = smeta.education

smeta.SEX = smeta.sex

smeta.BFEEDDUR = map(eachrow(smeta)) do row
    nolonger = row."BreastfeedingDone::noLongerFeedBreastmilkAge"
    bfperc = row."BreastfeedingStill::breastFedPercent"
    if !ismissing(nolonger)
        return nolonger
    elseif !ismissing(bfperc) && bfperc > 0
        return row.assessmentAgeMonths
    else
        return missing
    end
end

smeta.BFEED = map(eachrow(smeta)) do row
    ismissing(row.BFEEDDUR) ? missing : row.BFEEDDUR > 1.5
end


## for now, fill these with random variables

smeta.ABX = fill(false, nrow(smeta)) # TODO: add correct value


smeta.STOOLAGE = map(eachrow(smeta)) do row
    (row.ageMonths / 12 * 52 + (40 - row.GESTAGE)) / 52 * 365 # correctedAge + weeks early in days
end

smeta.AGE = map(x-> ismissing(x) ? missing : floor(Int, x), smeta.AGE)
smeta.STOOLAGE = map(x-> ismissing(x) ? missing : floor(Int, x), smeta.STOOLAGE)

smeta.AGEGROUP = map(smeta.STOOLAGE) do age
    age = age / 12
    ismissing(age) && return missing
    age < 0.5 && return 1
    1 < age < 2 && return 2
    age >= 2 && return 4
    return 3
end

rename!(smeta, "birthType" => "BMODE")

select!(smeta, r"^(subject|omni|cogAssessment|MARRIED|PAROUS|MAGE|FAGE|SMOKE|GESTAGE|AGE|SCORE|MED|SEX|BFEEDDUR|BFEED|ABX|STOOLAGE|AGEGROUP|BMODE)")
unique!(smeta, "omni")
rename!(smeta, "omni"=>"sample_base")

joined = unique(leftjoin(smeta, select(wide, Not(["sample", "file", "read_depth"])); on="sample_base"), "sample_base")
subset!(joined, "subject"=> ByRow(!ismissing), "SCORE"=> ByRow(!ismissing))

## write files

CSV.write(scratchfiles("hannah_wppsi.csv"), joined)


## Plot

using CairoMakie

fig = Figure()
ax = Axis(fig[1,1])

hist(collect(skipmissing(meta.STOOLAGE)))