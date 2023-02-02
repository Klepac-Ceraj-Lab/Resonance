# ---
# title = "VKC Lab data generation"
# author = "Kevin Bonham, PhD"
# ---
# 
# # VCK Lab data generation
#
# Here, I'm conforming the taxonomic profiles from the VCK lab ECHO cohort
# to fit the format needed to plug in to Hannah's script.
#
# Source code for functions can be found in the module defined in
# `src/ECHOSRS2.jl`.
# 

@info VERSION # julia version
using Pkg
Pkg.status() # installed packages and versions
@info pwd()

using CSV
using Dates
using Chain
using Resonance

function asv_df(taxdict, tidydf)
    tidydf = copy(tidydf)
    tidydf.taxon = [taxdict[t] for t in tidydf.taxon] 
    asvtab = unstack(unique(tidydf, [:sample, :taxon]), :taxon, :abundance)
    asvtab[!, r"taxon\d+"] = coalesce.(asvtab[!, r"taxon\d+"], 0.0)
    return asvtab
end

tidytax = Resonance.load_raw_metaphlan()
mdata = Resonance.load_raw_metadata()


filter!("taxon"=> (t-> occursin("s__", t)), tidytax) # only include species

taxtab, taxdict = taxa_df(tidytax.taxon)

##

asvtab = asv_df(taxdict, tidytax)
asvtab.sample = map(s-> replace(replace(s, r"_S\d+"=>""), "-"=>"_"), asvtab.sample)

## metadata

using Distributions
using ResonanceMicrobiome

samp = ResonanceMicrobiome.airtable_metadata()
rename!(samp, "subject"=>"studyID")

tp = CSV.read("data/Timepoint_centric_041421.csv", DataFrame)
filter!(:timepoint=> t-> floor(Int, t) == t, tp) # why the F do we still have 2.5 in here?
tp.timepoint = Int.(tp.timepoint)

meta = leftjoin(samp, tp, on=["studyID", "timepoint"], makeunique=true)

subj = CSV.read("data/Subject_Specific_050721.csv", DataFrame)
meta = leftjoin(meta, subj, on="studyID")

redcap = CSV.read("data/SRS_Redcap_07192021.csv", DataFrame)
filter!(:studyID=>!ismissing, redcap)

subj.MARRIED = map(eachrow(subj)) do row
    old = row."Marital_Status_old"
    mom = row."Marital_Status_Mother"
    ismissing(old) && ismissing(mom) && return missing
    !ismissing(old) && old == "Married" && return true
    !ismissing(mom) && occursin(r"^Married", mom) && return true
    return false
end

subj.PAROUS = map(subj."BasicFamilyAndChild::siblingInStudy") do sib
    ismissing(sib) && return false
    return occursin("yes", lowercase(sib))
end

subj.kidbday = [ismissing(bday) ? missing : Date(bday, "m/d/y") for bday in subj."BasicFamilyAndChild::childBirthday"]
# subj.mombday = [ismissing(bday) ? missing : Date(bday, "m/d/y") for bday in subj."BasicFamilyAndChild::motherBirthday"]
subj.dadbday = [ismissing(bday) ? missing : Date(bday, "m/d/y") for bday in subj."BasicFamilyAndChild::fatherBirthday"]

subj.mombday = map(collect(subj."BasicFamilyAndChild::motherBirthday")) do bday
    ismissing(bday) && return missing
    try
        Date(bday, "m/d/y")
    catch e
        missing
    end
end

subj.MAGE = map(eachrow(subj)) do row
    mage = row.kidbday - row.mombday
    ismissing(mage) && return missing
    return mage.value / 365
end

subj.FAGE = map(eachrow(subj)) do row
    fage = row.kidbday - row.dadbday
    ismissing(fage) && return missing
    return fage.value / 365
end

subj.SMOKE = [ismissing(smoke) ? missing : smoke=="Yes" for smoke in subj."Prenatal::smoke"]

subj.GESTAGE = subj."NewbornInfo::childGestationalPeriodWeeks"

select!(subj, r"^(studyID|mother_HHS_Education|BreastfeedingDone|MARRIED|PAROUS|MAGE|FAGE|SMOKE|GESTAGE)")


rename!(subj, "studyID" => "subject")
meta = leftjoin(meta, subj, on="subject")

tp = CSV.read("data/Timepoint_centric_041421.csv", DataFrame)
select!(tp, r"^(studyID|timepoint|BreastfeedingStill|assessmentAge)")
rename!(tp, "studyID" => "subject")
meta = leftjoin(meta, tp, on=["subject", "timepoint"])

unique!(meta, "sample")
# filter!("SRS2_T"=> !ismissing, meta)
filter!("correctedAgeDays" => !ismissing, meta)

sort!(meta.correctedAgeDays)

meta.AGE = map(row-> (row.assessmentAgeMonths / 12 + row.assessmentAgeDays / 365) * 365, eachrow(meta))
meta = @chain meta begin
    groupby("subject")
    transform(["AGE", "SRS2_T"] => ( (age, srs2) -> [
        all(ismissing, srs2) ? (missing, missing) : 
               ismissing(sc) ? (age[findfirst(!ismissing, srs2)], 
                                srs2[findfirst(!ismissing, srs2)]) : 
               (age[i], srs2[i]) for (i, sc) in enumerate(srs2)
        ] ) => [:AGE, :SCORE])
end

meta.MED = map(collect(meta.mother_HHS)) do hhs
    ismissing(hhs) && return missing
    hhs >= 6 # higher ed or not
end

meta.SEX = map(collect(meta.childSex)) do s
    ismissing(s) && return missing
    lowercase(s) == "male" && return true
    lowercase(s) == "female" && return false
    error("unknown sex $s")
end

meta.BFEEDDUR = map(eachrow(meta)) do row
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

meta.BFEED = map(eachrow(meta)) do row
    ismissing(row.BFEEDDUR) ? missing : row.BFEEDDUR > 1.5
end


## for now, fill these with random variables

meta.ABX = fill(false, nrow(meta)) # TODO: add correct value


meta.STOOLAGE = map(eachrow(meta)) do row
    (row.correctedAgeDays / 365 * 52 + (40 - row.GESTAGE)) / 52 * 365 # correctedAge + weeks early in days
end

meta.AGE = map(x-> ismissing(x) ? missing : floor(Int, x), meta.AGE)
meta.STOOLAGE = map(x-> ismissing(x) ? missing : floor(Int, x), meta.STOOLAGE)

meta.AGEGROUP = map(meta.STOOLAGE) do age
    age = age / 12
    ismissing(age) && return missing
    age < 0.5 && return 1
    1 < age < 2 && return 2
    age >= 2 && return 4
    return 3
end

rename!(meta, "birthType" => "BMODE")



filter!("sample"=> (s-> s in asvtab.sample), meta) # remove samples that we don't have metadata for
filter!("sample"=> (s-> s in meta.sample), asvtab) # remoce samples that we don't have profile for
asvtab = leftjoin(select(meta, "sample"), asvtab, on="sample") # get rows in same order

## write files

CSV.write("data/VKC_TAXTAB.csv", taxtab)
CSV.write("data/VKC_ASVTAB.csv", asvtab)
CSV.write("data/VCK_METATAB.csv", select(meta, [:sample, :SCORE, :BMODE, :AGE, :MED, :MARRIED, :PAROUS, :MAGE, :FAGE, :SMOKE, :ABX, :SEX, :GESTAGE, :STOOLAGE, :BFEEDDUR, :BFEED]))


## Plot

using CairoMakie

fig = Figure()
ax = Axis(fig[1,1])

hist(collect(skipmissing(meta.STOOLAGE)))