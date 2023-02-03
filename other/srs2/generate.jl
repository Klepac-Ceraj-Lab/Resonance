# ---
# title = "VKC Lab data generation"
# author = "Kevin Bonham, PhD"
# ---
# 
# # VCK Lab data generation
#
# Here, I'm conforming the taxonomic profiles from the VCK lab ECHO cohort
# to fit the format needed to plug in to Hannah's script.


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

taxa = Resonance.load_raw_metaphlan()
mdata = Resonance.load_raw_metadata()

species = filter(f-> taxrank(f) == :species, taxa)


## metadata

mdata.MARRIED = map(eachrow(mdata)) do row
    old = row."Marital_Status_old"
    mom = row."Marital_Status_Mother"
    ismissing(old) && ismissing(mom) && return missing
    !ismissing(old) && old == "Married" && return true
    !ismissing(mom) && occursin(r"^Married", mom) && return true
    return false
end

mdata.PAROUS = map(mdata."BasicFamilyAndChild::siblingInStudy") do sib
    ismissing(sib) && return false
    return occursin("yes", lowercase(sib))
end

mdata.kidbday = mdata."BasicFamilyAndChild::childBirthday"
mdata.mombday = mdata."BasicFamilyAndChild::motherBirthday"
mdata.dadbday = mdata."BasicFamilyAndChild::fatherBirthday"

mdata.MAGE = map(eachrow(mdata)) do row
    mage = row.kidbday - row.mombday
    ismissing(mage) && return missing
    return mage.value / 365
end

mdata.FAGE = map(eachrow(mdata)) do row
    fage = row.kidbday - row.dadbday
    ismissing(fage) && return missing
    return fage.value / 365
end

mdata.SMOKE = [ismissing(smoke) ? missing : smoke=="Yes" for smoke in mdata."Prenatal::smoke"]

mdata.GESTAGE = mdata."NewbornInfo::childGestationalPeriodWeeks"

mdata.SRS2_T = map(t-> coalesce(t...), zip(mdata."SchoolageSRS::schoolAgeTotalTScore", mdata."PreschoolSRS::preschoolTotalTScore"))
subset!(mdata, "SRS2_T"=> ByRow(!ismissing))
subset!(mdata, "ageMonths"=> ByRow(!ismissing))
subset!(mdata, "omni"=> ByRow(!ismissing))
unique!(mdata, "omni")

sort!(mdata.ageMonths)

transform!(mdata, "ageMonths"=> ByRow(x-> x / 12) =>"AGE")
@chain mdata begin
    groupby("subject")
    transform!(["AGE", "SRS2_T"] => ( (age, srs2) -> [
        all(ismissing, srs2) ? (missing, missing) : 
               ismissing(sc) ? (age[findfirst(!ismissing, srs2)], 
                                srs2[findfirst(!ismissing, srs2)]) : 
               (age[i], srs2[i]) for (i, sc) in enumerate(srs2)
        ] ) => [:AGE, :SCORE])
end

mdata.MED = map(collect(mdata.mother_HHS_Education)) do hhs
    ismissing(hhs) && return missing
    hhs >= 6 # higher ed or not
end

mdata.SEX = map(collect(mdata.sex)) do s
    ismissing(s) && return missing
    s == "Male" && return true
    s == "Female" && return false
    error("unknown sex $s")
end

mdata.BFEEDDUR = map(eachrow(mdata)) do row
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

mdata.BFEED = map(eachrow(mdata)) do row
    ismissing(row.BFEEDDUR) ? missing : row.BFEEDDUR > 1.5
end


## for now, fill these with random variables

mdata.ABX .= fill(false, nrow(mdata)) # TODO: add correct value


mdata.STOOLAGE = map(eachrow(mdata)) do row
    (row.ageMonths / 12 * 52 + (40 - row.GESTAGE)) / 52 * 365 # correctedAge + weeks early in days
end

mdata.AGE = map(x-> ismissing(x) ? missing : floor(Int, x), mdata.AGE)
mdata.STOOLAGE = map(x-> ismissing(x) ? missing : floor(Int, x), mdata.STOOLAGE)

mdata.AGEGROUP = map(mdata.STOOLAGE) do age
    age = age / 12
    ismissing(age) && return missing
    age < 0.5 && return 1
    1 < age < 2 && return 2
    age >= 2 && return 4
    return 3
end

rename!(mdata, "birthType" => "BMODE")

function generate_taxtab()
    taxa = Set(String[])
    for t in filter(f-> contains(basename(f), "profile"), readdir(analysisfiles("metaphlan"), join=true))
        df = CSV.read(t, DataFrame; comment="#", header=["taxon", "ncbi", "abundance", "other"])
        union!(taxa, subset(df, "taxon"=>ByRow(tax-> contains(tax, "s__") && !contains(tax, "t__"))).taxon)
    end
    df = DataFrame([NamedTuple(k=> v for (v, k) in zip(split(replace(taxstring, r"\[|\]"=>""), "|"), [:kingdom, :phylum, :class, :order, :family, :genus, :species])) for taxstring in taxa])
    sort!(df, :species)
    df.taxon = replace.(df.species, "s__"=>"")
    select(df, Cols(:taxon, :))
end

taxtab = generate_taxtab()
spec_subset = let sub = species[:, [s in mdata.omni for s in get(species, :sample_base)]]
    sub = sub[vec(prevalence(sub)) .> 0, :]
end
asvtab = comm2wide(spec_subset)

subset!(mdata, "omni"=> ByRow(o-> o in get(spec_subset, :sample_base)))# remove samples that we don't have profile for
unique!(mdata, "subject")
mdata.sample = mdata.omni
subset!(taxtab, "taxon"=> ByRow(t-> !(t in featurenames(spec_subset))))
## write files

CSV.write("/brewster/kevin/resonance_data/exports/VKC_TAXTAB.csv", taxtab)
CSV.write("/brewster/kevin/resonance_data/exports/VKC_ASVTAB.csv",  comm2wide(species))
CSV.write("/brewster/kevin/resonance_data/exports/VCK_METATAB.csv", select(mdata, [:sample, :subject, :SCORE, :BMODE, :AGE, :MED, :MARRIED, :PAROUS, :MAGE, :FAGE, :SMOKE, :ABX, :SEX, :GESTAGE, :STOOLAGE, :BFEEDDUR, :BFEED]))
