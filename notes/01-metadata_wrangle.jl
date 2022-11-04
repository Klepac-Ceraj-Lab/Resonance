# # Metadata wrangling

# ## Get sample-specific metadata from airtable database
#
# This section requires the presence of input data,
# the location of which should be encoded in environmental variables:
#
# - `"AIRTABLE_KEY"`: an airtable API key
# - `"SCRATCH_SPACE"`: The location for non-persistant data. Eg `/brewster/kevin/resonance_scratch`
# - `"ANALYSIS_FILES"`: The location of biobakery analysis files (metaphlan / humann profiles). Eg `/grace/sequencing/processed/mgx`
# - `"DATA_FILES"`: The location of non-biobakery inputs. Eg. `/brewster/kevin/resonance_data`
#
# I set these variables in a `.envrc` file that can be loaded by `direnv`
#
# ```sh
# layout julia
#
# export AIRTABLE_KEY="<secret>"
# PATH_add "/home/kevin/.julia/conda/3/envs/BiobakeryUtils/bin"
#
# export SCRATCH_SPACE="/brewster/kevin/resonance_scratch"
# export ANALYSIS_FILES="/grace/sequencing/processed/mgx"
# export DATA_FILES="/brewster/kevin/resonance_data"
# ```

using Resonance

perform_airtable_write = false

samplemeta = airtable_metadata() # having set ENV["AIRTABLE_KEY"]

#-
# Now to add important metadata to samples

@chain samplemeta begin
    groupby([:subject, :timepoint])
    transform!(:Metabolomics_batch => (x-> any(!ismissing, x)) => :has_metabolomics)
    sort!([:subject, :timepoint])
end


fmp_samples = DataFrame(XLSX.readtable(datafiles("resonance_fmp", "Sample_Centric_071322.xlsx"), "Sheet1", infer_eltypes=true))
rename!(fmp_samples, Dict(:studyID=>:subject, :collectionNum=> :timepoint))
subset!(fmp_samples, :subject => ByRow(s-> s in samplemeta.subject),
                     :timepoint => ByRow(!ismissing)
)

samplemeta = leftjoin(samplemeta, select(fmp_samples, [:subject, :timepoint, :collectionDate]), on=[:subject, :timepoint])

unique!(samplemeta)

# ## Metadata stored in filemaker pro database
#
# Read in the different tables as `DataFrame`s,
# then normalize certain columns.
# First, subject-specific data

fmp_subject = DataFrame(XLSX.readtable(datafiles("resonance_fmp", "Subject_Centric_071322.xlsx"), "Sheet1", infer_eltypes=true))
rename!(fmp_subject, Dict(:studyID=>:subject))
subset!(fmp_subject, :subject => ByRow(s-> s in samplemeta.subject))

# Then, timepoint-specific data

fmp_timepoint = DataFrame(XLSX.readtable(datafiles("resonance_fmp", "Timepoint_Centric_071322.xlsx"), "Sheet1", infer_eltypes=true))
rename!(fmp_timepoint, Dict(:studyID=>:subject))

fmp_timepoint = let cs_perc = CSV.read(datafiles("cogscore_percentiles.csv"), DataFrame)
    leftjoin(fmp_timepoint, select(cs_perc, "subject", "timepoint", "cogScorePercentile"), on=["subject", "timepoint"])
end

# and COVID-specific samples

fmp_covid = DataFrame(XLSX.readtable(datafiles("resonance_fmp", "COVID_Fecal_040122.xlsx"), "Sheet1", infer_eltypes=true))
rename!(fmp_covid, Dict(:studyID=>:subject))
fmp_covid = leftjoin(fmp_covid, fmp_subject, on=:subject)

# ## Getting data joined together
#
# The next stemp is to merge all of the tables.

subset!(fmp_timepoint, "subject" => ByRow(!ismissing))
fmp_alltp = leftjoin(fmp_timepoint, fmp_subject, on=[:subject])


codebreastfeeding!(fmp_alltp)

count(!ismissing, fmp_alltp.bfcalculated)

# ## Dealing with ages

fmp_alltp.ageMonths = map(eachrow(fmp_alltp)) do row
    if ismissing(row.scanAgeMonths)
        (row.assessmentAgeDays ./ 365 .* 12) .+ row.assessmentAgeMonths
    else
        (row.scanAgeDays ./ 365 .* 12) .+ row.scanAgeMonths
    end
end

# Add info about brain data (just if it's there)

brain = XLSX.readxlsx(datafiles("brain", "AllResonanceVols_Resonance_only.xlsx"))
volumes = let 
    tb = brain["RESONANCE"]["A2:CW1181"]

    df = DataFrame(tb[2:end, :], string.(vec(tb[1, :])))
    contains(brain["RESONANCE"]["CZ2"], "CSF") || error("CSF col missing")
    contains(brain["RESONANCE"]["DA2"], "GM") || error("GM col missing")
    contains(brain["RESONANCE"]["DB2"], "WM") || error("GM col missing")

    df."CSF (ignore)" = vec(brain["RESONANCE"]["CZ3:CZ1181"])
    df."Gray-matter" = vec(brain["RESONANCE"]["DA3:DA1181"])
    df."White-matter" = vec(brain["RESONANCE"]["DB3:DB1181"])

    rename!(df, "missing" => "Sex")
    for col in ["Cohort", "ID", "Session", "Sex"]
        df[!, col] = String.(df[!, col])
    end
    df.AgeInDays = Float64[ismissing(x) ? missing : x for x in df.AgeInDays]

    for (cidx, row) in zip(vec(brain["Tissue Labels"]["B2:B97"]), eachrow(brain["Tissue Labels"]["C2:F97"]))
        lab = join(strip.(coalesce.(row, ""), Ref(['“', '”'])), ' ')
        lab = replace(lab, " missing"=> "", r" $"=> "", " "=>"-")
        lab = strip(lab, '-')
        cidx isa Number && rename!(df, string(cidx) => lab)
        df[!, lab] = Float64[ismissing(x) ? 0. : x for x in df[!, lab]]
    end


    df.subject = map(df.ID) do id
        m = match(r"^sub-BAMBAM(\d+)$", id)
        isnothing(m) && error("$id is wrong format")
        parse(Int, m.captures[1])
    end

    df.timepoint = map(df.Session) do s
        m = match(r"^ses-(\d+)$", s)
        isnothing(m) && error("$s is wrong format")
        parse(Int, m.captures[1])
    end
    df
end

volumes.has_segmentation .= true

## Validation

@assert all(volumes.Cohort .== "RESONANCE")

fmp_alltp = leftjoin(fmp_alltp, select(volumes, [:subject, :timepoint, :AgeInDays, :Sex, :has_segmentation]), on=[:subject, :timepoint])
fmp_alltp.has_segmentation = [ismissing(x) ? false : x for x in fmp_alltp.has_segmentation]


## transform!(fmp_alltp, :sid_old => ByRow(id-> ismissing(id) ? id : replace(id, r"_(\d+)F_"=>s"_\1E_")) => :sid_old_etoh)
## etoh_map = Dict((old=>new for (old, new) in zip(samplemeta.sid_old, samplemeta.sample)))
## transform!(fmp_alltp, :sid_old_etoh => ByRow(id-> (ismissing(id) || !haskey(etoh_map, id)) ? missing : etoh_map[id]) => :sample_etoh)

sort!(fmp_alltp, [:subject, :timepoint])

# ## Getting values for presence of data types
#
# We want to plot set intersections for having various kinds of data.
# In some cases, additional wrangling is necessary

fmp_alltp.race = map(fmp_alltp."Merge_Dem_Child_Race") do r
    ismissing(r) && return missing
    r == "Unknown" && return missing
    r == "Decline to Answer" && return missing
    contains(r, "\n") && return "Mixed"
    r ∈ ("Mixed", "Mixed Race") && return "Mixed"
    r ∈ ("Other Asian", "Asian ") && return "Asian"
    return r
end

## unique!(fmp_alltp, [:subject, :timepoint])

fmp_alltp.has_race = .!ismissing.(fmp_alltp."race")

#-


subj = groupby(fmp_alltp, :subject)
transform!(subj, nrow => :n_samples; ungroup=false)

transform!(subj, AsTable(r"subject|breast"i) => (s -> begin
    if any(!ismissing, s.breastFedPercent)
        return fill(true, length(s[1]))
    else
        return fill(false, length(s[1]))
    end
    fill(0, length(s[1]))
end) => :has_bfperc_subj; ungroup=false)

transform!(subj, AsTable(r"subject|breast"i) => (s -> begin
    if any(!ismissing, s.breastFedPercent)
        return fill(true, length(s[1]))
    else
        return fill(false, length(s[1]))
    end
    fill(0, length(s[1]))
end) => :has_bfperc_subj; ungroup=false)

## fmp_alltp = DataFrames.transform(subj, AsTable([:timepoint, :sample]) => (s -> begin
##     has_stool = .!ismissing.(s.sample)
##     has_prevstool = fill(false, length(s[1]))
##     for i in 1:length(s[1])
##         i == 1 && continue
##         any(!ismissing, s.sample[1:i-1]) && (has_prevstool[i] = true)
##     end
##     (; has_stool, has_prevstool)
## end)=> [:has_stool, :has_prevstool])

fmp_alltp.has_everbreast = .!ismissing.(fmp_alltp."everBreastFed")
fmp_alltp.has_bfperc = .!ismissing.(fmp_alltp."breastFedPercent")



for n in names(fmp_alltp)
    fmp_alltp[!, n] = map(fmp_alltp[!, n]) do v
        v isa AbstractString ? replace(v, r"\n"=>"___") : v
    end
end

CSV.write(datafiles("wrangled", "timepoints.csv"), fmp_alltp)
CSV.write(datafiles("wrangled", "omnisamples.csv"), subset(samplemeta, "Fecal_EtOH" => ByRow(==("F")), "Mgx_batch"=> ByRow(!ismissing)))
CSV.write(datafiles("wrangled", "etohsamples.csv"), subset(samplemeta, "Fecal_EtOH" => ByRow(==("E")), "Metabolomics_batch"=> ByRow(!ismissing)))
CSV.write(datafiles("wrangled", "brain.csv"), select(volumes, Cols(:subject, :timepoint, Not(["Cohort", "ID", "Session"]))))
CSV.write(datafiles("wrangled", "covid.csv"), fmp_covid)

# ## Uploading to Airtable
#
# Some data is useful to feed back up to airtable.
# This requires an airtable API key that has write access.

if perform_airtable_write

    using Airtable

    base = AirBase("appSWOVVdqAi5aT5u")
    tab = AirTable("Samples", base)
    joinedsamples = leftjoin(select(samplemeta, [:airtable_id, :subject, :timepoint, :sample, :childAgeMonths]), 
                             select(fmp_alltp, [:subject, :timepoint, :ageMonths, :has_segmentation]), 
                             on=[:subject, :timepoint]
    )
    
    agepatches = AirRecord[]
    brainpatches = AirRecord[]
    
    for row in eachrow(joinedsamples)
        if !ismissing(row[:ageMonths]) && ismissing(row[:childAgeMonths])
            push!(agepatches, AirRecord(row[:airtable_id], tab, (; childAgeMonths=row[:ageMonths])))
        end
        if coalesce(row[:has_segmentation], false)
            push!(brainpatches, AirRecord(row[:airtable_id], tab, (; concurrent_scan=row[:has_segmentation])))
        end
    end
    
    !isempty(agepatches) && Airtable.patch!(tab, agepatches)
    
    unique!(brainpatches)
    !isempty(brainpatches) && Airtable.patch!(tab, brainpatches)

end