# functions for going from raw data to datatypes


function load_raw_metadata(;
    subject_centric_path = datafiles("resonance_fmp", "Subject_Centric_071322.xlsx"),
    sample_centric_path  = datafiles("resonance_fmp", "Sample_Centric_071322.xlsx"),
    timepoint_centric_path = datafiles("resonance_fmp", "Timepoint_Centric_071322.xlsx")
    )

    samplemeta = @chain airtable_metadata() begin
        groupby(["subject", "timepoint"])
        transform("Metabolomics_batch" => (x-> any(!ismissing, x)) => "has_metabolomics")
        sort(["subject", "timepoint"])
        unique("sample")
    end

    fmp_samples = @chain DataFrame(XLSX.readtable(sample_centric_path, "Sheet1", infer_eltypes=true)) begin
        rename!(Dict("studyID"=> "subject", "collectionNum"=> "timepoint", "SampleID"=>"sample"))
        subset!("subject"   => ByRow(s-> s in samplemeta.subject),
                "timepoint" => ByRow(!ismissing))
    end

    leftjoin!(samplemeta, select(fmp_samples, ["sample", "collectionDate"]), on="sample")
    @assert nrow(samplemeta) == nrow(unique(samplemeta, "sample"))

    fmp_subject = @chain DataFrame(XLSX.readtable(subject_centric_path, "Sheet1", infer_eltypes=true)) begin
        rename(Dict("studyID"=> "subject"))
        subset("subject" => ByRow(s-> s in samplemeta.subject))
    end 

    fmp_timepoint = DataFrame(XLSX.readtable(timepoint_centric_path, "Sheet1", infer_eltypes=true))
    rename!(fmp_timepoint, Dict("studyID"=> "subject"))
    cs_perc = CSV.read(datafiles("cogscore_percentiles.csv"), DataFrame)
    leftjoin!(fmp_timepoint, select(cs_perc, "subject", "timepoint", "cogScorePercentile"), on=["subject", "timepoint"])

    fmp_alltp = leftjoin(fmp_timepoint, fmp_subject, on=["subject"])

    fmp_alltp.ageMonths = map(eachrow(fmp_alltp)) do row
        if ismissing(row.scanAgeMonths)
            (row.assessmentAgeDays ./ 365 .* 12) .+ row.assessmentAgeMonths
        else
            (row.scanAgeDays ./ 365 .* 12) .+ row.scanAgeMonths
        end
    end

    @assert nrow(unique(fmp_alltp, ["subject", "timepoint"])) == nrow(fmp_alltp)

    DataFrames.transform!(groupby(fmp_alltp, :subject), :mother_HHS_Education => (r->coalesce(r...)) => :hhs)
    fmp_alltp.education = let
        ed = categorical(fmp_alltp.hhs; levels=[-8 , 2:7...], ordered=true)
        ed = recode(ed,
        -8 => missing,
        2 => "Junior high school",
        3 => "Some high school",
        4 => "High school grad",
        5 => "Some college",
        6 => "College grad",
        7 => "Grad/professional school")
        ed
    end

    DataFrames.transform!(groupby(fmp_alltp, :subject), :race => (r->coalesce(r...)) => :race)

    fmp_alltp.race = map(fmp_alltp."Merge_Dem_Child_Race") do r
        ismissing(r) && return missing
        r == "Unknown" && return missing
        r == "Decline to Answer" && return missing
        contains(r, "\n") && return "Mixed"
        r ∈ ("Mixed", "Mixed Race") && return "Mixed"
        r ∈ ("Other Asian", "Asian ") && return "Asian"
        return r
    end

    fmp_alltp.race = let
        race = categorical(fmp_alltp.race; ordered=true)
        race = recode(race, 
            "American Indian or Alaska Native"=> "Indiginous",
            "Some other race"                 => "Other",
            "Asian Indian"                    => "Asian",
            "Black or African American"       => "Black",
            missing                           => "Unknown"
        )
        droplevels!(race)
        levels!(race, ["White","Black","Asian","Mixed","Other","Unknown"])
        race
    end

    omni = @chain samplemeta begin
        subset("Fecal_EtOH" => ByRow(==("F")), "Mgx_batch"=> ByRow(!ismissing))
        select("subject", "timepoint", "sample")
        unique(["subject", "timepoint"])
        rename("sample"=> "omni")
    end
    etoh = @chain samplemeta begin
        subset("Fecal_EtOH" => ByRow(==("E")), "Metabolomics_batch"=> ByRow(!ismissing))
        select("subject", "timepoint", "sample")
        unique(["subject", "timepoint"])
        rename("sample"=> "etoh")
    end

    leftjoin!(fmp_alltp, omni; on=["subject", "timepoint"])
    leftjoin!(fmp_alltp, etoh; on=["subject", "timepoint"])
    return fmp_alltp
end


function load_raw_metaphlan()

end

# using Resonance

# metaphlan_files = filter(readdir(analysisfiles("metaphlan"), join=true)) do f
#     !contains(f, r"FE\d{5}") && contains(f, "profile.tsv")
# end

# #- 

