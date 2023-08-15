# functions for going from raw data to datatypes


function load_raw_metadata(;
    subject_centric_path = datafiles("Participant_Centric_081023.xlsx"),
    sample_centric_path  = datafiles("Sample_Centric_081023.xlsx"),
    timepoint_centric_path = datafiles("Timepoint_Centric_081023.xlsx"))

    samplemeta = @chain airtable_metadata() begin
        groupby(["subject", "timepoint"])
        transform("Metabolomics_batch" => (x-> any(!ismissing, x)) => "has_metabolomics")
        sort(["subject", "timepoint"])
        unique("sample")
    end

    fmp_samples = @chain DataFrame(XLSX.readtable(sample_centric_path, "Sheet1", infer_eltypes=true)) begin
        rename!(Dict("studyID"=> "subject", "collectionNum"=> "timepoint", "SampleID"=>"sample"))
        subset!("sample"    => ByRow(!ismissing),
                "subject"   => ByRow(s-> !ismissing(s) && s in samplemeta.subject),
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

    fmp_alltp = leftjoin(fmp_timepoint, fmp_subject, on=["subject"])

    transform!(fmp_alltp, AsTable(["ageMonths", "assessmentAgeDays", "assessmentAgeMonths", "scanAgeDays", "scanAgeMonths"]) => ByRow(row-> begin
        am, aad, aam, sad, sam = row
        calcage = coalesce(aad / 365 * 12 + aam, sad / 365 * 12 + sam)

        return coalesce(am, calcage)
        end) => "ageMonths"
    )

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

    DataFrames.transform!(groupby(fmp_alltp, :subject), 
                            :Merge_Dem_Child_Race => (r-> coalesce(r...)) => :Merge_Dem_Child_Race)

    fmp_alltp.race = let
        race = categorical(fmp_alltp.Merge_Dem_Child_Race)
        race[findall(r-> contains(string(r), '\n'), race)] .= "Mixed"        
        race = recode(race, 
            "Decline to Answer"                => missing,
            "Unknown"                          => missing,
            "Mixed Race"                       => "Mixed",
            "American Indian or Alaska Native" => "Indiginous",
            "Some other race"                  => "Other",
            "Other Asian"                      => "Asian",
            "Asian Indian"                     => "Asian",
            "Asian "                           => "Asian",
            "Black or African American"        => "Black",
        )
        droplevels!(race)
        levels!(race, ["White","Black","Asian","Indiginous", "Mixed","Other","Unknown"])
        race
    end
    fmp_alltp.sex = categorical(fmp_alltp.childGender)


    omni = @chain samplemeta begin
        subset("Fecal_EtOH" => ByRow(==("F")), "Mgx_batch"=> ByRow(!ismissing))
        select("subject", "timepoint", "sample"=>"omni", "collectionDate"=>"omni_collectionDate")
        unique(["subject", "timepoint"])
    end
    etoh = @chain samplemeta begin
        subset("Fecal_EtOH" => ByRow(==("E")), "Metabolomics_batch"=> ByRow(!ismissing))
        select("subject", "timepoint", "sample"=>"etoh", "collectionDate"=>"etoh_collectionDate")
        unique(["subject", "timepoint"])
    end

    leftjoin!(fmp_alltp, omni; on=["subject", "timepoint"])
    leftjoin!(fmp_alltp, etoh; on=["subject", "timepoint"])
    return fmp_alltp
end


function load_raw_metaphlan()
    df = DataFrame(file = filter(f-> contains(f, r"FG\d+_S\d+_profile"), readdir(analysisfiles("metaphlan"), join=true)))
    df.sample = map(s-> replace(s, "_profile.tsv"=> ""), basename.(df.file))
    df.sample_base = map(s-> replace(s, r"_S\d+"=>""), df.sample)

    knead = load(ReadCounts())
    taxa = metaphlan_profiles(df.file; samples = df.sample)
    set!(taxa, df)
    set!(taxa, select(knead, "sample_uid"=>"sample", AsTable(["final pair1", "final pair2"])=> ByRow(row-> row[1]+row[2]) =>"read_depth"))
    taxa
end


function load_raw_humann(; kind="genefamilies", overwrite=false, names=false, stratified=false, sample_filter=nothing)
    fname = joinpath(scratchfiles("genefunctions"), "$kind.arrow")
    @debug "writing to $fname"
    (!isfile(fname) || overwrite) && write_gfs_arrow(; kind, names, stratified, sample_filter)
    read_gfs_arrow(; kind)
end

function load_raw_brain()
    brain = XLSX.readxlsx(datafiles("AllResonanceVols_Resonance_only.xlsx"))
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
    return select(volumes, Cols(:subject, :timepoint, Not(["Cohort", "ID", "Session", "CSF (ignore)"])))
end