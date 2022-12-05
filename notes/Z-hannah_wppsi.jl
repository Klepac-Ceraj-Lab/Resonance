using Resonance

metaphlan_files = filter(readdir(analysisfiles("metaphlan"), join=true)) do f
    !contains(f, r"FE\d{5}") && contains(f, "profile.tsv")
end

#- 

samplemeta = airtable_metadata() # having set ENV["AIRTABLE_KEY"]

#-

@chain samplemeta begin
    groupby([:subject, :timepoint])
    transform!(:Metabolomics_batch => (x-> any(!ismissing, x)) => :has_metabolomics)
    sort!([:subject, :timepoint])
end


fmp_samples = DataFrame(XLSX.readtable(datafiles("resonance_fmp", "Sample_Centric_071322.xlsx"),
                        "Sheet1", infer_eltypes=true))
rename!(fmp_samples, Dict("studyID"=>"subject", "collectionNum"=> "timepoint", "SampleID"=>"sample"))
subset!(fmp_samples, :subject => ByRow(s-> s in samplemeta.subject),
                     :timepoint => ByRow(!ismissing)
)

samplemeta = leftjoin(samplemeta, select(fmp_samples, ["sample", "collectionDate"]), on="sample")

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

fmp_timepoint = DataFrame(XLSX.readtable(datafiles("resonance_fmp", "Timepoint_Centric_071322.xlsx"), 
        "Sheet1", infer_eltypes=true))
rename!(fmp_timepoint, Dict(:studyID=>:subject))

fmp_timepoint = let cs_perc = CSV.read(datafiles("cogscore_percentiles.csv"), DataFrame)
    leftjoin(fmp_timepoint, select(cs_perc, "subject", "timepoint", "cogScorePercentile"), on=["subject", "timepoint"])
end

fmp_alltp = leftjoin(fmp_timepoint, fmp_subject, on=[:subject])

fmp_alltp.ageMonths = map(eachrow(fmp_alltp)) do row
    if ismissing(row.scanAgeMonths)
        (row.assessmentAgeDays ./ 365 .* 12) .+ row.assessmentAgeMonths
    else
        (row.scanAgeDays ./ 365 .* 12) .+ row.scanAgeMonths
    end
end

@assert nrow(unique(fmp_alltp, ["subject", "timepoint"])) == nrow(fmp_alltp)

final = leftjoin(select(samplemeta, "sample", "subject", "timepoint", "DOC"=>"collectionDate", "Mgx_batch"),
                         fmp_alltp, on=["subject", "timepoint"]
)
subset!(final, "sample"=> ByRow(s->startswith(s, "FG")))
sort!(final, "sample")

count(eachrow(unique(final, ["subject", "timepoint"]))) do row
    !ismissing(row.Mgx_batch) && !ismissing(row.ageMonths) && !ismissing(row.cogScore) &&
    2.5*12 <= row.ageMonths <= 7.5*12
end

count(eachrow(unique(final, ["subject", "timepoint"]))) do row
    !ismissing(row.ageMonths) && !ismissing(row.cogScore) &&
    2.5*12 <= row.ageMonths <= 7.5*12
end

cm_assessments = let asmt = final[map(eachrow(final)) do row
        !ismissing(row.ageMonths) && !ismissing(row.cogScore) &&
        2.5*12 <= row.ageMonths <= 7.5*12
    end, :].cogAssessment

    Dict(a=> count(t-> t == a, asmt) for a in unique(asmt))
end


## Candance

taxa = metaphlan_profiles(metaphlan_files, :species)
samp_all = unique(s-> replace(s, r"_S\d{1,2}_profile"=>""), samplenames(taxa))
taxa = taxa[:, samp_all]


taxa = CommunityProfile(abundances(taxa), features(taxa),
        [MicrobiomeSample(replace(s, r"_S\d{1,2}_profile"=>"")) for s in samplenames(taxa)])

set!(taxa, final)
taxa = taxa[:, .!ismissing.(get(taxa, :subject))]

mdt = DataFrame(
    sample = samplenames(taxa),
    subject = get(taxa, :subject),
    timepoint = get(taxa, :timepoint),
    ageMonths = get(taxa, :ageMonths),
    sex = get(taxa, :childGender)
    collectionDate = get(taxa, :collectionDate)    
)

for f in features(taxa)
    mdt[!, name(f)] = vec(abundances(taxa[f, :]))
end

subjects = Set(parse.(Int, readlines("data/candace_ids.txt")))

CSV.write("data/candace_wide.csv", subset(mdt, "subject" => ByRow(s-> s in subjects)))