# #Metadata wrangling

# ## Get sample-specific metadata from airtable database

using Resonance

samplemeta = airtable_metadata() # having set ENV["AIRTABLE_KEY"]

# Now to add important metadata to samples

samples = MicrobiomeSample[]

let fields = [
    "timepoint",
    "subject",
    "sid_old",
    "CovidCollectionNumber",
    "Mgx_batch",
    "MaternalID"]

    for row in eachrow(samplemeta)
        startswith(row.sample, "FE") && continue # skip ethanol samples
        s = MicrobiomeSample(row.sample)
        set!(s, NamedTuple(Symbol(k) => v for (k, v) in pairs(row[Not(:sample)]) if !ismissing(v)))
        push!(samples, s)
    end
end

# ## Metadata stored in filemaker pro database
#
# Read in the different tables as `DataFrame`s,
# then normalize certain columns.
# First, samples

fmp_sample = DataFrame(XLSX.readtable("data/Sample_Centric_10252021.xlsx", "Sheet1", infer_eltypes=true)...)
rename!(fmp_sample, Dict(:SampleID=> :sample, :studyID=>:subject, :collectionNum=> :timepoint))


# Then, subject-specific data

fmp_subject = DataFrame(XLSX.readtable("data/Subject_Centric_10252021.xlsx", "Sheet1", infer_eltypes=true)...)
rename!(fmp_subject, Dict(:studyID=>:subject))

# Then, timepoint-specific data

fmp_timepoint = DataFrame(XLSX.readtable("data/Timepoint_Centric_10252021.xlsx", "Sheet1", infer_eltypes=true)...)
rename!(fmp_timepoint, Dict(:studyID=>:subject))

# and COVID-specific samples

fmp_covid = DataFrame(XLSX.readtable("data/COVID_Fecal_10252021.xlsx", "Sheet1", infer_eltypes=true)...)
rename!(fmp_covid, Dict(:studyID=>:subject))

# ## Getting data joined together

fmp_timed = outerjoin(fmp_timepoint,
                      subset(fmp_sample, :timepoint=> ByRow(!ismissing)), # filter out covid samples
                      on=[:subject, :timepoint])

fmp_alltp = leftjoin(fmp_timed, fmp_subject, on=[:subject])

# Add info about brain data (just if it's there)

brain = let 
    fcleft = brain_ingest("data/freesurfer_curvature_leftHemi_oct2021.csv")
    fcright = brain_ingest("data/freesurfer_curvature_rightHemi_oct2021.csv")
    ftleft = brain_ingest("data/freesurfer_thickness_leftHemi_oct2021.csv")
    ftright = brain_ingest("data/freesurfer_thickness_rightHemi_oct2021.csv")
    seg = brain_ingest("/home/kevin/Repos/Resonance/data/segmentationVolumeMeasurements_oct2021.csv")

    outerjoin(fcleft, fcright, ftleft, ftright, seg, on=[:subject, :timepoint], makeunique=true)
end

