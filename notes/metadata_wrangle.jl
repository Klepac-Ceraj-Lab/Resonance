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

using XLSX
using DataFrames

fmp_sample = DataFrame(XLSX.readtable("data/Sample_Centric_10252021.xlsx", "Sheet1")...)
fmp_subject = DataFrame(XLSX.readtable("data/Subject_Centric_10252021.xlsx", "Sheet1")...)
fmp_timepoint = DataFrame(XLSX.readtable("data/Timepoint_Centric_10252021.xlsx", "Sheet1")...)

fmp_covid = DataFrame(XLSX.readtable("data/COVID_Fecal_10252021.xlsx", "Sheet1")...)

