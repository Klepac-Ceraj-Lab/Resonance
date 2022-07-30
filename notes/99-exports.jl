using Resonance
using Chain
using Arrow

#-

metadata = let df = CSV.read("data/timepoints_final.csv", DataFrame)
    df.has_segmentation = map(s-> s == "true", df.has_segmentation)
    df[df.has_segmentation .| .!ismissing.(df.cogScore) .| .!ismissing.(df.omni)  .| .!ismissing.(df.etoh), :]
end

CSV.write("data/timepoint_metadata.csv", select(metadata,
    "subject", "timepoint", # Identifiers
    "ageMonths", "race", "ed" => "maternalEd", # Demographics
    "assessmentDate", "scanDate", # Dates 
    "cogScore", "omni", "etoh", "has_segmentation" # Measurements
))

CSV.write("data/brain_measures.csv", select(metadata, 
    "subject", "timepoint",
    Resonance.brainmeta...
))

keepomni = Set(skipmissing(metadata.omni))

#-

species = CSV.read(datafiles("species.csv"), DataFrame)
select!(species, "features", names(species)[map(n-> n in keepomni, names(species))]...) # Missing 11 samples?!?!
CSV.write("data/species.csv", species)


#-

allfiles = String[]

@info "getting files"
for (root, dirs, files) in walkdir(ENV["ANALYSIS_FILES"])
    filter!(f-> contains(splitext(f)[1], Regex(string("genefamilies", '$'))) && !contains(f, r"^FE\d+"), files)
    append!(allfiles, joinpath.(Ref(root), files))
end

unique!(allfiles) do f
    first(split(basename(f), '_')) |> String
end
filter(f-> String(first(split(basename(f), '_'))) in keepomni, allfiles)

if !isempty(keepsamples)
    filter!(f-> String(first(split(basename(f), '_'))) in keepsamples, allfiles)
    length(allfiles) > 0 || throw(ErrorException("No sample intersection found"))
end

#-

Resonance.write_gfs_arrow("data"; keepsamples = keepomni)
Resonance.write_gfs_arrow("data"; kind=keepsamples = keepomni)