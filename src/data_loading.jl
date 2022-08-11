abstract type Dataset end

struct Metadata <: Dataset end
struct TaxonomicProfiles <: Dataset end
struct UnirefProfiles <: Dataset end
struct KOProfiles <: Dataset end
struct ECProfiles <: Dataset end
struct MetabolicProfiles <: Dataset end
struct Neuroimaging <: Dataset end

load(ds::Dataset; kwargs...) = MethodError("load has not been implemented for $(typeof(ds))")

load(::Metadata) = CSV.read(datafiles("exports", "timepoint_metadata.csv"), DataFrame;
    types = [
        Int64,                   # subject
        Int64,                   # timepoint
        Union{Missing, Float64}, # ageMonths
        Union{Missing, String},  # race
        Union{Missing, String},  # maternalEd
        Union{Missing, Date},    # assessmentDate
        Union{Missing, Date},    # scanDate
        Union{Missing, Float64}, # cogScore
        Union{Missing, String},  # omni
        Union{Missing, String},  # etoh
        Bool,                     # has_segmentation
        Union{Missing, Float64}  # read_depth
    ]        
)

function load(::TaxonomicProfiles; timepoint_metadata = load(Metadata()))
    tbl = Arrow.Table(datafiles("exports", "taxa.arrow"))
    mat = sparse(tbl.fidx, tbl.sidx, tbl.abundance)
    fs = [taxon(last(split(line, "|"))) for line in eachline(datafiles("exports", "taxa_features.txt"))]
    ss = [MicrobiomeSample(line) for line in eachline(datafiles("exports", "samples.txt"))]
    comm = CommunityProfile(mat, fs, ss)
    insert!(comm, timepoint_metadata; namecol=:omni)
    return comm
end

function load(::UnirefProfiles; timepoint_metadata = load(Metadata()))
    tbl = Arrow.Table(datafiles("exports", "genefamilies.arrow"))
    mat = sparse(tbl.fidx, tbl.sidx, tbl.value)
    fs = [genefunction(line) for line in eachline(datafiles("exports", "genefamilies_features.txt"))]
    ss = [MicrobiomeSample(line) for line in eachline(datafiles("exports", "samples.txt"))]
    comm = CommunityProfile(mat, fs, ss)
    insert!(comm, timepoint_metadata; namecol=:omni)
    return comm
end

function load(::KOProfiles; timepoint_metadata = load(Metadata()))
    tbl = Arrow.Table(datafiles("exports", "kos.arrow"))
    mat = sparse(tbl.fidx, tbl.sidx, tbl.value)
    fs = [genefunction(line) for line in eachline(datafiles("exports", "kos_features.txt"))]
    ss = [MicrobiomeSample(line) for line in eachline(datafiles("exports", "samples.txt"))]
    comm = CommunityProfile(mat, fs, ss)
    insert!(comm, timepoint_metadata; namecol=:omni)
    return comm
end

function load(::ECProfiles; timepoint_metadata = load(Metadata()))
    tbl = Arrow.Table(datafiles("exports", "ecs.arrow"))
    mat = sparse(tbl.fidx, tbl.sidx, tbl.value)
    fs = [genefunction(line) for line in eachline(datafiles("exports", "ecs_features.txt"))]
    ss = [MicrobiomeSample(line) for line in eachline(datafiles("exports", "samples.txt"))]
    comm = CommunityProfile(mat, fs, ss)
    insert!(comm, timepoint_metadata; namecol=:omni)
    return comm
end

function load(::MetabolicProfiles; timepoint_metadata = load(Metadata()))
    df = CSV.read(datafiles("exports", "metabolites.csv"), DataFrame)

    ms = [Metabolite(row[:uid], row[:Metabolite], row[:MZ], row[:RT]) for row in eachrow(df)]
    comm = CommunityProfile(Matrix(df[!, r"^FE"]), ms, MicrobiomeSample.(names(df, r"^FE")))
    set!(comm, timepoint_metadata; namecol=:etoh)
    return comm
end