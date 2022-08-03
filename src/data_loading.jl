abstract type Dataset end

struct Metadata <: Dataset end
struct TaxonomicProfiles <: Dataset end
struct UnirefProfiles <: Dataset end
struct KOProfiles <: Dataset end
struct ECProfiles <: Dataset end
struct MetabolicProfiles <: Dataset end


load(::Metadata) = CSV.read(datafiles("exports", "timepoint_metadata.csv"), DataFrame)

function load(::TaxonomicProfiles; timepoint_metadata = load(Metadata()))
    tbl = Arrow.Table(datafiles("exports", "taxa.arrow"))
    mat = sparse(tbl.fidx, tbl.sidx, tbl.value)
    fs = [taxon(line) for line in eachline(datafiles("taxa_features.txt"))]
    ss = [MicrobiomeSample(line) for line in eachline(joinpath(scratch, "$(kind)_samples.txt"))]
    comm = CommunityProfile(mat, fs, ss)
    insert!(comm, timepoint_metadata)
    return comm
end