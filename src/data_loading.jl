# Microbiome.jl compat with neuroimaging

abstract type Hemisphere end

struct Left <: Hemisphere end
struct Right <: Hemisphere end

function Hemisphere(h::Symbol)
    h == :left && return Left()
    h == :right && return Right()
    throw(ArgumentError("Hemisphere must be left or right"))
end

Base.String(::Left) = "left"
Base.String(::Right) = "right"


struct BrainVolume <: Microbiome.AbstractFeature
    name::String
    hemisphere::Union{Missing, Hemisphere}
    
    BrainVolume(s::AbstractString, ::Missing) = new(s, missing)
    BrainVolume(s::AbstractString, h::Symbol) = new(s, Hemisphere(h))
end

Microbiome.name(bv::BrainVolume) = bv.name
hemisphere(bv::BrainVolume) = bv.hemisphere
hashemisphere(bv::BrainVolume) = !ismissing(hemisphere(bv))

function Base.String(bv::BrainVolume)
    if !hashemisphere(bv)
        return name(bv)
    else
        return join((String(hemisphere(bv)), name(bv)), "-")
    end
end

function brainvolume(s::String)
    if startswith(s, r"(left|right)-")
        parts = split(s, '-')
        BrainVolume(join(parts[2:end], '-'), Symbol(first(parts)))
    else
        BrainVolume(s, missing)
    end
end

@testset "Brain Volumes" begin
    bv1 = brainvolume("left-lateral-ventricle")
    bv2 = brainvolume("right-pericalcarine")
    bv3 = brainvolume("Brain-stem")

    @test hemisphere(bv1) isa Left
    @test hemisphere(bv2) isa Right
    @test ismissing(hemisphere(bv3))
    @test hashemisphere(bv1)
    @test hashemisphere(bv2)
    @test !hashemisphere(bv3)
    @test name(bv1) == "lateral-ventricle"
    @test name(bv2) == "pericalcarine"
    @test name(bv3) == "Brain-stem"
    @test String(bv1) == "left-lateral-ventricle"
    @test String(bv2) == "right-pericalcarine"
    @test String(bv3) == "Brain-stem"
end

abstract type Dataset end

struct Metadata <: Dataset end
struct TaxonomicProfiles <: Dataset end
struct UnirefProfiles <: Dataset end
struct KOProfiles <: Dataset end
struct ECProfiles <: Dataset end
# struct MetabolicProfiles <: Dataset end
struct Neuroimaging <: Dataset end

# function map_string_timepoint(concatstring)::Tuple{Int64, Int64}
#     letter_to_timepoint = Dict( [ x => y for (x,y) in zip('a':'l',1:14)] ) 
#     isletter(concatstring[end]) ? (return(parse(Int64, concatstring[1:end-1]), letter_to_timepoint[concatstring[end]])) : (return(parse(Int64, concatstring), 1))
# end # end function

# function split_subject_timepoint(stringsvector)
#     mappings = map(map_string_timepoint, stringsvector)
#     subjects = map(x -> x[1], mappings)
#     timepoints = map(x -> x[2], mappings)
#     return subjects, timepoints
# end # end function

load(ds::Dataset; kwargs...) = MethodError("load has not been implemented for $(typeof(ds))")

function load(::Metadata)
    Setup.datadownload(Setup.Timepoints(); inputdir=inputfiles())
    df = CSV.read(inputfiles("timepoint_metadata.csv"), DataFrame;
        types = [
            Int64,                   # subject
            Int64,                   # timepoint
            Union{Missing, Float64}, # ageMonths
            String,                  # sex
            String,                  # race
            Union{Missing, String},  # education
            Union{Missing, Date},    # assessmentDate
            Union{Missing, Date},    # scanDate
            Union{Missing, Float64}, # cogScore
            Union{Missing, Float64}, # cogScorePercentile
            Union{Missing, String},  # omni
            Union{Missing, String},  # etoh
            Bool,                    # has_segmentation
            Union{Missing, Float64}  # read_depth
        ]        
    )

    df.sex = categorical(df.sex)
    df.race = categorical(df.race)
    df.education = categorical(df.education; ordered = true, levels = ["Junior high school", "Some high school", "High school grad", "Some college", "College grad", "Grad/professional school"])
    return df
end

function load(::TaxonomicProfiles; timepoint_metadata = load(Metadata()))
    Setup.datadownload(Setup.Taxa(); inputdir=inputfiles())
    tbl = Arrow.Table(inputfiles("taxa.arrow"))
    mat = sparse(tbl.fidx, tbl.sidx, tbl.abundance)
    fs = [taxon(last(split(line, "|"))) for line in eachline(inputfiles("taxa_features.txt"))]
    ss = [MicrobiomeSample(line) for line in eachline(inputfiles("samples.txt"))]
    comm = CommunityProfile(mat, fs, ss)
    insert!(comm, timepoint_metadata; namecol=:omni)
    return comm
end

function load(::UnirefProfiles; timepoint_metadata = load(Metadata()))
    Setup.datadownload(Setup.Unirefs(); inputdir=inputfiles())
    tbl = Arrow.Table(inputfiles("genefamilies.arrow"))
    mat = sparse(tbl.fidx, tbl.sidx, tbl.value)
    fs = [genefunction(line) for line in eachline(inputfiles("genefamilies_features.txt"))]
    ss = [MicrobiomeSample(line) for line in eachline(inputfiles("samples.txt"))]
    comm = CommunityProfile(mat, fs, ss)
    insert!(comm, timepoint_metadata; namecol=:omni)
    return comm
end

function load(::KOProfiles; timepoint_metadata = load(Metadata()))
    Setup.datadownload(Setup.KOs(); inputdir=inputfiles())
    tbl = Arrow.Table(inputfiles("kos.arrow"))
    mat = sparse(tbl.fidx, tbl.sidx, tbl.value)
    fs = [genefunction(line) for line in eachline(inputfiles("kos_features.txt"))]
    ss = [MicrobiomeSample(line) for line in eachline(inputfiles("samples.txt"))]
    comm = CommunityProfile(mat, fs, ss)
    insert!(comm, timepoint_metadata; namecol=:omni)
    return comm
end

function load(::ECProfiles; timepoint_metadata = load(Metadata()))
    Setup.datadownload(Setup.ECs(); inputdir=inputfiles())
    tbl = Arrow.Table(inputfiles("ecs.arrow"))
    mat = sparse(tbl.fidx, tbl.sidx, tbl.value)
    fs = [genefunction(line) for line in eachline(inputfiles("ecs_features.txt"))]
    ss = [MicrobiomeSample(line) for line in eachline(inputfiles("samples.txt"))]
    comm = CommunityProfile(mat, fs, ss)
    insert!(comm, timepoint_metadata; namecol=:omni)
    return comm
end

# function load(::MetabolicProfiles; timepoint_metadata = load(Metadata()))
#     df = CSV.read(inputfiles("metabolites.csv"), DataFrame)
#
#     ms = [Metabolite(row[:uid], row[:Metabolite], row[:MZ], row[:RT]) for row in eachrow(df)]
#     comm = CommunityProfile(Matrix(df[!, r"^FE"]), ms, MicrobiomeSample.(names(df, r"^FE")))
#     set!(comm, timepoint_metadata; namecol=:etoh)
#     return comm
# end

function load(::Neuroimaging; timepoint_metadata = load(Metadata()), samplefield = "omni")
    Setup.datadownload(Setup.Neuro(); inputdir=inputfiles())
    df = CSV.read(inputfiles("brain_normalized.csv"), DataFrame)
    mat = Matrix(select(df, Not(["subject", "timepoint", "Sex", "AgeInDays"])))' |> collect

    feats = brainvolume.(names(df, Not(["subject", "timepoint", "Sex", "AgeInDays"])))
    
    grp = groupby(timepoint_metadata, ["subject", "timepoint"])

    samps = map(eachrow(df)) do row
        g = get(grp, (; subject=row.subject, timepoint=row.timepoint), nothing)
        sidx = isnothing(g) ? nothing : findfirst(!ismissing, g[!, samplefield])
        if isnothing(sidx)
            s = MicrobiomeSample("sub$(row.subject)_tp$(row.timepoint)")
            set!(s, :subject, row.subject)
            set!(s, :timepoint, row.timepoint)
            set!(s, :ageMonths, row.AgeInDays / 365 * 12)
            set!(s, :sex, row.Sex)
            set!(s, :hassample, false)
        else
            s = MicrobiomeSample(g[sidx, samplefield])
            set!(s, :hassample, true)
        end
        return s
    end

    comm = CommunityProfile(mat, feats, samps)
    set!(comm, timepoint_metadata; namecol=Symbol(samplefield))
    return comm
end