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
struct ReadCounts <: Dataset end
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

load(ds::Dataset; kwargs...) = throw(MethodError("load has not been implemented for $(typeof(ds))"))

function load(::Metadata)
    Setup.datadownload(Setup.Timepoints(); inputdir=inputfiles())
    df = CSV.read(inputfiles("complete_filtered_dataset.csv"), DataFrame;
    #     types = [
    #         Int64,                   # subject
    #         Int64,                   # timepoint
    #         Union{Missing, Float64}, # ageMonths
    #         String,                  # sex
    #         String,                  # race
    #         Union{Missing, String},  # education
    #         Union{Missing, Date},    # date
    #         Union{Missing, Float64}, # cogScore
    #         String,                  # sample
    #         String,                  # sample_base
    #         Union{Missing, Float64}, # read_depth
    #         Bool,                    # filter_00to120 
    #         Bool,                    # filter_00to06
    #         Bool,                    # filter_18to120 
    #     ]        
        stringtype = String
    )
    DataFrames.select!(subset!(df, "has_concurrent_stool_cog"=> identity), 
        "subject",
        "timepoint",
        "ageMonths",
        "sex",
        "education",
        "race",
        "cogScore",
        r"filter_",
        "seqid"=> "sample",
        "omni",
        r"^Mullen.+Composite$"
    )

    df.sex = categorical(df.sex)
    df.race = categorical(df.race)
    df.education = categorical(df.education; ordered = true, levels = ["Junior high school", "Some high school", "High school grad", "Some college", "College grad", "Grad/professional school"])
    return df
end

function load(::TaxonomicProfiles; timepoint_metadata = load(Metadata()))
    Setup.datadownload(Setup.Taxa(); inputdir=inputfiles())
    comm = read_arrow(inputfiles("taxa.arrow"); featurefunc = taxon)
    insert!(comm, timepoint_metadata; namecol=:sample)
    return comm[:, timepoint_metadata.sample]
end

function load(::UnirefProfiles; timepoint_metadata = load(Metadata()))
    Setup.datadownload(Setup.Unirefs(); inputdir=inputfiles())
    comm = read_arrow(inputfiles("genefamilies.arrow"); featurefunc = genefunction)
    insert!(comm, timepoint_metadata; namecol=:sample)
    return comm[:, timepoint_metadata.sample]
end

function load(::KOProfiles; timepoint_metadata = load(Metadata()))
    Setup.datadownload(Setup.KOs(); inputdir=inputfiles())
    comm = read_arrow(inputfiles("kos.arrow"); featurefunc = genefunction)
    insert!(comm, timepoint_metadata; namecol=:sample)
    return comm[:, timepoint_metadata.sample]
end

function load(::ECProfiles; timepoint_metadata = load(Metadata()))
    Setup.datadownload(Setup.ECs(); inputdir=inputfiles())
    comm = read_arrow(inputfiles("ecs.arrow"); featurefunc = genefunction)
    insert!(comm, timepoint_metadata; namecol=:sample)
    return comm[:, timepoint_metadata.sample]
end

# function load(::MetabolicProfiles; timepoint_metadata = load(Metadata()))
#     df = CSV.read(inputfiles("metabolites.csv"), DataFrame)
#
#     ms = [Metabolite(row[:uid], row[:Metabolite], row[:MZ], row[:RT]) for row in eachrow(df)]
#     comm = CommunityProfile(Matrix(df[!, r"^FE"]), ms, MicrobiomeSample.(names(df, r"^FE")))
#     set!(comm, timepoint_metadata; namecol=:etoh)
#     return comm
# end

function load(::Neuroimaging; timepoint_metadata = load(Metadata()), samplefield = "sample")
    Setup.datadownload(Setup.Neuro(); inputdir=inputfiles())
    df = CSV.read(inputfiles("brain_normalized.csv"), DataFrame)
    fact = df."White-matter" .+ df."Gray-matter"
    for f in names(df, Not(["subject", "timepoint"]))
        df[!, f] ./= fact
    end

    mat = Matrix(select(df, Not(["subject", "timepoint"])))' |> collect
    feats = brainvolume.(names(df, Not(["subject", "timepoint"])))
    
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

function load(::ReadCounts; srcdir=analysisfiles("kneaddata"), readcountfile="read_counts.csv", overwrite=false)
    destfile = joinpath(srcdir, readcountfile)
    if isfile(destfile) && !overwrite
        return CSV.read(destfile, DataFrame; types=(i, col) -> contains(string(col), "sample") ? String : Float64)
    else
        run(Cmd([
            "kneaddata_read_count_table",
            "--input", srcdir,
            "--output", destfile
        ]))

        df = CSV.read(destfile, DataFrame; missingstring="NA")
        for col in names(df, r"Homo_sapiens")
            df[!, col] = map(p-> coalesce(p...), zip(df[!, col], df[!, replace(col, "Homo_sapiens"=>"hg37dec_v0.1")]))
        end
        df."decontaminated Homo_sapiens pair1" = map(x-> ismissing(x) ? missing :  x isa Real ? x : parse(Float64, x), df."decontaminated Homo_sapiens pair1")
        df."decontaminated Homo_sapiens orphan2" = map(x-> ismissing(x) ? missing :  x isa Real ? x : parse(Float64, x), df."decontaminated Homo_sapiens orphan2")
        for row in eachrow(df)
            ismissing(row."final pair1") || continue
            for i in (1,2)
                ismissing(row["decontaminated Homo_sapiens pair$i"]) && (row["final pair$i"] = row["trimmed pair$i"])
                ismissing(row["decontaminated Homo_sapiens orphan$i"]) && (row["final orphan$i"] = row["trimmed pair$i"])
            end
        end

        df = select(df, Not(r"hg37"))
        # subset!(df, "Sample"=> ByRow(s-> contains(s, r"^FG\d+")))
        df.sample_uid = map(s-> replace(s, r"_kneaddata"=>""), df.Sample)
        df.sample = map(s-> replace(s, r"_S\d+"=>""), df.sample_uid)
        select!(df, Cols("sample", "sample_uid", r".+[12]"))
        CSV.write(destfile, df)
        return df
    end
end

"""
Given 2 vectors of (:subject, :timepoint) tuples,
find the overlap, and put them in the same order.
"""
function stp_overlap(ds1, ds2; lt = (x,y)-> x[1] == y[1] ? x[2] < y[2] : x[1] < y[1]) 
    srt1 = sortperm(ds1; lt)
    srt2 = sortperm(ds2; lt)
    
    keep1 = _keep_unique(ds1, ds2, srt1)
    keep2 = _keep_unique(ds2, ds1, srt2)

    return srt1[keep1], srt2[keep2]
end

function _keep_unique(ds1, ds2, srt)
    recip = Set(ds2)
    keep = Bool[]
    used = Set(eltype(ds1)[])
    for d in ds1[srt]
        if d ∈ used
            push!(keep, false)
        else
            push!(used, d)
            push!(keep, d ∈ recip)
        end
    end
    return keep
end