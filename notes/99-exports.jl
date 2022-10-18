using Resonance
using Chain
using Arrow
using Dictionaries
using Statistics

#-

mdata = let df = CSV.read("data/timepoints_final.csv", DataFrame)
    df.has_segmentation = map(s-> s == "true", df.has_segmentation)
    df[df.has_segmentation .| .!ismissing.(df.cogScore) .| .!ismissing.(df.omni)  .| .!ismissing.(df.etoh), :]
end

sort!(mdata, ["subject", "timepoint"])

knead = Resonance.load_knead()
knead.omni = replace.(knead.Sample, r"_S\d+_kneaddata"=>"")
let keep = Set(skipmissing(mdata.omni))
    subset!(knead, :omni=> ByRow(s-> s in keep))
end
unique!(knead, :omni)
transform!(knead, AsTable(r"final") => ByRow(sum) => "read_depth")
leftjoin!(mdata, select(knead, "omni", "read_depth"); on="omni", matchmissing=:equal)

CSV.write(datafiles("exports", "timepoint_metadata.csv"), select(mdata,
    "subject", "timepoint", # Identifiers
    "ageMonths", "sex", "race", "education", # Demographics
    "assessmentDate", "scanDate", # Dates 
    "cogScore", "cogScorePercentile", "omni", "etoh", "has_segmentation", "read_depth" # Measurements
))


metabolites = CSV.read(datafiles("wrangled", "metabolites.csv"), DataFrame)

#-

@info "Unique subjects:" N = mdata.subject |> unique |> length
@info "Unique stool samples:" N = count(!ismissing, mdata.omni)
@info "Unique metabolomics samples:" N = count(!ismissing, mdata.etoh)
@info "Ages (in months)" min = minimum(skipmissing(mdata.ageMonths)) max = maximum(skipmissing(mdata.ageMonths)) mean = mean(skipmissing(mdata.ageMonths)) median = median(skipmissing(mdata.ageMonths))

#

keepomni = Set(skipmissing(mdata.omni))
dropetoh = setdiff(names(metabolites, r"^FE"), Set(skipmissing(mdata.etoh)))

#-

CSV.write(datafiles("exports", "metabolites.csv"), select(metabolites, Not(dropetoh)))

#-

allsamples = sort(collect(skipmissing(mdata.omni)))

allfiles = String[]
for (root, dirs, files) in walkdir(ENV["ANALYSIS_FILES"])
    filter!(f-> contains(f, r"(FG\d+)_S(\d+_)"), files)
    filter!(f-> String(first(split(basename(f), '_'))) in keepomni, files)
    unique!(f-> match(r"(FG\d+)_S\d+_(.+)", f).captures, files)
    append!(allfiles, joinpath.(root, files))
end

taxafiles = filter(f-> occursin("_profile.tsv", f), allfiles)
gffiles = filter(f-> occursin("_genefamilies.tsv", f), allfiles)
kofiles = filter(f-> occursin("_kos_rename.tsv", f), allfiles)
ecfiles = filter(f-> occursin("_ecs_rename.tsv", f), allfiles)


@assert allsamples == sort([String(first(split(basename(f), '_'))) for f in taxafiles]) ==
                      sort([String(first(split(basename(f), '_'))) for f in gffiles]) ==
                      sort([String(first(split(basename(f), '_'))) for f in kofiles]) ==
                      sort([String(first(split(basename(f), '_'))) for f in ecfiles])

sample_index = Dictionary((String(s) for s in allsamples), eachindex(allsamples))

#-

isdir(datafiles("exports")) || mkdir(datafiles("exports"))

open(datafiles("exports", "samples.txt"), "w") do io
    for s in allsamples
        println(io, s)
    end
end

# ## Write Taxonomic Profiles

let
    features = Set(String[])
    scratch = get(ENV, "SCRATCH_SPACE", "./scratch")
    isdir(scratch) || mkpath(scratch)

    tmp = tempname(scratch)

    @info "writing temporary arrow file for taxa"

    tbls = Tables.partitioner(taxafiles) do f
        samplename = replace(splitext(basename(f))[1], r"_S\d+_profile" => "")
        @info samplename

        sdf = CSV.read(f, DataFrame; header = ["feature", "NCBI_taxid", "abundance", "additional_species"], skipto=5)
        sdf.sample .= samplename
        sdf.sidx .= sample_index[samplename]
        select!(sdf, ["feature", "sample", "abundance", "sidx"])

        union!(features, sdf.feature)

        sdf
    end

    open(tmp, "w") do io
        Arrow.write(io, tbls)
    end

    @info "Making feature dictionary"
    fdic = Dictionary((f for f in features), 1:length(features))

    @info "building new table"
    df = DataFrame(Arrow.Table(tmp))
    df.fidx = [fdic[f] for f in df.feature]

    @info "Writing table"
    Arrow.write(datafiles("exports", "taxa.arrow"), df)
    open(datafiles("exports", "taxa_features.txt"), "w") do io
        for f in keys(fdic)
            println(io, f)
        end
    end
end

# ## Write Gene families (UniRef90) Profiles

let
    features = Set(String[])
    scratch = get(ENV, "SCRATCH_SPACE", "./scratch")
    isdir(scratch) || mkpath(scratch)

    tmp = tempname(scratch)

    @info "writing temporary arrow file for gene families"

    tbls = Tables.partitioner(gffiles) do f
        samplename = replace(splitext(basename(f))[1], r"_S\d+_genefamilies" => "")

        sdf = CSV.read(f, DataFrame; header=["feature", "value"], skipto=2)
        sdf.sample .= samplename
        sdf.sidx .= sample_index[samplename]

        union!(features, sdf.feature)

        sdf
    end

    open(tmp, "w") do io
        Arrow.write(io, tbls)
    end

    @info "Making feature dictionary"
    fdic = Dictionary((f for f in features), 1:length(features))

    @info "building new table"
    df = DataFrame(Arrow.Table(tmp))
    df.fidx = [fdic[f] for f in df.feature]

    @info "Writing table"
    Arrow.write(datafiles("exports", "genefamilies.arrow"), df)
    open(datafiles("exports", "genefamilies_features.txt"), "w") do io
        for f in keys(fdic)
            println(io, f)
        end
    end
end

# ## Write Kegg Orthology (KO) Profiles

let
    features = Set(String[])
    scratch = get(ENV, "SCRATCH_SPACE", "./scratch")
    isdir(scratch) || mkpath(scratch)

    tmp = tempname(scratch)

    @info "writing temporary arrow file for KOs"

    tbls = Tables.partitioner(kofiles) do f
        samplename = replace(splitext(basename(f))[1], r"_S\d+_kos_rename" => "")

        sdf = CSV.read(f, DataFrame; header=["feature", "value"], skipto=2)
        sdf.sample .= samplename
        sdf.sidx .= sample_index[samplename]

        union!(features, sdf.feature)

        sdf
    end

    open(tmp, "w") do io
        Arrow.write(io, tbls)
    end

    @info "Making feature dictionary"
    fdic = Dictionary((f for f in features), 1:length(features))

    @info "building new table"
    df = DataFrame(Arrow.Table(tmp))
    df.fidx = [fdic[f] for f in df.feature]

    @info "Writing table"
    Arrow.write(datafiles("exports", "kos.arrow"), df)
    open(datafiles("exports", "kos_features.txt"), "w") do io
        for f in keys(fdic)
            println(io, f)
        end
    end
end

# ## Write Level 4 Enzyme Class (EC) Profiles

let
    features = Set(String[])
    scratch = get(ENV, "SCRATCH_SPACE", "./scratch")
    isdir(scratch) || mkpath(scratch)

    tmp = tempname(scratch)

    @info "writing temporary arrow file for ecs"

    tbls = Tables.partitioner(ecfiles) do f
        samplename = replace(splitext(basename(f))[1], r"_S\d+_ecs_rename" => "")

        sdf = CSV.read(f, DataFrame; header=["feature", "value"], skipto=2)
        sdf.sample .= samplename
        sdf.sidx .= sample_index[samplename]

        union!(features, sdf.feature)

        sdf
    end

    open(tmp, "w") do io
        Arrow.write(io, tbls)
    end

    @info "Making feature dictionary"
    fdic = Dictionary((f for f in features), 1:length(features))

    @info "building new table"
    df = DataFrame(Arrow.Table(tmp))
    df.fidx = [fdic[f] for f in df.feature]

    @info "Writing table"
    Arrow.write(datafiles("exports", "ecs.arrow"), df)
    open(datafiles("exports", "ecs_features.txt"), "w") do io
        for f in keys(fdic)
            println(io, f)
        end
    end
end

