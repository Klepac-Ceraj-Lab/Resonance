using Microbiome
using CSV
using DataFrames
using ECHOSRS2
using SparseArrays
using ProgressLogging
using Dictionaries
using CodecZlib
using MultipleTesting
using Statistics
using HypothesisTests
using CairoMakie

outdir = "/dartfs-hpc/rc/lab/M/MRKepistor7/collab/KevinBonham/ResultsFiles"
figsdir = joinpath(outdir, "figures")
isdir(figsdir) || mkdir(figsdir)

meta = CSV.read("data/NHBCS.DATA.MI.csv", DataFrame)
meta.BFEED = map(bf-> bf=="NA" ? missing : parse(Int, bf), meta.BFEED)
gfs = DataFrame()

files = readdir("/dartfs-hpc/rc/lab/M/MRKepistor7/collab/KevinBonham/SourceFiles/gene_families/raw_output/", join=true)
filter!(f-> occursin(r"_genefamilies\.tsv\.gz", f), files)

for f in files
    s = first(split(basename(f), "_"))
    row = findfirst(==(s), meta.sfdid)
    isnothing(row) && continue

    sample = MicrobiomeSample(s)
    for col in propertynames(meta)[2:end]
        set!(sample, col, meta[row, col])
    end

    gf = CSV.read(GzipDecompressorStream(open(f)), DataFrame)
    rename!(gf, ["feature", "abundance"])
    filter!("feature" => (x-> !occursin(r"\|[gu]", x)), gf)

    gf.sample = fill(sample, nrow(gf))
    gf.feature = GeneFunction.(gf.feature)
    append!(gfs, gf)
end

cp = let
    features = unique(gfs.feature)
    samples = unique(gfs.sample)
    fdict = Dictionary(name.(features), eachindex(features))
    sdict = Dictionary(name.(samples), eachindex(samples))

    mat = spzeros(length(features), length(samples))
    @progress for row in eachrow(gfs)
        mat[fdict[name(row.feature)], sdict[name(row.sample)]] = row.abundance
    end
    CommunityProfile(mat, features, samples)
end

cpmeta = DataFrame(get(cp))
accessory = findall(p-> 0.1 < p < 0.9, vec(prevalence(cp)))

unirefaccessory = cp[accessory, :]

unirefnames = map(u-> match.(r"UniRef90_(\w+)",u).captures[1], featurenames(unirefaccessory))
neuroactive = getneuroactive(unirefnames)

allneuroactive = union([neuroactive[k] for k in keys(neuroactive)]...)
metadatums = [:SCORE,
              :RESID,
              :BFEED]


##

sixweeks = filter(:REALTIME=> ==("6W"), meta)

mdcors = Dict(m=>Float64[] for m in metadatums)

fig = Figure()

for md in metadatums[2:end]
    filt = map(!ismissing, sixweeks[!,md])
    cors = cor(sixweeks[filt, md], abundances(unirefaccessory[:, sixweeks.sfdid])[:,filt], dims=2)'
    mdcors[md] = filter(!isnan, cors)
end

##


fig = Figure()

h = mn = mdn = nothing

for (i, md) in enumerate(metadatums)
    row = i > 2 ? 2 : 1
    i = (i + 1) % 2 + 1
    ax = Axis(fig[i,row], title="Gene function correlations with $md")
    h = hist!(ax, mdcors[md])
    mn = vlines!(ax, [mean(mdcors[md])])
    mdn = vlines!(ax, [median(mdcors[md])], linestyle=:dash)
end



leg = Legend(fig[2,2], [h, mn, mdn], ["frequency", "mean", "median"]; tellwidth=false)
save(joinpath(figsdir, "gf_all_correlations.pdf"), fig)
fig