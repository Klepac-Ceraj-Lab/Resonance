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

cpmeta = DataFrame(metadata(cp))
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
allfsea = DataFrame(
            geneset   = String[],
            metadatum = String[],
            median    = Float64[],
            U         = Float64[],
            f         = Float64[],
            pvalue    = Float64[],
            cors      = Vector{Float64}[],
            ext       = Tuple[])

mdcors = Dict(m=>Float64[] for m in metadatums)


for md in metadatums
    @info "Working on $md"
    filt = map(!ismissing, sixweeks[!,md])
    cors = cor(sixweeks[filt, md], abundances(unirefaccessory[:, sixweeks.sfdid])[:,filt], dims=2)'
    srtcors = sort(cors)
    mdcors[md] = filter(!isnan, cors)
    
    for (key, pos) in pairs(neuroactive)
        allcors = filter(!isnan, cors[pos])
        notcors = filter(!isnan, cors[Not(pos)])
        length(allcors) < 4 && continue
        @info "    $key"
        mwu = MannWhitneyUTest(allcors, notcors)
        u = mwu.U
        f = u / (mwu.nx * mwu.ny)
        m = median(allcors)
        p = pvalue(mwu)
        push!(allfsea, (geneset=key, metadatum=String(md), median=m, U=u, f=f, pvalue=p, cors=allcors, ext=extrema(mdcors[md])))
    end
end

allfsea.qvalue = adjust(allfsea.pvalue, BenjaminiHochberg())
allfsea.allmed = [median(mdcors[Symbol(md)]) for md in allfsea.metadatum]
sort!(allfsea, :qvalue)
CSV.write(joinpath(outdir, "6w_fsea.csv"), allfsea[!, Not(:cors)])

##

for md in metadatums
    @info "Plotting $md"
    df = filter(row-> row.qvalue < 0.2 && row.metadatum == string(md), allfsea)
    @info size(df)
    fig = Figure(resolution = (900, 1200))
    
    sq = ceil(Int, sqrt(nrow(df)))
    coln = 1
    rown = 1

    for row in eachrow(df)
        ax = Axis(fig[rown,coln], yticks=-0.3:0.1:0.3, yminorticks=IntervalsBetween(2), yminorticksvisible = true)
        hidexdecorations!(ax)
        clr = row.median < row.allmed ? :purple : :orange

        hlines!(ax, row.cors, color=clr)
        hlines!(ax, [row.allmed], color=:black, linestyle=:dash)
        hlines!(ax, [row.ext...], color=:black, linestyle=:dot)
        coln += 1
        Label(fig[rown,coln], string(row.geneset, " | f = ", round(row.f, digits=3), " | q = ", round(row.qvalue, digits=4)), rotation=-pi/2, tellheight=false)
        coln += 1
        if coln > sq * 2
            rown += 1
            coln = 1
        end
    end


    leg = Legend(fig[:,end+1],
        [
            LineElement(color = :purple, linestyle = nothing),
            LineElement(color = :orange, linestyle = nothing),
            LineElement(color = :black, linestyle = :dash),
            LineElement(color = :black, linestyle = :dot),
        ], ["negative assoc.", "positive assoc.", "median (all genes)", "extrema"]
    )
    supertitle = Label(fig[0, 1:end-1], "6W FSEA for $md", textsize = 20)
    save(joinpath(figsdir, "6W_fsea_$md.pdf"), fig)
    fig
end
fig
##


##

twelvemonths = filter(:REALTIME=> ==("12M"), meta)
allfsea = DataFrame(
            geneset   = String[],
            metadatum = String[],
            median    = Float64[],
            U         = Float64[],
            f         = Float64[],
            pvalue    = Float64[],
            cors      = Vector{Float64}[],
            ext       = Tuple[])

mdcors = Dict(m=>Float64[] for m in metadatums)


for md in metadatums
    @info "Working on $md"
    filt = map(!ismissing, twelvemonths[!,md])
    cors = cor(twelvemonths[filt, md], abundances(unirefaccessory[:, twelvemonths.sfdid])[:,filt], dims=2)'
    srtcors = sort(cors)
    mdcors[md] = filter(!isnan, cors)
    
    for (key, pos) in pairs(neuroactive)
        allcors = filter(!isnan, cors[pos])
        notcors = filter(!isnan, cors[Not(pos)])
        length(allcors) < 4 && continue
        @info "    $key"
        mwu = MannWhitneyUTest(allcors, notcors)
        u = mwu.U
        f = u / (mwu.nx * mwu.ny)
        m = median(allcors)
        p = pvalue(mwu)
        push!(allfsea, (geneset=key, metadatum=String(md), median=m, U=u, f=f, pvalue=p, cors=allcors, ext=extrema(mdcors[md])))
    end
end


allfsea.qvalue = adjust(allfsea.pvalue, BenjaminiHochberg())
sort!(allfsea, :qvalue)
allfsea.allmed = [median(mdcors[Symbol(md)]) for md in allfsea.metadatum]

CSV.write(joinpath(outdir, "12m_fsea.csv"), allfsea[!, Not(:cors)])

##

for md in metadatums
    @info "Plotting $md"
    df = filter(row-> row.qvalue < 0.2 && row.metadatum == string(md), allfsea)
    nrow(df) == 0 && continue
    fig = Figure(resolution = (900, 1200))
    
    sq = ceil(Int, sqrt(nrow(df)))
    coln = 1
    rown = 1

    for row in eachrow(df)
        ax = Axis(fig[rown,coln], yticks=-0.3:0.1:0.3, yminorticks=IntervalsBetween(2), yminorticksvisible = true)
        hidexdecorations!(ax)
        clr = row.median < row.allmed ? :purple : :orange

        hlines!(ax, row.cors, color=clr)
        hlines!(ax, [row.allmed], color=:black, linestyle=:dash)
        hlines!(ax, [row.ext...], color=:black, linestyle=:dot)
        coln += 1
        Label(fig[rown,coln], string(row.geneset, " | f = ", round(row.f, digits=3), " | q = ", round(row.qvalue, digits=4)), rotation=-pi/2, tellheight=false)
        coln += 1
        if coln > sq * 2
            rown += 1
            coln = 1
        end
    end


    leg = Legend(fig[:,end+1],
        [
            LineElement(color = :purple, linestyle = nothing),
            LineElement(color = :orange, linestyle = nothing),
            LineElement(color = :black, linestyle = :dash),
            LineElement(color = :black, linestyle = :dot),
        ], ["negative assoc.", "positive assoc.", "median (all genes)", "extrema"]
    )
    supertitle = Label(fig[0, 1:end-1], "12M FSEA for $md", textsize = 20)
    save(joinpath(figsdir, "12M_fsea_$md.pdf"), fig)
    fig
end
fig