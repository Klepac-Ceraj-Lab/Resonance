using Microbiome
using CSV
using DataFrames
using ECHOSRS2
using SparseArrays
using ProgressLogging
using Dictionaries
using MultipleTesting
using Statistics
using HypothesisTests
using CairoMakie

outdir = "/grace/facstaff/kevin/srs2/"
figsdir = joinpath(outdir, "figures")
isdir(figsdir) || mkpath(figsdir)

meta = CSV.read("data/VKC.DATA.MI.csv", DataFrame)

gfs = DataFrame()

for batch in 1:14
    b = lpad(batch, 3, '0')
    @info "working on batch$b"
    for (root, dirs, files) in walkdir("/grace/echo/analysis/biobakery3/batch$b/")
        occursin("humann", root) || continue
        occursin("main", root) || continue
        filter!(f-> occursin(r"_genefamilies\.tsv", f), files)
        
        for f in joinpath.(Ref(root), files)
            s = first(basename(f), 11)
            s = replace(s, "-"=>"_")
            row = findfirst(==(s), meta.MBL_ID)
            isnothing(row) && continue

            sample = MicrobiomeSample(s)
            set!(sample, :batch, batch)
            for col in propertynames(meta)[2:end]
                set!(sample, col, meta[row, col])
            end

            gf = CSV.read(f, DataFrame)
            rename!(gf, ["feature", "abundance"])
            gf.sample = fill(sample, nrow(gf))
            

            filter!("feature" => (x-> !occursin(r"\|[gu]", x)), gf)
            gf.feature = GeneFunction.(gf.feature)
            append!(gfs, gf)
        end
    end
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
              :RESID]

##



youngkid = filter(:REALTIME=> (a-> a == "6W" || a == "4M"), meta)
allfsea = DataFrame(
            geneset   = String[],
            metadatum = String[],
            median    = Float64[],
            pvalue    = Float64[],
            cors      = Vector{Float64}[],
            ext       = Tuple[])

mdcors = Dict(m=>Float64[] for m in metadatums)


genesets = Set([
    "Propionate degradation",
    "Propionate synthesis",
    "Acetate synthesis",
    "Butyrate synthesis",
])


for md in metadatums
    @info "Working on $md"
    filt = map(!ismissing, youngkid[!,md])
    cors = cor(youngkid[filt, md], abundances(unirefaccessory[:, youngkid.MBL_ID])[:,filt], dims=2)'
    srtcors = sort(cors)
    mdcors[md] = filter(!isnan, cors)
    
    for (key, pos) in pairs(neuroactive)
        in(key, genesets) || continue
        allcors = filter(!isnan, cors[pos])
        notcors = filter(!isnan, cors[Not(pos)])
        length(allcors) < 4 && continue
        @info "    $key"
        mwu = MannWhitneyUTest(allcors, notcors)
        m = median(allcors)
        p = pvalue(mwu)
        push!(allfsea, (geneset=key, metadatum=String(md), median=m, pvalue=p, cors=allcors, ext=extrema(mdcors[md])))
    end
end

allfsea.qvalue = adjust(allfsea.pvalue, BenjaminiHochberg())
allfsea.allmed = [median(mdcors[Symbol(md)]) for md in allfsea.metadatum]
sort!(allfsea, :qvalue)
CSV.write(joinpath(outdir, "young_kids_fsea.csv"), allfsea[!, Not(:cors)])

##

fig = Figure(resolution = (1600, 600))
coln = 1
for row in eachrow(allfsea)
    ax = Axis(fig[1,coln], title=row.metadatum, yticks=-0.6:0.2:0.8, yminorticks=IntervalsBetween(2), yminorticksvisible = true)
    hidexdecorations!(ax)
    clr = row.median < row.allmed ? :purple : :orange

    hlines!(ax, row.cors, color=clr)
    hlines!(ax, [row.allmed], color=:black, linestyle=:dash)
    hlines!(ax, [row.ext...], color=:black, linestyle=:dot)
    coln += 1
    Label(fig[1,coln], string(row.geneset, " | p = ", round(row.pvalue, digits=4)), rotation=-pi/2, tellheight=false)
    coln += 1
end

leg = Legend(fig[1,coln],
    [
        LineElement(color = :purple, linestyle = nothing),
        LineElement(color = :orange, linestyle = nothing),
        LineElement(color = :black, linestyle = :dash),
        LineElement(color = :black, linestyle = :dot),
    ], ["negative assoc.", "positive assoc.", "median (all genes)", "extrema"]
)

fig
##

save(joinpath(figsdir, "young_kids_fsea.pdf"), fig)

##



oldkid = filter(:REALTIME=> (a-> !(a == "6W" || a == "4M" || a == "6M")), meta)
allfsea = DataFrame(
            geneset   = String[],
            metadatum = String[],
            median    = Float64[],
            pvalue    = Float64[],
            cors      = Vector{Float64}[],
            ext       = Tuple[])

mdcors = Dict(m=>Float64[] for m in metadatums)


genesets = Set([
    "Menaquinone synthesis",
    "Isovaleric acid synthesis",
])


for md in metadatums
    @info "Working on $md"
    filt = map(!ismissing, oldkid[!,md])
    cors = cor(oldkid[filt, md], abundances(unirefaccessory[:, oldkid.MBL_ID])[:,filt], dims=2)'
    srtcors = sort(cors)
    mdcors[md] = filter(!isnan, cors)
    
    for (key, pos) in pairs(neuroactive)
        if !in(key, genesets)
            @info "skipping $key"
            continue
        end
        allcors = filter(!isnan, cors[pos])
        notcors = filter(!isnan, cors[Not(pos)])
        length(allcors) < 4 && continue
        @info "    $key"
        mwu = MannWhitneyUTest(allcors, notcors)
        m = median(allcors)
        p = pvalue(mwu)
        push!(allfsea, (geneset=key, metadatum=String(md), median=m, pvalue=p, cors=allcors, ext=extrema(mdcors[md])))
    end
end

allfsea.qvalue = adjust(allfsea.pvalue, BenjaminiHochberg())
allfsea.allmed = [median(mdcors[Symbol(md)]) for md in allfsea.metadatum]
sort!(allfsea, :qvalue)
CSV.write(joinpath(outdir, "old_kids_fsea.csv"), allfsea[!, Not(:cors)])

##

fig = Figure(resolution = (800, 600))
coln = 1
for row in eachrow(allfsea)
    ax = Axis(fig[1,coln], title=row.metadatum, yticks=-0.6:0.2:0.8, yminorticks=IntervalsBetween(2), yminorticksvisible = true)
    hidexdecorations!(ax)
    clr = row.median < row.allmed ? :purple : :orange

    hlines!(ax, row.cors, color=clr)
    hlines!(ax, [row.allmed], color=:black, linestyle=:dash)
    hlines!(ax, [row.ext...], color=:black, linestyle=:dot)
    coln += 1
    Label(fig[1,coln], string(row.geneset, " | p = ", round(row.pvalue, digits=4)), rotation=-pi/2, tellheight=false)
    coln += 1
end

leg = Legend(fig[1,coln],
    [
        LineElement(color = :purple, linestyle = nothing),
        LineElement(color = :orange, linestyle = nothing),
        LineElement(color = :black, linestyle = :dash),
        LineElement(color = :black, linestyle = :dot),
    ], ["negative assoc.", "positive assoc.", "median (all genes)", "extrema"]
)

fig
##

save(joinpath(figsdir, "old_kids_fsea.pdf"), fig)

##