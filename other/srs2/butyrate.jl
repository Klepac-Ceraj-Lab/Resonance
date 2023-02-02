using Microbiome
using Microbiome.Dictionaries
using BiobakeryUtils
using ResonanceMicrobiome
using ECHOSRS2
using CairoMakie
using CSV
using CSV.Tables
using CodecZlib

##

neuro = get_neuroactive_kos()
# butdeg = neuro["Butyrate degradation"]
butsyn = neuro["Butyrate synthesis"]

kos2uniref = Dictionary{String, Vector{String}}()
for line in eachline(GzipDecompressorStream(open("data/map_ko_uniref90.txt.gz")))
    line = split(line, '\t')
    insert!(kos2uniref, line[1], map(x-> String(match(r"UniRef90_(\w+)", x).captures[1]), line[2:end]))
end

butuniref = Set(Iterators.flatten(kos2uniref[s] for s in butsyn))

##

mss = map(CSV.File("data/NHBCS.DATA.MI.csv")) do row
    MicrobiomeSample(row.sfdid, Dict(k=> v for (k, v) in pairs(row) if k != :sfdid) |> dictionary)
end

##

files = readdir("/dartfs-hpc/rc/lab/M/MRKepistor7/collab/KevinBonham/SourceFiles/gene_families/raw_output/", join=true)
filter!(f-> occursin(r"_genefamilies\.tsv\.gz", f), files)

profs = CommunityProfile[]

for f in files
    @info f
    samp = first(split(basename(f), "_"))
    sample = findfirst(s-> name(s) == samp, mss)
    isnothing(sample) && continue

    sample = mss[sample]
    gfs = GeneFunction[]
    abunds = Float64[]

    for row in CSV.File(GzipDecompressorStream(open(f)))
        m = match(r"UniRef90_(\w+)", row[1])
        isnothing(m) && continue
        string(m.captures[1]) ∈ butuniref || continue

        genestring = replace(row[1], r"g__\w+\.s__"=> "s__")
        push!(gfs, genefunction(genestring))
        push!(abunds, row[2])
    end
    
    push!(profs, CommunityProfile(reshape(abunds, length(abunds), 1), gfs, [sample]))
end 

##

comm = commjoin(profs...)
commspec = filter(hastaxon, comm)

##

utaxa = unique(name.(taxon.(features(commspec))))
totals = map(utaxa) do t
    c = filter(f-> name(taxon(f)) == t, commspec)
    sum(featuretotals(c))
end

srt = sortperm(totals)

##

extrema(totals)

fig = Figure()
ax = Axis(fig[1,1], title="Total Contributions", xticks=(1:length(totals), utaxa[srt]), 
            xticklabelrotation=π/2,
            ylims = (1, 1.5e6),
            yscale=Makie.pseudolog10,
            yminorticksvisible=true,
            yminorticks=IntervalsBetween(10),
            yticks = 10 .^ collect(0:6))
barplot!(ax, 1:length(totals), totals[srt])

fig

##

ecoli = vec(sampletotals(filter(f-> contains(name(taxon(f)), "coli"), commspec)))
findall(>(0), ecoli)