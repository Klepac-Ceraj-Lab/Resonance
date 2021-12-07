using Resonance
using CairoMakie

md = CSV.read("data/wrangled.csv", DataFrame)

##

reads = CSV.read("data/read_counts.csv", DataFrame)
(fig, ax, hs) = hist(reads.count, axis=(
    xlabel="Total reads (single)",
    ylabel="# of samples"
))

save("slides/assets/kneaddatacounts.png", fig)

##

@rsubset!(md, !ismissing(:subject),
                    !ismissing(:sample),
                    )
@rsubset!(md, !startswith(:sample, "C"),
                    !startswith(:sample, "z"),
                    )

spec = metaphlan_profiles(filter(f-> contains(f, "profile") && any(s-> contains(f, s), md.sample), 
    readdir("/grace/echo/analysis/biobakery3/links/metaphlan", join=true)), :species)

sns = map(samplenames(spec)) do s
    replace(s, r"_S\d+_profile"=>"")
end

idx = unique([findfirst(==(x), sns) for x in sns])

spec = CommunityProfile(abundances(spec)[:, idx], features(spec), map(samplenames(spec)[idx]) do s
    MicrobiomeSample(replace(s, r"_S\d+_profile"=>""))
end)

set!(spec, md)

for sampl in samples(spec)
    sn = name(sampl)
    set!(spec, sn, :richness, count(!=(0), abundances(spec[:, sampl])))
end


(fig, ax, hs) = hist([get(spec, sn, :richness) for sn in samplenames(spec)],
    axis=(
        xlabel="Richness (# of taxa)",
        ylabel="# of samples"
    ))

save("slides/assets/echo_richness.png", fig)
fig

##

fig = Figure()
ax = Axis(fig[1,1], ylabel="# of samples", xlabel = "Richness (# of taxa)")

idx = [any(col-> get(spec, s, col, false), [:age0to3mo, :age3to6mo, :age6to12mo]) for s in samplenames(spec)]
hist!(ax, [get(spec, sn, :richness, missing) for sn in samplenames(spec)][Not(idx)], label="Age > 12mo")
hist!(ax, [get(spec, sn, :richness, missing) for sn in samplenames(spec)][idx], label="Age < 12mo")
axislegend(ax)

save("slides/assets/echo_richness_byage.png", fig)
fig

## 

metabs = CSV.read("data/metabolites.csv", DataFrame)

aw

##

func = humann_profile.(filter(f-> contains(f, "genefamilies.tsv") && any(s-> contains(f, s), md.sample), 
    readdir("/grace/echo/analysis/biobakery3/links/humann/genefamilies", join=true)), stratified=false)
func = commjoin(func...)

sns = map(samplenames(func)) do s
    replace(s, r"_S\d+_genefamilies"=>"")
end

idx = unique([findfirst(==(x), sns) for x in sns])

func = CommunityProfile(abundances(func)[:, idx], features(func), map(samplenames(func)[idx]) do s
    MicrobiomeSample(replace(s, r"_S\d+_genefamilies"=>""))
end)

set!(func, md)

for sampl in samples(func)
    sn = name(sampl)
    set!(func, sn, :richness, count(!=(0), abundances(func[:, sampl])))
end

##

(fig, ax, hs) = hist([get(func, sn, :richness) for sn in samplenames(func)],
    axis=(
        xlabel="Richness (# of gene functions)",
        ylabel="# of samples"
    ))

save("slides/assets/echo_richness_func.png", fig)
fig

##

fig = Figure()
ax = Axis(fig[1,1], ylabel="# of samples", xlabel = "Richness (# of gene functions)")

idx = [any(col-> get(func, s, col, false), [:age0to3mo, :age3to6mo, :age6to12mo]) for s in samplenames(func)]
hist!(ax, [get(func, sn, :richness, missing) for sn in samplenames(func)][Not(idx)], label="Age > 12mo")
hist!(ax, [get(func, sn, :richness, missing) for sn in samplenames(func)][idx], label="Age < 12mo")
axislegend(ax)

save("slides/assets/echo_richness__func_byage.png", fig)
fig