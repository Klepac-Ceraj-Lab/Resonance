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
    readdir(joinpath(ENV["ANALYSIS_FILES"], "metaphlan"), join=true)), :species)

sns = map(samplenames(spec)) do s
    replace(s, r"_S\d+_profile"=>"")
end

idx = unique([findfirst(==(x), sns) for x in sns])

spec = CommunityProfile(abundances(spec)[:, idx], features(spec), MicrobiomeSample.(sns[idx]))

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

yidx = [any(col-> get(spec, s, col, false), [:age0to3mo]) for s in samplenames(spec)]
oidx = [any(col-> get(spec, s, col, false), [:age12moplus]) for s in samplenames(spec)]

ratios = map(featurenames(spec)) do feature
    log(mean(abundances(spec[feature,oidx])) / mean(abundances(spec[feature, yidx])))
end

onlyyoung = count(isnan, ratios)
onlyold = count(isinf, ratios)

filter!(x-> !isinf(x) && !isnan(x), ratios)

#

fig = Figure(resolution=(1200, 800))
yngax = Axis(fig[1,1], title="Only in <3 mo", ylabel="# of species")
bothax = Axis(fig[1, 2:3], title="Present in both", xlabel="Log ratio abundances")
oldax = Axis(fig[1,4], title="Only in >12 mo")

barplot!(yngax, [1], [onlyyoung])
hist!(bothax, ratios)
barplot!(oldax, [1], [onlyold])

hidexdecorations!.([yngax, oldax])
bothax.ygridvisible = true
oldax.ygridvisible = true

tightlimits!(yngax, Bottom())
tightlimits!(bothax, Left(), Right(), Bottom())
tightlimits!(oldax, Bottom())

linkyaxes!(yngax, bothax, oldax)
save("slides/assets/age_taxon_ratios.png", fig)
fig
## 

metabs = CSV.read("data/metabolites.csv", DataFrame)



##

humann_join(joinpath(ENV["ANALYSIS_FILES"], "humann", "genefamilies"), 
            joinpath(ENV["ANALYSIS_FILES"], "humann", "all_genefamilies.tsv"); 
            file_name="genefamilies")

func = humann_profiles(joinpath(ENV["ANALYSIS_FILES"], "humann", "all_genefamilies.tsv"), stratified=false)

sns = map(samplenames(func)) do s
    replace(s, r"_S\d+_genefamilies"=>"")
end

idx = unique([findfirst(==(x), sns) for x in sns])

func = CommunityProfile(abundances(func)[:, idx], features(func), sns[idx])

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