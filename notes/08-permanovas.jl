using Resonance
using CairoMakie

omni, etoh, tps, complete_brain, metabolites, species = startup()
genes = Resonance.read_gfs_arrow()
ecs = Resonance.read_gfs_arrow(kind="ecs_rename")
ecs = filter(!hastaxon, ecs)
kos = Resonance.read_gfs_arrow(kind="kos_rename")
kos = filter(!hastaxon, kos)


unique!(tps, ["subject", "timepoint"])

set!(species, leftjoin(select(omni, ["subject", "timepoint", "sample"]), tps, on=[:subject, :timepoint]))
set!(genes, leftjoin(select(omni, ["subject", "timepoint", "sample"]), tps, on=[:subject, :timepoint]))
set!(metabolites, leftjoin(select(etoh, ["subject", "timepoint", "sample"]), tps, on=[:subject, :timepoint]))

#-

specu6mo = species[:, map(a-> !ismissing(a) && a < 6, get(species, :ageMonths))]
specu6mo = specu6mo[vec(prevalence(specu6mo)) .> 0, :]

genesu6mo = genes[:, map(a-> !ismissing(a) && a < 6, get(genes, :ageMonths))]
genesu6mo = genesu6mo[vec(prevalence(genesu6mo)) .> 0, :]

metabu6mo = metabolites[:, map(a-> !ismissing(a) && a < 6, get(metabolites, :ageMonths))]
metabu6mo = metabu6mo[vec(prevalence(metabu6mo)) .> 0, :]


#-

fig, ax, hm = plot_permanovas([specu6mo, genesu6mo, metabu6mo], [:ageMonths, :cogScore, :rce, :ed];
    commlabels = ["species", "genes", "metabolites"],
    mdlabels = ["age", "cogScore", "race", "education"]
)

Label(fig[0,1], "Kids under 6 mo"; tellwidth=false)
save("figures/kids_under6mo_permanovas.png", fig)

fig

#-

speco1y = species[:, map(a-> !ismissing(a) && a > 12, get(species, :ageMonths))]
speco1y = speco1y[:, unique(DataFrame(metadata(speco1y)), [:subject, :timepoint]).sample]
speco1y = speco1y[vec(prevalence(speco1y)) .> 0, :]

geneso1y = genes[:, map(a-> !ismissing(a) && a > 12, get(genes, :ageMonths))]
geneso1y = geneso1y[:, unique(DataFrame(metadata(geneso1y)), [:subject, :timepoint]).sample]
geneso1y = geneso1y[vec(prevalence(geneso1y)) .> 0, :]


#-

fig, ax, hm = plot_permanovas([speco1y, geneso1y], [:ageMonths, :cogScore, :rce, :ed];
    commlabels = ["species", "genes"],
    mdlabels = ["age", "cogScore", "race", "education"]
)

Label(fig[0,1], "Kids over 1 year"; tellwidth=false)
save("figures/kids_over1y_permanovas.png", fig)

fig

#-

speco1y = species[:, map(a-> !ismissing(a) && a > 12, get(species, :ageMonths))]
speco1y = speco1y[:, unique(DataFrame(metadata(speco1y)), [:subject, :timepoint]).sample]
speco1y = speco1y[vec(prevalence(speco1y)) .> 0, :]

geneso1y = genes[:, map(a-> !ismissing(a) && a > 12, get(genes, :ageMonths))]
geneso1y = geneso1y[:, unique(DataFrame(metadata(geneso1y)), [:subject, :timepoint]).sample]
geneso1y = geneso1y[vec(prevalence(geneso1y)) .> 0, :]


#-

fig, ax, hm = plot_permanovas([speco1y, geneso1y], [:ageMonths, :cogScore, :rce, :ed];
    commlabels = ["species", "genes"],
    mdlabels = ["age", "cogScore", "race", "education"]
)

Label(fig[0,1], "Kids over 1 year"; tellwidth=false)
save("figures/kids_over1y_permanovas.png", fig)

fig

#-

specallkids = species[:, map(a-> !ismissing(a), get(species, :ageMonths))]
specallkids = specallkids[:, unique(DataFrame(metadata(specallkids)), ["subject", "timepoint"]).sample]
specallkids = specallkids[vec(prevalence(specallkids)) .> 0, :]

genesallkids = genes[:, map(a-> !ismissing(a), get(genes, :ageMonths))]
genesallkids = genesallkids[:, unique(DataFrame(metadata(genesallkids)), ["subject", "timepoint"]).sample]
genesallkids = genesallkids[vec(prevalence(genesallkids)) .> 0, :]

metaballkids = metabolites[:, map(a-> !ismissing(a), get(metabolites, :ageMonths))]
metaballkids = metaballkids[:, unique(DataFrame(metadata(metaballkids)), ["subject", "timepoint"]).sample]
metaballkids = metaballkids[vec(prevalence(metaballkids)) .> 0, :]


#-

fig, ax, hm = plot_permanovas([specallkids, genesallkids, metaballkids], [:ageMonths, :cogScore, :rce, :ed];
    commlabels = ["species", "genes", "metabolites"],
    mdlabels = ["age", "cogScore", "race", "education"]
)

Label(fig[0,1], "All kids"; tellwidth=false)
save("figures/kids_all_permanovas.png", fig)

fig

get(specallkids, :ageMonths) |> extrema # age range in months
nsamples(specallkids)
