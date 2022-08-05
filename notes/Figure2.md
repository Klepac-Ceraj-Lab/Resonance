# Figure 2 - Per-feature, cross-sectional tests on cognitive function scores

First, load packages that will be used throughout this notebook.

```julia
using Resonance
using CairoMakie # for plotting
using GLM
```

Then, we'll load in the different data sources.

```julia
mtdt = Resonance.load(Metadata())

taxa = Resonance.load(TaxonomicProfiles(); timepoint_metadata = mtdt)
species = filter(t-> taxrank(t) == :species, taxa)

kos = Resonance.load(KOProfiles(); timepoint_metadata = mtdt)
kos = filter(!hastaxon, kos)

```

## Data filtering

Because there is a major shift in microbial composition
upon the introduction of solid foods,
we are going to split the datasets into stools collected prior to 6 months old
(most kids are on liquid diets, breast milk and/or formula)
and over 12 months old (most kids are eating at least some solid foods).

```julia
specmtdt = select(DataFrame(metadata(species)), ["subject", "timepoint", "ageMonths", "cogScore"])
specmtdt.sample = samplenames(species)
sort!(specmtdt, ["subject", "timepoint"])

specu6 = subset(specmtdt, "ageMonths" => ByRow(<(6)), "cogScore"=> ByRow(!ismissing))
unique!(specu6, "subject")

speco12 = subset(specmtdt, "ageMonths" => ByRow(>(12)), "cogScore"=> ByRow(!ismissing))
unique!(speco12, "subject")


komtdt = select(DataFrame(metadata(kos)), ["subject", "timepoint", "ageMonths", "cogScore"])
komtdt.sample = samplenames(kos)
sort!(komtdt, ["subject", "timepoint"])

kou6 = subset(komtdt, "ageMonths" => ByRow(<(6)), "cogScore"=> ByRow(!ismissing))
unique!(kou6, "subject")

koo12 = subset(komtdt, "ageMonths" => ByRow(>(12)), "cogScore"=> ByRow(!ismissing))
unique!(koo12, "subject")
```

## Adding features

```julia
for f in features(species)
    colu6 = vec(abundances(species[f, specu6.sample]))
    count(>(0), colu6) / length(colu6) > 0.1 && (specu6[!, name(f)] = colu6)

    colo12 = vec(abundances(species[f, speco12.sample]))
    count(>(0), colo12) / length(colo12) > 0.1 && (speco12[!, name(f)] = colo12)
end

for f in features(kos)
    colu6 = vec(abundances(kos[f, kou6.sample]))
    count(>(0), colu6) / length(colu6) > 0.1 && (kou6[!, name(f)] = colu6)

    colo12 = vec(abundances(kos[f, koo12.sample]))
    count(>(0), colo12) / length(colo12) > 0.1 && (koo12[!, name(f)] = colo12)
end

```

## GLMs




Then, we will start constructing the figure.
See the [Makie documentation](https://makie.juliaplots.org/stable/tutorials/layout-tutorial/) for more information.
