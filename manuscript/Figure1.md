# Figure 1 - Cohort Characteristics and Data Summaries

First, load packages that will be used throughout this notebook.

```julia
using Resonance
using FileIO
using CairoMakie # for plotting
using Distances
using MultivariateStats
using CategoricalArrays
```

Then, we'll load in the different data sources.

```julia
mdata = Resonance.load(Metadata())
subset!(mdata, "ageMonths"=> ByRow(<(120)))

taxa = Resonance.load(TaxonomicProfiles(); timepoint_metadata = mdata)
taxa = taxa[:, findall(!ismissing, get(taxa, :ageMonths))]
species = filter(t-> taxrank(t) == :species, taxa)

unirefs = Resonance.load(UnirefProfiles(); timepoint_metadata = mdata) # this can take a bit
unirefs = unirefs[:, findall(!ismissing, get(unirefs, :ageMonths))]
unirefs = filter(!hastaxon, unirefs) # don't use species stratification for summaries

ecs = Resonance.load(ECProfiles(); timepoint_metadata = mdata)
ecs = ecs[:, findall(!ismissing, get(ecs, :ageMonths))]
ecs = filter(!hastaxon, ecs)

kos = Resonance.load(KOProfiles(); timepoint_metadata = mdata)
kos = kos[:, findall(!ismissing, get(kos, :ageMonths))]
kos = filter(!hastaxon, kos)

# metabolites = Resonance.load(MetabolicProfiles(); timepoint_metadata = mdata)
brain = Resonance.load(Neuroimaging(), timepoint_metadata = mdata)


@assert all(samplenames(species) .== samplenames(unirefs))
@assert all(samplenames(species) .== samplenames(ecs))
@assert all(samplenames(species) .== samplenames(kos))

uidx = let md = DataFrame(Microbiome.metadata(taxa))
    smp = Set(unique(md, :subject).sample)
    findall(s-> s in smp, samplenames(taxa))
end

buidx = let md = DataFrame(Microbiome.metadata(brain))
    grp = groupby(md, :subject)
    smp = map(keys(grp)) do g
        idx = findfirst(grp[g].hassample)
        return isnothing(idx) ? first(grp[g].sample) : grp[g].sample[idx]
    end
    findall(s-> s in smp, samplenames(brain))
end
```

Both PERMANOVA and Mantel tests
require beta diversity distance metrics.
We'll calculate them once for each profile rather
rather than separately for each test.

```julia
spedm = CSV.read(outputfiles("spedm.csv"), DataFrame) |> Matrix
unidm = CSV.read(outputfiles("unidm.csv"), DataFrame) |> Matrix
ecsdm = CSV.read(outputfiles("ecsdm.csv"), DataFrame) |> Matrix
kosdm = CSV.read(outputfiles("kosdm.csv"), DataFrame) |> Matrix
# metdm = CSV.read(outputfiles("metdm.csv"), DataFrame) |> Matrix
brndm = CSV.read(outputfiles("brndm.csv"), DataFrame) |> Matrix
```


Then, we will start constructing the figure.
See the [Makie documentation](https://makie.juliaplots.org/stable/tutorials/layout-tutorial/) for more information.


```julia
figure = Figure(; resolution = (2000, 1200))
A = GridLayout(figure[1,1]; alignmode=Outside())
BC = GridLayout(figure[2,1:2]; alignmode=Outside())
DEF = GridLayout(figure[1,2]; alignmode=Outside());
```


## Summaries

### 1A - Cohort diagram and ages

For now, the graphic here is a place-holder, but we'll have something like this for 1A.

```julia
A_img = Axis(A[2,1:2]; aspect = DataAspect(), alignmode=Inside())
A_histleft = Axis(A[1,1]; ylabel="Samples (N)", alignmode=Mixed(; left=45))
A_histright = Axis(A[1,2]; xlabel = "Age (years)", xticks=2:2:16, alignmode=Inside())
hidedecorations!(A_img)
hidespines!(A_img)

image!(A_img, rotr90(load("figures/fig1a_placeholder.png")))

hist!(A_histleft, filter(<=(24), mdata.ageMonths); color = :darkgray)
hist!(A_histright, filter(>(24), mdata.ageMonths) ./ 12; color = :darkgray, bins=8)

rowgap!(A, -30)
rowsize!(A, 1, Relative(1/5))
colsize!(A, 2, Relative(1/3))
linkyaxes!(A_histleft, A_histright)

colsize!(figure.layout, 1, Relative(3/7))
A_histleft.xticks = ([0,3,6,12], ["birth", "3m", "6m", "12m"])
figure
```

### 1B-C: Omnibus tests

#### PERMANOVA

This is a permutation test of variance

```julia
B = GridLayout(BC[1,1])
# Ba = Axis(B[1:2,1]; alignmode=Outside())
Bb = Axis(B[1,1]; title="Under 6mo")
Bc = Axis(B[1,2]; title="Over 18mo")
```

```julia
# plot_permanovas!(Ba, CSV.read(outputfiles("permanovas_all.csv"), DataFrame))
plot_permanovas!(Bb, CSV.read(outputfiles("permanovas_u6mo.csv"), DataFrame))
hideydecorations!(Bc)
plot_permanovas!(Bc, CSV.read(outputfiles("permanovas_o18mo.csv"), DataFrame))

colsize!(BC, 1, Relative(1/2))
Label(BC[0,1], "PERMANOVAs")
```


### 1C - Mantel tests

```julia
C = GridLayout(BC[1,2])

# Ca = Axis(C[1:2, 1]; alignmode=Outside())
Cb = Axis(C[1,1]; title="Under 6mo")
Cc = Axis(C[1,2]; title="Over 18mo")
hideydecorations!(Cc)

# plot_mantel!(Ca, CSV.read(outputfiles("mantel_all.csv"), DataFrame))
plot_mantel!(Cb, CSV.read(outputfiles("mantel_u6.csv"), DataFrame))
plot_mantel!(Cc, CSV.read(outputfiles("mantel_o18.csv"), DataFrame))

Label(BC[0,2], "Mantel"; tellwidth=false)
```

## Ordinations

```julia
D = Axis(DEF[1,1]; title = "Taxa")
spepco = fit(MDS, spedm; distances=true)
sc = plot_pcoa!(D, spepco; color=get(species, :ageMonths))

E = GridLayout(DEF[2,1])
Ea = Axis(E[1,1]; title = "Bacteroidetes")
Eb = Axis(E[1,2]; title = "Firmicutes")
# Ec = Axis(E[1,3]; title = "Actinobacteria")

plot_pcoa!(Ea, spepco; color=vec(abundances(filter(t-> taxrank(t) == :phylum, taxa)[r"Bacteroidetes", :])), colormap=:Purples)
plot_pcoa!(Eb, spepco; color=vec(abundances(filter(t-> taxrank(t) == :phylum, taxa)[r"Firmicutes", :])), colormap=:Purples)
hideydecorations!(Eb)
# pco = plot_pcoa!(Ec, spepco; color=vec(abundances(filter(t-> taxrank(t) == :phylum, taxa)[r"Actinobacteria", :])), colormap=:Purples)
F = Axis(DEF[1,2]; title = "Functions")

unipco = fit(MDS, unidm; distances=true)
plot_pcoa!(F, unipco; color=get(unirefs, :ageMonths))

G = Axis(DEF[2,2]; title = "Neuroimaging")

brnpco = fit(MDS, brndm; distances=true)
plot_pcoa!(G, brnpco; color=get(brain, :ageMonths))

Colorbar(DEF[1:2, 3], sc; label="Age (months)", flipaxis=true)
```

```julia

Label(A[1, 1, TopLeft()], "A",
        textsize = 26,
        font = "Open Sans Bold",
        padding = (0, 5, 5, 0),
        halign = :right
)
Label(DEF[1, 1, TopLeft()], "B",
        textsize = 26,
        font = "Open Sans Bold",
        padding = (0, 5, 5, 0),
        halign = :right
)
Label(DEF[2, 1, TopLeft()], "C",
        textsize = 26,
        font = "Open Sans Bold",
        padding = (0, 5, 5, 0),
        halign = :right
)
Label(DEF[1, 2, TopLeft()], "D",
        textsize = 26,
        font = "Open Sans Bold",
        padding = (0, 5, 5, 0),
        halign = :right
)
Label(DEF[2, 2, TopLeft()], "E",
        textsize = 26,
        font = "Open Sans Bold",
        padding = (0, 5, 5, 0),
        halign = :right
)

Label(B[1, 1, TopLeft()], "F",
        textsize = 26,
        font = "Open Sans Bold",
        padding = (0, 5, 5, 0),
        halign = :right
)
Label(C[1, 1, TopLeft()], "G",
        textsize = 26,
        font = "Open Sans Bold",
        padding = (0, 5, 5, 0),
        halign = :right
)

save(figurefiles("Figure1.svg"), figure)
save(figurefiles("Figure1.png"), figure)
figure
```

![](figures/Figure1.png)

## Supplement

### More Taxa PCoA

```julia
fig = Figure()
ax = Axis(fig[1,1], title = "Bacteroidetes")
ax2 = Axis(fig[1,2], title = "Prevotella")
ax3 = Axis(fig[2,1], title = "Bacteroides")
ax4 = Axis(fig[2,2], title = "Alistipes")
plot_pcoa!(ax, spepco; color=vec(abundances(filter(t-> taxrank(t) == :phylum, taxa)[r"Bacteroidetes", :])),
        colormap=:Purples,
        strokecolor=:black,
        strokewidth=1
)
plot_pcoa!(ax2, spepco; color=vec(abundances(filter(t-> taxrank(t) == :genus, taxa)[r"Prevotella", :])),
        colormap=:Purples,
        strokecolor=:black,
        strokewidth=1
)
plot_pcoa!(ax3, spepco; color=vec(abundances(filter(t-> taxrank(t) == :genus, taxa)[r"Bacteroides", :])),
        colormap=:Purples,
        strokecolor=:black,
        strokewidth=1
)
plot_pcoa!(ax4, spepco; color=vec(abundances(filter(t-> taxrank(t) == :genus, taxa)[r"Alistipes", :])),
        colormap=:Purples,
        strokecolor=:black,
        strokewidth=1
)
fig
```

### Brain PCoA

```julia
fig = Figure()
brain_pco = fit(MDS, brndm; distances=true)
ax = Axis(fig[1,1])
plot_pcoa!(ax, brain_pco; color=brain.AgeInDays)
fig
```

### Multivariate permanovas

Maternal education & race.

```julia
using PERMANOVA

df = DataFrame(education = get(species, :education)[uidx], 
               race       = get(species, :race)[uidx])
oidx = findall(.!ismissing.(df.education) .& .!ismissing.(df.race))
df = df[oidx, :]

permanova(df, spedm[uidx[oidx], uidx[oidx]], @formula(1 ~ education + race), 1000)

```