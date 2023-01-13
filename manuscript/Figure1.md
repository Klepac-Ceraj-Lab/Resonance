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

species = filter(t-> taxrank(t) == :species, taxa)
unirefs = Resonance.load(UnirefProfiles(); timepoint_metadata = mdata) # this can take a bit
ecs = Resonance.load(ECProfiles(); timepoint_metadata = mdata)
kos = Resonance.load(KOProfiles(); timepoint_metadata = mdata)
# metabolites = Resonance.load(MetabolicProfiles(); timepoint_metadata = mdata)
brain = Resonance.load(Neuroimaging(), timepoint_metadata = mdata)

@assert all(samplenames(species) .== samplenames(unirefs) .== samplenames(ecs) .== samplenames(kos))
```

Both PERMANOVA and Mantel tests
require beta diversity distance metrics.
We'll calculate them once for each profile rather
rather than separately for each test.

```julia
spedm = CSV.read(scratchfiles("spedm.csv"), DataFrame) |> Matrix
unidm = CSV.read(scratchfiles("unidm.csv"), DataFrame) |> Matrix
ecsdm = CSV.read(scratchfiles("ecsdm.csv"), DataFrame) |> Matrix
kosdm = CSV.read(scratchfiles("kosdm.csv"), DataFrame) |> Matrix
# metdm = CSV.read(scratchfiles("metdm.csv"), DataFrame) |> Matrix
brndm = CSV.read(scratchfiles("brndm.csv"), DataFrame) |> Matrix
```


Then, we will start constructing the figure.
See the [Makie documentation](https://makie.juliaplots.org/stable/tutorials/layout-tutorial/) for more information.


```julia
figure = Figure(; resolution = (2000, 1200))
A = GridLayout(figure[1,1]; alignmode=Outside())
BCDE = GridLayout(figure[1,2]; alignmode=Outside());
FG = GridLayout(figure[2,1:2]; alignmode=Outside())
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

# image!(A_img, rotr90(load("figures/fig1a_placeholder.png")))

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


### Ordinations

```julia
B = Axis(BCDE[1,1]; title = "Taxa")
spepco = fit(MDS, spedm; distances=true)
sc = plot_pcoa!(B, spepco; color=get(species, :ageMonths))

C = GridLayout(BCDE[2,1])
Ca = Axis(E[1,1]; title = "Bacteroidetes")
Cb = Axis(E[1,2]; title = "Firmicutes")
# Cc = Axis(E[1,3]; title = "Actinobacteria")

# plot_pcoa!(Ca, spepco; color=vec(abundances(filter(t-> taxrank(t) == :phylum, taxa)[r"Bacteroidetes", :])), colormap=:Purples)
# plot_pcoa!(Cb, spepco; color=vec(abundances(filter(t-> taxrank(t) == :phylum, taxa)[r"Firmicutes", :])), colormap=:Purples)
hideydecorations!(Cb)
# pco = plot_pcoa!(Ec, spepco; color=vec(abundances(filter(t-> taxrank(t) == :phylum, taxa)[r"Actinobacteria", :])), colormap=:Purples)
D = Axis(BCDE[1,2]; title = "Functions")

unipco = fit(MDS, unidm; distances=true)
plot_pcoa!(D, unipco; color=get(unirefs, :ageMonths))

E = Axis(BCDE[2,2]; title = "Neuroimaging")

brnpco = fit(MDS, brndm; distances=true)
plot_pcoa!(E, brnpco; color=get(brain, :ageMonths))

Colorbar(BCDE[1:2, 3], sc; label="Age (months)", flipaxis=true)

```


### 1F-G: Omnibus tests

#### 1F - PERMANOVA

This is a permutation test of variance

```julia
F = GridLayout(FG[1,1])
# Ba = Axis(B[1:2,1]; alignmode=Outside())
Fa = Axis(F[1,1]; title="Under 6mo")
Fb = Axis(F[1,2]; title="Over 18mo")
```

```julia
# plot_permanovas!(Ba, CSV.read(scratchfiles("permanovas_all.csv"), DataFrame))
plot_permanovas!(Fa, CSV.read(scratchfiles("permanovas_00to06.csv"), DataFrame))
hideydecorations!(Fb)
plot_permanovas!(Fb, CSV.read(scratchfiles("permanovas_18to120.csv"), DataFrame))

colsize!(FG, 1, Relative(1/2))
Label(FG[0,1], "PERMANOVAs")
```


#### 1G - Mantel tests

```julia
G = GridLayout(FG[1,2])

# Ca = Axis(C[1:2, 1]; alignmode=Outside())
Ga = Axis(G[1,1]; title="Under 6mo")
Gb = Axis(G[1,2]; title="Over 18mo")
hideydecorations!(Gb)

# plot_mantel!(Ca, CSV.read(scratchfiles("mantel_all.csv"), DataFrame))
plot_mantel!(Ga, CSV.read(scratchfiles("mantel_00to06.csv"), DataFrame))
plot_mantel!(Gb, CSV.read(scratchfiles("mantel_18to120.csv"), DataFrame))

Label(FG[0,2], "Mantel"; tellwidth=false)
```

### Labels & saving

```julia

Label(A[1, 1, TopLeft()], "A",
        textsize = 26,
        font = "Open Sans Bold",
        padding = (0, 5, 5, 0),
        halign = :right
)
Label(BCDE[1, 1, TopLeft()], "B",
        textsize = 26,
        font = "Open Sans Bold",
        padding = (0, 5, 5, 0),
        halign = :right
)
Label(BCDE[2, 1, TopLeft()], "C",
        textsize = 26,
        font = "Open Sans Bold",
        padding = (0, 5, 5, 0),
        halign = :right
)
Label(BCDE[1, 2, TopLeft()], "D",
        textsize = 26,
        font = "Open Sans Bold",
        padding = (0, 5, 5, 0),
        halign = :right
)
Label(BCDE[2, 2, TopLeft()], "E",
        textsize = 26,
        font = "Open Sans Bold",
        padding = (0, 5, 5, 0),
        halign = :right
)

Label(F[1, 1, TopLeft()], "F",
        textsize = 26,
        font = "Open Sans Bold",
        padding = (0, 5, 5, 0),
        halign = :right
)
Label(G[1, 1, TopLeft()], "G",
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