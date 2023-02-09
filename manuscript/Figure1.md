# Figure 1 - Cohort Characteristics and Data Summaries

First, load packages that will be used throughout this notebook.

```julia
using Resonance
using FileIO
using CairoMakie # for plotting
using ColorSchemes
using Distances
using MultivariateStats
using CategoricalArrays
```

Then, we'll load in the different data sources.

```julia
mdata = Resonance.load(Metadata())

taxa = Resonance.load(TaxonomicProfiles(); timepoint_metadata = mdata) # this can take a bit
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

divr = shannon(species) |> vec
```


Then, we will start constructing the figure.
See the [Makie documentation](https://makie.juliaplots.org/stable/tutorials/layout-tutorial/) for more information.


```julia
figure = Figure(; resolution = (2000, 1200))
A = GridLayout(figure[1,1])
BCDE = GridLayout(figure[1,2])
F = GridLayout(figure[1,3])
GH = GridLayout(figure[2,1:3])
```


## Summaries

### 1A - Cohort diagram and ages

For now, the graphic here is a place-holder, but we'll have something like this for 1A.

```julia
A_img = Axis(A[2,1:3]; aspect = DataAspect(), alignmode=Inside())
A_histleft = Axis(A[1,2]; ylabel="Samples (N)")
A_histright = Axis(A[1,3]; xlabel = "Age (years)", xticks=2:2:16, yticklabelsvisible=false, yticksvisible=false)
hidedecorations!(A_img)
hidespines!(A_img)

image!(A_img, rotr90(load("manuscript/assets/Figure1A.jpg")))

hist!(A_histleft, filter(<=(24), mdata.ageMonths); color = :darkgray)
hist!(A_histright, filter(>(24), mdata.ageMonths) ./ 12; color = :darkgray, bins=8)

rowgap!(A, -30)
rowsize!(A, 1, Relative(1/5))
colsize!(A, 2, Relative(1/3))
linkyaxes!(A_histleft, A_histright)

A_histleft.xticks = ([3,6,12], ["", "", ""])
```


### Ordinations

```julia
spepco = fit(MDS, spedm; distances=true)

B = GridLayout(BCDE[1,1])
Ba = Axis(B[1,1]; ylabel= "Age (months)", xlabel = mdsaxis(spepco, 1), 
                  title = "Taxa", yminorticksvisible = true, yticks=0:24:120, yminorticks = IntervalsBetween(2),
                  alignmode = Mixed(; left=-15))
sc1 = scatter!(Ba, Resonance.loadings(spepco, 1), get(species, :ageMonths);
        color = divr, colormap=:plasma)
hlines!(Ba, [6, 18]; linestyle=:dash, color=:darkgray)
Colorbar(B[1, 2], sc1; label="Shannon Diversity", flipaxis=true)

C = Axis(BCDE[2,1]; title = "Taxa", alignmode = Mixed(; left=-15))
sc2 = plot_pcoa!(C, spepco; color=get(species, :ageMonths))


# C = GridLayout(BCDE[2,1])
# Ca = Axis(C[1,1]; title = "Bacteroidetes")
# Cb = Axis(C[1,2]; title = "Firmicutes")
# Cc = Axis(E[1,3]; title = "Actinobacteria")

# plot_pcoa!(Ca, spepco; color=vec(abundances(filter(t-> taxrank(t) == :phylum, taxa)[r"Bacteroidetes", :])), colormap=:Purples)
# plot_pcoa!(Cb, spepco; color=vec(abundances(filter(t-> taxrank(t) == :phylum, taxa)[r"Firmicutes", :])), colormap=:Purples)
# hideydecorations!(Cb)
# pco = plot_pcoa!(Ec, spepco; color=vec(abundances(filter(t-> taxrank(t) == :phylum, taxa)[r"Actinobacteria", :])), colormap=:Purples)
D = Axis(BCDE[1,2]; title = "Functions", alignmode = Mixed(; left=-15))

unipco = fit(MDS, unidm; distances=true)
plot_pcoa!(D, unipco; color=get(unirefs, :ageMonths))

E = Axis(BCDE[2,2]; title = "Neuroimaging", alignmode = Mixed(; left=-15))

brnpco = fit(MDS, brndm; distances=true)
plot_pcoa!(E, brnpco; color=get(brain, :ageMonths))

Colorbar(BCDE[1:2, 3], sc2; label="Age (months)", flipaxis=true, ticks=0:24:120)
```

### 1F - cogscores

```julia
fax1 = Axis(F[1,1]; xlabel = "Age (months)", ylabel = "cogScore", xticks=(4:4:24), limits=((2,24), nothing), alignmode = Mixed(; left=-15))
fax2 = Axis(F[1,2]; xlabel = "Age (years)", xticks = (4:2:12), limits=((2,nothing), nothing))
hideydecorations!(fax2)
linkyaxes!(fax1, fax2)

let
    u2y = findall(p-> !ismissing(p[2]) && p[1] <= 24, collect(zip(mdata.ageMonths, mdata.cogScore)))
    o2y = findall(p-> !ismissing(p[2]) && p[1] > 24, collect(zip(mdata.ageMonths, mdata.cogScore)))

    cs = ColorSchemes.colorschemes[:Set2_7]
    function colorage(age)
        age <= 36 && return cs[1] # mullen
        age <= 60 && return cs[2] # WPPSI
        return cs[3] # WISC
    end
    ages = mdata.ageMonths[u2y]
    scatter!(fax1, ages, mdata.cogScore[u2y]; color=colorage.(ages))

    ages = mdata.ageMonths[o2y]
    scatter!(fax2, ages ./ 12, mdata.cogScore[o2y]; color=colorage.(ages))
    
    vlines!(fax1, [6, 18]; linestyle=:dash, color=:gray)
    # TODO: add colors for training/test sets in later models

    Legend(F[2, 1:2], [MarkerElement(; marker=:circle, color=c) for c in cs[1:3]],
                      ["Mullen", "WPPSI", "WISC"], "Assessment";
                      orientation = :horizontal, tellheight=true, tellwidth=false
    )
end

colsize!(figure.layout, 1, Relative(1/3))
colsize!(figure.layout, 3, Relative(1/4))
colsize!(A, 1, Fixed(12))
colsize!(A, 2, Relative(1/2))
colgap!(F, Fixed(4))
figure
```

### 1G-H: Omnibus tests

#### 1G - PERMANOVA

This is a permutation test of variance

```julia
G = GridLayout(GH[1,1])
# Ba = Axis(B[1:2,1]; alignmode=Outside())
Ga = Axis(G[1,1]; title="Under 6mo", alignmode = Mixed(; left=-15))
Gb = Axis(G[1,2]; title="Over 18mo")
```

```julia
# plot_permanovas!(Ba, CSV.read(scratchfiles("permanovas_all.csv"), DataFrame))
plot_permanovas!(Ga, CSV.read(scratchfiles("permanovas_00to06.csv"), DataFrame))
hideydecorations!(Gb)
plot_permanovas!(Gb, CSV.read(scratchfiles("permanovas_18to120.csv"), DataFrame))

colsize!(GH, 1, Relative(1/2))
Ga.alignmode = Mixed(; left=0)
colsize!(G, 1, Relative(4/7))
```


#### 1H - Mantel tests

```julia
H = GridLayout(GH[1,2])

# Ca = Axis(C[1:2, 1]; alignmode=Outside())
Ha = Axis(H[1,1]; title="Under 6mo", alignmode = Mixed(; left=-15))
Hb = Axis(H[1,2]; title="Over 18mo")
hideydecorations!(Hb)

# plot_mantel!(Ca, CSV.read(scratchfiles("mantel_all.csv"), DataFrame))
plot_mantel!(Ha, CSV.read(scratchfiles("mantel_00to06.csv"), DataFrame))
plot_mantel!(Hb, CSV.read(scratchfiles("mantel_18to120.csv"), DataFrame))

# Label(GH[0,2], "Mantel"; tellwidth=false)
```

#### 1H - cogScores

### Labels & saving

```julia

Label(A[1, 1, TopLeft()], "A",
        fontsize = 26,
        font = "Open Sans Bold",
        halign = :right
)
Label(BCDE[1, 1, TopLeft()], "B",
        fontsize = 26,
        font = "Open Sans Bold",
        halign = :right
)
Label(BCDE[2, 1, TopLeft()], "C",
        fontsize = 26,
        font = "Open Sans Bold",
        halign = :right
)
Label(BCDE[1, 2, TopLeft()], "D",
        fontsize = 26,
        font = "Open Sans Bold",
        halign = :right
)
Label(BCDE[2, 2, TopLeft()], "E",
        fontsize = 26,
        font = "Open Sans Bold",
        halign = :right
)

Label(F[1, 1, TopLeft()], "F",
        fontsize = 26,
        font = "Open Sans Bold",
        halign = :right
)
Label(G[1, 1, TopLeft()], "G",
        fontsize = 26,
        font = "Open Sans Bold",
        halign = :right
)
Label(H[1, 1, TopLeft()], "H",
        fontsize = 26,
        font = "Open Sans Bold",
        halign = :right
)

save(figurefiles("Figure1.svg"), figure)
save("manuscript/assets/Figure1.png", figure)
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