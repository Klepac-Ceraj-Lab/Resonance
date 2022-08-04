# Figure 1 - Cohort Characteristics and Data Summaries

First, load packages that will be used throughout this notebook.

```julia
using Resonance
using FileIO
using CairoMakie # for plotting
```

Then, we'll load in the different data sources.

```julia
mtdt = Resonance.load(Metadata())

taxa = Resonance.load(TaxonomicProfiles(); timepoint_metadata = mtdt)
species = filter(t-> taxrank(t) == :species, taxa)

unirefs = Resonance.load(UnirefProfiles(); timepoint_metadata = mtdt) # this can take a bit
unirefs = filter(!hastaxon, unirefs)

ecs = Resonance.load(ECProfiles(); timepoint_metadata = mtdt)
ecs = filter(!hastaxon, ecs)

kos = Resonance.load(KOProfiles(); timepoint_metadata = mtdt)
kos = filter(!hastaxon, kos)

metabolites = Resonance.load(MetabolicProfiles(); timepoint_metadata = mtdt)
```



Then, we will start constructing the figure.
See the [Makie documentation](https://makie.juliaplots.org/stable/tutorials/layout-tutorial/) for more information.


```julia
figure = Figure(; resolution = (850, 1100));
A = GridLayout(figure[1,1])
BC = GridLayout(figure[2,1])
DEF = GridLayout(figure[1:2,2])
```


## Summaries



### 1A - Cohort diagram and ages

For now, the graphic here is a place-holder, but we'll have something like this for 1A.

```julia

A_histleft = Axis(A[1,1]; xlabel = "Age (months)", ylabel="Samples (N)", alignmode=Outside()) # TODO: set xticks to match cartoon
A_histright = Axis(A[1,2]; xlabel = "Age (years)", xticks=2:2:16, alignmode=Outside())
A_img = Axis(A[2,1:2]; aspect = DataAspect())
hidedecorations!(A_img)
hidespines!(A_img)

image!(A_img, rotr90(load("figures/fig1a_placeholder.png")))

hist!(A_histleft, filter(<=(24), mtdt.ageMonths); color = :darkgray)
hist!(A_histright, filter(>(24), mtdt.ageMonths) ./ 12; color = :darkgray, bins=8)

rowgap!(A, -20)
rowsize!(A, 1, Relative(1/4))
colsize!(A, 2, Relative(1/3))
linkyaxes!(A_histleft, A_histright)

figure
```

### 1B - PERMANOVA

```julia
B = Axis(BC[1,1]; alignmode=Outside())

commlabels = ["taxa", #="UniRef90s",=# "ECs", "KOs", "metabolites"]
mdlabels = ["Cog. score", "Age", "Race", "Maternal Edu."]

perms = let permout = outputfiles("permanovas_all.csv")
    if isfile(permout)
        p = CSV.read(permout, DataFrame)
    else
        p = permanovas([species, #=unirefs,=# ecs, kos, metabolites],
                    [:cogScore, :ageMonths, :race, :maternalEd];
                    commlabels, mdlabels
        )
        CSV.write(permout, p)
    end
    p
end

CSV.write("output/permanovas_all.csv", p)

plot_permanovas!(B, perms)
figure
```


### 1C - Mantel tests

```julia
C = Axis(BC[1,2]; alignmode=Outside())

mdf = let mantout = outputfiles("mantel_all.csv")
    if isfile(mantout)
        m = CSV.read(mantout, DataFrame)
    else
        m = mantel([species, #=unirefs,=# ecs, kos, metabolites]; commlabels)
        CSV.write(mantout, m)
    end
    m
end

plot_mantel!(C, mdf)
figure
```

## Ordinations

```julia
D = Axis(DEF[1,1])

plot_pcoa!(D, pcoa(species); color=get(species, :ageMonths))

E = Axis(DEF[2,1])

plot_pcoa!(E, pcoa(unirefs); color=get(species, :ageMonths))

F = Axis(DEF[3,1])

plot_pcoa!(F, pcoa(metabolites); color=get(species, :ageMonths))

save(figurefiles("Figure1.svg"), figure)
save(figurefiles("Figure1.png"), figure)
figure
```

