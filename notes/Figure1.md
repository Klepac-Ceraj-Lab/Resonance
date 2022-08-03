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
metabolites = Resonance.load(MetabolicProfiles(); timepoint_metadata = mtdt)
```



Then, we will start constructing the figure.
See the [Makie documentation](https://makie.juliaplots.org/stable/tutorials/layout-tutorial/) for more information.


## Figure 1


```julia
figure = Figure(; resolution = (850, 1100));
A = GridLayout(figure[1,1])
BC = GridLayout(figure[2,1])

```

### 1A - Cohort diagram and ages

For now, the graphic here is a place-holder, but we'll have something like this for 1A.

```julia

A_histleft = Axis(A[1,1]; xlabel = "Age (months)", ylabel="Samples (N)") # TODO: set xticks to match cartoon
A_histright = Axis(A[1,2]; xlabel = "Age (years)", xticks=2:2:16)
A_img = Axis(A[2,1:2]; aspect = DataAspect())
hidedecorations!(A_img)
hidespines!(A_img)

image!(A_img, rotr90(load("figures/fig1a_placeholder.png")))

hist!(A_histleft, filter(<=(24), mtdt.ageMonths); color = :darkgray)
hist!(A_histright, filter(>(24), mtdt.ageMonths) ./ 12; color = :darkgray, bins=8)

rowgap!(A, -20)
rowsize!(A, 1, Relative(1/3))
colsize!(A, 2, Relative(1/3))
linkyaxes!(A_histleft, A_histright)

figure
```

### 1B - PERMANOVA

```julia
B = Axis(BC[1,1])

commlabels = ["taxa", "genes", "metabolites"]
mdlabels = ["Cog. score", "Age", "Race", "Maternal Edu."]

perms = let permout = outputfiles("permanovas_all.csv")
    if !isfile(permout)
        p = permanovas([species, unirefs, metabolites],
                        [:cogScore, :ageMonths, :race, :maternalEd];
                        commlabels, mdlabels)
        CSV.write(permout, p)
    else
        p = CSV.read(permout, DataFrame)
    end
    p
end


Resonance.plot_permanovas!(B, perms)
figure
```

### 1C - Mantel tests

```julia




```