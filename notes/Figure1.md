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
unirefs = Resonance.load(UnirefProfiles(); timepoint_metadata = mtdt) # this can take a bit
metabolites = Resonance.load(MetabolicProfiles(); timepoint_metadata = mtdt)
```



Then, we will start constructing the figure.
See the [Makie documentation](https://makie.juliaplots.org/stable/tutorials/layout-tutorial/) for more information.

```julia
figure1 = Figure(; resolution = (850, 1100));
```

For now, this graphic is a place-holder, but we'll have something like this for 1A.

```julia
fig1a = GridLayout(figure1[1,1])

fig1a_hist = Axis(fig1a[1,1])
fig1a_img = Axis(fig1a[2,1]; aspect = DataAspect())

image!(fig1a_img, rotr90(load("figures/fig1a_placeholder.png")))

figure1

```




