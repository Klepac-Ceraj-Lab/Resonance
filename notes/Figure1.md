# Figure 1 - Cohort Characteristics and Data Summaries

First, load packages that will be used throughout this notebook.

```julia
using Resonance
using FileIO
using CairoMakie # for plotting
using Distance
```

Then, we'll load in the different data sources.

```julia
mtdt = Resonance.load(Metadata())

taxa = Resonance.load(TaxonomicProfiles(); timepoint_metadata = mtdt)
species = filter(t-> taxrank(t) == :species, taxa)

unirefs = Resonance.load(UnirefProfiles(); timepoint_metadata = mtdt) # this can take a bit
unirefs = filter(!hastaxon, unirefs) # don't use species stratification for summaries

ecs = Resonance.load(ECProfiles(); timepoint_metadata = mtdt)
ecs = filter(!hastaxon, ecs)

kos = Resonance.load(KOProfiles(); timepoint_metadata = mtdt)
kos = filter(!hastaxon, kos)

metabolites = Resonance.load(MetabolicProfiles(); timepoint_metadata = mtdt)

@assert all(samplenames(species) .== samplenames(unirefs))
@assert all(samplenames(species) .== samplenames(ecs))
@assert all(samplenames(species) .== samplenames(kos))
```



Then, we will start constructing the figure.
See the [Makie documentation](https://makie.juliaplots.org/stable/tutorials/layout-tutorial/) for more information.


```julia
figure = Figure(; resolution = (1200, 1200));
A = GridLayout(figure[1,1])
BC = GridLayout(figure[1:2,2])
DEF = GridLayout(figure[2,1])
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

### 1B-C: Omnibus tests

Both PERMANOVA and Mantel tests
require beta diversity distance metrics.
We'll calculate them once for each profile rather
rather than separately for each test.

```julia
spedm = braycurtis(species)
unidm = braycurtis(unirefs)
ecsdm = braycurtis(ecs)
kosdm = braycurtis(kos)
metdm = braycurtis(metabolites)
```

#### PERMANOVA

This is a permutation test of variance

```julia
B = Axis(BC[1,1]; alignmode=Outside())

commlabels = ["taxa", "UniRef90s", "ECs", "KOs"]
mdlabels = ["Cog. score", "Age", "Race", "Maternal Edu."]



perms = let permout = outputfiles("permanovas_all.csv")
    if isfile(permout)
        p = CSV.read(permout, DataFrame)
    else
        p = permanovas([spedm, unidm, ecsdm, kosdm], [
                                get(species, :cogScore), 
                                get(species, :ageMonths), 
                                get(species, :race), 
                                get(species, :maternalEd)
                    ]; commlabels, mdlabels
        )
        p2 = permanovas(metdm, [
                                get(metabolites, :cogScore), 
                                get(metabolites, :ageMonths), 
                                get(metabolites, :race), 
                                get(metabolites, :maternalEd)
                    ]; mdlabels
        )
        p2.label .= "metabolites"

        append!(p, p2)

        CSV.write(permout, p)
    end
    p
end

CSV.write("output/permanovas_all.csv", perms)

plot_permanovas!(B, perms)
figure
```


### 1C - Mantel tests

```julia
C = Axis(BC[2,1]; alignmode=Outside())

mdf = let mantout = outputfiles("mantel_all.csv")
    if isfile(mantout)
        mdf = CSV.read(mantout, DataFrame)
    else
        mdf = mantel([spedm, unidm, ecsdm, kosdm]; commlabels)

        (ol1, ol2) = stp_overlap(
                collect(zip(get(species, :subject), get(species, :timepoint))),
                collect(zip(get(metabolites, :subject), get(metabolites, :timepoint)))
        )
        m2 = DataFrame()
        for (i, dm1) in enumerate([spedm, unidm, ecsdm, kosdm])
            m, p = mantel(dm1[ol1, ol1], metdm[ol2, ol2])
            push!(m2, (; stat=m, pvalue=p, thing1=commlabels[i], thing2="metabolites"))
        end
        append!(mdf, m2)
        CSV.write(mantout, mdf)
    end
    mdf
end

plot_mantel!(C, mdf)
figure
```

## Ordinations

```julia
D = Axis(DEF[1,1])

spepco = fit(MDS, spedm)
plot_pcoa!(D, spepco; color=get(species, :ageMonths))

E = Axis(DEF[2,1])

unipco = fit(MDS, unidm)
plot_pcoa!(E, unipco; color=get(unirefs, :ageMonths))

F = Axis(DEF[3,1])

metpco = fit(MDS, metdm)
plot_pcoa!(F, metpco; color=get(metabolites, :ageMonths))

save(figurefiles("Figure1.svg"), figure)
save(figurefiles("Figure1.png"), figure)
figure
```

![](figures/Figure1.png)