# Figure 1 - Cohort Characteristics and Data Summaries

First, load packages that will be used throughout this notebook.

```julia
using Resonance
using FileIO
using CairoMakie # for plotting
using MultivariateStats
using CategoricalArrays
```

Then, we'll load in the different data sources.

```julia
mdata = Resonance.load(Metadata())
mdata.maternalEd = categorical(map(e-> e == "-8" ? missing : parse(Int, e), mdata.maternalEd); ordered=true)
    

taxa = Resonance.load(TaxonomicProfiles(); timepoint_metadata = mdata)
species = filter(t-> taxrank(t) == :species, taxa)

unirefs = Resonance.load(UnirefProfiles(); timepoint_metadata = mdata) # this can take a bit
unirefs = filter(!hastaxon, unirefs) # don't use species stratification for summaries

ecs = Resonance.load(ECProfiles(); timepoint_metadata = mdata)
ecs = filter(!hastaxon, ecs)

kos = Resonance.load(KOProfiles(); timepoint_metadata = mdata)
kos = filter(!hastaxon, kos)

metabolites = Resonance.load(MetabolicProfiles(); timepoint_metadata = mdata)
brain = Resonance.load(Neuroimaging())


@assert all(samplenames(species) .== samplenames(unirefs))
@assert all(samplenames(species) .== samplenames(ecs))
@assert all(samplenames(species) .== samplenames(kos))

uidx = let md = DataFrame(Microbiome.metadata(taxa))
    smp = Set(unique(md, :subject).sample)
    findall(s-> s in smp, samplenames(taxa))
end
```



Then, we will start constructing the figure.
See the [Makie documentation](https://makie.juliaplots.org/stable/tutorials/layout-tutorial/) for more information.


```julia
figure = Figure(; resolution = (1600, 1200))
A = GridLayout(figure[1,1])
BC = GridLayout(figure[2,1:2])
DEF = GridLayout(figure[1,2]);
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
A_histleft.xticks = ([0,3,6,12], ["", "", "", ""])
figure
```

### 1B-C: Omnibus tests

Both PERMANOVA and Mantel tests
require beta diversity distance metrics.
We'll calculate them once for each profile rather
rather than separately for each test.

```julia
isdefined(Main, :spedm) || (spedm = braycurtis(species))
isdefined(Main, :unidm) || (unidm = braycurtis(unirefs))
isdefined(Main, :ecsdm) || (ecsdm = braycurtis(ecs))
isdefined(Main, :kosdm) || (kosdm = braycurtis(kos))
isdefined(Main, :metdm) || (metdm = braycurtis(metabolites))
```

#### PERMANOVA

This is a permutation test of variance

```julia
B = GridLayout(BC[1,1])
Ba = Axis(B[1:2,1]; alignmode=Outside())

commlabels = ["taxa", "UniRef90s", "ECs", "KOs"]
mdlabels = ["Cog. score", "Age", "Race", "Maternal Edu."]


perms = let permout = outputfiles("permanovas_all.csv")
    if isfile(permout)
        p = CSV.read(permout, DataFrame)
    else
        p = permanovas([spedm[uidx, uidx], unidm[uidx, uidx], ecsdm[uidx, uidx], kosdm[uidx, uidx]], [
                                get(species, :cogScore)[uidx], 
                                get(species, :ageMonths)[uidx], 
                                get(species, :race)[uidx], 
                                get(species, :maternalEd)[uidx]
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

plot_permanovas!(Ba, perms)
```

```julia
Bb = Axis(B[1,2]; alignmode=Outside())
Bc = Axis(B[2,2]; alignmode=Outside())
Label(B[1,3], "Under 6mo"; tellwidth=true, tellheight=false, rotation=-π/2)
Label(B[2,3], "Over 18mo"; tellwidth=true, tellheight=false, rotation=-π/2)

perms = let permout = outputfiles("permanovas_u6mo.csv")
    if isfile(permout)
        p = CSV.read(permout, DataFrame)
    else
        idx = findall(<(6), get(species, :ageMonths))
        p = permanovas([spedm[idx, idx], unidm[idx, idx]], [
                                get(species, :cogScore)[idx], 
                                get(species, :ageMonths)[idx], 
                                get(species, :race)[idx], 
                                get(species, :maternalEd)[idx]
                    ]; commlabels=["taxa", "UniRef90s"], mdlabels
        )
        idx = findall(<(6), get(metabolites, :ageMonths))
        p2 = permanovas(metdm[idx,idx], [
                                get(metabolites, :cogScore)[idx], 
                                get(metabolites, :ageMonths)[idx], 
                                get(metabolites, :race)[idx], 
                                get(metabolites, :maternalEd)[idx]
                    ]; mdlabels
        )
        p2.label .= "metabolites"

        append!(p, p2)

        CSV.write(permout, p)
    end
    p
end

plot_permanovas!(Bb, perms)
perms = let permout = outputfiles("permanovas_o12mo.csv")
    if isfile(permout)
        p = CSV.read(permout, DataFrame)
    else
        idx = findall(>(18), get(species, :ageMonths))
        p = permanovas([spedm[idx, idx], unidm[idx, idx]], [
                                get(species, :cogScore)[idx], 
                                get(species, :ageMonths)[idx], 
                                get(species, :race)[idx], 
                                get(species, :maternalEd)[idx]
                    ]; commlabels=["taxa", "UniRef90s"], mdlabels
        )

        CSV.write(permout, p)
    end
    p
end

plot_permanovas!(Bc, perms)

colsize!(BC, 1, Relative(2/3))
Label(BC[0,1], "PERMANOVAs")
```


### 1C - Mantel tests

```julia
C = Axis(BC[1,2]; alignmode=Outside())
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
Label(BC[0,2], "Mantel"; tellwidth=false)

```

## Ordinations

```julia
D = Axis(DEF[1,1])

spepco = fit(MDS, spedm; distances=true)
sc = plot_pcoa!(D, spepco; color=get(species, :ageMonths))
Colorbar(DEF[1, 2], sc; label="Age (months)", flipaxis=true)

E = GridLayout(DEF[2,1])
Ea = Axis(E[1,1]; title = "Bacteroidetes")
Eb = Axis(E[1,2]; title = "Firmicutes")
Ec = Axis(E[1,3]; title = "Actinobacteria")

plot_pcoa!(Ea, spepco; color=vec(abundances(filter(t-> taxrank(t) == :phylum, taxa)[r"Bacteroidetes", :])), colormap=:Purples)
plot_pcoa!(Eb, spepco; color=vec(abundances(filter(t-> taxrank(t) == :phylum, taxa)[r"Firmicutes", :])), colormap=:Purples)
pco = plot_pcoa!(Ec, spepco; color=vec(abundances(filter(t-> taxrank(t) == :phylum, taxa)[r"Actinobacteria", :])), colormap=:Purples)

Colorbar(DEF[2, 2], pco; label="Relative abundance (%)", colormap = :Purples)

# F = Axis(DEF[3,1])

# metpco = fit(MDS, metdm)
# plot_pcoa!(F, metpco; color=get(metabolites, :ageMonths))


save(figurefiles("Figure1.svg"), figure)
save(figurefiles("Figure1.png"), figure)
figure
```

![](figures/Figure1.png)

## Supplement

### Multivariate permanovas

Maternal education & race.

```julia
using PERMANOVA

df = DataFrame(maternalEd = get(species, :maternalEd)[uidx], 
               race       = get(species, :race)[uidx])
oidx = findall(.!ismissing.(df.maternalEd) .& .!ismissing.(df.race))
df = df[oidx, :]

permanova(df, spedm[uidx[oidx], uidx[oidx]], @formula(1 ~ maternalEd + race), 1000)

```