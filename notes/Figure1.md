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

taxa = Resonance.load(TaxonomicProfiles(); timepoint_metadata = mdata)
species = filter(t-> taxrank(t) == :species, taxa)

unirefs = Resonance.load(UnirefProfiles(); timepoint_metadata = mdata) # this can take a bit
unirefs = filter(!hastaxon, unirefs) # don't use species stratification for summaries

ecs = Resonance.load(ECProfiles(); timepoint_metadata = mdata)
ecs = filter(!hastaxon, ecs)

kos = Resonance.load(KOProfiles(); timepoint_metadata = mdata)
kos = filter(!hastaxon, kos)

metabolites = Resonance.load(MetabolicProfiles(); timepoint_metadata = mdata)
brain = Resonance.load(Neuroimg(), timepoint_metadata = mdata)


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



Then, we will start constructing the figure.
See the [Makie documentation](https://makie.juliaplots.org/stable/tutorials/layout-tutorial/) for more information.


```julia
figure = Figure(; resolution = (2000, 1200))
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
isdefined(Main, :spedm) || (spedm = Microbiome.braycurtis(species))
isdefined(Main, :unidm) || (unidm = Microbiome.braycurtis(unirefs))
isdefined(Main, :ecsdm) || (ecsdm = Microbiome.braycurtis(ecs))
isdefined(Main, :kosdm) || (kosdm = Microbiome.braycurtis(kos))
isdefined(Main, :metdm) || (metdm = Microbiome.braycurtis(metabolites))
isdefined(Main, :brndm) || (brndm = Microbiome.braycurtis(brain))
```

#### PERMANOVA

This is a permutation test of variance

```julia
B = GridLayout(BC[1,1])
Ba = Axis(B[1:2,1]; alignmode=Outside())

commlabels = ["taxa", "UniRef90s", "ECs", "KOs"]
mdlabels = ["Cog. score", "Age", "Sex", "Maternal Edu."]


perms = let permout = outputfiles("permanovas_all.csv")
    if isfile(permout)
        p = CSV.read(permout, DataFrame)
    else
        p = permanovas([spedm[uidx, uidx], unidm[uidx, uidx], ecsdm[uidx, uidx], kosdm[uidx, uidx]], [
                                get(species, :cogScore)[uidx], 
                                get(species, :ageMonths)[uidx], 
                                get(species, :sex)[uidx], 
                                get(species, :education)[uidx]
                    ]; commlabels, mdlabels
        )
        p2 = permanovas(metdm, [
                                get(metabolites, :cogScore), 
                                get(metabolites, :ageMonths), 
                                get(metabolites, :sex), 
                                get(metabolites, :education)
                    ]; mdlabels
        )
        p2.label .= "metabolites"
        append!(p, p2)

        p3 = permanovas(brndm[buidx, buidx], [
                                get(brain, :cogScore)[buidx], 
                                get(brain, :ageMonths)[buidx], 
                                get(brain, :sex)[buidx], 
                                get(brain, :education)[buidx]
                    ]; mdlabels
        )
        p3.label .= "neuroimg"
        append!(p, p3)

        CSV.write(permout, p)
    end
    p
end

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
                                get(species, :sex)[idx], 
                                get(species, :education)[idx]
                    ]; commlabels=["taxa", "UniRef90s"], mdlabels
        )
        idx = findall(<(6), get(metabolites, :ageMonths))
        p2 = permanovas(metdm[idx,idx], [
                                get(metabolites, :cogScore)[idx], 
                                get(metabolites, :ageMonths)[idx], 
                                get(metabolites, :sex)[idx], 
                                get(metabolites, :education)[idx]
                    ]; mdlabels
        )
        p2.label .= "metabolites"

        append!(p, p2)

        idx = findall(<(6), get(brain, :ageMonths))
        p3 = permanovas(brndm[idx, idx], [
                                get(brain, :cogScore)[idx], 
                                get(brain, :ageMonths)[idx], 
                                get(brain, :sex)[idx], 
                                get(brain, :education)[idx]
                    ]; mdlabels
        )
        p3.label .= "neuroimg"
        append!(p, p3)

        CSV.write(permout, p)
    end
    p
end

plot_permanovas!(Bb, perms)

perms = let permout = outputfiles("permanovas_o18mo.csv")
    if isfile(permout)
        p = CSV.read(permout, DataFrame)
    else
        idx = findall(>(18), get(species, :ageMonths))
        p = permanovas([spedm[idx, idx], unidm[idx, idx]], [
                                get(species, :cogScore)[idx], 
                                get(species, :ageMonths)[idx], 
                                get(species, :sex)[idx], 
                                get(species, :education)[idx]
                    ]; commlabels=["taxa", "UniRef90s"], mdlabels
        )
        bidx = findall(>(18), get(brain, :ageMonths))
        p3 = permanovas(brndm[bidx, bidx], [
                                get(brain, :cogScore)[bidx], 
                                get(brain, :ageMonths)[bidx], 
                                get(brain, :sex)[bidx], 
                                get(brain, :education)[bidx]
                    ]; mdlabels
        )
        p3.label .= "neuroimg"
        append!(p, p3)
        CSV.write(permout, p)
    end
    p
end

plot_permanovas!(Bc, perms)

colsize!(BC, 1, Relative(1/2))
Label(BC[0,1], "PERMANOVAs")
```


### 1C - Mantel tests

```julia
C = GridLayout(BC[1,2])

Ca = Axis(C[1:2, 1]; alignmode=Outside())
Cb = Axis(C[1,2]; alignmode=Outside())
Cc = Axis(C[2,2]; alignmode=Outside())

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

        (ol3, ol4) = stp_overlap(
                collect(zip(get(species, :subject), get(species, :timepoint))),
                collect(zip(get(brain, :subject), get(brain, :timepoint)))
        )
        m3 = DataFrame()
        for (i, dm1) in enumerate([spedm, unidm, ecsdm, kosdm])
            m, p = mantel(dm1[ol3, ol3], brndm[ol4, ol4])
            push!(m3, (; stat=m, pvalue=p, thing1=commlabels[i], thing2="neuroimg"))
        end
        append!(mdf, m3)
        
        (ol5, ol6) = stp_overlap(
                collect(zip(get(metabolites, :subject), get(metabolites, :timepoint))),
                collect(zip(get(brain, :subject), get(brain, :timepoint))),
        )

        m, p = mantel(metdm[ol5, ol5], brndm[ol6, ol6])
        push!(mdf, (; stat=m, pvalue=p, thing1="metabolites", thing2="neuroimg"))        
        
        CSV.write(mantout, mdf)
    end
    mdf
end

mdfu6 = let mantout = outputfiles("mantel_u6.csv")
    if isfile(mantout)
        mdf = CSV.read(mantout, DataFrame)
    else
        speidx = get(species, :ageMonths) .< 6
        metidx = get(metabolites, :ageMonths) .< 6
        brnidx = get(brain, :ageMonths) .< 6


        speu6dm = spedm[speidx, speidx]
        uniu6dm = unidm[speidx, speidx]
        ecsu6dm = ecsdm[speidx, speidx]
        kosu6dm = kosdm[speidx, speidx]
        metu6dm = metdm[metidx, metidx]
        brnu6dm = brndm[brnidx, brnidx]
        

        mdf = mantel([speu6dm, uniu6dm, ecsu6dm, kosu6dm]; commlabels)

        (ol1, ol2) = stp_overlap(
                collect(zip(get(species, :subject)[speidx],
                            get(species, :timepoint)[speidx])
                        ),
                collect(zip(get(metabolites, :subject)[metidx],
                            get(metabolites, :timepoint)[metidx])
                        )
        )



        m2 = DataFrame()
        for (i, dm1) in enumerate([speu6dm, uniu6dm, ecsu6dm, kosu6dm])
            m, p = mantel(dm1[ol1, ol1], metu6dm[ol2, ol2])
            push!(m2, (; stat=m, pvalue=p, thing1=commlabels[i], thing2="metabolites"))
        end
        append!(mdf, m2)


        (ol3, ol4) = stp_overlap(
                collect(zip(get(species, :subject)[speidx],
                            get(species, :timepoint)[speidx])
                        ),
                collect(zip(get(brain, :subject)[brnidx],
                            get(brain, :timepoint)[brnidx])
                        )
        )
        m3 = DataFrame()
        for (i, dm1) in enumerate([speu6dm, uniu6dm, ecsu6dm, kosu6dm])
            m, p = mantel(dm1[ol3, ol3], brnu6dm[ol4, ol4])
            push!(m3, (; stat=m, pvalue=p, thing1=commlabels[i], thing2="neuroimg"))
        end
        append!(mdf, m3)
        
        (ol5, ol6) = stp_overlap(
                collect(zip(get(metabolites, :subject)[metidx],
                            get(metabolites, :timepoint)[metidx])
                        ),
                collect(zip(get(brain, :subject)[brnidx],
                            get(brain, :timepoint)[brnidx])
                        ),
        )

        m, p = mantel(metu6dm[ol5, ol5], brnu6dm[ol6, ol6])
        push!(mdf, (; stat=m, pvalue=p, thing1="metabolites", thing2="neuroimg"))        
        
        CSV.write(mantout, mdf)
    end
    mdf
end

mdfo18 = let mantout = outputfiles("mantel_o18.csv")
    if isfile(mantout)
        mdf = CSV.read(mantout, DataFrame)
    else
        speidx = get(species, :ageMonths) .> 18
        brnidx = get(brain, :ageMonths) .> 18


        speo18dm = spedm[speidx, speidx]
        unio18dm = unidm[speidx, speidx]
        ecso18dm = ecsdm[speidx, speidx]
        koso18dm = kosdm[speidx, speidx]
        brno18dm = brndm[brnidx, brnidx]
        

        mdf = mantel([speo18dm, unio18dm, ecso18dm, koso18dm]; commlabels)


        (ol3, ol4) = stp_overlap(
                collect(zip(get(species, :subject)[speidx],
                            get(species, :timepoint)[speidx])
                        ),
                collect(zip(get(brain, :subject)[brnidx],
                            get(brain, :timepoint)[brnidx])
                        )
        )
        m3 = DataFrame()
        for (i, dm1) in enumerate([speo18dm, unio18dm, ecso18dm, koso18dm])
            m, p = mantel(dm1[ol3, ol3], brno18dm[ol4, ol4])
            push!(m3, (; stat=m, pvalue=p, thing1=commlabels[i], thing2="neuroimg"))
        end
        append!(mdf, m3)   
        
        CSV.write(mantout, mdf)
    end
    mdf
end

plot_mantel!(Ca, mdf)
plot_mantel!(Cb, mdfu6)
plot_mantel!(Cc, mdfo18)

Label(BC[0,2], "Mantel"; tellwidth=false)
Label(C[1,3], "Under 6mo"; tellwidth=true, tellheight=false, rotation=-π/2)
Label(C[2,3], "Over 18mo"; tellwidth=true, tellheight=false, rotation=-π/2)

```

## Ordinations

```julia
D = Axis(DEF[1,1]; title = "Taxa")

spepco = fit(MDS, spedm; distances=true)

sc = plot_pcoa!(D, spepco; color=get(species, :ageMonths))
Colorbar(DEF[1, 2], sc; label="Age (months)", flipaxis=true)

E = GridLayout(DEF[2,1])
Ea = Axis(E[1,1]; title = "Bacteroidetes")
Eb = Axis(E[1,2]; title = "Firmicutes")
# Ec = Axis(E[1,3]; title = "Actinobacteria")


plot_pcoa!(Ea, spepco; color=vec(abundances(filter(t-> taxrank(t) == :phylum, taxa)[r"Bacteroidetes", :])), colormap=:Purples)
plot_pcoa!(Eb, spepco; color=vec(abundances(filter(t-> taxrank(t) == :phylum, taxa)[r"Firmicutes", :])), colormap=:Purples)
# pco = plot_pcoa!(Ec, spepco; color=vec(abundances(filter(t-> taxrank(t) == :phylum, taxa)[r"Actinobacteria", :])), colormap=:Purples)

Colorbar(DEF[2, 2], pco; label="Relative abundance (%)")

F = Axis(DEF[1,3]; title = "Functions")

unipco = fit(MDS, unidm; distances=true)
plot_pcoa!(F, unipco; color=get(unirefs, :ageMonths))

G = Axis(DEF[2,3]; title = "Neuroimaging")

brnpco = fit(MDS, brndm; distances=true)
plot_pcoa!(G, brnpco; color=get(brain, :ageMonths))

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