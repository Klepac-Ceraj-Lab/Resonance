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
figure = Figure(; resolution = (1500, 900))
A = GridLayout(figure[1,1])
B = GridLayout(figure[2,1])
CDEF = GridLayout(figure[1:2,2])
GH = GridLayout(figure[3,1:2])
```


## Summaries

### 1A-B - Cohort diagram and ages, cogscores

For now, the graphic here is a place-holder, but we'll have something like this for 1A.

```julia
A_img = Axis(A[1,1]; aspect = DataAspect(), alignmode=Outside())
hidedecorations!(A_img)
hidespines!(A_img)

image!(A_img, rotr90(load("manuscript/assets/Figure1A.jpg")))

bax1 = Axis(B[1,1]; xlabel = "Age (months)", ylabel = "cogScore", xticks=(4:4:24), limits=((2,24), nothing), alignmode = Mixed(; left=-15))
bax2 = Axis(B[1,2]; xlabel = "Age (years)", xticks = (4:2:12), limits=((2,nothing), nothing))
hideydecorations!(bax2)
linkyaxes!(bax1, bax2)

let
    u2y = findall(p-> !ismissing(p[2]) && p[1] <= 24, collect(zip(mdata.ageMonths, mdata.cogScore)))
    o2y = findall(p-> !ismissing(p[2]) && p[1] > 24, collect(zip(mdata.ageMonths, mdata.cogScore)))

    cs = ColorSchemes.colorschemes[:tab20b][[1,7,15]]
    function colorage(age)
        age <= 36 && return cs[1] # mullen
        age <= 60 && return cs[2] # WPPSI
        return cs[3] # WISC
    end
    ages = mdata.ageMonths[u2y]
    scatter!(bax1, ages, mdata.cogScore[u2y]; color=colorage.(ages))

    ages = mdata.ageMonths[o2y]
    scatter!(bax2, ages ./ 12, mdata.cogScore[o2y]; color=colorage.(ages))
    
    vlines!(bax1, [6, 18]; linestyle=:dash, color=:gray)
    # TODO: add colors for training/test sets in later models

    Legend(B[1, 3], [MarkerElement(; marker=:circle, color=c) for c in cs[1:3]],
                      ["MSEL", "WPPSI", "WISC"], "Assessment";
    )
end


colgap!(B, Fixed(4))
```


### Ordinations

```julia
spepco = fit(MDS, spedm; distances=true)

C = GridLayout(CDEF[1,1])
Ca = Axis(C[1,1]; ylabel= "Age (months)", xlabel = mdsaxis(spepco, 1), 
                  title = "Taxa", yminorticksvisible = true, yticks=0:24:120, yminorticks = IntervalsBetween(2),
                  alignmode = Mixed(; left=-15))
sc1 = scatter!(Ca, Resonance.loadings(spepco, 1), get(species, :ageMonths);
        color = divr, colormap=:plasma)
hlines!(Ca, [6, 18]; linestyle=:dash, color=:darkgray)
Colorbar(C[1, 2], sc1; label="Shannon Diversity", flipaxis=true)

D = Axis(CDEF[2,1]; title = "Taxa", alignmode = Mixed(; left=-15))
sc2 = plot_pcoa!(D, spepco; color=get(species, :ageMonths))

E = Axis(CDEF[1,2]; title = "Functions", alignmode = Mixed(; left=-15))

unipco = fit(MDS, unidm; distances=true)
plot_pcoa!(E, unipco; color=get(unirefs, :ageMonths))

F = Axis(CDEF[2,2]; title = "Neuroimaging", alignmode = Mixed(; left=-15))

brnpco = fit(MDS, brndm; distances=true)
plot_pcoa!(F, brnpco; color=get(brain, :ageMonths))

Colorbar(CDEF[1:2, 3], sc2; label="Age (months)", flipaxis=true, ticks=0:24:120)
```

### 1G-H: Omnibus tests

#### 1G - PERMANOVA

This is a permutation test of variance

```julia
G = GridLayout(GH[1,1])
# Ba = Axis(B[1:2,1]; alignmode=Outside())
Ga = Axis(G[1,1]; title="Under 6mo", alignmode = Mixed(; left=-15))
Gb = Axis(G[1,2]; title="Over 18mo")

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
        fontsize = 20,
        font = "Open Sans Bold",
        halign = :right
)
Label(B[1, 1, TopLeft()], "B",
        fontsize = 20,
        font = "Open Sans Bold",
        halign = :right
)
Label(CDEF[1, 1, TopLeft()], "C",
        fontsize = 20,
        font = "Open Sans Bold",
        halign = :right
)
Label(CDEF[2, 1, TopLeft()], "D",
        fontsize = 20,
        font = "Open Sans Bold",
        halign = :right
)
Label(CDEF[1, 2, TopLeft()], "E",
        fontsize = 20,
        font = "Open Sans Bold",
        halign = :right
)
Label(CDEF[2, 2, TopLeft()], "F",
        fontsize = 20,
        font = "Open Sans Bold",
        halign = :right
)

Label(G[1, 1, TopLeft()], "G",
        fontsize = 20,
        font = "Open Sans Bold",
        halign = :right
)
Label(H[1, 1, TopLeft()], "H",
        fontsize = 20,
        font = "Open Sans Bold",
        halign = :right
)

 
colsize!(figure.layout, 1, Relative(1/3))
figure
```

```julia
aw
```

![](figures/Figure1.png)

## Supplement

### Sample collection histograms

```julia

let fig = Figure()
        histleft = Axis(fig[1,1];  xlabel = "Age (months)", ylabel="Samples (N)", xticks=(3:3:24))
        histright = Axis(fig[1,2]; xlabel = "Age (years)", xticks=2:2:10, yticklabelsvisible=false, yticksvisible=false)
        hist!(histleft, filter(<=(24), mdata.ageMonths); color = :darkgray)
        hist!(histright, filter(>(24), mdata.ageMonths) ./ 12; color = :darkgray, bins=8)

        linkyaxes!(histleft, histright)
        
        tightlimits!(histleft, Bottom())
        tightlimits!(histright, Bottom())
        save(figurefiles("Supp_Figure1.svg"), fig)
        save("manuscript/assets/Supp_Figure1.png", fig)
        fig
end
```

### More Taxa PCoA

```julia
fulltax = Resonance.load_raw_metaphlan()[:, samplenames(taxa)]

let
    fig = Figure()
    ax = Axis(fig[1,1], title = "Bacteroidetes")
    ax2 = Axis(fig[1,2], title = "Prevotella")
    ax3 = Axis(fig[2,1], title = "Firmicutes")
    ax4 = Axis(fig[2,2], title = "Bifidobacterium")


    plot_pcoa!(ax, spepco; color=vec(abundances(filter(t-> taxrank(t) == :phylum, fulltax)[r"Bacteroidetes$", :])),
            colormap=:Purples,
            strokecolor=:black,
            strokewidth=1,
            colorrange=(0,100)
    )
    plot_pcoa!(ax2, spepco; color=vec(abundances(filter(t-> taxrank(t) == :genus, fulltax)[r"Prevotella$", :])),
            colormap=:Purples,
            strokecolor=:black,
            strokewidth=1,
            colorrange=(0,100)
    )
    plot_pcoa!(ax3, spepco; color=vec(abundances(filter(t-> taxrank(t) == :phylum, fulltax)[r"Firmicutes$", :])),
            colormap=:Purples,
            strokecolor=:black,
            strokewidth=1,
            colorrange=(0,100)
    )
    sc = plot_pcoa!(ax4, spepco; color=vec(abundances(filter(t-> taxrank(t) == :genus, fulltax)[r"Bifidobacterium$", :])),
            colormap=:Purples,
            strokecolor=:black,
            strokewidth=1,
            colorrange=(0,100)
    )
    Colorbar(fig[1:2, 3], sc; label="Relative abundance (%)")
    save(figurefiles("Supp_Figure2.svg"), fig)
    save("manuscript/assets/Supp_Figure2.png", fig)
    fig
end
```

### Cog function / COVID


```julia
let
    fig = Figure()
    ax = Axis(fig[1,1]; ylabel = "cogScore", xlabel = "age group (months)", xticks=(1:5, ["0-6", "6-12", "12-18", "18-24", "> 24"]))
    df = DataFrame(age = get(unirefs, :ageMonths), score = get(unirefs, :cogScore), date = get(unirefs, :date))
    subset!(df, AsTable(["age", "score", "date"]) => ByRow(row-> all(!ismissing, values(row))))
    df.grp = categorical(map(df.age) do a
        a < 6 && return "0-6"
        a < 12 && return "6-12"
        a < 18 && return "12-18"
        a < 24 && return "18-24"
        return "> 24"
    end; ordered=true, levels = ["0-6", "6-12", "12-18", "18-24", "> 24"])

    grp = groupby(df, "grp")
    transform!(grp, "grp"=> ByRow(levelcode) => "x", "score" => (x-> x .< mean(x)) => "low")
    scatter!(ax, df.x .+ rand(Normal(0, 0.05), size(df, 1)) .+ [x < Date("2020-03-01") ? -0.15 : 0.15 for x in df.date], df.score; 
            color = [x < Date("2020-03-01") ? (:dodgerblue, 0.3) : (:orangered, 0.3) for x in df.date])
    
    
    Legend(fig[2,1], [MarkerElement(; color = :dodgerblue, marker=:circle), MarkerElement(; color = :orangered, marker=:circle)],
                ["Pre-covid", "Post-covid"]; orientation=:horizontal, tellheight=true, tellwidth=false, framevisible=false)
    save(figurefiles("Supp_Figure3.svg"), fig)
    save("manuscript/assets/Supp_Figure3.png", fig)
    fig 
end

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

maximum([
        vec(abundances(filter(t-> taxrank(t) == :phylum, fulltax)[r"Bacteroidetes$", :]));
        vec(abundances(filter(t-> taxrank(t) == :genus, fulltax)[r"Prevotella$", :]));
        vec(abundances(filter(t-> taxrank(t) == :phylum, fulltax)[r"Firmicutes$", :]));
        vec(abundances(filter(t-> taxrank(t) == :genus, fulltax)[r"Bifidobacterium$", :]))])

```

