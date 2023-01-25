# Figure2 - functional analysis

```julia
using Resonance
using CairoMakie
using Statistics
using HypothesisTests
using MultipleTesting
using KernelDensity
using MultivariateStats
using Distributions
using CategoricalArrays
using ThreadsX
using ColorSchemes
using GLM
using Dates
using FileIO
using JLD2
```

## Data Loading

```julia
mdata = Resonance.load(Metadata())
taxa = Resonance.load(TaxonomicProfiles(); timepoint_metadata = mdata)
taxdf = comm2wide(taxa)
```


## Load FSEA


```julia
speclms_00to120 = subset(CSV.read(tablefiles("lms_species_00to120.csv"), DataFrame), "kind"=> ByRow(==("cogScore")))
speclms_00to06 = subset(CSV.read(tablefiles("lms_species_00to06.csv"), DataFrame), "kind"=> ByRow(==("cogScore")))
speclms_18to120 = subset(CSV.read(tablefiles("lms_species_18to120.csv"), DataFrame), "kind"=> ByRow(==("cogScore")))
speclms = subset(CSV.read(tablefiles("lms_species_18to120.csv"), DataFrame), "kind"=> ByRow(==("cogScore")))
speclms_pa = CSV.read(tablefiles("lms_species_18to120_pa.csv"), DataFrame)
fsdf_00to120 = CSV.read(scratchfiles("figure2", "fsea_consolidated_00to120.csv"), DataFrame)
fsdf2_00to120 = CSV.read(scratchfiles("figure2", "fsea_all_00to120.csv"), DataFrame)
fsdf_00to06 = CSV.read(scratchfiles("figure2", "fsea_consolidated_00to06.csv"), DataFrame)
fsdf2_00to06 = CSV.read(scratchfiles("figure2", "fsea_all_00to06.csv"), DataFrame)
fsdf_18to120 = CSV.read(scratchfiles("figure2", "fsea_consolidated_18to120.csv"), DataFrame)
fsdf2_18to120 = CSV.read(scratchfiles("figure2", "fsea_all_18to120.csv"), DataFrame)

cors_00to120 = CSV.read(tablefiles("lms_unirefs_00to120.csv"), DataFrame)
cors_00to06 = CSV.read(tablefiles("lms_unirefs_00to06.csv"), DataFrame)
cors_18to120 = CSV.read(tablefiles("lms_unirefs_18to120.csv"), DataFrame)
neuroactive_00to120 = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), cors_00to120.species); consolidate=false)
neuroactive_00to06 = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), cors_00to06.species); consolidate=false)
neuroactive_18to120 = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), cors_18to120.species); consolidate=false)
```


```julia
sigs = subset(speclms, "qvalue"=> ByRow(<(0.2)))
for i in eachindex(sigs.species)
    sp = sigs.species[i]
    prev_00to120 = prevalence(taxa[Regex(sp), mdata.filter_00to120])
    prev_00to06 = prevalence(taxa[Regex(sp), mdata.filter_00to06])
    prev_18to120 = prevalence(taxa[Regex(sp), mdata.filter_18to120])
    meanab_00to120 = mean(abundances(taxa[Regex(sp), mdata.filter_00to120]))
    meanab_00to06 = mean(abundances(taxa[Regex(sp), mdata.filter_00to06]))
    meanab_18to120 = mean(abundances(taxa[Regex(sp), mdata.filter_18to120]))

    
end
```

## Plotting

```julia
figure = Figure(resolution=(1200, 900))

AB = GridLayout(figure[1,1]; alignmode=Outside())
CDEF = GridLayout(figure[1:2,2]; alignmode=Outside())
C = GridLayout(CDEF[1,1])
D = GridLayout(CDEF[2,1])
E = GridLayout(CDEF[3,1])
F = GridLayout(CDEF[4,1])
G = GridLayout(figure[3,1:2]; alignmode=Outside())
```


```julia
aax1 = Axis(AB[1,1]; title = "under 6mo", ylabel=L"$-log_2(P)$", xlabel = "Coef.")
aax2 = Axis(AB[1,2]; title = "over 18mo", ylabel=L"$-log_2(P)$", xlabel = "Coef.")
scatter!(aax1, speclms_00to06."Coef.", -1 .* log2.(speclms_00to06.qvalue), color = map(q-> q < 0.2 ? :orange : :dodgerblue, speclms_00to06.qvalue))
scatter!(aax2, speclms_18to120."Coef.", -1 .* log2.(speclms_18to120.qvalue), color = map(q-> q < 0.2 ? :orange : :dodgerblue, speclms_18to120.qvalue))

bax1 = Axis(AB[1,2]; )



```


```julia

let
    gs = "Acetate synthesis II"
    panel = C
    ixs = neuroactive_00to06[gs]
    cs = filter(!isnan, cors_00to06.t[ixs])
    acs = filter(!isnan, cors_00to06.t[Not(ixs)])

    Resonance.plot_fsea!(panel, cs, acs;
        label = replace(gs, "degradation"=> "degr.", "synthesis"=> "synth.", " (vitamin K2)"=> ""),
        ylabel = "")
    Label(panel[1,2], "Under 6m"; tellheight=false, tellwidth=true)


    gs = "Glutamate synthesis II"
    panel = D
    ixs = neuroactive_18to120[gs]
    cs = filter(!isnan, cors_18to120.t[ixs])
    acs = filter(!isnan, cors_18to120.t[Not(ixs)])

    Resonance.plot_fsea!(panel, cs, acs;
        label = replace(gs, "degradation"=> "degr.", "synthesis"=> "synth.", " (vitamin K2)"=> ""),
        ylabel = "", xticks = -0.2:0.1:0.0)
    Label(panel[1,2], "Over 18m"; tellheight=false, tellwidth=true)

    gs = "Propionate degradation I"
    panel = E
    ixs = neuroactive_00to06[gs]
    cs = filter(!isnan, cors_00to06.t[ixs])
    acs = filter(!isnan, cors_00to06.t[Not(ixs)])

    Resonance.plot_fsea!(panel, cs, acs;
        label = replace(gs, "degradation"=> "degr.", "synthesis"=> "synth.", " (vitamin K2)"=> ""),
        ylabel = "", xticks = -0.7:0.2:0.0)
    Label(panel[1,2], "Under 6m"; tellheight=false, tellwidth=true)

    gs = "Propionate degradation I"
    panel = F
    ixs = neuroactive_18to120[gs]
    cs = filter(!isnan, cors_18to120.t[ixs])
    acs = filter(!isnan, cors_18to120.t[Not(ixs)])

    Resonance.plot_fsea!(panel, cs, acs;
        label = replace(gs, "degradation"=> "degr.", "synthesis"=> "synth.", " (vitamin K2)"=> ""),
        ylabel = "")
    Label(panel[1,2], "Over 18m"; tellheight=false, tellwidth=true)
end

```


```julia
let
    genesets = union(subset(fsdf2_00to120, "qvalue"=> ByRow(<(0.2)), "cortest"=> ByRow(==("cogScore"))).geneset,
                     subset(fsdf2_00to06, "qvalue"=> ByRow(<(0.2)), "cortest"=> ByRow(==("cogScore"))).geneset,
                     subset(fsdf2_18to120, "qvalue"=> ByRow(<(0.2)), "cortest"=> ByRow(==("cogScore"))).geneset
    )
  
    df = sort(subset(fsdf2_00to120, "geneset"=> ByRow(gs-> gs in genesets), "cortest"=>ByRow(==("cogScore"))), :geneset; rev=true)
    ax = Axis(G[1,1]; yticks = (1:nrow(df), replace.(df.geneset, r" \(.+\)" => "", "synthesis"=>"syn.", "degradation"=>"deg.")), 
                xlabel="T stat", title="All ages")
    m = median(filter(x-> !isnan(x) && x < 7, cors_00to120.t))
    colors = ColorSchemes.colorschemes[:RdBu_7]

    for (i, row) in enumerate(eachrow(df))
        sign = row.enrichment < 0 ? "neg" : "pos"
        c = row.qvalue > 0.2 ? :gray : 
            row.qvalue > 0.05 ? (sign == "neg" ? colors[3] : colors[5]) :
            row.qvalue > 0.01 ? (sign == "neg" ? colors[2] : colors[6]) :
            sign == "neg" ? colors[1] : colors[7]

        y = filter(x-> !isnan(x) && x < 7, cors_00to120.t[neuroactive_00to120[row.geneset]])
        scatter!(ax, y, rand(Normal(0, 0.1), length(y)) .+ i; color=(c,0.3), strokecolor=:gray, strokewidth=0.5)
        row.qvalue < 0.2 && lines!(ax, fill(median(y), 2), [i-0.4, i+0.4]; color = c, linewidth=2)
    end
    vlines!(ax, m; linestyle=:dash, color=:darkgray)

    ####

    df = sort(subset(fsdf2_00to06, "geneset"=> ByRow(gs-> gs in genesets), "cortest"=>ByRow(==("cogScore"))), :geneset; rev=true)
    ax = Axis(G[1,2]; yticks = (1:nrow(df), replace.(df.geneset, r" \(.+\)" => "", "synthesis"=>"syn.", "degradation"=>"deg.")), 
                xlabel="T stat", title="Under 6mo")
    hideydecorations!(ax, grid=false)

    m = median(filter(x-> !isnan(x) && x < 7, cors_00to06.t))
    colors = ColorSchemes.colorschemes[:RdBu_7]

    for (i, row) in enumerate(eachrow(df))
        sign = row.enrichment < 0 ? "neg" : "pos"
        c = row.qvalue > 0.2 ? :gray : 
            row.qvalue > 0.05 ? (sign == "neg" ? colors[3] : colors[5]) :
            row.qvalue > 0.01 ? (sign == "neg" ? colors[2] : colors[6]) :
            sign == "neg" ? colors[1] : colors[7]

        y = filter(x-> !isnan(x) && x < 7, cors_00to06.t[neuroactive_00to06[row.geneset]])
        scatter!(ax, y, rand(Normal(0, 0.1), length(y)) .+ i; color=(c,0.3), strokecolor=:gray, strokewidth=0.5)
        row.qvalue < 0.2 && lines!(ax, fill(median(y), 2), [i-0.4, i+0.4]; color = c, linewidth=2)
    end
    vlines!(ax, m; linestyle=:dash, color=:darkgray)

    ####

    df = sort(subset(fsdf2_18to120, "geneset"=> ByRow(gs-> gs in genesets), "cortest"=>ByRow(==("cogScore"))), :geneset; rev=true)
    ax = Axis(G[1,3]; yticks = (1:nrow(df), replace.(df.geneset, r" \(.+?\)" => "", "synthesis"=>"syn.", "degradation"=>"deg.")), 
                xlabel="T stat", title="over 18")
    hideydecorations!(ax, grid=false)
    m = median(filter(x-> !isnan(x) && x < 7, cors_18to120.t))
    colors = ColorSchemes.colorschemes[:RdBu_7]

    for (i, row) in enumerate(eachrow(df))
        sign = row.enrichment < 0 ? "neg" : "pos"
        c = row.qvalue > 0.2 ? :gray : 
            row.qvalue > 0.05 ? (sign == "neg" ? colors[3] : colors[5]) :
            row.qvalue > 0.01 ? (sign == "neg" ? colors[2] : colors[6]) :
            sign == "neg" ? colors[1] : colors[7]

        y = filter(x-> !isnan(x) && x < 7, cors_18to120.t[neuroactive_18to120[row.geneset]])
        scatter!(ax, y, rand(Normal(0, 0.1), length(y)) .+ i; color=(c,0.3), strokecolor=:gray, strokewidth=0.5)
        row.qvalue < 0.2 && lines!(ax, fill(median(y), 2), [i-0.4, i+0.4]; color = c, linewidth=2)
    end
    vlines!(ax, m; linestyle=:dash, color=:darkgray)

    Legend(G[1,4], [MarkerElement(; color = (c, 0.3),
                                    marker=:circle,
                                    strokecolor=:gray,
                                    strokewidth=0.5) for c in colors[[1:3..., 5:7...]]],
                   ["(-) q < 0.01", "(-) q < 0.05", "(-) q < 0.2", 
                    "(+) q < 0.20", "(+) q < 0.05", "(+) p < 0.01"])
end

rowsize!(figure.layout, 3, Relative(2/5))
colgap!(G, Fixed(4))
```

```julia
for (label, layout) in zip(["A", "B"], [A, B])
    Label(layout[1, 1, TopLeft()], label,
        textsize = 26,
        font = "Open Sans Bold",
        padding = (0, 50, 5, 0),
        halign = :right)
end
for (label, layout) in zip(["C", "D", "E", "F"], [C, D, E, F])
    Label(layout[1, 1, TopLeft()], label,
        textsize = 26,
        font = "Open Sans Bold",
        padding = (0, 40, 5, 0),
        halign = :right)
end

Label(G[1, 1, TopLeft()], "G",
        textsize = 26,
        font = "Open Sans Bold",
        padding = (0, 180, 5, 0),
        halign = :right
)

# colgap!(figure.layout, 2, -35)
save(figurefiles("Figure2.svg"), figure)
save(figurefiles("Figure2.png"), figure)
figure
```


## Supplement

### Other FSEA Plots

```julia
fig = Figure(resolution=(800, 2000))
A = GridLayout(fig[1,1])
B = GridLayout(fig[2,1])
C = GridLayout(fig[3,1])
D = GridLayout(fig[4,1])

gs = "GABA synthesis"
ixs = neuroactive[gs]
cs = filter(!isnan, cors[ixs])
acs = filter(!isnan, cors[Not(ixs)])

Resonance.plot_fsea!(A, cs, acs; label=gs)

gs = "GABA synthesis I"
ixs = neuroactive_full[gs]
cs = filter(!isnan, cors[ixs])
acs = filter(!isnan, cors[Not(ixs)])

Resonance.plot_fsea!(B, cs, acs; label=gs)

gs = "GABA synthesis II"
ixs = neuroactive_full[gs]
cs = filter(!isnan, cors[ixs])
acs = filter(!isnan, cors[Not(ixs)])

Resonance.plot_fsea!(C, cs, acs; label=gs)


gs = "GABA synthesis III"
ixs = neuroactive_full[gs]
cs = filter(!isnan, cors[ixs])
acs = filter(!isnan, cors[Not(ixs)])

Resonance.plot_fsea!(D, cs, acs; label=gs)
fig
```

```julia
fig = Figure(resolution=(800, 2000))
A = GridLayout(fig[1,1])
B = GridLayout(fig[2,1])
C = GridLayout(fig[3,1])
D = GridLayout(fig[4,1])

gs = "Propionate synthesis"
ixs = neuroactive[gs]
cs = filter(!isnan, cors[ixs])
acs = filter(!isnan, cors[Not(ixs)])

Resonance.plot_fsea!(A, cs, acs; label=gs)

gs = "Propionate synthesis I"
ixs = neuroactive_full[gs]
cs = filter(!isnan, cors[ixs])
acs = filter(!isnan, cors[Not(ixs)])

Resonance.plot_fsea!(B, cs, acs; label=gs)

gs = "Propionate synthesis II"
ixs = neuroactive_full[gs]
cs = filter(!isnan, cors[ixs])
acs = filter(!isnan, cors[Not(ixs)])

Resonance.plot_fsea!(C, cs, acs; label=gs)


gs = "Propionate synthesis III"
ixs = neuroactive_full[gs]
cs = filter(!isnan, cors[ixs])
acs = filter(!isnan, cors[Not(ixs)])

Resonance.plot_fsea!(D, cs, acs; label=gs)
fig
```


```julia
fig = Figure(resolution=(800, 1600))
A = GridLayout(fig[1,1])
B = GridLayout(fig[2,1])
C = GridLayout(fig[3,1])

gs = "Isovaleric acid synthesis"
ixs = neuroactive[gs]
cs = filter(!isnan, cors[ixs])
acs = filter(!isnan, cors[Not(ixs)])

Resonance.plot_fsea!(A, cs, acs; label=gs)

gs = "Isovaleric acid synthesis I (KADH pathway)"
ixs = neuroactive_full[gs]
cs = filter(!isnan, cors[ixs])
acs = filter(!isnan, cors[Not(ixs)])

Resonance.plot_fsea!(B, cs, acs; label=gs)

gs = "Isovaleric acid synthesis II (KADC pathway)"
ixs = neuroactive_full[gs]
cs = filter(!isnan, cors[ixs])
acs = filter(!isnan, cors[Not(ixs)])

Resonance.plot_fsea!(C, cs, acs; label=gs)

fig

figure = Figure()
A = GridLayout(figure[1,1])
B = GridLayout(figure[2,1])
CDEF = GridLayout(figure[1,2])
GHIJ = GridLayout(figure[2,2])

Axis(A[1,1])
Legend(A[2,1], [MarkerElement(; color=c, marker=:rect) for c in (:red,:blue,:orange,:teal)], ["thing$i" for i in 1:4]; orientation=:horizontal, tellheight=true, tellwidth=false)
Axis(B[1,1])
Legend(B[2,1], [MarkerElement(; color=c, marker=:rect) for c in (:red,:blue,:orange,:teal)], ["thing$i" for i in 1:4]; orientation=:horizontal, tellheight=true, tellwidth=false)

for (i, j) in zip((1,2,1,2), (1,1,2,2))
    Axis(CDEF[i,j]; title = "$i and $j")
end
for (i, j) in zip((1,2,1,2), (1,1,2,2))
    Axis(GHIJ[i,j]; title = "$i and $j")
end

colsize!(figure.layout, 1, Relative(1/3))
figure
```
