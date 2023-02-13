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
speclms_00to120 = subset(CSV.read(tablefiles("lms_species_00to120.csv"), DataFrame), "t"=> ByRow(!isnan))
speclms_00to06 = subset(CSV.read(tablefiles("lms_species_00to06.csv"), DataFrame), "t"=> ByRow(!isnan))
speclms_18to120 = subset(CSV.read(tablefiles("lms_species_18to120.csv"), DataFrame), "t"=> ByRow(!isnan))

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
lms_mat = let 
    df = vcat(speclms_00to120, speclms_00to06, speclms_18to120)
    df.group = [fill("00to120", nrow(speclms_00to120)); fill("00to06", nrow(speclms_00to06)); fill("18to120", nrow(speclms_18to120))]
    sigs = subset(df, "qvalue"=> ByRow(<(0.2)))
    sig_feats = unique(sigs.feature)
    gdf = groupby(df, ["feature", "group"])
    DataFrame(ThreadsX.map(eachindex(sig_feats)) do i
        ft = sig_feats[i]
        q_00to120 = haskey(gdf, (; feature=ft, group="00to120")) ? only(gdf[(; feature=ft, group="00to120")].qvalue) : 1.0
        q_00to06 = haskey(gdf, (; feature=ft, group="00to06")) ? only(gdf[(; feature=ft, group="00to06")].qvalue) : 1.0
        q_18to120 = haskey(gdf, (; feature=ft, group="18to120")) ? only(gdf[(; feature=ft, group="18to120")].qvalue) : 1.0
        corr_00to120 = cor(taxdf[taxdf.filter_00to120, "cogScore"], taxdf[taxdf.filter_00to120, ft])
        corr_00to06 = cor(taxdf[taxdf.filter_00to06, "cogScore"], taxdf[taxdf.filter_00to06, ft])
        corr_18to120 = cor(taxdf[taxdf.filter_18to120, "cogScore"], taxdf[taxdf.filter_18to120, ft])
        (; feature = ft, 
           prev_00to120 = only(prevalence(taxa[Regex(ft), mdata.filter_00to120])),
           meanab_00to120 = mean(abundances(taxa[Regex(ft), mdata.filter_00to120])),
           corr_00to120,
           q_00to120,
           prev_00to06 = only(prevalence(taxa[Regex(ft), mdata.filter_00to06])),
           meanab_00to06 = mean(abundances(taxa[Regex(ft), mdata.filter_00to06])),
           corr_00to06,
           q_00to06, 
           prev_18to120 = only(prevalence(taxa[Regex(ft), mdata.filter_18to120])),
           meanab_18to120 = mean(abundances(taxa[Regex(ft), mdata.filter_18to120])),
           corr_18to120,
           q_18to120,
        )
   end)
end
```

## Plotting

```julia
figure = Figure(resolution=(1200, 900))

AB = GridLayout(figure[1,1])
CDEF = GridLayout(figure[2,1])
C = GridLayout(CDEF[1, 1])
D = GridLayout(CDEF[1, 2])
E = GridLayout(CDEF[1, 3])
F = GridLayout(CDEF[1, 4])
G = GridLayout(figure[3,1], alignmode=Mixed(; left=0))
```

```julia
aax = Axis(AB[1,1]; title = "over 18mo", ylabel=L"$-log_2(P)$", xlabel = "Coef.",
                xticks = ([-1.5e-3, 0.0, 1.5e-3], ["-1.5e-3", "0.0", "1.5e-3"]),
                limits = ((-1.5e-3, 1.5e-3), nothing),
                xminorticksvisible=true, xminorticks = IntervalsBetween(3))

scatter!(aax, speclms_18to120.coef, -1 .* log2.(speclms_18to120.qvalue);
    color = map(q-> q < 0.2 ? :orange : :dodgerblue, speclms_18to120.qvalue))
```
```julia
B = GridLayout(AB[1, 2:3])

bax1 = let groups = ["0-120m", "0-6m", "18-120m"]
    Axis(B[1, 1]; xticks = (1.5:length(groups) + 0.5, groups), 
                  yticks = (1.5:nrow(lms_mat) + 0.5, replace.(lms_mat.feature, "_"=>" ")),
                  yticklabelfont = "TeX Gyre Heros Makie Italic",
                  yticklabelsize = 14,
                  xticklabelsize = 12
                  )
end
bax2 = let groups = ["0-120m", "0-6m", "18-120m"]
    Axis(B[1, 2]; xticks = (1.5:length(groups) + 0.5, groups), 
                  xticklabelsize = 12)
end
bax2.yticklabelsvisible = false
bax2.yticksvisible = false
bax3 = let groups = ["0-120m", "0-6m", "18-120m"]
    Axis(B[1, 3]; xticks = (1.5:length(groups) + 0.5, groups), 
                  xticklabelsize = 12)
end
bax3.yticksvisible = false
bax3.yticklabelsvisible = false

for ax in (bax1, bax2, bax3)
    tightlimits!(ax)
end

let clrs = vcat((lms_mat[:, "prev_$group"] for group in ("00to120", "00to06", "18to120"))...)
    poly!(bax1, [Rect(j, i, 1, 1) for j in 1:3 for i in eachindex(lms_mat.feature)]; 
        color = clrs, colormap = :viridis)
    Colorbar(B[2,1]; colormap=:viridis, colorrange=extrema(clrs),
                    vertical = false, flipaxis = false,
                    label="Prevalence",
                    )
end

let clrs = log.(vcat((lms_mat[:, "meanab_$group"] for group in ("00to120", "00to06", "18to120"))...))
    poly!(bax2, [Rect(j, i, 1, 1) for j in 1:3 for i in eachindex(lms_mat.feature)]; 
        color = clrs, colormap = :magma)
    Colorbar(B[2,2]; colormap=:magma, colorrange=extrema(clrs),
                    vertical = false, flipaxis = false,
                    label="Mean abundance (log)",
                    )
end

let clrs = [isnan(x) ? 0.0 : x for x in vcat((lms_mat[:, "corr_$group"] for group in ("00to120", "00to06", "18to120"))...)]
    
    poly!(bax3, [Rect(j, i, 1, 1) for j in 1:3 for i in eachindex(lms_mat.feature)]; 
        color = clrs, colorrange=(-0.3, 0.3),
        colormap = :RdBu)
    
    
    Colorbar(B[2,3]; colormap=:RdBu,
                    vertical = false, flipaxis = false,
                    label="Correlation",
                    limits=(-0.3, 0.3))

    qs = vcat((lms_mat[:, "q_$group"] for group in ("00to120", "00to06", "18to120"))...)
    annotations!(bax3, map(x-> x < 0.05 ? "**" : x < 0.2 ? "*" : "", qs), [Point2f(j + 0.5, i + 0.2) for j in 1:3 for i in eachindex(lms_mat.feature)];
    color = [abs(x) > 0.2 ? :white : :black for x in clrs],
    align=(:center, :center))
end


rowgap!(B, Fixed(4))
```


```julia
let
    gs = "Acetate synthesis II"
    panel = C
    ixs = neuroactive_00to06[gs]
    cs = filter(!isnan, cors_00to06.t[ixs])
    acs = filter(!isnan, cors_00to06.t[Not(ixs)])

    (_, cax, _) = Resonance.plot_fsea!(panel, cs, acs;
        label = replace(gs, "degradation"=> "degr.", "synthesis"=> "synth.", " (vitamin K2)"=> ""),
        ylabel = "")
    Label(panel[3,1], "Under 6m"; tellheight=true, tellwidth=false)
    rowgap!(panel, 2, Fixed(4))
    
    gs = "Propionate degradation I"
    panel = D
    ixs = neuroactive_00to06[gs]
    cs = filter(!isnan, cors_00to06.t[ixs])
    acs = filter(!isnan, cors_00to06.t[Not(ixs)])

    (_, dax, _) = Resonance.plot_fsea!(panel, cs, acs;
        label = replace(gs, "degradation"=> "degr.", "synthesis"=> "synth.", " (vitamin K2)"=> ""),
        ylabel = "", xticks = -0.7:0.2:0.0)
    Label(panel[3,1], "Under 6m"; tellheight=true, tellwidth=false)
    rowgap!(panel, 2, Fixed(4))
    
    gs = "Glutamate synthesis II"
    panel = E
    ixs = neuroactive_18to120[gs]
    cs = filter(!isnan, cors_18to120.t[ixs])
    acs = filter(!isnan, cors_18to120.t[Not(ixs)])

    (_, eax, _) =  Resonance.plot_fsea!(panel, cs, acs;
        label = replace(gs, "degradation"=> "degr.", "synthesis"=> "synth.", " (vitamin K2)"=> ""),
        ylabel = "", xticks = -0.2:0.1:0.0)
    Label(panel[3,1], "Over 18m"; tellheight=true, tellwidth=false)
    rowgap!(panel, 2, Fixed(4))


    gs = "Propionate degradation I"
    panel = F
    ixs = neuroactive_18to120[gs]
    cs = filter(!isnan, cors_18to120.t[ixs])
    acs = filter(!isnan, cors_18to120.t[Not(ixs)])

    (_, fax, _) = Resonance.plot_fsea!(panel, cs, acs;
        label = replace(gs, "degradation"=> "degr.", "synthesis"=> "synth.", " (vitamin K2)"=> ""),
        ylabel = "")
    Label(panel[3,1], "Over 18m"; tellheight=true, tellwidth=false)
    rowgap!(panel, 2, Fixed(4))
    linkyaxes!(cax, dax, eax, fax)
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
                xlabel="T stat", title="All ages", alignmode=Mixed(; left=-15))
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

# rowsize!(figure.layout, 3, Relative(2/5))
colgap!(G, Fixed(4))
```

```julia
Label(AB[1, 1, TopLeft()], "A",
    fontsize = 26,
    font = "Open Sans Bold",
    padding = (0, 30, 5, 0),
    halign = :right)
Label(AB[1, 2, TopLeft()], "B",
    fontsize = 26,
    font = "Open Sans Bold",
    padding = (0, 150, 5, 0),
    halign = :right)

for (label, layout) in zip(["C", "D", "E", "F"], [C, D, E, F])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 26,
        font = "Open Sans Bold",
        padding = (0, 30, 5, 0),
        halign = :right)
end

Label(figure[3, 1, TopLeft()], "G";
             fontsize = 26,
             font = "Open Sans Bold",
             padding = (0, 30, 5, 0),
             halign = :right
       )

rowsize!(figure.layout, 1, Relative(2/5))
colsize!(G, 1, Relative(3/8))
save("manuscript/assets/Figure2.png", figure)
save(figurefiles("Figure2.svg"), figure)
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
