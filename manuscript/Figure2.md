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
cors_00to120              = JLD2.load(scratchfiles("figure2", "figure2_data.jld2"), "cors_00to120")
cors_00to06               = JLD2.load(scratchfiles("figure2", "figure2_data.jld2"), "cors_00to06")
cors_18to120              = JLD2.load(scratchfiles("figure2", "figure2_data.jld2"), "cors_18to120")
cors_00to120_age          = JLD2.load(scratchfiles("figure2", "figure2_data.jld2"), "cors_00to120_age")
cors_00to06_age           = JLD2.load(scratchfiles("figure2", "figure2_data.jld2"), "cors_00to06_age")
cors_18to120_age          = JLD2.load(scratchfiles("figure2", "figure2_data.jld2"), "cors_18to120_age")
neuroactive_00to120       = JLD2.load(scratchfiles("figure2", "figure2_data.jld2"), "neuroactive_00to120")
neuroactive_00to06        = JLD2.load(scratchfiles("figure2", "figure2_data.jld2"), "neuroactive_00to06")
neuroactive_00to120       = JLD2.load(scratchfiles("figure2", "figure2_data.jld2"), "neuroactive_00to120")
neuroactive_full_00to120  = JLD2.load(scratchfiles("figure2", "figure2_data.jld2"), "neuroactive_full_00to120")
neuroactive_full_00to06   = JLD2.load(scratchfiles("figure2", "figure2_data.jld2"), "neuroactive_full_00to06")
neuroactive_full_00to120  = JLD2.load(scratchfiles("figure2", "figure2_data.jld2"), "neuroactive_full_00to120")
unimdata                  = JLD2.load(scratchfiles("figure2", "figure2_data.jld2"), "unimdata")
```

## Load FSEA


```julia
speclms = CSV.read(tablefiles("lms_species_18to120.csv"), DataFrame)
speclms_pa = CSV.read(tablefiles("lms_species_18to120_pa.csv"), DataFrame)
fsdf_00to120 = CSV.read(scratchfiles("figure2", "fsea_consolidated_00to120.csv"), DataFrame)
fsdf2_00to120 = CSV.read(scratchfiles("figure2", "fsea_all_00to120.csv"), DataFrame)
fsdf_00to06 = CSV.read(scratchfiles("figure2", "fsea_consolidated_00to06.csv"), DataFrame)
fsdf2_00to06 = CSV.read(scratchfiles("figure2", "fsea_all_00to06.csv"), DataFrame)
fsdf_18to120 = CSV.read(scratchfiles("figure2", "fsea_consolidated_18to120.csv"), DataFrame)
fsdf2_18to120 = CSV.read(scratchfiles("figure2", "fsea_all_18to120.csv"), DataFrame)
```


## Plotting

```julia
figure = Figure(resolution=(900, 900))

A = GridLayout(figure[1,1]; alignmode=Outside())
B = GridLayout(figure[2,1]; alignmode=Outside())
CDEF = GridLayout(figure[1:2,2]; alignmode=Outside())
C = GridLayout(CDEF[1,1])
D = GridLayout(CDEF[2,1])
E = GridLayout(CDEF[3,1])
F = GridLayout(CDEF[4,1])
G = GridLayout(figure[3,1:2]; alignmode=Outside())
```


```julia


colgap!(A, Fixed(4))

bax = Axis(B[1,1]; ylabel = "cogScore", xlabel = "age group (months)", xticks=(1:5, ["0-6", "6-12", "12-18", "18-24", "> 24"]))
let
    df = DataFrame(age = unimdata.ageMonths, score = unimdata.cogScore, date = unimdata.date)
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
    scatter!(bax, df.x .+ rand(Normal(0, 0.05), size(df, 1)) .+ [x < Date("2020-03-01") ? -0.15 : 0.15 for x in df.date], df.score; 
            color = [x < Date("2020-03-01") ? (:dodgerblue, 0.3) : (:orangered, 0.3) for x in df.date])
end

Legend(B[1,2], [MarkerElement(; color = :dodgerblue, marker=:circle), MarkerElement(; color = :orangered, marker=:circle)],
               ["Pre-covid", "Post-covid"]
)
colgap!(B, Fixed(4))
```



```julia

let
    gs = "Acetate synthesis II"
    panel = C
    ixs = neuroactive_full_00to06[gs]
    cs = filter(!isnan, cors_00to06[ixs])
    acs = filter(!isnan, cors_00to06[Not(ixs)])

    Resonance.plot_fsea!(panel, cs, acs;
        label = replace(gs, "degradation"=> "degr.", "synthesis"=> "synth.", " (vitamin K2)"=> ""),
        ylabel = "")

end

let
    gs = "Menaquinone synthesis (vitamin K2) I"
    panel = D
    ixs = neuroactive_full_00to06[gs]
    cs = filter(!isnan, cors_00to06[ixs])
    acs = filter(!isnan, cors_00to06[Not(ixs)])

    Resonance.plot_fsea!(panel, cs, acs;
        label = replace(gs, "degradation"=> "degr.", "synthesis"=> "synth.", " (vitamin K2)"=> ""),
        ylabel = "", xticks = -0.2:0.1:0.0)
end

let
    gs = "Propionate degradation I"
    panel = E
    ixs = neuroactive_full_00to06[gs]
    cs = filter(!isnan, cors_00to06[ixs])
    acs = filter(!isnan, cors_00to06[Not(ixs)])

    Resonance.plot_fsea!(panel, cs, acs;
        label = replace(gs, "degradation"=> "degr.", "synthesis"=> "synth.", " (vitamin K2)"=> ""),
        ylabel = "", xticks = -0.7:0.2:0.0)
end

let
    gs = "Propionate synthesis I"
    panel = F
    ixs = neuroactive_full_00to06[gs]
    cs = filter(!isnan, cors_00to06[ixs])
    acs = filter(!isnan, cors_00to06[Not(ixs)])

    Resonance.plot_fsea!(panel, cs, acs;
        label = replace(gs, "degradation"=> "degr.", "synthesis"=> "synth.", " (vitamin K2)"=> ""),
        ylabel = "")
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
                xlabel="correlation", title="All ages")
    m = median(cors_00to120)
    colors = ColorSchemes.colorschemes[:RdBu_7]

    for (i, row) in enumerate(eachrow(df))
        sign = row.enrichment < 0 ? "neg" : "pos"
        c = row.qvalue > 0.2 ? :gray : 
            row.qvalue > 0.05 ? (sign == "neg" ? colors[3] : colors[5]) :
            row.qvalue > 0.01 ? (sign == "neg" ? colors[2] : colors[6]) :
            sign == "neg" ? colors[1] : colors[7]

        y = filter(!isnan, cors_00to120[neuroactive_full_00to120[row.geneset]])
        scatter!(ax, y, rand(Normal(0, 0.1), length(y)) .+ i; color=(c,0.3), strokecolor=:gray, strokewidth=0.5)
        row.qvalue < 0.2 && lines!(ax, fill(median(y), 2), [i-0.4, i+0.4]; color = c, linewidth=2)
    end
    vlines!(ax, m; linestyle=:dash, color=:darkgray)

    ####

    df = sort(subset(fsdf2_00to06, "geneset"=> ByRow(gs-> gs in genesets), "cortest"=>ByRow(==("cogScore"))), :geneset; rev=true)
    ax = Axis(G[1,2]; yticks = (1:nrow(df), replace.(df.geneset, r" \(.+\)" => "", "synthesis"=>"syn.", "degradation"=>"deg.")), 
                xlabel="correlation", title="Under 6mo")
    hideydecorations!(ax, grid=false)

    m = median(cors_00to06)
    colors = ColorSchemes.colorschemes[:RdBu_7]

    for (i, row) in enumerate(eachrow(df))
        sign = row.enrichment < 0 ? "neg" : "pos"
        c = row.qvalue > 0.2 ? :gray : 
            row.qvalue > 0.05 ? (sign == "neg" ? colors[3] : colors[5]) :
            row.qvalue > 0.01 ? (sign == "neg" ? colors[2] : colors[6]) :
            sign == "neg" ? colors[1] : colors[7]

        y = filter(!isnan, cors_00to06[neuroactive_full_00to06[row.geneset]])
        scatter!(ax, y, rand(Normal(0, 0.1), length(y)) .+ i; color=(c,0.3), strokecolor=:gray, strokewidth=0.5)
        row.qvalue < 0.2 && lines!(ax, fill(median(y), 2), [i-0.4, i+0.4]; color = c, linewidth=2)
    end
    vlines!(ax, m; linestyle=:dash, color=:darkgray)

    ####

    df = sort(subset(fsdf2_18to120, "geneset"=> ByRow(gs-> gs in genesets), "cortest"=>ByRow(==("cogScore"))), :geneset; rev=true)
    ax = Axis(G[1,3]; yticks = (1:nrow(df), replace.(df.geneset, r" \(.+?\)" => "", "synthesis"=>"syn.", "degradation"=>"deg.")), 
                xlabel="correlation", title="over 18")
    hideydecorations!(ax, grid=false)
    m = median(cors_18to120)
    colors = ColorSchemes.colorschemes[:RdBu_7]

    for (i, row) in enumerate(eachrow(df))
        @info i
        sign = row.enrichment < 0 ? "neg" : "pos"
        c = row.qvalue > 0.2 ? :gray : 
            row.qvalue > 0.05 ? (sign == "neg" ? colors[3] : colors[5]) :
            row.qvalue > 0.01 ? (sign == "neg" ? colors[2] : colors[6]) :
            sign == "neg" ? colors[1] : colors[7]

        y = filter(!isnan, cors_18to120[neuroactive_full_18to120[row.geneset]])
        scatter!(ax, y, rand(Normal(0, 0.1), length(y)) .+ i; color=(c,0.3), strokecolor=:gray, strokewidth=0.5)
        row.qvalue < 0.2 && lines!(ax, fill(median(y), 2), [i-0.4, i+0.4]; color = c, linewidth=2)
    end
    vlines!(ax, m; linestyle=:dash, color=:darkgray)

    Legend(G[1,4], [MarkerElement(; color = (c, 0.3),
                                    marker=:circle,
                                    strokecolor=:gray,
                                    strokewidth=0.5) for c in colors[[1:3..., 5:7...]]],
                   ["(-) q < 0.01", "(-) q < 0.05", "(-) q < 0.2", 
                    "(+) q < 0.01", "(+) q < 0.05", "(+) p < 0.2"])
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
