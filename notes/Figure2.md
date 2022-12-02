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
allcors              = JLD2.load(outputfiles("figure2_data.jld2"), "allcors")
u6cors               = JLD2.load(outputfiles("figure2_data.jld2"), "u6cors")
o18cors              = JLD2.load(outputfiles("figure2_data.jld2"), "o18cors")
allcors_age          = JLD2.load(outputfiles("figure2_data.jld2"), "allcors_age")
u6cors_age           = JLD2.load(outputfiles("figure2_data.jld2"), "u6cors_age")
o18cors_age          = JLD2.load(outputfiles("figure2_data.jld2"), "o18cors_age")
all_neuroactive      = JLD2.load(outputfiles("figure2_data.jld2"), "all_neuroactive")
u6_neuroactive       = JLD2.load(outputfiles("figure2_data.jld2"), "u6_neuroactive")
o18_neuroactive      = JLD2.load(outputfiles("figure2_data.jld2"), "o18_neuroactive")
all_neuroactive_full = JLD2.load(outputfiles("figure2_data.jld2"), "all_neuroactive_full")
u6_neuroactive_full  = JLD2.load(outputfiles("figure2_data.jld2"), "u6_neuroactive_full")
o18_neuroactive_full = JLD2.load(outputfiles("figure2_data.jld2"), "o18_neuroactive_full")
unimdata             = JLD2.load(outputfiles("figure2_data.jld2"), "unimdata")
```

## Load FSEA


```julia
allfsdf = CSV.read(outputfiles("fsea_all_consolidated.csv"), DataFrame)
allfsdf2 = CSV.read(outputfiles("fsea_all.csv"), DataFrame)
u6fsdf = CSV.read(outputfiles("fsea_u6_consolidated.csv"), DataFrame)
u6fsdf2 = CSV.read(outputfiles("fsea_u6.csv"), DataFrame)
o18fsdf = CSV.read(outputfiles("fsea_o18_consolidated.csv"), DataFrame)
o18fsdf2 = CSV.read(outputfiles("fsea_o18.csv"), DataFrame)
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
aax1 = Axis(A[1,1]; xlabel = "Age (months)", ylabel = "cogScore", xticks=(4:4:24))
aax2 = Axis(A[1,2]; xlabel = "Age (years)")
hideydecorations!(aax2)
linkyaxes!(aax1, aax2)

let
    u2y = findall(p-> !ismissing(p[2]) && p[1] <= 24, collect(zip(unimdata.ageMonths, unimdata.cogScore)))
    o2y = findall(p-> !ismissing(p[2]) && p[1] > 24, collect(zip(unimdata.ageMonths, unimdata.cogScore)))

    cs = ColorSchemes.colorschemes[:Set2_7]
    function colorage(age)
        age <= 36 && return cs[1] # mullen
        age <= 60 && return cs[2] # WPPSI
        return cs[3] # WISC
    end
    ages = unimdata.ageMonths[u2y]
    scatter!(aax1, ages, unimdata.cogScore[u2y]; color=colorage.(ages))

    ages = unimdata.ageMonths[o2y]
    scatter!(aax2, ages ./ 12, unimdata.cogScore[o2y]; color=colorage.(ages))
    
    vlines!(aax1, [6, 12, 18]; linestyle=:dash, color=:gray)
    # TODO: add colors for training/test sets in later models

    Legend(A[1, 3], [MarkerElement(; marker=:circle, color=c) for c in cs[1:3]],
                      ["Mullen", "WPPSI", "WISC"], "Assessment"; 
    )
end

colgap!(A, Fixed(4))

bax = Axis(B[1,1]; ylabel = "cogScore", xlabel = "age group (months)", xticks=(1:5, ["0-6", "6-12", "12-18", "18-24", "> 24"]))
let
    df = DataFrame(age = unimdata.ageMonths, score = unimdata.cogScore, date = unimdata.assessmentDate)
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
    ixs = u6_neuroactive_full[gs]
    cs = filter(!isnan, u6cors[ixs])
    acs = filter(!isnan, u6cors[Not(ixs)])

    Resonance.plot_fsea!(panel, cs, acs;
        label = replace(gs, "degradation"=> "degr.", "synthesis"=> "synth.", " (vitamin K2)"=> ""),
        ylabel = "")

end

let
    gs = "Menaquinone synthesis (vitamin K2) I"
    panel = D
    ixs = u6_neuroactive_full[gs]
    cs = filter(!isnan, u6cors[ixs])
    acs = filter(!isnan, u6cors[Not(ixs)])

    Resonance.plot_fsea!(panel, cs, acs;
        label = replace(gs, "degradation"=> "degr.", "synthesis"=> "synth.", " (vitamin K2)"=> ""),
        ylabel = "", xticks = -0.2:0.1:0.0)
end

let
    gs = "Propionate degradation I"
    panel = E
    ixs = u6_neuroactive_full[gs]
    cs = filter(!isnan, u6cors[ixs])
    acs = filter(!isnan, u6cors[Not(ixs)])

    Resonance.plot_fsea!(panel, cs, acs;
        label = replace(gs, "degradation"=> "degr.", "synthesis"=> "synth.", " (vitamin K2)"=> ""),
        ylabel = "", xticks = -0.7:0.2:0.0)
end

let
    gs = "Propionate synthesis I"
    panel = F
    ixs = u6_neuroactive_full[gs]
    cs = filter(!isnan, u6cors[ixs])
    acs = filter(!isnan, u6cors[Not(ixs)])

    Resonance.plot_fsea!(panel, cs, acs;
        label = replace(gs, "degradation"=> "degr.", "synthesis"=> "synth.", " (vitamin K2)"=> ""),
        ylabel = "")
end

```


```julia
let
    genesets = union(subset(allfsdf2, "qvalue"=> ByRow(<(0.2)), "cortest"=> ByRow(==("cogScorePercentile"))).geneset,
                     subset(u6fsdf2, "qvalue"=> ByRow(<(0.2)), "cortest"=> ByRow(==("cogScorePercentile"))).geneset,
                     subset(o18fsdf2, "qvalue"=> ByRow(<(0.2)), "cortest"=> ByRow(==("cogScorePercentile"))).geneset
    )
  
    df = sort(subset(allfsdf2, "geneset"=> ByRow(gs-> gs in genesets), "cortest"=>ByRow(==("cogScorePercentile"))), :geneset; rev=true)
    ax = Axis(G[1,1]; yticks = (1:nrow(df), replace.(df.geneset, r" \(.+\)" => "", "synthesis"=>"syn.", "degradation"=>"deg.")), 
                xlabel="correlation", title="All ages")
    m = median(allcors)
    colors = ColorSchemes.colorschemes[:RdBu_7]

    for (i, row) in enumerate(eachrow(df))
        sign = row.enrichment < 0 ? "neg" : "pos"
        c = row.qvalue > 0.2 ? :gray : 
            row.qvalue > 0.05 ? (sign == "neg" ? colors[3] : colors[5]) :
            row.qvalue > 0.01 ? (sign == "neg" ? colors[2] : colors[6]) :
            sign == "neg" ? colors[1] : colors[7]

        y = filter(!isnan, allcors[all_neuroactive_full[row.geneset]])
        scatter!(ax, y, rand(Normal(0, 0.1), length(y)) .+ i; color=(c,0.3), strokecolor=:gray, strokewidth=0.5)
        row.qvalue < 0.2 && lines!(ax, fill(median(y), 2), [i-0.4, i+0.4]; color = c, linewidth=2)
    end
    vlines!(ax, m; linestyle=:dash, color=:darkgray)

    ####

    df = sort(subset(u6fsdf2, "geneset"=> ByRow(gs-> gs in genesets), "cortest"=>ByRow(==("cogScorePercentile"))), :geneset; rev=true)
    ax = Axis(G[1,2]; yticks = (1:nrow(df), replace.(df.geneset, r" \(.+\)" => "", "synthesis"=>"syn.", "degradation"=>"deg.")), 
                xlabel="correlation", title="Under 6mo")
    hideydecorations!(ax, grid=false)

    m = median(u6cors)
    colors = ColorSchemes.colorschemes[:RdBu_7]

    for (i, row) in enumerate(eachrow(df))
        sign = row.enrichment < 0 ? "neg" : "pos"
        c = row.qvalue > 0.2 ? :gray : 
            row.qvalue > 0.05 ? (sign == "neg" ? colors[3] : colors[5]) :
            row.qvalue > 0.01 ? (sign == "neg" ? colors[2] : colors[6]) :
            sign == "neg" ? colors[1] : colors[7]

        y = filter(!isnan, u6cors[u6_neuroactive_full[row.geneset]])
        scatter!(ax, y, rand(Normal(0, 0.1), length(y)) .+ i; color=(c,0.3), strokecolor=:gray, strokewidth=0.5)
        row.qvalue < 0.2 && lines!(ax, fill(median(y), 2), [i-0.4, i+0.4]; color = c, linewidth=2)
    end
    vlines!(ax, m; linestyle=:dash, color=:darkgray)

    ####

    df = sort(subset(o18fsdf2, "geneset"=> ByRow(gs-> gs in genesets), "cortest"=>ByRow(==("cogScorePercentile"))), :geneset; rev=true)
    ax = Axis(G[1,3]; yticks = (1:nrow(df), replace.(df.geneset, r" \(.+?\)" => "", "synthesis"=>"syn.", "degradation"=>"deg.")), 
                xlabel="correlation", title="over 18")
    hideydecorations!(ax, grid=false)
    m = median(o18cors)
    colors = ColorSchemes.colorschemes[:RdBu_7]

    for (i, row) in enumerate(eachrow(df))
        sign = row.enrichment < 0 ? "neg" : "pos"
        c = row.qvalue > 0.2 ? :gray : 
            row.qvalue > 0.05 ? (sign == "neg" ? colors[3] : colors[5]) :
            row.qvalue > 0.01 ? (sign == "neg" ? colors[2] : colors[6]) :
            sign == "neg" ? colors[1] : colors[7]

        y = filter(!isnan, o18cors[o18_neuroactive_full[row.geneset]])
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
save(figurefiles("Figure2.svg"), fig)
save(figurefiles("Figure2.png"), fig)
figure
```


## Supplement

### Comparing normalized to un-normalized cogScores

```julia
fig = Figure()
A = GridLayout(fig[1,1])
B = GridLayout(fig[2,1])

aax1 = Axis(A[1,1]; xlabel = "Age (months)", ylabel = "cogScore", xticks=(4:4:24))
aax2 = Axis(A[1,2]; xlabel = "Age (years)")
hideydecorations!(aax2)
linkyaxes!(aax1, aax2)

aax3 = Axis(A[1,3]; xlabel = "Age (months)", ylabel = "cogScorePercentile", xticks=(4:4:24))
aax4 = Axis(A[1,4]; xlabel = "Age (years)")

let
    u2y = findall(p-> !ismissing(p[2]) && p[1] <= 24, collect(zip(unimdata.ageMonths, unimdata.cogScore)))
    o2y = findall(p-> !ismissing(p[2]) && p[1] > 24, collect(zip(unimdata.ageMonths, unimdata.cogScore)))

    cs = ColorSchemes.colorschemes[:Set2_7]
    function colorage(age)
        age <= 36 && return cs[1] # mullen
        age <= 60 && return cs[2] # WPPSI
        return cs[3] # WISC
    end
    ages = unimdata.ageMonths[u2y]
    scatter!(aax1, ages, unimdata.cogScore[u2y]; color=colorage.(ages))
    scatter!(aax3, ages, unimdata.cogScorePercentile[u2y]; color=colorage.(ages))

    ages = unimdata.ageMonths[o2y]
    scatter!(aax2, ages ./ 12, unimdata.cogScore[o2y]; color=colorage.(ages))
    scatter!(aax4, ages ./ 12, unimdata.cogScorePercentile[o2y]; color=colorage.(ages))
    
    vlines!(aax1, [6, 12, 18]; linestyle=:dash, color=:gray)
    vlines!(aax3, [6, 12, 18]; linestyle=:dash, color=:gray)
    # TODO: add colors for training/test sets in later models

    Legend(A[2, 1:4], [MarkerElement(; marker=:circle, color=c) for c in cs[1:3]],
                      ["Mullen", "WPPSI", "WISC"]; orientation=:horizontal, tellheight=true, tellwidth=false
    )
end



bax1 = Axis(B[1,1]; ylabel = "cogScore", xlabel = "age group (months)", xticks=(1:5, ["0-6", "6-12", "12-18", "18-24", "> 24"]))
let
    df = DataFrame(age = unimdata.ageMonths, score = unimdata.cogScore, date = unimdata.assessmentDate)
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
    scatter!(bax1, df.x .+ rand(Normal(0, 0.05), size(df, 1)) .+ [x < Date("2020-03-01") ? -0.15 : 0.15 for x in df.date], df.score; 
            color = [x < Date("2020-03-01") ? (:dodgerblue, 0.3) : (:orangered, 0.3) for x in df.date])
end

bax2 = Axis(B[1,2]; ylabel = "cogScorePercentile", xlabel = "age group (months)", xticks=(1:5, ["0-6", "6-12", "12-18", "18-24", "> 24"]))
let
    df = DataFrame(age = unimdata.ageMonths, score = unimdata.cogScorePercentile, date = unimdata.assessmentDate)
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
    scatter!(bax2, df.x .+ rand(Normal(0, 0.05), size(df, 1)) .+ [x < Date("2020-03-01") ? -0.15 : 0.15 for x in df.date], df.score; 
            color = [x < Date("2020-03-01") ? (:dodgerblue, 0.3) : (:orangered, 0.3) for x in df.date])
end

Legend(B[2,1:2], [MarkerElement(; color = :dodgerblue, marker=:circle), 
                  MarkerElement(; color = :orangered, marker=:circle)],
            ["Pre-covid", "Post-Covid"]; orientation=:horizontal, tellheight=true, tellwidth=false
)


save(figurefiles("Supp_Figure2.svg"), fig)
save(figurefiles("Supp_Figure2.png"), fig)
fig
```

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
