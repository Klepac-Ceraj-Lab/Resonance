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
```

## Data Loading

```julia
mdata = Resonance.load(Metadata())

unirefs = Resonance.load(UnirefProfiles(); timepoint_metadata = mdata) # this can take a bit
unirefs = unirefs[:, map(s-> !ismissing(s) && s < 120, get(unirefs, :ageMonths))]
unistrat =  filter(f-> hastaxon(f) || name(f) == "UNMAPPED", unirefs)
unirefs = filter(!hastaxon, unirefs) # don't use species stratification for summaries

metabolites = Resonance.load(MetabolicProfiles(); timepoint_metadata=mdata)
metabolites = metabolites[:, [!ismissing(a) && a < 14 for a in get(metabolites, :ageMonths)]]
isdefined(Main, :metdm) || (metdm = braycurtis(metabolites))
metpco = fit(MDS, metdm; distances=true)

brain = Resonance.load(Neuroimaging(); timepoint_metadata=mdata)

brain_roi = [
    "right-lateral-occipital",
    "left-lateral-occipital",
    "right-inferior-parietal",
    "left-inferior-parietal",
    "right-middle-temporal",
    "left-middle-temporal",
    "right-cerebellum-white-matter",
    "left-cerebellum-white-matter",
    "right-thalamus-proper",
    "left-thalamus-proper",
]
brainmeta = let
    brainsub = brain[:, startswith.(samplenames(brain), "FG")]
    df = DataFrame(:sample => samplenames(brainsub), (Symbol(reg) => vec(abundances(brainsub[reg, :])) for reg in brain_roi)...)
end

set!(unirefs, brainmeta)
```

## Calculate correlations

```julia
unimdata = DataFrame(Microbiome.metadata(unirefs))
allages = unique(subset(unimdata, :cogScorePercentile => ByRow(!ismissing)), :subject)
u6 = unique(subset(unimdata, :ageMonths => ByRow(<(6)), :cogScorePercentile => ByRow(!ismissing)), :subject)
o18 = unique(subset(unimdata, :ageMonths => ByRow(>(18)), :cogScorePercentile => ByRow(!ismissing)), :subject)

allcomm = let keepuni = vec(prevalence(unirefs[:, allages.sample]) .> 0)
    unirefs[keepuni, allages.sample]
end

u6comm = let keepuni = vec(prevalence(unirefs[:, u6.sample]) .> 0)
    unirefs[keepuni, u6.sample]
end

o18comm = let keepuni = vec(prevalence(unirefs[:, o18.sample]) .> 0)
    unirefs[keepuni, o18.sample]
end

allcors = vec(cor(get(allcomm, :cogScorePercentile), abundances(allcomm), dims=2))

u6cors = vec(cor(get(u6comm, :cogScorePercentile), abundances(u6comm), dims=2))
o18cors = vec(cor(get(o18comm, :cogScorePercentile), abundances(o18comm), dims=2))

allcors_age = vec(cor(get(allcomm, :ageMonths), abundances(allcomm), dims=2))
u6cors_age = vec(cor(get(u6comm, :ageMonths), abundances(u6comm), dims=2))
o18cors_age = vec(cor(get(o18comm, :ageMonths), abundances(o18comm), dims=2))
```

```julia
all_neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(allcomm)))
all_neuroactive_full = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(allcomm)); consolidate=false)
u6_neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(u6comm)))
u6_neuroactive_full = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(u6comm)); consolidate=false)
o18_neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(o18comm)))
o18_neuroactive_full = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(o18comm)); consolidate=false)
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
DEFG = GridLayout(figure[1:3,2]; alignmode=Outside())
C = GridLayout(DEFG[1,1])
D = GridLayout(DEFG[2,1])
E = GridLayout(DEFG[3,1])
F = GridLayout(DEFG[4,1])
G = GridLayout(figure[4,1:2]; alignmode=Outside())
```


```julia
aax1 = Axis(A[1,1]; xlabel = "Age (months)", ylabel = "cogScore", xticks=(4:4:24))
aax2 = Axis(A[1,2]; xlabel = "Age (years)")
hideydecorations!(aax2)
linkyaxes!(aax1, aax2)
colgap!(A, Fixed(4))
figure
let
    u2y = findall(p-> !ismissing(p[2]) && p[1] <= 24, collect(zip(get(unirefs, :ageMonths), get(unirefs, :cogScore))))
    o2y = findall(p-> !ismissing(p[2]) && p[1] > 24, collect(zip(get(unirefs, :ageMonths), get(unirefs, :cogScore))))

    scatter!(aax1, get(unirefs, :ageMonths)[u2y], get(unirefs, :cogScore)[u2y])
    scatter!(aax2, get(unirefs, :ageMonths)[o2y] ./ 12, get(unirefs, :cogScore)[o2y])
    
    vlines!(aax1, [6, 12, 18]; linestyle=:dash, color=:gray)
    # TODO: add colors for training/test sets in later models
end

bax = Axis(B[1,1]; ylabel = "cogScore", xlabel = "age group (months)", xticks=(1:5, ["0-6", "6-12", "12-18", "18-24", "> 24"]))
let
    df = DataFrame(age = get(unirefs, :ageMonths), score = get(unirefs, :cogScore), date = get(unirefs, :assessmentDate))
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

Legend(figure[3,1], [MarkerElement(; color = :dodgerblue, marker=:circle), MarkerElement(; color = :orangered, marker=:circle)],
               ["Pre-covid", "Post-covid"]; orientation=:horizontal, tellheight=true, tellwidth=false, framevisible=false)
rowgap!(figure.layout, 2, Fixed(0))
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

rowsize!(figure.layout, 4, Relative(2/5))
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


### Brain


```julia
fsdfbrain = let
    if isfile(outputfiles("fsea_brain.csv"))
        tmp = CSV.read(outputfiles("fsea_brain.csv"), DataFrame)
    else
        tmp = DataFrame()
        has_brain = ThreadsX.map(!ismissing, get(unirefs, Symbol(first(brain_roi))))
        keep = vec(prevalence(unirefs[:, has_brain]) .> 0)
        age = get(unirefs, :ageMonths)[has_brain]
        sex = categorical(get(unirefs, :sex)[has_brain])        
        neuroactive_full = Resonance.getneuroactive(
            map(f-> replace(f, "UniRef90_"=>""), featurenames(unirefs[keep, has_brain]));
            consolidate=false
        )
        mat = abundances(unirefs[keep, has_brain])
        hmin = map(eachrow(mat)) do row
           row = filter(>(0), row)
           isempty(row) ? 0 : minimum(row) / 2
        end
        
        mat .+= hmin

        ubr = CommunityProfile(log2.(mat), features(unirefs)[keep], samples(unirefs)[has_brain])

        for roi in brain_roi
            @info roi
            cors = ThreadsX.map(features(ubr)) do f
                df = DataFrame((;
                    brain = get(ubr, Symbol(roi)),
                    feature = vec(abundances(ubr[f, :])),
                    age, sex))
                mod = lm(@formula(brain ~ feature + age + sex), df; dropcollinear=false)
                DataFrame(coeftable(mod))."Coef."[2]
            end

            tmp2 = DataFrame(
                ThreadsX.map(collect(keys(neuroactive_full))) do gs
                    ixs = neuroactive_full[gs]
                    isempty(ixs) && return (; geneset = gs, U = NaN, median = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

                    cs = filter(!isnan, cors[ixs])
                    isempty(cs) && return (; geneset = gs, U = NaN, median = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

                    acs = filter(!isnan, cors[Not(ixs)])
                    mwu = MannWhitneyUTest(cs, acs)

                    return (; geneset = gs, U = mwu.U, median = mwu.median, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
                end
            )

            subset!(tmp2, :pvalue=> ByRow(!isnan))
            tmp2.roi .= roi
            append!(tmp, tmp2)
        end
        grp = groupby(tmp, :roi)
        DataFrames.transform!(grp, "pvalue" => (p-> adjust(collect(p), BenjaminiHochberg())) => "qvalue")
        sort!(tmp, :qvalue)
        # CSV.write(outputfiles("fsea_brain.csv"), tmp)
    end
    tmp
end

```


## Other

```julia
mfeats = commonname.(features(metabolites))
pyrdidx = only(ThreadsX.findall(f-> !ismissing(f) && contains(f, r"pyridoxine"i), mfeats))
gabaidx = only(ThreadsX.findall(f-> !ismissing(f) && contains(f, r"gamma-aminobutyric"i), mfeats))
glutidx = ThreadsX.findall(f-> !ismissing(f) && contains(f, r"^glutamic"i), mfeats)
butyidx = only(ThreadsX.findall(f-> !ismissing(f) && contains(f, r"^buty"i), mfeats))
propidx = only(ThreadsX.findall(f-> !ismissing(f) && contains(f, r"^propio"i), mfeats))
isovidx = only(ThreadsX.findall(f-> !ismissing(f) && contains(f, r"^isoval"i), mfeats))

pyrd = vec(abundances(metabolites[pyrdidx, :]))
gaba = vec(abundances(metabolites[gabaidx, :]))
glut1 = vec(abundances(metabolites[glutidx[1], :]))
glut2 = vec(abundances(metabolites[glutidx[2], :]))
buty = vec(abundances(metabolites[butyidx, :]))
prop = vec(abundances(metabolites[propidx, :]))
isov = vec(abundances(metabolites[isovidx, :]))

```

```julia

Ga = Axis(G[1,1])
Gb = Axis(G[1,2]; xlabel = L"log_e(GABA)", ylabel = L"log_e(Glutamate)")

pcom = plot_pcoa!(Ga, metpco; color=get(metabolites, :ageMonths))
sc = scatter!(Gb, log.(gaba), log.((glut1 .+ glut2) ./ 2); color = get(metabolites, :ageMonths))
Colorbar(G[1, 3], sc; label = "Age (Months)")

```

## Supplement

### Comparing normalized to un-normalized cogScores

```julia
fig = Figure()

aax = Axis(fig[1,1]; ylabel = "cogScore", xlabel = "age group (months)", xticks=(1:5, ["0-6", "6-12", "12-18", "18-24", "> 24"]))
let
    df = DataFrame(age = get(unirefs, :ageMonths), score = get(unirefs, :cogScore), date = get(unirefs, :assessmentDate))
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
    scatter!(aax, df.x .+ rand(Normal(0, 0.05), size(df, 1)) .+ [x < Date("2020-03-01") ? -0.15 : 0.15 for x in df.date], df.score; 
            color = [x < Date("2020-03-01") ? (:dodgerblue, 0.3) : (:orangered, 0.3) for x in df.date])
end

bax = Axis(fig[1,2]; ylabel = "cogScorePercentile", xlabel = "age group (months)", xticks=(1:5, ["0-6", "6-12", "12-18", "18-24", "> 24"]))
let
    df = DataFrame(age = get(unirefs, :ageMonths), score = get(unirefs, :cogScorePercentile), date = get(unirefs, :assessmentDate))
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

Legend(fig[2,1:2], [MarkerElement(; color = :dodgerblue, marker=:circle), 
                  MarkerElement(; color = :orangered, marker=:circle)],
            ["Pre-covid", "Post-Covid"]; orientation=:horizontal, tellheight=true, tellwidth=false
)

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
```
### Metabolites

```julia
bmoi = [ # metabolites of interest
    "pyridoxine", # vit B6
    "gaba",
    "glutamate",
    "acetate",
    "propionate",
    "butyrate"
]

mfeats = commonname.(features(metabolites))
pyrdidx = only(ThreadsX.findall(f-> !ismissing(f) && contains(f, r"pyridoxine"i), mfeats))
gabaidx = only(ThreadsX.findall(f-> !ismissing(f) && contains(f, r"gamma-aminobutyric"i), mfeats))
glutidx = ThreadsX.findall(f-> !ismissing(f) && contains(f, r"^glutamic"i), mfeats)
butyidx = only(ThreadsX.findall(f-> !ismissing(f) && contains(f, r"^buty"i), mfeats))
propidx = only(ThreadsX.findall(f-> !ismissing(f) && contains(f, r"^propio"i), mfeats))
isovidx = only(ThreadsX.findall(f-> !ismissing(f) && contains(f, r"^isoval"i), mfeats))

pyrd = vec(abundances(metabolites[pyrdidx, :]))
gaba = vec(abundances(metabolites[gabaidx, :]))
glut1 = vec(abundances(metabolites[glutidx[1], :]))
glut2 = vec(abundances(metabolites[glutidx[2], :]))
buty = vec(abundances(metabolites[butyidx, :]))
prop = vec(abundances(metabolites[propidx, :]))
isov = vec(abundances(metabolites[isovidx, :]))

```

```julia

Ga = Axis(G[1,1])
Gb = Axis(G[1,2]; xlabel = L"log_e(GABA)", ylabel = L"log_e(Glutamate)")

pcom = plot_pcoa!(Ga, metpco; color=get(metabolites, :ageMonths))
sc = scatter!(Gb, log.(gaba), log.((glut1 .+ glut2) ./ 2); color = get(metabolites, :ageMonths))
Colorbar(G[1, 3], sc; label = "Age (Months)")

save(figurefiles("Figure2.svg"), figure)
save(figurefiles("Figure2.png"), figure)
figure

```


```julia
fig = Figure()
ax1 = Axis(fig[1,1]; xlabel = L"log_e(pyridoxine)", ylabel = L"log_e(GABA)")
ax2 = Axis(fig[1,2]; xlabel = L"log_e(pyridoxine)", ylabel = L"log_e(Glutamate)")
ax3 = Axis(fig[2,1]; xlabel = L"log_e(GABA)", ylabel = L"log_e(Glutamate)")
ax4 = Axis(fig[2,2]; xlabel = L"log_e(GABA[1])", ylabel = L"log_e(Glutamate[2])")

sc1 = scatter!(ax1, log.(pyrd), log.(gaba); color = get(metabolites, :ageMonths))
sc2 = scatter!(ax2, log.(pyrd), log.((glut1 .+ glut2) ./ 2); color = get(metabolites, :ageMonths))

sc4 = scatter!(ax4, log.(glut1), log.(glut2); color = get(metabolites, :ageMonths))


```

## GABA and GABA genes

```julia
mtbsubtp = collect(zip(get(metabolites, :subject), get(metabolites, :timepoint)))
uniol = findall(stp -> stp in mtbsubtp, collect(zip(get(unirefs, :subject), get(unirefs, :timepoint))))
ur = unirefs[:, uniol]
ursubtp = collect(zip(get(ur, :subject), get(ur, :timepoint)))
mtbol = findall(stp -> stp in ursubtp, collect(zip(get(metabolites, :subject), get(metabolites, :timepoint))))

cors = vec(cor(glut1[mtbol] .+ glut2[mtbol], abundances(ur), dims=2))

ixs = neuroactive_full["Glutamate degradation II"]
cs = filter(!isnan, cors[ixs])
acs = filter(!isnan, cors[Not(ixs)])
mwu = MannWhitneyUTest(cs, acs)

Resonance.plot_fsea(cs, acs; label="Glutamate synthesis II")

```

```julia
fig = Figure()
ax = Axis(fig[1,1])

scatter!(ax, vec(abundances(ur[ixs[10], :])), log.((glut1[mtbol] .+ glut2[mtbol]) ./ 2))
fig
```


## Gene / metabolite correspondance

```julia
gs = "GABA synthesis I"
ixs = neuroactive_full[gs]

```

