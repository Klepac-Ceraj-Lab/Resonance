# Figure3 - functional analysis

```julia
using Resonance
using CairoMakie
using Statistics
using HypothesisTests
using MultipleTesting
using KernelDensity
using ThreadsX
```

## Data Loading

```julia
mdata = Resonance.load(Metadata())

unirefs = Resonance.load(UnirefProfiles(); timepoint_metadata = mdata) # this can take a bit
unirefs = filter(!hastaxon, unirefs) # don't use species stratification for summaries
```

## Calculate correlations

```julia

unimdata = DataFrame(metadata(unirefs))
allages = unique(subset(unimdata, :cogScore => ByRow(!ismissing)), :subject)
u6 = unique(subset(unimdata, :ageMonths => ByRow(<(6)), :cogScore => ByRow(!ismissing)), :subject)
o18 = unique(subset(unimdata, :ageMonths => ByRow(>(18)), :cogScore => ByRow(!ismissing)), :subject)

keepuni = vec(prevalence(unirefs[:, allages.sample]) .> 0)

cors = vec(cor(allages.cogScore, abundances(unirefs[keepuni, allages.sample]), dims=2))
neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(unirefs[keepuni, :])))
neuroactive_full = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(unirefs[keepuni, :])); consolidate=false)
```

```julia
fsdf = DataFrame(
    ThreadsX.map(collect(keys(neuroactive))) do gs
        ixs = neuroactive[gs]
        isempty(ixs) && return (; geneset = gs, U = NaN, median = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, cors[ixs])
        isempty(cs) && return (; geneset = gs, U = NaN, median = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, cors[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)

        return (; geneset = gs, U = mwu.U, median = mwu.median, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end
)

subset!(fsdf, :pvalue=> ByRow(!isnan))
fsdf.qvalue = adjust(fsdf.pvalue, BenjaminiHochberg())
sort!(fsdf, :qvalue)
CSV.write(outputfiles("fsea_consolidated.csv"), fsdf)
```

```julia
fsdf2 = DataFrame(
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

subset!(fsdf2, :pvalue=> ByRow(!isnan))
fsdf2.qvalue = adjust(fsdf2.pvalue, BenjaminiHochberg())
sort!(fsdf2, :qvalue)
CSV.write(outputfiles("fsea_all.csv"), fsdf2)
```

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


```julia
k = kde(hcat(vec(prevalence(unirefs[keepuni, allages.sample])), cors))
```

```julia
figure = Figure(resolution=(1200, 1200))

A = Axis(figure[1,1]; xlabel="prevalence", ylabel="correlation")
BCDE = GridLayout(figure[1:2,2])
B = GridLayout(BCDE[1,1])
C = GridLayout(BCDE[2,1])
D = GridLayout(BCDE[3,1])
E = Axis(BCDE[4,1])
```


```julia
heatmap!(A,k.x,k.y, log2.(k.density .+ minimum(k.density)/2))


for (gs, panel) in zip(("Propionate degradation I", "Glutamate degradation I", "GABA synthesis I"), (B,C,D))
    ixs = neuroactive_full[gs]
    cs = filter(!isnan, cors[ixs])
    acs = filter(!isnan, cors[Not(ixs)])

    Resonance.plot_fsea!(panel, cs, acs; label=gs)
end

Resonance.plot_corrband!(E, cors)

figure
```

```julia
gs = "Glutamate degradation I"
ixs = neuroactive_full[gs]
cs = filter(!isnan, cors[ixs])
acs = filter(!isnan, cors[Not(ixs)])

Resonance.plot_fsea!(C, cs, acs; label=gs, bandres=0)
```

```julia
gs = 
ixs = neuroactive_full[gs]
cs = filter(!isnan, cors[ixs])
acs = filter(!isnan, cors[Not(ixs)])

Resonance.plot_fsea!(D, cs, acs; label=gs)
figure
```

