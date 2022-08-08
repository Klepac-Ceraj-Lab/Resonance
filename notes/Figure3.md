# Figure3 - functional analysis

```julia
using Resonance
using CairoMakie
using Statistics
using HypothesisTests
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
CSV.write(outputfiles("fsea_all.csv"), fsdf)

pretty_table(first(fsdf, 10); backend = Val(:latex))
```

```julia
figure = Figure(resolution=(1200, 800))

A = Axis(figure[1,1]; xlabel="prevalence", ylabel="correlation")
B = GridLayout(figure[1,2])
C = GridLayout(figure[2,1])
D = GridLayout(figure[2,2])
```

```julia

```

```julia
gs = "Menaquinone synthesis"
ixs = neuroactive[gs]
cs = filter(!isnan, cors[ixs])
acs = filter(!isnan, cors[Not(ixs)])

Resonance.plot_fsea(cs, acs; label=gs)
```

```julia
gs = "Propionate degradation"
ixs = neuroactive[gs]
cs = filter(!isnan, cors[ixs])
acs = filter(!isnan, cors[Not(ixs)])

Resonance.plot_fsea!(B, cs, acs; label=gs)
```

```julia
gs = "Glutamate degradation"
ixs = neuroactive[gs]
cs = filter(!isnan, cors[ixs])
acs = filter(!isnan, cors[Not(ixs)])

Resonance.plot_fsea!(C, cs, acs; label=gs)
```

```julia
gs = "Tryptophan synthesis"
ixs = neuroactive[gs]
cs = filter(!isnan, cors[ixs])
acs = filter(!isnan, cors[Not(ixs)])

Resonance.plot_fsea(cs, acs; label=gs)
```

```julia
gs = "GABA synthesis"
ixs = neuroactive[gs]
cs = filter(!isnan, cors[ixs])
acs = filter(!isnan, cors[Not(ixs)])

Resonance.plot_fsea!(D, cs, acs; label=gs)
figure
```

```julia
gs = "Acetate synthesis"
ixs = neuroactive[gs]
cs = filter(!isnan, cors[ixs])
acs = filter(!isnan, cors[Not(ixs)])

Resonance.plot_fsea(cs, acs; label=gs)
```

```julia
gs = "Quinolinic acid degradation"
ixs = neuroactive[gs]
cs = filter(!isnan, cors[ixs])
acs = filter(!isnan, cors[Not(ixs)])

Resonance.plot_fsea(cs, acs; label=gs)
```

```julia
gs = "S-Adenosylmethionine synthesis"
ixs = neuroactive[gs]
cs = filter(!isnan, cors[ixs])
acs = filter(!isnan, cors[Not(ixs)])

Resonance.plot_fsea(cs, acs; label=gs)
```

```julia
gs = "ClpB"
ixs = neuroactive[gs]
cs = filter(!isnan, cors[ixs])
acs = filter(!isnan, cors[Not(ixs)])

Resonance.plot_fsea(cs, acs; label=gs)
```

