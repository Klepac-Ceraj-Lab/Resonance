# Figure3 - functional analysis

```julia
using Resonance
using CairoMakie
using Statistics
using HypothesisTests
using MultipleTesting
using KernelDensity
using MultivariateStats
using ThreadsX
```

## Data Loading

```julia
mdata = Resonance.load(Metadata())

unirefs = Resonance.load(UnirefProfiles(); timepoint_metadata = mdata) # this can take a bit
unirefs = filter(!hastaxon, unirefs) # don't use species stratification for summaries
isdefined(Main, :unidm) || (unidm = braycurtis(unirefs))
unipco = fit(MDS, unidm; distances=true)

metabolites = Resonance.load(MetabolicProfiles(); timepoint_metadata=mdata)
metabolites = metabolites[:, [!ismissing(a) && a < 14 for a in get(metabolites, :ageMonths)]]
isdefined(Main, :metdm) || (metdm = braycurtis(metabolites))
metpco = fit(MDS, metdm; distances=true)
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
fsdf = let 
    if isfile(outputfiles("fsea_consolidated.csv"))
        tmp = CSV.read(outputfiles("fsea_consolidated.csv"), DataFrame)
    else
        tmp = DataFrame(ThreadsX.map(collect(keys(neuroactive))) do gs
            ixs = neuroactive[gs]
            isempty(ixs) && return (; geneset = gs, U = NaN, median = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

            cs = filter(!isnan, cors[ixs])
            isempty(cs) && return (; geneset = gs, U = NaN, median = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

            acs = filter(!isnan, cors[Not(ixs)])
            mwu = MannWhitneyUTest(cs, acs)

            return (; geneset = gs, U = mwu.U, median = mwu.median, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
        end)

        subset!(tmp, :pvalue=> ByRow(!isnan))
        tmp.qvalue = adjust(tmp.pvalue, BenjaminiHochberg())
        sort!(tmp, :qvalue)
        CSV.write(outputfiles("fsea_consolidated.csv"), tmp)
    end
    tmp
end

```

```julia
fsdf2 = let
    if isfile(outputfiles("fsea_all.csv"))
        tmp = CSV.read(outputfiles("fsea_all.csv"), DataFrame)
    else
        tmp = DataFrame(
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

        subset!(tmp, :pvalue=> ByRow(!isnan))
        tmp.qvalue = adjust(tmp.pvalue, BenjaminiHochberg())
        sort!(tmp, :qvalue)
        CSV.write(outputfiles("fsea_all.csv"), tmp)
    end
    tmp
end
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
figure = Figure(resolution=(900, 900))

A = GridLayout(figure[1,1])
B = GridLayout(figure[2,1])
CDEF = GridLayout(figure[1:2,2])
C = GridLayout(CDEF[1,1])
D = GridLayout(CDEF[2,1])
E = GridLayout(CDEF[3,1])
F = Axis(CDEF[4,1])
G = GridLayout(figure[3,1:2])
```


```julia
aax = Axis(A[1,1])
pco = plot_pcoa!(aax, unipco; color=get(unirefs, :ageMonths))
Colorbar(A[1,2], pco; label="Age (months)")
ax, hm = heatmap(B[1,1],k.x,k.y, k.density .^ (1/4); colormap=:plasma, axis = (;xlabel = "prevalence", ylabel="correlation"))
Colorbar(B[1,2], hm; label=L"$Density^{\frac{1}{4}}$")

for (gs, panel) in zip(("Propionate degradation I", "Glutamate degradation I", "GABA synthesis I"), (C,D,E))
    ixs = neuroactive_full[gs]
    cs = filter(!isnan, cors[ixs])
    acs = filter(!isnan, cors[Not(ixs)])

    Resonance.plot_fsea!(panel, cs, acs; label=replace(gs, "degradation"=> "degr.", "synthesis"=> "synth."))
end

Resonance.plot_corrband!(F, cors)

rowsize!(CDEF, 4, Relative(1/8))
```

## Metabolites

```julia
bmoi = [ # metabolites of interest
    "pyridoxine", # vit B6
    "gaba",
    "glutamate"
]
mfeats = commonname.(features(metabolites))
pyrdidx = only(ThreadsX.findall(f-> !ismissing(f) && contains(f, r"pyridoxine"i), mfeats))
gabaidx = only(ThreadsX.findall(f-> !ismissing(f) && contains(f, r"gamma-aminobutyric"i), mfeats))
glutidx = ThreadsX.findall(f-> !ismissing(f) && contains(f, r"^glutamic"i), mfeats)

pyrd = vec(abundances(metabolites[pyrdidx, :]))
gaba = vec(abundances(metabolites[gabaidx, :]))
glut1 = vec(abundances(metabolites[glutidx[1], :]))
glut2 = vec(abundances(metabolites[glutidx[2], :]))

```

```julia

Ga = Axis(G[1,1])
Gb = Axis(G[1,2]; xlabel = L"log_e(GABA)", ylabel = L"log_e(Glutamate)")

pcom = plot_pcoa!(Ga, metpco; color=get(metabolites, :ageMonths))
sc = scatter!(Gb, log.(gaba), log.((glut1 .+ glut2) ./ 2); color = get(metabolites, :ageMonths))
Colorbar(G[1, 3], sc; label = "Age (Months)")

save(figurefiles("Figure3.svg"), figure)
save(figurefiles("Figure3.png"), figure)
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
