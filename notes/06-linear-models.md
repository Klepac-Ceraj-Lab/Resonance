# Figure 2 - Per-feature, cross-sectional tests on cognitive function scores

First, load packages that will be used throughout this notebook.

```julia
using Resonance
using CairoMakie # for plotting
using GLM
using Statistics
using MultipleTesting
using CategoricalArrays
using ThreadsX
```

Then, we'll load in the different data sources.

```julia
mdata = Resonance.load(Metadata())
mdata.education


let (lq, uq) = [0.25, 0.75]
    mdata.quartile = categorical(map(mdata.cogScorePercentile) do cs
        ismissing(cs) && return missing
        cs <= lq && return "lower"
        cs >= uq && return "upper"
        return "middle"
    end; levels = ["lower", "middle", "upper"], ordered = true)
end

# https://github.com/JuliaStats/GLM.jl/issues/431
mdata.read_depth ./= 1e6
```

```julia

taxa = Resonance.load(TaxonomicProfiles(); timepoint_metadata = mdata)
species = filter(t-> taxrank(t) == :species, taxa)

kos = Resonance.load(KOProfiles(); timepoint_metadata = mdata)
kos = filter(!hastaxon, kos)

```

## Data filtering

Because there is a major shift in microbial composition
upon the introduction of solid foods,
we are going to split the datasets into stools collected prior to 6 months old
(most kids are on liquid diets, breast milk and/or formula)
and over 18 months old (most kids are eating solid foods).

```julia
specmdata = select(DataFrame(Microbiome.metadata(species)),
                ["subject", "timepoint", "ageMonths", "cogScorePercentile", "quartile", "read_depth", "education"]
)

specmdata.sample = samplenames(species)
sort!(specmdata, ["subject", "timepoint"])

specu6 = subset(specmdata, "ageMonths" => ByRow(<(6)), "cogScorePercentile"=> ByRow(!ismissing))
unique!(specu6, "subject")

speco18 = subset(specmdata, "ageMonths" => ByRow(>(18)), "cogScorePercentile"=> ByRow(!ismissing))
unique!(speco18, "subject")

komdata = select(DataFrame(Microbiome.metadata(kos)),
                ["subject", "timepoint", "ageMonths", "cogScorePercentile", "quartile", "read_depth", "education"]
)
komdata.sample = samplenames(kos)
sort!(komdata, ["subject", "timepoint"])

kou6 = subset(komdata, "ageMonths" => ByRow(<(6)), "cogScorePercentile"=> ByRow(!ismissing))
unique!(kou6, "subject")

koo18 = subset(komdata, "ageMonths" => ByRow(>(18)), "cogScorePercentile"=> ByRow(!ismissing))
unique!(koo18, "subject")
```

## Adding features

```julia
for f in features(species)
    colu6 = vec(abundances(species[f, specu6.sample]))
    specu6[!, name(f)] = colu6

    colo18 = vec(abundances(species[f, speco18.sample]))
    speco18[!, name(f)] = colo18
end

for f in features(kos)
    colu6 = vec(abundances(kos[f, kou6.sample]))
    kou6[!, name(f)] = colu6

    colo18 = vec(abundances(kos[f, koo18.sample]))
    koo18[!, name(f)] = colo18
end

```

## Running the models

### Kids under 6mo, species

```julia
specu6_lmresults = DataFrame()

for spc in names(specu6, Not(["subject", "timepoint", "ageMonths", "cogScorePercentile", "quartile", "sample", "read_depth", "education"]))
    count(>(0), specu6[!, spc]) / size(specu6, 1) > 0.1 || continue
    
    @info spc
    over0 = specu6[!, spc] .> 0
    ab = collect(specu6[over0, spc] .+ (minimum(filter(>(0), specu6[!, spc])) / 2)) # add half-minimum non-zerovalue

    df = specu6[over0, ["ageMonths", "cogScorePercentile", "quartile", "read_depth", "education"]]
    df.bug = log2.(ab)

    mod = lm(@formula(bug ~ cogScorePercentile + ageMonths + read_depth + education), df; dropcollinear=false)
    ct = DataFrame(coeftable(mod))
    ct.species .= spc
    ct.kind .= "cogScorePercentile"
    append!(specu6_lmresults, ct)
    subset!(df, :quartile=> ByRow(q-> q in ("upper", "lower")))
    droplevels!(df.quartile)
    length(unique(df.quartile)) > 1 || continue 

    try
        mod = lm(@formula(bug ~ quartile + ageMonths + read_depth + education), df; dropcollinear=false)
    catch e
        @info "$spc failed due to $e"
        continue
    end
    ct = DataFrame(coeftable(mod))
    ct.species .= spc
    ct.kind .= "quartile"
    append!(specu6_lmresults, ct)

end    
select!(specu6_lmresults, Cols(:species, :Name, :))
rename!(specu6_lmresults, "Pr(>|t|)"=>"pvalue");

```



```julia
@chain specu6_lmresults begin
    subset!(:Name => ByRow(x->
        !any(y-> contains(x, y), 
            ("(Intercept)", "ageMonths", "read_depth", "education", r"education")
            )
        )
    )

    groupby(:kind)
    transform!(:pvalue => (col-> adjust(collect(col), BenjaminiHochberg())) => :qvalue)
    sort!(:qvalue)
end

CSV.write(outputfiles("lms_u6mo_species.csv"), specu6_lmresults)
specu6_lmresults[:, Not(["Lower 95%", "Upper 95%"])]
```

### Kids over 18mo, species


```julia
speco18_lmresults = DataFrame()


for spc in names(speco18, Not(["subject", "timepoint", "ageMonths", "cogScorePercentile", "quartile", "sample", "read_depth", "education"]))
    count(>(0), speco18[!, spc]) / size(speco18, 1) > 0.1 || continue
    
    @info spc
    over0 = speco18[!, spc] .> 0
    ab = collect(speco18[over0, spc] .+ (minimum(filter(>(0), speco18[!, spc])) / 2)) # add half-minimum non-zerovalue

    df = speco18[over0, ["ageMonths", "cogScorePercentile", "quartile", "read_depth", "education"]]
    df.bug = log2.(ab)

    mod = lm(@formula(bug ~ cogScorePercentile + ageMonths + read_depth + education), df; dropcollinear=false)
    ct = DataFrame(coeftable(mod))
    ct.species .= spc
    ct.kind .= "cogScorePercentile"
    append!(speco18_lmresults, ct)
    subset!(df, :quartile=> ByRow(q-> q in ("upper", "lower")))
    droplevels!(df.quartile)
    length(unique(df.quartile)) > 1 || continue 

    try
        mod = lm(@formula(bug ~ quartile + ageMonths + read_depth + education), df; dropcollinear=false)
    catch e
        @info "$spc failed due to $e"
        continue
    end
    ct = DataFrame(coeftable(mod))
    ct.species .= spc
    ct.kind .= "quartile"
    append!(speco18_lmresults, ct)

end    
select!(speco18_lmresults, Cols(:species, :Name, :))
rename!(speco18_lmresults, "Pr(>|t|)"=>"pvalue");

```



```julia
@chain speco18_lmresults begin
    subset!(:Name => ByRow(x->
        !any(y-> contains(x, y), 
            ("(Intercept)", "ageMonths", "read_depth", "education", r"education")
            )
        )
    )

    groupby(:kind)
    transform!(:pvalue => (col-> adjust(collect(col), BenjaminiHochberg())) => :qvalue)
    sort!(:qvalue)
end

CSV.write(outputfiles("lms_o18mo_species.csv"), speco18_lmresults)
speco18_lmresults[:, Not(["Lower 95%", "Upper 95%"])]
```


### Kids under 6mo, ko

```julia
kou6_lmresults = DataFrame()

for ko in names(kou6, Not(["subject", "timepoint", "ageMonths", "cogScorePercentile", "quartile", "sample", "read_depth", "education"]))
    count(>(0), kou6[!, ko]) / size(kou6, 1) > 0.1 || continue
    
    @info ko

    ab = kou6[!, ko] .+ (minimum(filter(>(0), kou6[!, ko])) / 2) # add half-minimum non-zerovalue

    df = kou6[:, ["ageMonths", "cogScorePercentile", "quartile", "read_depth", "education"]]
    df.func = log2.(ab)

    mod = lm(@formula(func ~ cogScorePercentile + ageMonths + read_depth + education), df)
    ct = DataFrame(coeftable(mod))
    ct.ko .= ko
    ct.kind .= "cogScorePercentile"
    append!(kou6_lmresults, ct)

    subset!(df, :quartile=> ByRow(q-> q in ("upper", "lower")))
    droplevels!(df.quartile)

    mod = lm(@formula(func ~ quartile + ageMonths + read_depth + education), df)
    ct = DataFrame(coeftable(mod))
    ct.ko .= ko
    ct.kind .= "quartile"
    append!(kou6_lmresults, ct)

end    
select!(kou6_lmresults, Cols(:ko, :Name, :))
rename!(kou6_lmresults, "Pr(>|t|)"=>"pvalue");

```



```julia
@chain kou6_lmresults begin
    subset!(:Name => ByRow(x->
        !any(y-> contains(x, y), 
            ("(Intercept)", "ageMonths", "read_depth", "education", r"education")
            )
        )
    )

    groupby(:kind)
    transform!(:pvalue => (col-> adjust(collect(col), BenjaminiHochberg())) => :qvalue)
    sort!(:qvalue)
end

CSV.write(outputfiles("lms_u6mo_ko.csv"), kou6_lmresults)
kou6_lmresults[:, Not(["Lower 95%", "Upper 95%"])]
```

### Kids over 18mo, ko


```julia
koo18_lmresults = DataFrame()

for ko in names(koo18, Not(["subject", "timepoint", "ageMonths", "cogScorePercentile", "quartile", "sample", "read_depth", "education"]))
    count(>(0), koo18[!, ko]) / size(koo18, 1) > 0.1 || continue
    
    @info ko

    ab = koo18[!, ko] .+ (minimum(filter(>(0), koo18[!, ko])) / 2) # add half-minimum non-zerovalue

    df = koo18[:, ["ageMonths", "cogScorePercentile", "quartile", "read_depth", "education"]]
    df.func = log2.(ab)

    mod = lm(@formula(func ~ cogScorePercentile + ageMonths + read_depth + education), df)
    ct = DataFrame(coeftable(mod))
    ct.ko .= ko
    ct.kind .= "cogScorePercentile"
    append!(koo18_lmresults, ct)

    subset!(df, :quartile=> ByRow(q-> q in ("upper", "lower")))
    droplevels!(df.quartile)

    mod = lm(@formula(func ~ quartile + ageMonths + read_depth + education), df)
    ct = DataFrame(coeftable(mod))
    ct.ko .= ko
    ct.kind .= "quartile"
    append!(koo18_lmresults, ct)

end    
select!(koo18_lmresults, Cols(:ko, :Name, :))
rename!(koo18_lmresults, "Pr(>|t|)"=>"pvalue");

```



```julia
@chain koo18_lmresults begin
    subset!(:Name => ByRow(x->
        !any(y-> contains(x, y), 
            ("(Intercept)", "ageMonths", "read_depth", "education", r"education")
            )
        )
    )

    groupby(:kind)
    transform!(:pvalue => (col-> adjust(collect(col), BenjaminiHochberg())) => :qvalue)
    sort!(:qvalue)
end

CSV.write(outputfiles("lms_o18mo_ko.csv"), koo18_lmresults)
koo18_lmresults[:, Not(["Lower 95%", "Upper 95%"])]
```

## Plots

Then, we will start constructing the figure.
See the [Makie documentation](https://makie.juliaplots.org/stable/tutorials/layout-tutorial/) for more information.

```julia
figure = Figure(; resolution = (1200, 800))

Alo = GridLayout(figure[1,1])

A = Axis(Alo[1,1]; xlabel="Age (months)", ylabel="Cognitive function score")

scatter!(A, mdata.ageMonths, mdata.cogScorePercentile;
    color = map(mdata.omni) do s
        if ismissing(s)
            return (:gray, 0.4)
        elseif s in specu6.sample
            return :purple
        elseif s in speco18.sample
            return :dodgerblue
        else
            return (:seagreen, 0.4)
        end
    end
)

Legend(Alo[1,2], 
    [MarkerElement(; color, marker=:circle) for color in (:purple, :dodgerblue, (:seagreen, 0.4), (:gray, 0.4))],
    ["Under 6mo", "over 18mo", "not included", "no stool"])

figure
```

```julia
using RCall
spc = "Romboutsia_ilealis"

df = speco18[:, ["ageMonths", "cogScorePercentile", "read_depth", "subject", "education"]]

halfmin = minimum(filter(>(0), speco18[:, spc])) / 2
df.bug = collect(log2.((speco18[:, spc] .+ halfmin)))

mod1 = lm(@formula(cogScorePercentile ~ ageMonths + read_depth), df)
mod2 = lm(@formula(cogScorePercentile ~ bug + ageMonths + read_depth), df)

reval("library('lme4')")
@rput df

reval("mod3 <- lm(cogScorePercentile ~ 1 + bug + ageMonths + read_depth, df)")
reval("summary(mod3)")

```

```julia
cc = completecases(speco18, ["ageMonths", "cogScorePercentile", "read_depth", "subject"])

scatter(predict(mod1), df[cc, "cogScorePercentile"])
abline!(0,1)

scatter!(predict(mod2), df[cc, "cogScorePercentile"]; color=:orange)
current_figure()

cor(predict(mod1), df[cc, "cogScorePercentile"])
cor(predict(mod2), df[cc, "cogScorePercentile"])

df2 = copy(df)
df2.bug .= 0

scatter!(predict(mod2, df2[cc, :]), df[cc, "cogScorePercentile"]; color=:purple)
current_figure()
```