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

mdata.maternalEd = categorical([ismissing(m) || m == "-8" ? missing : parse(Int, m) for m in mdata.maternalEd]; ordered=true)
droplevels!(mdata.maternalEd)

let (lq, uq) = quantile(collect(skipmissing(mdata.cogScore)), [0.25, 0.75])
    mdata.quartile = categorical(map(mdata.cogScore) do cs
        ismissing(cs) && return missing
        cs <= lq && return "lower"
        cs >= uq && return "upper"
        return "middle"
    end; levels = ["lower", "middle", "upper"], ordered = true)
end


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
specmdata = select(DataFrame(metadata(species)),
                ["subject", "timepoint", "ageMonths", "cogScore", "quartile", "read_depth", "maternalEd"]
)

specmdata.sample = samplenames(species)
sort!(specmdata, ["subject", "timepoint"])

specu6 = subset(specmdata, "ageMonths" => ByRow(<(6)), "cogScore"=> ByRow(!ismissing))
unique!(specu6, "subject")

speco18 = subset(specmdata, "ageMonths" => ByRow(>(18)), "cogScore"=> ByRow(!ismissing))
unique!(speco18, "subject")


komdata = select(DataFrame(metadata(kos)),
                ["subject", "timepoint", "ageMonths", "cogScore", "quartile", "read_depth", "maternalEd"]
)
komdata.sample = samplenames(kos)
sort!(komdata, ["subject", "timepoint"])

kou6 = subset(komdata, "ageMonths" => ByRow(<(6)), "cogScore"=> ByRow(!ismissing))
unique!(kou6, "subject")

koo18 = subset(komdata, "ageMonths" => ByRow(>(18)), "cogScore"=> ByRow(!ismissing))
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

## MaAsLin

### Export for MaAsLin2

```julia
isdir(outputfiles("maaslin")) || mkdir(outputfiles("maaslin"))

let df = select(mdata, "omni"=>"sample", "ageMonths", "cogScore", "quartile", "read_depth", "maternalEd")
    subset!(df, "sample"=> ByRow(!ismissing))
    unique!(df, "sample")
    CSV.write(outputfiles("maaslin", "metadata.tsv"), df; delim='\t')
end

let df = select(specu6, Cols("sample", Not(["subject", "timepoint", "ageMonths", "cogScore", "quartile", "read_depth", "maternalEd"])))
    CSV.write(outputfiles("maaslin", "specu6.tsv"), df; delim='\t')
end

let df = select(speco18, Cols("sample", Not(["subject", "timepoint", "ageMonths", "cogScore", "quartile", "read_depth", "maternalEd"])))
    CSV.write(outputfiles("maaslin", "speco18.tsv"), df; delim='\t')
end

let df = select(kou6, Cols("sample", Not(["subject", "timepoint", "ageMonths", "cogScore", "quartile", "read_depth", "maternalEd"])))
    CSV.write(outputfiles("maaslin", "kou6.tsv"), df; delim='\t')
end

let df = select(koo18, Cols("sample", Not(["subject", "timepoint", "ageMonths", "cogScore", "quartile", "read_depth", "maternalEd"])))
    CSV.write(outputfiles("maaslin", "koo18.tsv"), df; delim='\t')
end
```

### Run MaAsLin

First, set MaAsLin2 CLI executable in `PATH`


```julia

for (df, label) in zip((specu6, speco18, kou6, koo18), ("specu6", "speco18", "kou6", "koo18"))
    c = Cmd([
        "Maaslin2.R",
        "-f", "ageMonths,quartile,read_depth,maternalEd",
        "--standardize=FALSE",
        outputfiles("maaslin", "$label.tsv"),
        outputfiles("maaslin", "metadata.tsv"),
        outputfiles("maaslin", string(label, "_cogScore")),
    ])

    run(c)
end


for (df, label) in zip((specu6, speco18, kou6, koo18), ("specu6", "speco18", "kou6", "koo18"))
    c = Cmd([
        "Maaslin2.R",
        "-f", "ageMonths,quartile,read_depth,maternalEd",
        "--standardize=FALSE", "-d", "quartile,middle",
        outputfiles("maaslin", "$label.tsv"),
        outputfiles("maaslin", "metadata.tsv"),
        outputfiles("maaslin", string(label, "_quartile")),
    ])

    run(c)
end

```

### Read MaAsLin

```julia
specu6_maaslin = CSV.read(outputfiles("maaslin", "specu6", "all_results.tsv"), DataFrame; delim='\t')
subset!(specu6_maaslin, :metadata=>ByRow(==("cogScore")))
specu6_maaslin.qvalue = adjust(specu6_maaslin.pval, BenjaminiHochberg())

speco18_maaslin = CSV.read(outputfiles("maaslin", "speco18", "all_results.tsv"), DataFrame; delim='\t')
subset!(speco18_maaslin, :metadata=>ByRow(==("cogScore")))
speco18_maaslin.qvalue = adjust(speco18_maaslin.pval, BenjaminiHochberg())

kou6_maaslin = CSV.read(outputfiles("maaslin", "kou6", "all_results.tsv"), DataFrame; delim='\t')
subset!(kou6_maaslin, :metadata=>ByRow(==("cogScore")))
kou6_maaslin.qvalue = adjust(collect(kou6_maaslin.pval), BenjaminiHochberg())

koo18_maaslin = CSV.read(outputfiles("maaslin", "koo18", "all_results.tsv"), DataFrame; delim='\t')
subset!(koo18_maaslin, :metadata=>ByRow(==("cogScore")))
koo18_maaslin.qvalue = adjust(collect(koo18_maaslin.pval), BenjaminiHochberg())
```


## GLM.jl


### Kids under 6mo, species

```julia
specu6_lmresults = DataFrame()

for spc in names(specu6, Not(["subject", "timepoint", "ageMonths", "cogScore", "sample", "read_depth", "maternalEd"]))
    count(>(0), specu6[!, spc]) / size(specu6, 1) > 0.1 || continue
    
    @info spc

    ab = specu6[!, spc] .+ (minimum(filter(>(0), specu6[!, spc])) / 2) # add half-minimum non-zerovalue

    df = specu6[:, ["ageMonths", "cogScore", "read_depth", "maternalEd"]]
    df.bug = log2.(ab)

    mod = lm(@formula(cogScore ~ bug + ageMonths + read_depth + maternalEd), df)
    ct = DataFrame(coeftable(mod))
    ct.species .= spc
    append!(specu6_lmresults, ct)
end    
select!(specu6_lmresults, Cols(:species, :Name, :))
rename!(specu6_lmresults, "Pr(>|t|)"=>"pvalue");

```



```julia
@chain specu6_lmresults begin
    subset!(:Name => ByRow(x->
        !any(y-> contains(x, y), 
            ("(Intercept)", "ageMonths", "read_depth", "maternalEd", r"maternalEd")
            )
        )
    )
    transform!(:pvalue => (col-> adjust(collect(col), BenjaminiHochberg())) => :qvalue)
    sort!(:qvalue)
end

CSV.write(outputfiles("lms_u6mo_species.csv"), specu6_lmresults)
specu6_lmresults
```

### Kids over 18mo, species


```julia
speco18_lmresults = DataFrame()

for spc in names(speco18, Not(["subject", "timepoint", "ageMonths", "cogScore", "sample", "read_depth", "maternalEd"]))
    count(>(0), speco18[!, spc]) / size(speco18, 1) > 0.1 || continue
    @info spc

    ab = speco18[!, spc] .+ (minimum(filter(>(0), speco18[!, spc])) / 2) # add half-minimum non-zerovalue

    df = speco18[:, ["ageMonths", "cogScore", "read_depth", "maternalEd"]]
    df.bug = log2.(ab)

    mod = lm(@formula(cogScore ~ bug + ageMonths + read_depth + maternalEd), df)
    ct = DataFrame(coeftable(mod))
    ct.species .= spc
    append!(speco18_lmresults, ct)
end    
select!(speco18_lmresults, Cols(:species, :Name, :))
rename!(speco18_lmresults, "Pr(>|t|)"=>"pvalue");

```



```julia
@chain speco18_lmresults begin
    subset!(:Name => ByRow(x->
        !any(y-> contains(x, y), 
            ("(Intercept)", "ageMonths", "read_depth", "maternalEd", r"maternalEd")
            )
        )
    )
    transform!(:pvalue => (col-> adjust(collect(col), BenjaminiHochberg())) => :qvalue)
    sort!(:qvalue)
end

CSV.write(outputfiles("lms_o18mo_species.csv"), speco18_lmresults)
speco18_lmresults
```


### Kids under 6mo, KOs

```julia
kou6_lmresults = DataFrame()

for ko in names(kou6, Not(["subject", "timepoint", "ageMonths", "cogScore", "sample", "read_depth", "maternalEd"]))
    count(>(0), kou6[!, ko]) / size(kou6, 1) > 0.1 || continue
    @info ko

    ab = kou6[!, ko] .+ (minimum(filter(>(0), kou6[!, ko])) / 2) # add half-minimum non-zerovalue

    df = kou6[:, ["ageMonths", "cogScore", "read_depth", "maternalEd"]]
    df.bug = log2.(ab)

    mod = lm(@formula(cogScore ~ bug + ageMonths + read_depth + maternalEd), df)
    ct = DataFrame(coeftable(mod))
    ct.ko .= ko
    append!(kou6_lmresults, ct)
end    
select!(kou6_lmresults, Cols(:ko, :Name, :))
rename!(kou6_lmresults, "Pr(>|t|)"=>"pvalue");

```



```julia
@chain kou6_lmresults begin
    subset!(:Name => ByRow(x->
        !any(y-> contains(x, y), 
            ("(Intercept)", "ageMonths", "read_depth", "maternalEd", r"maternalEd")
            )
        )
    )
    transform!(:pvalue => (col-> adjust(collect(col), BenjaminiHochberg())) => :qvalue)
    sort!(:qvalue)
end

CSV.write(outputfiles("lms_u6mo_ko.csv"), kou6_lmresults)
kou6_lmresults
```

```julia
kou6_ecs = DataFrame(filter(!isnothing, ThreadsX.map(subset(kou6_lmresults, :qvalue=> ByRow(<(0.1))).ko) do ko
    m = match(r"\[EC\:(\d+\.\d+\.\d+\.\d+)\]", ko)
    if isnothing(m)
        return nothing
    else
        (l1, l2, l3, l4) = parse.(Int, split(m.captures[1], "."))
        return (; l1, l2, l3, l4)
    end
end))

sort!(kou6_ecs)

```


### Kids over 18mo, KOs


```julia
koo18_lmresults = DataFrame()

for ko in names(koo18, Not(["subject", "timepoint", "ageMonths", "cogScore", "sample", "read_depth", "maternalEd"]))
    count(>(0), koo18[!, ko]) / size(koo18, 1) > 0.1 || continue
    @info ko

    ab = koo18[!, ko] .+ (minimum(filter(>(0), koo18[!, ko])) / 2) # add half-minimum non-zerovalue

    df = koo18[:, ["ageMonths", "cogScore", "read_depth", "maternalEd"]]
    df.bug = log2.(ab)

    mod = lm(@formula(cogScore ~ bug + ageMonths + read_depth + maternalEd), df)
    ct = DataFrame(coeftable(mod))
    ct.ko .= ko
    append!(koo18_lmresults, ct)
end    
select!(koo18_lmresults, Cols(:ko, :Name, :))
rename!(koo18_lmresults, "Pr(>|t|)"=>"pvalue");

```



```julia
@chain koo18_lmresults begin
    subset!(:Name => ByRow(x->
        !any(y-> contains(x, y), 
            ("(Intercept)", "ageMonths", "read_depth", "maternalEd", r"maternalEd")
            )
        )
    )
    transform!(:pvalue => (col-> adjust(collect(col), BenjaminiHochberg())) => :qvalue)
    sort!(:qvalue)
end

CSV.write(outputfiles("lms_o18mo_ko.csv"), koo18_lmresults)
koo18_lmresults
```

## Plots

Then, we will start constructing the figure.
See the [Makie documentation](https://makie.juliaplots.org/stable/tutorials/layout-tutorial/) for more information.

```julia
figure = Figure(; resolution = (1200, 800))

Alo = GridLayout(figure[1,1])

A = Axis(Alo[1,1]; xlabel="Age (months)", ylabel="Cognitive function score")

scatter!(A, mdata.ageMonths, mdata.cogScore;
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
# for sp in subset(speco18_lmresults, :qvalue => ByRow(<(0.1))).species

spc = "Eubacterium_eligens"

df = speco18[:, ["ageMonths", "cogScore", "read_depth", "subject", "maternalEd"]]
@show size(df)  

halfmin = minimum(filter(>(0), speco18[:, spc])) / 2
df.bug = log2.((speco18[:, spc] .+ halfmin))

mod1 = lm(@formula(cogScore ~ ageMonths + read_depth + maternalEd), df)
mod2 = lm(@formula(cogScore ~ bug + ageMonths + read_depth + maternalEd), df)
# end

cc = completecases(speco18, ["ageMonths", "cogScore", "read_depth", "subject", "maternalEd"])

scatter(predict(mod1), df[cc, "cogScore"])
abline!(0,1)

scatter!(predict(mod2), df[cc, "cogScore"]; color=:orange)
current_figure()

cor(predict(mod1), df[cc, "cogScore"])
cor(predict(mod2), df[cc, "cogScore"])

df2 = copy(df)
df2.bug .= 0

scatter!(predict(mod2, df2[cc, :]), df[cc, "cogScore"]; color=:purple)
current_figure()
```