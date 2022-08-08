# Figure 2 - Per-feature, cross-sectional tests on cognitive function scores

First, load packages that will be used throughout this notebook.

```julia
using Resonance
using CairoMakie # for plotting
using GLM
using Statistics
using MultipleTesting
using CategoricalArrays
```

Then, we'll load in the different data sources.

```julia
mdata = Resonance.load(Metadata())

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
                ["subject", "timepoint", "ageMonths", "cogScore", "read_depth", "maternalEd"]
)

specmdata.maternalEd = categorical([x == "-8" ? missing : x for x in specmdata.maternalEd]; ordered = true)

specmdata.sample = samplenames(species)
sort!(specmdata, ["subject", "timepoint"])

specu6 = subset(specmdata, "ageMonths" => ByRow(<(6)), "cogScore"=> ByRow(!ismissing))
unique!(specu6, "subject")

speco18 = subset(specmdata, "ageMonths" => ByRow(>(18)), "cogScore"=> ByRow(!ismissing))
unique!(speco18, "subject")


komdata = select(DataFrame(metadata(kos)),
                ["subject", "timepoint", "ageMonths", "cogScore", "read_depth"]
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
    count(>(0), colu6) / length(colu6) > 0.1 && (specu6[!, name(f)] = colu6)

    colo18 = vec(abundances(species[f, speco18.sample]))
    count(>(0), colo18) / length(colo18) > 0.1 && (speco18[!, name(f)] = colo18)
end

for f in features(kos)
    colu6 = vec(abundances(kos[f, kou6.sample]))
    count(>(0), colu6) / length(colu6) > 0.1 && (kou6[!, name(f)] = colu6)

    colo18 = vec(abundances(kos[f, koo18.sample]))
    count(>(0), colo18) / length(colo18) > 0.1 && (koo18[!, name(f)] = colo18)
end

```

## GLMs

### Export for MaAsLin2

```julia
isdir(outputfiles("maaslin")) || mkdir(outputfiles("maaslin"))

let df = select(speco18, Cols("sample", Not(["subject", "timepoint", "ageMonths", "cogScore", "read_depth", "maternalEd"])))
    CSV.write(outputfiles("maaslin", "bugs.tsv"), df; delim='\t')
end

let df = select(speco18, "sample", "ageMonths", "cogScore", "read_depth", "maternalEd")
    CSV.write(outputfiles("maaslin", "metadata.tsv"), df; delim='\t')
end

```

### Kids under 6mo, species

```julia
specu6_lmresults = DataFrame()

for spc in names(specu6, Not(["subject", "timepoint", "ageMonths", "cogScore", "sample", "read_depth", "maternalEd"]))
    @info spc

    ab = specu6[!, spc]
    over0 = ab .> 0

    df = specu6[over0, ["ageMonths", "cogScore", "read_depth", "maternalEd"]]
    df.bug = log.(specu6[over0, spc] ./ 100)

    mod = lm(@formula(bug ~ ageMonths + cogScore + read_depth + maternalEd), df)
    ct = DataFrame(coeftable(mod))
    ct.species .= spc
    append!(specu6_lmresults, ct)
end    
select!(specu6_lmresults, Cols(:species, :Name, :))
rename!(specu6_lmresults, "Pr(>|t|)"=>"pvalue");

```



```julia
@chain specu6_lmresults begin
    subset!(:Name=>ByRow(x-> !in(x, ("(Intercept)", "ageMonths", "read_depth"))))
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
    @info spc

    over0 = speco18[!, spc] .> 0

    df = speco18[over0, ["ageMonths", "cogScore", "read_depth", "subject", "maternalEd"]]
    df.bug = log.(speco18[over0, spc] ./ 100)

    mod = lm(@formula(cogScore ~ bug  + ageMonths + read_depth + maternalEd), df)

    ct = DataFrame(coeftable(mod))
    ct.species .= spc
    append!(speco18_lmresults, ct)
end    

select!(speco18_lmresults, Cols(:species, :Name, :))
rename!(speco18_lmresults, "Pr(>|t|)"=>"pvalue");

```



```julia
@chain speco18_lmresults begin
    subset!(:Name => ByRow(x-> !any(y-> contains(x, y), ("(Intercept)", "ageMonths", "read_depth", "maternalEd", r"maternalEd"))))
    transform!(:pvalue => (col-> adjust(collect(col), BenjaminiHochberg())) => :qvalue)
    sort!(:qvalue)
end

CSV.write(outputfiles("lms_o18mo_species.csv"), speco18_lmresults)
speco18_lmresults
```


### Kids under 6mo, KOs

```julia
kou6_lmresults = DataFrame()

for ko in names(kou6, Not(["subject", "timepoint", "ageMonths", "cogScore", "sample", "read_depth"]))
    @info ko

    ab = kou6[!, ko]
    over0 = ab .> 0

    df = kou6[over0, ["ageMonths", "cogScore", "read_depth"]]
    df.bug = log.(kou6[over0, ko] ./ 100)

    mod = lm(@formula(bug ~ ageMonths + cogScore + read_depth), df)
    ct = DataFrame(coeftable(mod))
    ct.ko .= ko
    append!(kou6_lmresults, ct)
end    
select!(kou6_lmresults, Cols(:ko, :Name, :))
rename!(kou6_lmresults, "Pr(>|t|)"=>"pvalue");

```



```julia
@chain kou6_lmresults begin
    subset!(:Name=>ByRow(x-> !in(x, ("(Intercept)", "ageMonths", "read_depth"))))
    transform!(:pvalue => (col-> adjust(collect(col), BenjaminiHochberg())) => :qvalue)
    sort!(:qvalue)
end

CSV.write(outputfiles("lms_u6mo_ko.csv"), kou6_lmresults)
kou6_lmresults
```

### Kids over 18mo, KOs


```julia
koo18_lmresults = DataFrame()

for ko in names(koo18, Not(["subject", "timepoint", "ageMonths", "cogScore", "sample", "read_depth"]))
    @info ko

    ab = koo18[!, ko]
    over0 = ab .> 0

    df = koo18[over0, ["ageMonths", "cogScore", "read_depth", "subject"]]
    df.bug = log.(koo18[over0, ko] ./ 100)

    mod = lm(@formula(bug ~ ageMonths + cogScore + read_depth), df)
    ct = DataFrame(coeftable(mod))
    ct.ko .= ko
    append!(koo18_lmresults, ct)
end    

select!(koo18_lmresults, Cols(:ko, :Name, :))
rename!(koo18_lmresults, "Pr(>|t|)"=>"pvalue");

```



```julia
@chain koo18_lmresults begin
    subset!(:Name => ByRow(x-> !in(x, ("(Intercept)", "ageMonths", "read_depth"))))
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

spc = "Faecalibacterium_prausnitzii"

over0 = findall(row-> !ismissing(row.maternalEd) && row[spc] > 0, eachrow(speco18))

df = speco18[over0, ["ageMonths", "cogScore", "read_depth", "subject", "maternalEd"]]
@show size(df)
df.bug = log.(speco18[over0, spc] ./ 100)

mod1 = lm(@formula(cogScore ~ ageMonths + read_depth + maternalEd), df)
mod2 = lm(@formula(cogScore ~ bug + ageMonths + read_depth + maternalEd), df)
# end

scatter(predict(mod1), speco18[over0, "cogScore"])
lines!([40, 140], [40, 140])

scatter!(predict(mod2), speco18[over0, "cogScore"]; color=:orange)
current_figure()

cor(predict(mod1), speco18[over0, "cogScore"])
cor(predict(mod2), speco18[over0, "cogScore"])

df2 = copy(df)
df2.bug .= 0

scatter!(predict(mod2, df2), speco18[over0, "cogScore"]; color=:purple)
current_figure()
```