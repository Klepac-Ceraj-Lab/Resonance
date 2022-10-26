# Tables to be included in the manuscript

## Overview

Using `PrettyTables.jl` to spit out LaTeX tables
or `DataFrames.jl` to generate CSVs.

```julia
using Resonance
using CategoricalArrays
using Statistics
using Chain
using DataFrames.PrettyTables

isdir(tablefiles()) || mkpath(tablefiles())
```

## Table 1 - Cohort overview

```julia
mdata = Resonance.load(Metadata())
taxa = Resonance.load(TaxonomicProfiles())
subset!(mdata,
    :omni=> ByRow(o-> !ismissing(o) && o in samplenames(taxa)),
    :cogScore=> ByRow(!ismissing),
    :ageMonths=> ByRow(<(120)),
    :education => ByRow(!ismissing)
)

mdata.ageGroup = map(mdata.ageMonths) do a
    ismissing(a) && return missing
    a <=6 && return "under 6mo"
    a > 18 && return "over 8mo"
    return "interim"
end

demographics = DataFrame(
    group = [["N subjects"]; fill("samples", 3); fill("Age", 3); fill("Sex", 2); fill("Race", 6); fill("SES", 6) ],
    subgroup = ["", # N subjects
                "1", "2", ">2",
                "min", "max", "median", # Age
                "F", "M", # Sex
                "White", "Black", "Asian", "Mixed", "Other", "unknown/declined", # Race
                levels(mdata.education)... # SES
    ]
)

function count_perc(pred, vec)
    c = count(pred, vec)
    l = length(vec)
    return "$c ($(round(c / l * 100; digits=1))%)"
end

let gdf = groupby(mdata, :subject)
    demographics.all = Any[
        keys(gdf) |> length,
        let scounts = combine(gdf, :omni=> (o-> count(!ismissing, o)) => :samples).samples
            (count_perc(==(1), scounts), count_perc(==(2), scounts), count_perc(>(2), scounts))
        end...,
        extrema(combine(gdf, :ageMonths=>identity=> :ageMonths).ageMonths)...,
        median(combine(gdf, :ageMonths=>identity=> :ageMonths).ageMonths),
        let sex = combine(gdf, :sex=> first => :sex).sex
            (count_perc(==(s), sex) for s in levels(sex))
        end...,
        
        let race = combine(gdf, :race=> first => :race).race
            (count_perc(==(r), race) for r in ["White", "Black", "Asian", "Mixed", "Other", "Unknown"])
        end...,
        let education = combine(gdf, :education=> first => :education).education
            (count_perc(x-> !ismissing(x) && x == e, education) for e in levels(education))
        end...,
    ]
end

let gdf = groupby(subset(mdata, :ageMonths => ByRow(<(6))), :subject)
    demographics.under6 = Any[
        keys(gdf) |> length,
        let scounts = combine(gdf, :omni=> (o-> count(!ismissing, o)) => :samples).samples
            (count_perc(==(1), scounts), count_perc(==(2), scounts), count_perc(>(2), scounts))
        end...,
        extrema(combine(gdf, :ageMonths=>identity=> :ageMonths).ageMonths)...,
        median(combine(gdf, :ageMonths=>identity=> :ageMonths).ageMonths),
        let sex = combine(gdf, :sex=> first => :sex).sex
            (count_perc(==(s), sex) for s in levels(sex))
        end...,
        
        let race = combine(gdf, :race=> first => :race).race
            (count_perc(==(r), race) for r in ["White", "Black", "Asian", "Mixed", "Other", "Unknown"])
        end...,
        let education = combine(gdf, :education=> first => :education).education
            (count_perc(x-> !ismissing(x) && x == e, education) for e in levels(education))
        end...,
    ]
end

let gdf = groupby(subset(mdata, :ageMonths => ByRow(>(18))), :subject)
    demographics.over18 = Any[
        keys(gdf) |> length,
        let scounts = combine(gdf, :omni=> (o-> count(!ismissing, o)) => :samples).samples
            (count_perc(==(1), scounts), count_perc(==(2), scounts), count_perc(>(2), scounts))
        end...,
        extrema(combine(gdf, :ageMonths=>identity=> :ageMonths).ageMonths)...,
        median(combine(gdf, :ageMonths=>identity=> :ageMonths).ageMonths),
        let sex = combine(gdf, :sex=> first => :sex).sex
            (count_perc(==(s), sex) for s in levels(sex))
        end...,
        
        let race = combine(gdf, :race=> first => :race).race
            (count_perc(==(r), race) for r in ["White", "Black", "Asian", "Mixed", "Other", "Unknown"])
        end...,
        let education = combine(gdf, :education=> first => :education).education
            (count_perc(x-> !ismissing(x) && x == e, education) for e in levels(education))
        end...,
    ]
end

CSV.write(tablefiles("Table1.csv"), demographics)
pretty_table(demographics)
```


## Table 2 - FSEA results

```julia
fsea = let
    fsea_all = CSV.read(outputfiles("fsea_all.csv"), DataFrame)
    fsea_all.group .= "all"
    fsea_u6 = CSV.read(outputfiles("fsea_u6.csv"), DataFrame)
    fsea_u6.group .= "u6"
    fsea_o18 = CSV.read(outputfiles("fsea_o18.csv"), DataFrame)
    fsea_o18.group .= "o18"
    vcat(fsea_all, fsea_u6, fsea_o18)
end

fsea.group = categorical(fsea.group; levels=["all", "u6", "o18"])

@chain fsea begin
    subset!(:cortest=> ByRow(==("cogScorePercentile")))
    groupby(:geneset)
    transform!(:qvalue => (
        q-> any(<(0.2), q) ? fill(true, length(q)) : 
                             fill(false, length(q))) => :keep
    )

    subset!(:keep => identity)
    
    transform!(:geneset=> ByRow(
        g-> replace(g, r" \(.+?\)" => "", 
                        "synthesis"=>"syn.",
                        "degradation"=>"deg.")
        ) => :geneset
    )
    select!(["group", "geneset", "enrichment", "qvalue"])
    sort!(["geneset", "group"])
end

CSV.write(tablefiles("Table2.csv"), fsea)
h1 = Highlighter((data, i, j) -> data[i,4] < 0.2; bold=true, foreground=:blue)
pretty_table(fsea; highlighters=(h1,))
```