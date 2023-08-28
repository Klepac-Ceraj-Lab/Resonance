using Resonance
using Resonance.PrettyTables
using Dictionaries
using Chain
using CSV
using Statistics
const combine = DataFrames.combine

mdata = Resonance.load(Metadata())
levels!(mdata.race, ["White", "Black", "Indiginous", "Asian", "Mixed", "Other"])

groups = dictionary([
    "N subjects" => [""],
    "samples" => [
        "1",
        "2",
        ">2"
    ],
    "Age (months)" => [
        "min",
        "max",
        "median"
    ],
    "Sex" => [
        "F",
        "M"
    ],
    "Race" => [
        "White",
        "Black",
        "Indiginous",
        "Asian",
        "Mixed",
        "Other",
    ],
    "Maternal Ed." => [
        "Junior high school",
        "Some high school",
        "High school grad",
        "Some college",
        "College grad",
        "Grad/professional school"
    ]
])

Table1 = DataFrame(
    group = reduce(vcat, [fill(k, length(groups[k])) for k in keys(groups)]),
    subgroup = reduce(vcat, [groups[k] for k in keys(groups)])
)

Table1.all = let df = mdata; ufilt = df.filter_00to120
    nsub = [length(unique(df.subject))]
    
    ns = DataFrames.combine(groupby(df, :subject), "timepoint"=> length => "n sample")."n sample"
    ss = [count(==(1), ns), count(==(2), ns), count(>(2), ns)]
    ss = ["$s ($(round(s / length(ns) * 100; digits = 1))%)" for s in ss]

    ages = df.ageMonths[ufilt]
    ages = ["$(round(f(ages); digits=2))" for f in [minimum, maximum, median]]

    sexes = [count(si -> !ismissing(si) && si == s, df.sex[ufilt]) for s in levels(df.sex)]
    sexes = ["$s ($(round(s / sum(ufilt) * 100; digits = 1))%)" for s in sexes]

    races = [count(ri -> !ismissing(ri) && ri == r, df.race[ufilt]) for r in levels(df.race)]
    races = ["$r ($(round(r / sum(ufilt) * 100; digits = 1))%)" for r in races]


    eds = [count(ei -> !ismissing(ei) && ei == e, df.education[ufilt]) for e in levels(df.education)]
    eds = ["$e ($(round(e / sum(ufilt) * 100; digits = 1))%)" for e in eds]
    [nsub; ss; ages; sexes; races; eds]
end

Table1.under6mo = let df = subset(mdata, "ageMonths"=> ByRow(<=(6))); ufilt = df.filter_00to06
    nsub = [length(unique(df.subject))]
    
    ns = DataFrames.combine(groupby(df, :subject), "timepoint"=> length => "n sample")."n sample"
    ss = [count(==(1), ns), count(==(2), ns), count(>(2), ns)]
    ss = ["$s ($(round(s / length(ns) * 100; digits = 1))%)" for s in ss]

    ages = df.ageMonths[ufilt]
    ages = ["$(round(f(ages); digits=2))" for f in [minimum, maximum, median]]

    sexes = [count(si -> !ismissing(si) && si == s, df.sex[ufilt]) for s in levels(df.sex)]
    sexes = ["$s ($(round(s / sum(ufilt) * 100; digits = 1))%)" for s in sexes]

    races = [count(ri -> !ismissing(ri) && ri == r, df.race[ufilt]) for r in levels(df.race)]
    races = ["$r ($(round(r / sum(ufilt) * 100; digits = 1))%)" for r in races]


    eds = [count(ei -> !ismissing(ei) && ei == e, df.education[ufilt]) for e in levels(df.education)]
    eds = ["$e ($(round(e / sum(ufilt) * 100; digits = 1))%)" for e in eds]
    [nsub; ss; ages; sexes; races; eds]
end

Table1.over18mo = let df = subset(mdata, "ageMonths"=> ByRow(>(18))); ufilt = df.filter_18to120
    nsub = [length(unique(df.subject))]
    
    ns = DataFrames.combine(groupby(df, :subject), "timepoint"=> length => "n sample")."n sample"
    ss = [count(==(1), ns), count(==(2), ns), count(>(2), ns)]
    ss = ["$s ($(round(s / length(ns) * 100; digits = 1))%)" for s in ss]

    ages = df.ageMonths[ufilt]
    ages = ["$(round(f(ages); digits=2))" for f in [minimum, maximum, median]]

    sexes = [count(si -> !ismissing(si) && si == s, df.sex[ufilt]) for s in levels(df.sex)]
    sexes = ["$s ($(round(s / sum(ufilt) * 100; digits = 1))%)" for s in sexes]

    races = [count(ri -> !ismissing(ri) && ri == r, df.race[ufilt]) for r in levels(df.race)]
    races = ["$r ($(round(r / sum(ufilt) * 100; digits = 1))%)" for r in races]


    eds = [count(ei -> !ismissing(ei) && ei == e, df.education[ufilt]) for e in levels(df.education)]
    eds = ["$e ($(round(e / sum(ufilt) * 100; digits = 1))%)" for e in eds]
    [nsub; ss; ages; sexes; races; eds]
end

CSV.write(tablefiles("Table1.csv"), Table1)

pretty_table(Table1; backend = Val(:latex),
    label = "Table 1",
    nosubheader=true
)