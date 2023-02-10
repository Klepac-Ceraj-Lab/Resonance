using Resonance
using CategoricalArrays
using Statistics
using Chain
using DataFrames.PrettyTables

isdir(tablefiles()) || mkpath(tablefiles())

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

fseaout = DataFrame(geneset=unique(fsea.geneset))
for g in levels(fsea.group)
    df = select(subset(fsea, "group"=> ByRow(==(g))), "enrichment", "qvalue")
    rename!(df, ["$(g)_enrichment", "$(g)_qvalue"])
    fseaout = hcat(fseaout, df)
end


CSV.write(tablefiles("Table2.csv"), fseaout)
h1 = Highlighter((data, i, j) -> j in (3,5,7) && data[i,j] < 0.2; bold=true, foreground=:blue)
pretty_table(fseaout; highlighters=(h1,))