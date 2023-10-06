using Resonance
using CategoricalArrays
using Statistics
using Chain
using DataFrames.PrettyTables

isdir(tablefiles()) || mkpath(tablefiles())

fsea = let
    fsea_all = CSV.read(tablefiles("figure2", "fsea_all_00to120.csv"), DataFrame)
    fsea_all.group .= "00to120"
    fsea_u6 = CSV.read(tablefiles("figure2", "fsea_all_00to06.csv"), DataFrame)
    fsea_u6.group .= "00to06"
    fsea_o18 = CSV.read(tablefiles("figure2", "fsea_all_18to120.csv"), DataFrame)
    fsea_o18.group .= "18to120"
    vcat(fsea_all, fsea_u6, fsea_o18)
end

fsea.group = categorical(fsea.group; levels=["00to120", "00to06", "18to120"])

@chain fsea begin
    subset!(:cortest=> ByRow(==("cogScore")))
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
    df = select(leftjoin(select(fseaout, "geneset"), subset(fsea, "group"=> ByRow(==(g))); on="geneset"), "enrichment", "qvalue")
    rename!(df, ["$(g)_enrichment", "$(g)_qvalue"])
    fseaout = hcat(fseaout, df)
end
for col in names(fseaout, Not("geneset"))
    fseaout[!, col] = coalesce.(fseaout[!, col], NaN)
end

CSV.write(tablefiles("Table2.csv"), fseaout)
formatter = (v, i, j) -> v isa AbstractFloat ? round(v, digits = 3) : v
h1 = LatexHighlighter((data, i, j) -> !ismissing(j) && j in (3,5,7) && data[i,j] < 0.2 && data[i,j-1] < 0, ["color{blue}","textbf"])
h2 = LatexHighlighter((data, i, j) -> !ismissing(j) && j in (3,5,7) && data[i,j] < 0.2 && data[i,j-1] > 0, ["color{red}","textbf"])

pretty_table(fseaout; highlighters=(h1,h2), formatters=ft_round(3),
    header = (map(n-> replace(replace(n, r"_.+" => ""), "geneset"=>"Age group"), names(fseaout)),
              map(n-> replace(replace(n, r".+_enrichment" => "E.S."), r".+_qvalue"=>"q value"), names(fseaout))
    ), backend=Val(:latex)
    )