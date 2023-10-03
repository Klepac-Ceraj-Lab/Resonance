# Figure2 - functional analysis

```julia
using Resonance
using CairoMakie
using Statistics
using HypothesisTests
using MultipleTesting
using MultivariateStats
using Distributions
using CategoricalArrays
using ThreadsX
using ColorSchemes
using GLM
using Dates
using FileIO
using JLD2
```

## Data Loading

```julia
mdata = Resonance.load(Metadata())
seqs = subset(mdata, "sample"=> ByRow(!ismissing)) 
seqs.sample = [s for s in seqs.sample]
taxa = Resonance.load(TaxonomicProfiles(); timepoint_metadata = seqs)
taxdf = comm2wide(taxa)
```

## Load FSEA


```julia
speclms_00to120 = CSV.read(tablefiles("figure2", "lms_species_00to120.csv"), DataFrame)
speclms_00to06 = CSV.read(tablefiles("figure2", "lms_species_00to06.csv"), DataFrame)
speclms_18to120 = CSV.read(tablefiles("figure2", "lms_species_18to120.csv"), DataFrame)
subscale_lms = CSV.read(tablefiles("newfigures", "subscale_lms.csv"), DataFrame)
subscale_lms.kind[findall(k-> contains(k, "Expressivelan"), subscale_lms.kind)] .= "Mullen::mullen_ExpressiveLanguageT"


speclms_pa = CSV.read(tablefiles("figure2", "lms_species_18to120_pa.csv"), DataFrame)
fsdf_00to120 = CSV.read(tablefiles("figure2", "fsea_consolidated_00to120.csv"), DataFrame)
fsdf2_00to120 = CSV.read(tablefiles("figure2", "fsea_all_00to120.csv"), DataFrame)
fsdf_00to06 = CSV.read(tablefiles("figure2", "fsea_consolidated_00to06.csv"), DataFrame)
fsdf2_00to06 = CSV.read(tablefiles("figure2", "fsea_all_00to06.csv"), DataFrame)
fsdf_18to120 = CSV.read(tablefiles("figure2", "fsea_consolidated_18to120.csv"), DataFrame)
fsdf2_18to120 = CSV.read(tablefiles("figure2", "fsea_all_18to120.csv"), DataFrame)

cors_00to120 = CSV.read(tablefiles("figure2", "lms_unirefs_00to120.csv"), DataFrame)
cors_00to06 = CSV.read(tablefiles("figure2", "lms_unirefs_00to06.csv"), DataFrame)
cors_18to120 = CSV.read(tablefiles("figure2", "lms_unirefs_18to120.csv"), DataFrame)
neuroactive_00to120 = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), cors_00to120.feature); consolidate=false)
neuroactive_00to06 = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), cors_00to06.feature); consolidate=false)
neuroactive_18to120 = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), cors_18to120.feature); consolidate=false)
```


```julia
lms_mat = let 
    df = vcat(speclms_00to120, speclms_00to06, speclms_18to120)
    df.group = [fill("00to120", nrow(speclms_00to120)); fill("00to06", nrow(speclms_00to06)); fill("18to120", nrow(speclms_18to120))]
    sigs = subset(df, "qvalue"=> ByRow(<(0.2)))
    sort!(sigs, :qvalue)
    sig_feats = unique(sigs.feature)
    gdf = groupby(df, ["feature", "group"])
    DataFrame(ThreadsX.map(eachindex(sig_feats)) do i
        ft = sig_feats[i]
        q_00to120 = haskey(gdf, (; feature=ft, group="00to120")) ? only(gdf[(; feature=ft, group="00to120")].qvalue) : 1.0
        q_00to06 = haskey(gdf, (; feature=ft, group="00to06")) ? only(gdf[(; feature=ft, group="00to06")].qvalue) : 1.0
        q_18to120 = haskey(gdf, (; feature=ft, group="18to120")) ? only(gdf[(; feature=ft, group="18to120")].qvalue) : 1.0
        corr_00to120 = cor(taxdf[taxdf.filter_00to120, "cogScore"], taxdf[taxdf.filter_00to120, ft])
        corr_00to06 = cor(taxdf[taxdf.filter_00to06, "cogScore"], taxdf[taxdf.filter_00to06, ft])
        corr_18to120 = cor(taxdf[taxdf.filter_18to120, "cogScore"], taxdf[taxdf.filter_18to120, ft])
        (; feature = ft, 
           prev_00to120 = only(prevalence(taxa[Regex(ft), seqs.filter_00to120])),
           meanab_00to120 = mean(filter(>(0), abundances(taxa[Regex(ft), seqs.filter_00to120]))),
           corr_00to120,
           q_00to120,
           prev_00to06 = only(prevalence(taxa[Regex(ft), seqs.filter_00to06])),
           meanab_00to06 = mean(filter(>(0), abundances(taxa[Regex(ft), seqs.filter_00to06]))),
           corr_00to06,
           q_00to06, 
           prev_18to120 = only(prevalence(taxa[Regex(ft), seqs.filter_18to120])),
           meanab_18to120 = mean(filter(>(0), abundances(taxa[Regex(ft), seqs.filter_18to120]))),
           corr_18to120,
           q_18to120,
        )
   end)
end
```

```julia
subscale_mat = let
    df = @chain subscale_lms begin
           sort("species")
           groupby("species")
           subset("qvalue"=> (q-> any(<(0.2), q)))
           select("filter", "kind", "species", "qvalue", "t")
           groupby(["kind", "filter"])
    end
    newdf = DataFrame(species = sort(unique(DataFrames.combine(df, "species"=> identity=>"species").species)))
    for key in keys(df)
        col = string(replace(key[:kind], "Mullen::mullen_"=>"", "languageT"=> "LanguageT"), "_", key[:filter])
        leftjoin!(newdf, select(df[key], "species", "qvalue"=> "$(col)_q", "t"=> "$(col)_t"); on = "species")
    end
    for col in names(newdf, r"_q$")
        newdf[!, col] = coalesce.(newdf[!, col], 1.0)
    end
    for col in names(newdf, r"_t$")
        newdf[!, col] = coalesce.(newdf[!, col], 0.0)
    end
    newdf 
end

```

```julia
fig = Figure(;resolution = (1200, 500));
ax1 = Axis(fig[1,1]; xlabel = "t score", ylabel = "-log(q)", title="0 to 6m")
ax2 = Axis(fig[1,2]; xlabel = "t score", ylabel = "-log(q)", title="18 to 120m")
ax3 = Axis(fig[1,3]; xlabel = "t score", ylabel = "-log(q)", title="0 to 120m")

for (i, gdf) in enumerate(groupby(subset(subscale_lms, "filter"=> ByRow(==("00to06"))), "kind"))
    scheme = ColorSchemes.Paired_10
    c = [q > 0.2 ? (scheme[i], 0.2) : (scheme[i+1], 1) for q in gdf.qvalue]
    scatter!(ax1, gdf.t, -1 .* log.(gdf.qvalue); color = c)
end

for (i, gdf) in enumerate(groupby(subset(subscale_lms, "filter"=> ByRow(==("18to120"))), "kind"))
    scheme = ColorSchemes.Paired_10
    c = [q > 0.2 ? (scheme[i], 0.2) : (scheme[i+1], 1) for q in gdf.qvalue]
    scatter!(ax2, gdf.t, -1 .* log.(gdf.qvalue); color = c)
end

for (i, gdf) in enumerate(groupby(subset(subscale_lms, "filter"=> ByRow(==("00to120"))), "kind"))
    scheme = ColorSchemes.Paired_10
    c = [q > 0.2 ? (scheme[i], 0.2) : (scheme[i+1], 1) for q in gdf.qvalue]
    scatter!(ax3, gdf.t, -1 .* log.(gdf.qvalue); color = c)
end

Legend(fig[2, 1:3], [MarkerElement(; marker=:circle, color=ColorSchemes.Paired_10[i*2]) for i in 1:length(unique(subscale_lms.kind))],
                    [replace(k, "Mullen::mullen_"=>"") for k in unique(subscale_lms.kind)];
                    orientation=:horizontal,
                    nbanks = 2)

fig
```
  

## Plotting

```julia
figure = Figure(resolution=(1200, 1000));

AB = GridLayout(figure[1,1:3])
A = GridLayout(AB[1,1:2])
B = GridLayout(AB[1,3])
C = GridLayout(figure[2,1:2]; alignmode=Mixed(; left=0, bottom=0))
DE = GridLayout(figure[2,3])
D = GridLayout(DE[1,1])
E = GridLayout(DE[2,1])
F = GridLayout(figure[3,1:3]; alignmode=Mixed(; left=0))

#-

aax1 = Axis(A[1,1]; title = "over 18m", ylabel=L"$-log_2(P)$", xlabel = "coef.",
                xticks = ([-1.5e-3, 0.0, 1.5e-3], ["-1.5e-3", "0.0", "1.5e-3"]),
                limits = ((-1.5e-3, 1.5e-3), nothing),
                xminorticksvisible=true, xminorticks = IntervalsBetween(3))
aax2 = Axis(A[1,2]; title = "all ages", ylabel=L"$-log_2(P)$", xlabel = "coef.",
                xticks = ([-1.5e-3, 0.0, 1.5e-3], ["-1.5e-3", "0.0", "1.5e-3"]),
                limits = ((-1.5e-3, 1.5e-3), nothing),
                xminorticksvisible=true, xminorticks = IntervalsBetween(3))

scatter!(aax1, speclms_18to120.coef, -1 .* log2.(speclms_18to120.qvalue);
    color = map(q-> q < 0.2 ? ColorSchemes.tableau_10[3] : ColorSchemes.tableau_10[10], speclms_18to120.qvalue))
scatter!(aax2, speclms_00to120.coef, -1 .* log2.(speclms_00to120.qvalue);
    color = map(q-> q < 0.2 ? ColorSchemes.tableau_10[3] : ColorSchemes.tableau_10[10], speclms_00to120.qvalue))

#-

bax1 = Axis(B[1, 1];
    title = "over 18m",
    yticks = (1.5:nrow(lms_mat) + 0.5, format_species_labels(lms_mat.feature)),
    yticklabelfont = "TeX Gyre Heros Makie Italic",
    yticklabelsize = 14,
    xticklabelsize = 12
)

bax1.xticklabelsvisible = false
bax1.xticksvisible = false

bax2 = Axis(B[1, 2]; title = "all ages")

bax2.xticklabelsvisible = false
bax2.xticksvisible = false
bax2.yticklabelsvisible = false
bax2.yticksvisible = false

tightlimits!.((bax1, bax2))

let clrs = [isnan(x) ? 0.0 : x for x in lms_mat[:, "corr_18to120"]]
    
    poly!(bax1, [Rect(1, i, 1, 1) for i in eachindex(lms_mat.feature)]; 
        color = clrs, colorrange=(-0.3, 0.3),
        colormap = Reverse(:RdBu))
    
    
    qs = lms_mat[:, "q_18to120"]
    annotations!(bax1, map(x-> x < 0.05 ? "**" : x < 0.2 ? "*" : "", qs), [Point2f(1.5, i + 0.2) for i in eachindex(lms_mat.feature)];
    color = [abs(x) > 0.2 ? :white : :black for x in clrs],
    align=(:center, :center))
end
let clrs = [isnan(x) ? 0.0 : x for x in lms_mat[:, "corr_00to120"]]
    
    poly!(bax2, [Rect(1, i, 1, 1) for i in eachindex(lms_mat.feature)]; 
        color = clrs, colorrange=(-0.3, 0.3),
        colormap = Reverse(:RdBu))
    
    
    qs = lms_mat[:, "q_00to120"]
    annotations!(bax2, map(x-> x < 0.05 ? "**" : x < 0.2 ? "*" : "", qs), [Point2f(1.5, i + 0.2) for i in eachindex(lms_mat.feature)];
    color = [abs(x) > 0.2 ? :white : :black for x in clrs],
    align=(:center, :center))
end

Colorbar(B[1,3];
    colormap=Reverse(:RdBu),
    label="correlation",
    limits=(-0.3, 0.3)
)


colgap!(B, Fixed(4))

#-

scales = unique(replace.(names(subscale_mat, r"^[A-Z]"), r"_.+"=> ""))
scales_names = map(scales) do n
    m = match(r"([A-Z][a-z]+)([A-Z][a-z]+)T", n)
    isnothing(m) && error("$n doesn't match")
    return string(m[1], " ", replace(m[2], "Language"=> "Lang."))
end
cylabs = format_species_labels(subscale_mat.species)
filters = ["00to06", "18to120", "00to120"]

cax1 = Axis(C[1,1]; xlabel = "subscore", title="under 6m", 
    yticklabelfont = "TeX Gyre Heros Makie Italic",
    yticklabelsize = 14,
    yticks = (range(1.5; stop = size(subscale_mat, 1) + 0.5), cylabs),
    xticks = (range(1.5, length(scales)+0.5), lowercase.(scales_names)),
    xticklabelrotation = pi/4)

cax2 = Axis(C[1,2]; xlabel = "subscore", title="18 to 48m",
                    xticks = (range(1.5, length(scales)+0.5), lowercase.(scales_names)),
                    xticklabelrotation = pi/4)
cax2.yticksvisible = false
cax2.yticklabelsvisible = false

cax3 = Axis(C[1,3]; xlabel = "subscore", title="0 to 48m",
                    xticks = (range(1.5, length(scales)+0.5), lowercase.(scales_names)),
                    xticklabelrotation = pi/4)
cax3.yticksvisible = false
cax3.yticklabelsvisible = false

let
    clrs = vcat((subscale_mat[:, "$(scale)_00to06_t"] for scale in scales)...)
    qs = vcat((subscale_mat[:, "$(scale)_00to06_q"] for scale in scales)...)
    poly!(cax1, [Rect(j, i, 1, 1) for j in eachindex(scales) for i in eachindex(subscale_mat.species)]; 
        color = clrs, colormap = Reverse(:PuOr))

    annotations!(cax1, map(x-> x < 0.05 ? "**" : x < 0.2 ? "*" : "", qs), [Point2f(j + 0.5, i + 0.2) for j in eachindex(scales) for i in eachindex(subscale_mat.species)];
        color = [abs(x) > 0.2 ? :white : :black for x in clrs],
        align=(:center, :center))
end
let
    clrs = vcat((subscale_mat[:, "$(scale)_18to120_t"] for scale in scales)...)
    qs = vcat((subscale_mat[:, "$(scale)_18to120_q"] for scale in scales)...)
    poly!(cax2, [Rect(j, i, 1, 1) for j in eachindex(scales) for i in eachindex(subscale_mat.species)]; 
        color = clrs, colormap = Reverse(:PuOr))

    annotations!(cax2, map(x-> x < 0.05 ? "**" : x < 0.2 ? "*" : "", qs), [Point2f(j + 0.5, i + 0.2) for j in eachindex(scales) for i in eachindex(subscale_mat.species)];
        color = [abs(x) > 0.2 ? :white : :black for x in clrs],
        align=(:center, :center))
end
let
    clrs = vcat((subscale_mat[:, "$(scale)_00to120_t"] for scale in scales)...)
    qs = vcat((subscale_mat[:, "$(scale)_00to120_q"] for scale in scales)...)
    poly!(cax3, [Rect(j, i, 1, 1) for j in eachindex(scales) for i in eachindex(subscale_mat.species)]; 
        color = clrs, colormap = Reverse(:PuOr))

    annotations!(cax3, map(x-> x < 0.05 ? "**" : x < 0.2 ? "*" : "", qs), [Point2f(j + 0.5, i + 0.2) for j in eachindex(scales) for i in eachindex(subscale_mat.species)];
        color = [abs(x) > 0.2 ? :white : :black for x in clrs],
        align=(:center, :center))
end
tightlimits!.([cax1,cax2,cax3])

Colorbar(C[1,4];
    colormap=Reverse(:PuOr),
    label="T score",
    limits=round.(extrema(Matrix(select(subscale_mat, r"_t$"))); digits=1)
)
colgap!(C, Fixed(4))
#-

let
    gs = "Propionate degradation I"
    panel = D
    ixs = neuroactive_00to06[gs]
    cs = filter(!isnan, cors_00to06.stat[ixs])
    acs = filter(!isnan, cors_00to06.stat[Not(ixs)])

    (_, dax, _) = Resonance.plot_fsea!(panel, cs, acs;
        title = "Under 6m",
        xticks = -0.7:0.2:0.0)
    Label(panel[3,1], replace(lowercase(gs), "degradation"=> "degr.", "synthesis"=> "synth.", " (vitamin K2)"=> ""); tellheight=true, tellwidth=false)
    rowgap!(panel, 2, Fixed(4))
    gs = "Glutamate synthesis II"
    panel = E
    ixs = neuroactive_18to120[gs]
    cs = filter(!isnan, cors_18to120.stat[ixs])
    acs = filter(!isnan, cors_18to120.stat[Not(ixs)])

    (_, eax, _) = Resonance.plot_fsea!(panel, cs, acs;
        title = "Over 18m",
        )
    Label(panel[3,1], replace(lowercase(gs), "degradation"=> "degr.", "synthesis"=> "synth.", " (vitamin K2)"=> ""); tellheight=true, tellwidth=false)
    rowgap!(panel, 2, Fixed(4))
    linkyaxes!(dax, eax)
    ylims!(dax, -0.7, 0.5)
    ylims!(eax, -0.7, 0.5)
end

#-

let
    genesets = union(subset(fsdf2_00to06, "qvalue"=> ByRow(<(0.2)), "cortest"=> ByRow(==("cogScore"))).geneset,
                     subset(fsdf2_18to120, "qvalue"=> ByRow(<(0.2)), "cortest"=> ByRow(==("cogScore"))).geneset,
                     subset(fsdf2_00to120, "qvalue"=> ByRow(<(0.2)), "cortest"=> ByRow(==("cogScore"))).geneset
               )
    gsidx = Dict(gs=> i for (i, gs) in enumerate(genesets))
  
    df = sort(subset(fsdf2_00to06, "geneset"=> ByRow(gs-> gs in genesets), "cortest"=>ByRow(==("cogScore"))), :geneset; rev=true)
    ax1 = Axis(F[1,1]; yticks = (1:length(genesets), replace.(lowercase.(genesets), r" \(.+\)" => "", "synthesis"=>"syn.", "degradation"=>"deg.")), 
                xlabel="T stat", title="under 6m", alignmode=Mixed(; left=-15))
    m = median(filter(x-> !isnan(x) && x < 7, cors_00to06.stat))
    colors = ColorSchemes.colorschemes[:RdBu_7]

    for gs in genesets
        i = gsidx[gs]
        rowi = findfirst(==(gs), df.geneset)
        if isnothing(rowi)
            row = (; enrichment=0.0, qvalue=1.0, geneset=gs)
        else
            row = df[rowi,:]
        end 
        sign = row.enrichment < 0 ? "neg" : "pos"
        c = row.qvalue > 0.2 ? :white : 
            row.qvalue > 0.05 ? (sign == "pos" ? colors[3] : colors[5]) :
            row.qvalue > 0.01 ? (sign == "pos" ? colors[2] : colors[6]) :
            sign == "pos" ? colors[1] : colors[7]

        y = filter(x-> !isnan(x) && x < 7, cors_00to06.stat[neuroactive_00to06[row.geneset]])
        scatter!(ax1, y, rand(Normal(0, 0.1), length(y)) .+ i; color=(c,0.5), strokecolor=:gray, strokewidth=0.5)
        row.qvalue < 0.2 && lines!(ax1, fill(median(y), 2), [i-0.4, i+0.4]; color = c in colors[1:3] ? colors[1] : c in colors[5:7] ? colors[7] : :white , linewidth=2)
    end
    vlines!(ax1, m; linestyle=:dash, color=:darkgray)

    ####

    df = sort(subset(fsdf2_18to120, "geneset"=> ByRow(gs-> gs in genesets), "cortest"=>ByRow(==("cogScore"))), :geneset; rev=true)
    ax2 = Axis(F[1,2]; yticks = (1:length(genesets), replace.(genesets, r" \(.+\)" => "", "synthesis"=>"syn.", "degradation"=>"deg.")), 
                xlabel="T stat", title="over 18m")
    hideydecorations!(ax2, grid=false)

    m = median(filter(x-> !isnan(x) && x < 7, cors_18to120.stat))
    colors = ColorSchemes.colorschemes[:RdBu_7]

    for gs in genesets
        i = gsidx[gs]
        rowi = findfirst(==(gs), df.geneset)
        if isnothing(rowi)
            row = (; enrichment=0.0, qvalue=1.0, geneset=gs)
        else
            row = df[rowi, :] 
        end
        sign = row.enrichment < 0 ? "neg" : "pos"
        c = row.qvalue > 0.2 ? :white : 
            row.qvalue > 0.05 ? (sign == "pos" ? colors[3] : colors[5]) :
            row.qvalue > 0.01 ? (sign == "pos" ? colors[2] : colors[6]) :
            sign == "pos" ? colors[1] : colors[7]

        y = filter(x-> !isnan(x) && x < 7, cors_18to120.stat[neuroactive_18to120[row.geneset]])
        scatter!(ax2, y, rand(Normal(0, 0.1), length(y)) .+ i; color=(c,0.5), strokecolor=:gray, strokewidth=0.5)
        row.qvalue < 0.2 && lines!(ax2, fill(median(y), 2), [i-0.4, i+0.4]; color = c in colors[1:3] ? colors[1] : c in colors[5:7] ? colors[7] : :white , linewidth=2)
    end
    vlines!(ax2, m; linestyle=:dash, color=:darkgray)

    ####

    df = sort(subset(fsdf2_00to120, "geneset"=> ByRow(gs-> gs in genesets), "cortest"=>ByRow(==("cogScore"))), :geneset; rev=true)
    ax3 = Axis(F[1,3]; yticks = (1:length(genesets), replace.(genesets, r" \(.+?\)" => "", "synthesis"=>"syn.", "degradation"=>"deg.")), 
                xlabel="T stat", title="all ages")
    hideydecorations!(ax3, grid=false)
    m = median(filter(x-> !isnan(x) && x < 7, cors_00to120.stat))
    colors = ColorSchemes.colorschemes[:RdBu_7]

    for gs in genesets
        i = gsidx[gs]
        rowi = findfirst(==(gs), df.geneset)
        if isnothing(rowi)
            row = (; enrichment=0.0, qvalue=1.0, geneset=gs)
        else
            row = df[rowi, :]
        end
        sign = row.enrichment < 0 ? "neg" : "pos"
        c = row.qvalue > 0.2 ? :white : 
            row.qvalue > 0.05 ? (sign == "pos" ? colors[3] : colors[5]) :
            row.qvalue > 0.01 ? (sign == "pos" ? colors[2] : colors[6]) :
            sign == "pos" ? colors[1] : colors[7]

        y = filter(x-> !isnan(x) && x < 7, cors_00to120.stat[neuroactive_00to120[row.geneset]])
        scatter!(ax3, y, rand(Normal(0, 0.1), length(y)) .+ i; color=(c,0.5), strokecolor=:gray, strokewidth=0.5)
        row.qvalue < 0.2 && lines!(ax3, fill(median(y), 2), [i-0.4, i+0.4]; color = c in colors[1:3] ? colors[1] : c in colors[5:7] ? colors[7] : :white , linewidth=2)
    end
    vlines!(ax3, m; linestyle=:dash, color=:darkgray)
    linkyaxes!(ax1,ax2,ax3)

    Legend(F[1,4], [MarkerElement(; color = (c, 0.5),
                                    marker=:circle,
                                    strokecolor=:gray,
                                    strokewidth=0.5) for c in colors[[1:3..., 5:7...]]],
                   ["(+) q < 0.01", "(+) q < 0.05", "(+) q < 0.2", 
                    "(-) q < 0.20", "(-) q < 0.05", "(-) p < 0.01"])
end

# rowsize!(figure.layout, 3, Relative(2/5))
colgap!(F, Fixed(4))

#-

Label(AB[1, 1, TopLeft()], "A",
    fontsize = 26,
    font = "Open Sans Bold",
    padding = (0, 30, 5, 0),
    halign = :right)
Label(AB[1, 3, TopLeft()], "B",
    fontsize = 26,
    font = "Open Sans Bold",
    padding = (0, 150, 5, 0),
    halign = :right)

Label(C[1, 1, TopLeft()], "C",
    fontsize = 26,
    font = "Open Sans Bold",
    padding = (0, 0, 0, 0),
    tellwidth = false,
    halign = :right)

for (label, layout) in zip(["D", "E"], [D, E])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 26,
        font = "Open Sans Bold",
        padding = (0, 30, 5, 0),
        halign = :right)
end

Label(figure[3, 1, TopLeft()], "F";
             fontsize = 26,
             font = "Open Sans Bold",
             padding = (0, 30, 5, 0),
             halign = :right
       )

rowsize!(figure.layout, 1, Relative(1/5))
rowsize!(figure.layout, 3, Relative(1/5))
colsize!(AB, 3, Relative(1/4))
colsize!(F, 1, Relative(2/5))
figure
```

```julia
save("manuscript/assets/Figure2.png", figure)
save(figurefiles("Figure2.svg"), figure)
```


## Supplement

### Other FSEA Plots

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

figure = Figure()
A = GridLayout(figure[1,1])
B = GridLayout(figure[2,1])
CDEF = GridLayout(figure[1,2])
GHIJ = GridLayout(figure[2,2])

Axis(A[1,1])
Legend(A[2,1], [MarkerElement(; color=c, marker=:rect) for c in (:red,:blue,:orange,:teal)], ["thing$i" for i in 1:4]; orientation=:horizontal, tellheight=true, tellwidth=false)
Axis(B[1,1])
Legend(B[2,1], [MarkerElement(; color=c, marker=:rect) for c in (:red,:blue,:orange,:teal)], ["thing$i" for i in 1:4]; orientation=:horizontal, tellheight=true, tellwidth=false)

for (i, j) in zip((1,2,1,2), (1,1,2,2))
    Axis(CDEF[i,j]; title = "$i and $j")
end
for (i, j) in zip((1,2,1,2), (1,1,2,2))
    Axis(GHIJ[i,j]; title = "$i and $j")
end

colsize!(figure.layout, 1, Relative(1/3))
figure
```
