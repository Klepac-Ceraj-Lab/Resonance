using Resonance
using FileIO
using ThreadsX
using Statistics
using CairoMakie # for plotting
using ColorSchemes
using Distances
using MultivariateStats
using CategoricalArrays
using MLJ
using JLD2

#-

mdata = Resonance.load(Metadata())


gdf = groupby(mdata, "subject")
subset!(gdf, AsTable(r"filter_00to120") => (nt -> any(any(values(row)) for row in nt)); ungroup= false)
ys = Dict(row.subject => i for (i, row) in enumerate(
    eachrow(sort(DataFrames.combine(gdf, "ageMonths"=> minimum=>"min"), "min")))
)
mdata = subset(mdata, AsTable(r"filter_")=> ByRow(any))
seqs = subset(mdata, "sample"=> ByRow(!ismissing)) 
seqs.sample = [s for s in seqs.sample]
taxa = Resonance.load(TaxonomicProfiles(); timepoint_metadata = seqs)
taxdf = comm2wide(taxa)

#-

fig = Figure(; resolution=(400,660))

ax = Axis(fig[1,1]; xlabel="Age (months)", ylabel="Subject")
cs = [:darkorange4, :dodgerblue, :slateblue3, :gray70]

for sub in keys(gdf)
    y = Float64(ys[sub.subject])
    ages = gdf[sub].ageMonths
    has_stool = map(!ismissing, gdf[sub].sample)
    has_cog = map(!ismissing, gdf[sub].cogScore)
    c = map(zip(has_stool, has_cog)) do (s,c)
        s && c && return cs[3]
        s && return cs[1]
        c && return cs[2]
        return cs[4]
    end
    scatter!(ax, ages, fill(y, length(ages)); color=c, markersize=6)
    lines!(ax, [extrema(ages)...], fill(y, 2); color=:gray70, linestyle=:dash, linewidth=1)
end
tightlimits!(ax)
Legend(fig[2,1], [MarkerElement(; color=c, marker=:circle) for c in cs], ["Stool", "Cog", "Both", "Neither"];
    orientation=:horizontal, tellwidth=false, tellheight=true)
save("/home/kevin/Downloads/resonance-age-vs-subject.png", fig)
fig


#-

B = Figure(; resolution=(1200,600))

bax1 = Axis(B[1,1]; xlabel = "Age (months)", ylabel = "cogScore", xticks=(4:4:24), limits=((2,24), nothing))
bax2 = Axis(B[1,2]; xlabel = "Age (years)", xticks = (4:2:12), limits=((2,nothing), nothing))
hideydecorations!(bax2)
linkyaxes!(bax1, bax2)

bleg = let
    u2y = findall(p-> !ismissing(p[2]) && p[1] <= 24, collect(zip(mdata.ageMonths, mdata.cogScore)))
    o2y = findall(p-> !ismissing(p[2]) && p[1] > 24, collect(zip(mdata.ageMonths, mdata.cogScore)))

    cs = ColorSchemes.colorschemes[:tab20b][[1,7,15]]
    function colorage(age)
        age <= 36 && return cs[1] # mullen
        age <= 60 && return cs[2] # WPPSI
        return cs[3] # WISC
    end
    ages = mdata.ageMonths[u2y]
    scatter!(bax1, ages, mdata.cogScore[u2y]; color=colorage.(ages))

    ages = mdata.ageMonths[o2y]
    scatter!(bax2, ages ./ 12, mdata.cogScore[o2y]; color=colorage.(ages))
    
    vlines!(bax1, [6, 18]; linestyle=:dash, color=:gray)
    # TODO: add colors for training/test sets in later models

    Legend(B[1, 3], [MarkerElement(; marker=:circle, color=c) for c in cs[1:3]],
                      ["MSEL", "WPPSI", "WISC"];
    )
end

save("/home/kevin/Downloads/resonance-cog-vs-age.png", B)   
B

#-

speclms_00to120 = CSV.read(tablefiles("figure2", "lms_species_00to120.csv"), DataFrame)
speclms_00to06 = CSV.read(tablefiles("figure2", "lms_species_00to06.csv"), DataFrame)
speclms_18to120 = CSV.read(tablefiles("figure2", "lms_species_18to120.csv"), DataFrame)
subscale_lms = CSV.read(tablefiles("newfigures", "subscale_lms.csv"), DataFrame)
subscale_lms.kind[findall(k-> contains(k, "Expressivelan"), subscale_lms.kind)] .= "Mullen::mullen_ExpressiveLanguageT"

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

#-

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

#-


figure = Figure(resolution=(900, 500));

A = GridLayout(figure[1,1:2])
B = GridLayout(figure[1,3])
aax1 = Axis(A[1,1]; title = "over 18m", ylabel=L"$-log_2(P)$", xlabel = "coef.",
                xticks = ([-1.5e-3, 0.0, 1.5e-3], ["-1.5e-3", "0.0", "1.5e-3"]),
                limits = ((-1.5e-3, 1.5e-3), nothing),
                xminorticksvisible=true, xminorticks = IntervalsBetween(3))
aax2 = Axis(A[2,1]; title = "all ages", ylabel=L"$-log_2(P)$", xlabel = "coef.",
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
save("/home/kevin/Downloads/resonance-lms.png", figure) 
figure
#-

C = Figure(; resolution=(1200, 600))

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

cax2 = Axis(C[1,2]; xlabel = "subscore", title="18 to 36m",
                    xticks = (range(1.5, length(scales)+0.5), lowercase.(scales_names)),
                    xticklabelrotation = pi/4)
cax2.yticksvisible = false
cax2.yticklabelsvisible = false

cax3 = Axis(C[1,3]; xlabel = "subscore", title="0 to 36m",
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
colgap!(C.layout, Fixed(4))
save("/home/kevin/Downloads/resonance-subscale.png", C)
C

#-


RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree
# concurrent cogScore regression from taxonomic profiles
regression_currentCogScores_00to06mo_onlydemo = JLD2.load(
    modelfiles("regression_currentCogScores_00to06mo_onlydemo.jld"),
    "regression_currentCogScores_00to06mo_onlydemo"
);
regression_currentCogScores_00to06mo_onlytaxa = JLD2.load(
    modelfiles("regression_currentCogScores_00to06mo_onlytaxa.jld"),
    "regression_currentCogScores_00to06mo_onlytaxa"
);
regression_currentCogScores_00to06mo_demoplustaxa = JLD2.load(
    modelfiles("regression_currentCogScores_00to06mo_demoplustaxa.jld"),
    "regression_currentCogScores_00to06mo_demoplustaxa"
);
regression_currentCogScores_00to06mo_onlyecs = JLD2.load(
    modelfiles("regression_currentCogScores_00to06mo_onlyecs.jld"),
    "regression_currentCogScores_00to06mo_onlyecs"
);
regression_currentCogScores_00to06mo_demoplusecs = JLD2.load(
    modelfiles("regression_currentCogScores_00to06mo_demoplusecs.jld"),
    "regression_currentCogScores_00to06mo_demoplusecs"
);
regression_currentCogScores_18to120mo_onlydemo = JLD2.load(
    modelfiles("regression_currentCogScores_18to120mo_onlydemo.jld"),
    "regression_currentCogScores_18to120mo_onlydemo"
);
regression_currentCogScores_18to120mo_onlytaxa = JLD2.load(
    modelfiles("regression_currentCogScores_18to120mo_onlytaxa.jld"),
    "regression_currentCogScores_18to120mo_onlytaxa"
);
regression_currentCogScores_18to120mo_demoplustaxa = JLD2.load(
    modelfiles("regression_currentCogScores_18to120mo_demoplustaxa.jld"),
    "regression_currentCogScores_18to120mo_demoplustaxa"
);
regression_currentCogScores_18to120mo_onlyecs = JLD2.load(
    modelfiles("regression_currentCogScores_18to120mo_onlyecs.jld"),
    "regression_currentCogScores_18to120mo_onlyecs"
);
regression_currentCogScores_18to120mo_demoplusecs = JLD2.load(
    modelfiles("regression_currentCogScores_18to120mo_demoplusecs.jld"),
    "regression_currentCogScores_18to120mo_demoplusecs"
);
regression_futureCogScores_demoplustaxa = JLD2.load(
    modelfiles("regression_futureCogScores_demoplustaxa.jld"),
    "regression_futureCogScores_demoplustaxa"
);
regression_futureCogScores_onlytaxa = JLD2.load(
    modelfiles("regression_futureCogScores_onlytaxa.jld"),
    "regression_futureCogScores_onlytaxa"
);

#-


AB_subfig = Figure(resolution = (800, 400));
axA = Axis(
    AB_subfig[1, 1];
    xlabel = "-log(p) for LM coefficients",
    ylabel = "relative importance",
    title = "over 18mo",
)
t10 = ColorSchemes.tableau_10
let plot_colorset = [t10[1], t10[3], t10[7], (:white, 0.)]
# plot_comparative_lmvsrf_scatterplots!(axA, regression_currentCogScores_00to06mo_onlytaxa, tablefiles("figure2", "lms_species_00to06.csv"); plot_colorset = plot_colorset, strokewidth=1)
# plot_comparative_lmvsrf_scatterplots!(axB, regression_currentCogScores_18to120mo_onlytaxa, tablefiles("figure2", "lms_species_18to120.csv"); plot_colorset = plot_colorset, strokewidth=1)
    plot_comparative_lmvsrf_scatterplots!(axA, regression_currentCogScores_18to120mo_onlytaxa,
        tablefiles("figure2/lms_species_18to120.csv");
        plot_colorset, strokewidth=1
    ) 

    Legend(
        AB_subfig[1, 2],
        [
            MarkerElement(; marker=:circle, color=plot_colorset[1], strokewidth=1),
            MarkerElement(; marker=:circle, color=plot_colorset[2], strokewidth=1),
            MarkerElement(; marker=:circle, color=plot_colorset[3], strokewidth=1),
            MarkerElement(; marker=:circle, color=plot_colorset[4], strokewidth=1),
        ],
        [
            "> 60% ranked\nimportance",
            "q < 0.2 in LM",
            "Both",
            "None"
        ];
        tellheight = false,
        tellwidth = true,
    )
end


save("/home/kevin/Downloads/resonance-lm-vs-rf.png", AB_subfig)
AB_subfig

#-

axB = GridLayout(AB_subfig[2,1])

let plot_colorset = [t10[1], t10[6], t10[5], (:white, 0.)]
    Resonance.plot_comparative_rfvsrf_scatterplots!(axB,
        regression_currentCogScores_18to120mo_onlytaxa,
        regression_futureCogScores_onlytaxa;
        plot_colorset, strokewidth=1
    ) 

    Legend(
        AB_subfig[2, 2],
        [
            MarkerElement(; marker=:circle, color=plot_colorset[1], strokewidth=1),
            MarkerElement(; marker=:circle, color=plot_colorset[2], strokewidth=1),
            MarkerElement(; marker=:circle, color=plot_colorset[3], strokewidth=1),
            MarkerElement(; marker=:circle, color=plot_colorset[4], strokewidth=1),
        ],
        [
            "over 18mo",
            "future",
            "Both",
            "None"
        ];
        tellheight = false,
        tellwidth = true,
    )
end

#-

joined_importances_00to06 = compute_joined_importances(
    regression_currentCogScores_00to06mo_onlytaxa,
    regression_currentCogScores_00to06mo_demoplustaxa;
    imp_fun = weighted_hpimportances
)
joined_importances_18to120 = compute_joined_importances(
    regression_currentCogScores_18to120mo_onlytaxa,
    regression_currentCogScores_18to120mo_demoplustaxa;
    imp_fun = weighted_hpimportances
)

joined_importances_future = compute_joined_importances(
    regression_futureCogScores_onlytaxa,
    regression_futureCogScores_demoplustaxa;
    imp_fun = weighted_hpimportances
)


joined_importances_future = compute_joined_importances(
    regression_futureCogScores_onlytaxa,
    regression_futureCogScores_demoplustaxa;
    imp_fun = weighted_hpimportances
)
#-

nbars_toplot = 10

CD_subfig = Figure(resolution = (1200, 600));

axC = Axis(
    CD_subfig[1, 1];
    xlabel = "relative importance",
    yticks = (reverse(collect(1:nbars_toplot)), format_species_labels(joined_importances_18to120.variable[1:nbars_toplot])),
    ylabel = "Feature",
    title = "18 to 120 months",
    # yticklabelrotation= -pi/2
)

plot_comparativedemo_importance_barplots!(axC, joined_importances_18to120; n_rows = nbars_toplot)


Legend(
    CD_subfig[1, 1], [MarkerElement(; marker=:rect, color=:gray), MarkerElement(; marker=:star8, color=:red)], ["Microbiome alone", "Microbiome + demographics"];
    tellheight = false,
    tellwidth = false,
    margin = (10, 10, 10, 10),
    halign = :right, valign = :bottom, orientation = :vertical
)

save("/home/kevin/Downloads/resonance-importance-18to120.png", CD_subfig)

#-

