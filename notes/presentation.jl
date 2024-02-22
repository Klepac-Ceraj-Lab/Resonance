using ColorSchemes: grays
using Resonance
using FileIO
using ThreadsX
using Statistics
using CairoMakie # for plotting
using ColorSchemes
using Distances
using MultivariateStats
using CategoricalArrays
using Clustering
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

B = Figure(; resolution=(700, 350))

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
    scatter!(bax1, ages, mdata.cogScore[u2y]; color=[(colorage(age), 0.4) for age in ages])

    ages = mdata.ageMonths[o2y]
    scatter!(bax2, ages ./ 12, mdata.cogScore[o2y]; color=[(colorage(age), 0.4) for age in ages])
    
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

fig = Figure(; resolution=(1100,700))

axB = Axis(fig[1,2])

let plot_colorset = [t10[1], t10[6], t10[5], (:white, 0.)]
    Resonance.plot_comparative_rfvsrf_scatterplots!(axB,
        regression_currentCogScores_18to120mo_onlytaxa,
        regression_futureCogScores_onlytaxa;
        plot_colorset, strokewidth=1
    ) 

    Legend(
        fig[1, 3],
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



############
# Figure 4 #
############

# Composite cognitive scores and Mullen subscales
regression_currentCogScores_18to120mo_onlydemo = JLD2.load(
    modelfiles("regression_currentCogScores_18to120mo_onlydemo.jld")
    )["regression_currentCogScores_18to120mo_onlydemo"];
regression_currentCogScores_18to120mo_demoplustaxa = JLD2.load(
    modelfiles("regression_currentCogScores_18to120mo_demoplustaxa.jld")
    )["regression_currentCogScores_18to120mo_demoplustaxa"];
regression_currentCogScores_18to120mo_demoplusecs = JLD2.load(
    modelfiles("regression_currentCogScores_18to120mo_demoplusecs.jld")
    )["regression_currentCogScores_18to120mo_demoplusecs"];
regression_currentExpressiveLanguages_18to120mo_onlydemo = JLD2.load(
    modelfiles("regression_currentExpressiveLanguages_18to120mo_onlydemo.jld")
    )["regression_currentExpressiveLanguages_18to120mo_onlydemo"];
regression_currentExpressiveLanguages_18to120mo_demoplustaxa = JLD2.load(
    modelfiles("regression_currentExpressiveLanguages_18to120mo_demoplustaxa.jld")
    )["regression_currentExpressiveLanguages_18to120mo_demoplustaxa"];
regression_currentExpressiveLanguages_18to120mo_demoplusecs = JLD2.load(
    modelfiles("regression_currentExpressiveLanguages_18to120mo_demoplusecs.jld")
    )["regression_currentExpressiveLanguages_18to120mo_demoplusecs"];
regression_currentGrossMotors_18to120mo_onlydemo = JLD2.load(
    modelfiles("regression_currentGrossMotors_18to120mo_onlydemo.jld")
    )["regression_currentGrossMotors_18to120mo_onlydemo"];
regression_currentGrossMotors_18to120mo_demoplustaxa = JLD2.load(
    modelfiles("regression_currentGrossMotors_18to120mo_demoplustaxa.jld")
    )["regression_currentGrossMotors_18to120mo_demoplustaxa"];
regression_currentGrossMotors_18to120mo_demoplusecs = JLD2.load(
    modelfiles("regression_currentGrossMotors_18to120mo_demoplusecs.jld")
    )["regression_currentGrossMotors_18to120mo_demoplusecs"];
regression_currentVisualReceptions_18to120mo_onlydemo = JLD2.load(
    modelfiles("regression_currentVisualReceptions_18to120mo_onlydemo.jld")
    )["regression_currentVisualReceptions_18to120mo_onlydemo"];
regression_currentVisualReceptions_18to120mo_demoplustaxa = JLD2.load(
    modelfiles("regression_currentVisualReceptions_18to120mo_demoplustaxa.jld")
    )["regression_currentVisualReceptions_18to120mo_demoplustaxa"];
regression_currentVisualReceptions_18to120mo_demoplusecs = JLD2.load(
    modelfiles("regression_currentVisualReceptions_18to120mo_demoplusecs.jld")
    )["regression_currentVisualReceptions_18to120mo_demoplusecs"];

regression_currentCogScores_18to120mo_onlytaxa = JLD2.load(
    modelfiles("regression_currentCogScores_18to120mo_onlytaxa.jld")
    )["regression_currentCogScores_18to120mo_onlytaxa"];
regression_currentCogScores_18to120mo_onlyecs = JLD2.load(
    modelfiles("regression_currentCogScores_18to120mo_onlyecs.jld")
    )["regression_currentCogScores_18to120mo_onlyecs"];
regression_currentExpressiveLanguages_18to120mo_onlytaxa = JLD2.load(
    modelfiles("regression_currentExpressiveLanguages_18to120mo_onlytaxa.jld")
    )["regression_currentExpressiveLanguages_18to120mo_onlytaxa"];
regression_currentExpressiveLanguages_18to120mo_onlyecs = JLD2.load(
    modelfiles("regression_currentExpressiveLanguages_18to120mo_onlyecs.jld")
    )["regression_currentExpressiveLanguages_18to120mo_onlyecs"];
regression_currentGrossMotors_18to120mo_onlytaxa = JLD2.load(
    modelfiles("regression_currentGrossMotors_18to120mo_onlytaxa.jld")
    )["regression_currentGrossMotors_18to120mo_onlytaxa"];
regression_currentGrossMotors_18to120mo_onlyecs = JLD2.load(
    modelfiles("regression_currentGrossMotors_18to120mo_onlyecs.jld")
    )["regression_currentGrossMotors_18to120mo_onlyecs"];
regression_currentVisualReceptions_18to120mo_onlytaxa = JLD2.load(
    modelfiles("regression_currentVisualReceptions_18to120mo_onlytaxa.jld")
    )["regression_currentVisualReceptions_18to120mo_onlytaxa"];
regression_currentVisualReceptions_18to120mo_onlyecs = JLD2.load(
    modelfiles("regression_currentVisualReceptions_18to120mo_onlyecs.jld")
    )["regression_currentVisualReceptions_18to120mo_onlyecs"];
# concurrent brain regression from taxonomic profiles
brain_models = JLD2.load(modelfiles("brain_models.jld"))["brain_models"]
include("manuscript/Figure4-definitions.jl")

mean_brain_merits = reduce(
    vcat,
    [ DataFrame(:variable => ordered_brain_segments_list[i], :Train_Cor => mean(brain_models[ordered_brain_segments_list[i]].merits.Train_Cor), :Test_Cor => mean(brain_models[ordered_brain_segments_list[i]].merits.Test_Cor)) for i in eachindex(ordered_brain_segments_list) ]
)


weighted_brain_importances = dropmissing(reduce(
    (x, y) -> outerjoin(x, y, on = :variable, makeunique=true),
    [ rename!(weighted_hpimportances(j; normalize_importances=false), :weightedImportance => Symbol(split(j.name, '_')[2]))
      for (i, j) in brain_models ] ) )

relative_brain_importances = dropmissing(reduce(
    (x, y) -> outerjoin(x, y, on = :variable, makeunique=true),
    [ rename!(weighted_hpimportances(j; normalize_importances=true), :weightedImportance => Symbol(split(j.name, '_')[2]))
      for (i, j) in brain_models ] ))


function generate_cols(model_df, xs, grp, color)
    insertcols(DataFrames.combine(
        groupby(model_df.merits, :Hyperpar_Idx),
            "Test_Cor" => (x -> round(mean(x); digits=2))=> "Test_Cor_mean"), 1,
            "xs"=> xs, "grp"=> grp, "color"=> color
    )

end 
cm = [ColorSchemes.Spectral_10[i] for i in (10,2,4,8)]
panelA_plot_df = let
    vcat(
        generate_cols(regression_currentCogScores_18to120mo_demoplustaxa,           3,  1, cm[1] ),
        generate_cols(regression_currentCogScores_18to120mo_demoplusecs,            4,  1, cm[1] ),
        generate_cols(regression_currentExpressiveLanguages_18to120mo_demoplustaxa, 3,  2, cm[2] ),
        generate_cols(regression_currentExpressiveLanguages_18to120mo_demoplusecs,  4,  2, cm[2] ),
        generate_cols(regression_currentGrossMotors_18to120mo_demoplustaxa,         3,  3, cm[3] ),
        generate_cols(regression_currentGrossMotors_18to120mo_demoplusecs,          4,  3, cm[3] ),
        generate_cols(regression_currentVisualReceptions_18to120mo_demoplustaxa,    3,  4, cm[4] ),
        generate_cols(regression_currentVisualReceptions_18to120mo_demoplusecs,     4,  4, cm[4] ),
        generate_cols(regression_currentCogScores_18to120mo_onlytaxa,           1,  1, cm[1] ),
        generate_cols(regression_currentCogScores_18to120mo_onlyecs,            2,  1, cm[1] ),
        generate_cols(regression_currentExpressiveLanguages_18to120mo_onlytaxa, 1,  2, cm[2] ),
        generate_cols(regression_currentExpressiveLanguages_18to120mo_onlyecs,  2,  2, cm[2] ),
        generate_cols(regression_currentGrossMotors_18to120mo_onlytaxa,         1,  3, cm[3] ),
        generate_cols(regression_currentGrossMotors_18to120mo_onlyecs,          2,  3, cm[3] ),
        generate_cols(regression_currentVisualReceptions_18to120mo_onlytaxa,    1,  4, cm[4] ),
        generate_cols(regression_currentVisualReceptions_18to120mo_onlyecs,     2,  4, cm[4] ),
    )
end

composite_importances = rename(
    weighted_hpimportances(regression_currentCogScores_18to120mo_demoplustaxa;
        normalize_importances=true),
        :weightedImportance => :Composite
)
el_importances = rename(
    weighted_hpimportances(regression_currentExpressiveLanguages_18to120mo_demoplustaxa;
        normalize_importances=true),
        :weightedImportance => :ExpressiveLanguage
)
gm_importances = rename(
    weighted_hpimportances(regression_currentGrossMotors_18to120mo_demoplustaxa;
        normalize_importances=true),
        :weightedImportance => :GrossMotor
)
vr_importances = rename(
    weighted_hpimportances(regression_currentVisualReceptions_18to120mo_demoplustaxa;
        normalize_importances=true),
        :weightedImportance => :VisualReception
)

sidebyside_importances = reduce( (x, y) -> outerjoin(x, y, on = :variable; makeunique = true), [ composite_importances, el_importances, gm_importances, vr_importances ])
sidebyside_importances.avg = [ mean(skipmissing(collect(rr))) for rr in eachrow(sidebyside_importances[:, 2:end]) ]
sidebyside_importances.maxacross = [ maximum(skipmissing(collect(rr))) for rr in eachrow(sidebyside_importances[:, 2:end]) ]
sidebyside_importances.stdev = [ Statistics.std(skipmissing(collect(rr))) for rr in eachrow(sidebyside_importances[:, 2:end]) ]

important_bugs_composite = sort(dropmissing(sidebyside_importances, :Composite), :Composite; rev = true).variable[1:12]
important_bugs_ExpressiveLanguage = sort(dropmissing(sidebyside_importances, :ExpressiveLanguage), :ExpressiveLanguage; rev = true).variable[1:12]
important_bugs_GrossMotor = sort(dropmissing(sidebyside_importances, :GrossMotor), :GrossMotor; rev = true).variable[1:12]
important_bugs_VisualReception = sort(dropmissing(sidebyside_importances, :VisualReception), :VisualReception; rev = true).variable[1:12]

panelB_taxa = [
    "Bifidobacterium_pseudocatenulatum",
    "Faecalibacterium_prausnitzii",
    "Blautia_wexlerae",
    "Bifidobacterium_longum",
    "Roseburia_faecis",
    "Eubacterium_eligens",
    "Fusicatenibacter_saccharivorans",
    "Streptococcus_salivarius",
    "Ruminococcus_gnavus",
    "Anaerostipes_hadrus",
    "Intestinibacter_bartlettii",
    "Bacteroides_vulgatus",
    "Clostridium_innocuum",
    "Clostridium_symbiosum"
]

function comparison_table(scale_imp::DataFrame, colname::String, brain_imp::DataFrame, segment::String; nrank = 20)
    scl = sort(scale_imp, colname; rev = true)
    scl.scale_rank = 1:nrow(scl)
    brn = @chain brain_imp begin
        select(["variable", segment])
        sort(segment, rev = true)
        insertcols!( _ , 3, :brain_rank => 1:nrow( _ ))
    end
    return subset(innerjoin(scl, brn, on = :variable; makeunique=true), [:scale_rank, :brain_rank] => (x, y) -> (x .< nrank) .& (y .< nrank))
end

## EL and bilateral Thalamus
comparison_table(el_importances, "ExpressiveLanguage", relative_brain_importances, "left-thalamus-proper"; nrank = 20)
comparison_table(el_importances, "ExpressiveLanguage", relative_brain_importances, "right-thalamus-proper"; nrank = 20)

## GM and central opercular cortex
comparison_table(gm_importances, "GrossMotor", relative_brain_importances, "left-pars-opercularis"; nrank = 25)
comparison_table(gm_importances, "GrossMotor", relative_brain_importances, "right-pars-opercularis"; nrank = 25)
comparison_table(gm_importances, "GrossMotor", relative_brain_importances, "left-pars-orbitalis"; nrank = 25)
comparison_table(gm_importances, "GrossMotor", relative_brain_importances, "right-pars-orbitalis"; nrank = 25)
comparison_table(gm_importances, "GrossMotor", relative_brain_importances, "left-pars-triangularis"; nrank = 25)
comparison_table(gm_importances, "GrossMotor", relative_brain_importances, "right-pars-triangularis"; nrank = 25)
## Fusicatenibacter_saccharivorans

# Sean's paper; VR associated with occipital lobe white matter and visual cortex, cerebellum, and pulvinar nucleus of the thalamus associated with visual reception47,48
comparison_table(vr_importances, "VisualReception", relative_brain_importances, "left-thalamus-proper"; nrank = 30)
comparison_table(vr_importances, "VisualReception", relative_brain_importances, "right-thalamus-proper"; nrank = 30)
comparison_table(vr_importances, "VisualReception", relative_brain_importances, "cerebellar-vermal-lobules-I-V"; nrank = 30)
comparison_table(vr_importances, "VisualReception", relative_brain_importances, "cerebellar-vermal-lobules-VI-VII"; nrank = 30)
comparison_table(vr_importances, "VisualReception", relative_brain_importances, "cerebellar-vermal-lobules-VIII-X"; nrank = 30)
comparison_table(vr_importances, "VisualReception", relative_brain_importances, "left-cerebellum-white-matter"; nrank = 30)
comparison_table(vr_importances, "VisualReception", relative_brain_importances, "right-cerebellum-white-matter"; nrank = 30)

panelB_plot_df = @chain vcat(
    insertcols(
        weighted_hpimportances(
            regression_currentCogScores_18to120mo_demoplustaxa;
            normalize_importances=true
        ), 1, :grp => 1, :color => cm[1]),
    insertcols(
        weighted_hpimportances(
            regression_currentExpressiveLanguages_18to120mo_demoplustaxa;
            normalize_importances=true
        ), 1, :grp => 2, :color => cm[2]),
    insertcols(
        weighted_hpimportances(
            regression_currentGrossMotors_18to120mo_demoplustaxa;
            normalize_importances=true
        ), 1, :grp => 3, :color => cm[3]),
    insertcols(
        weighted_hpimportances(
            regression_currentVisualReceptions_18to120mo_demoplustaxa;
            normalize_importances=true
        ), 1, :grp => 4, :color => cm[4])
) begin
    subset(:variable => ( x -> x .∈ Ref(panelB_taxa) ))
    transform!(:variable => (x -> [ panelB_taxa_idxer[el] for el in x ]) => :xs; renamecols = false)
end

interesting_segments_idxes = mean_brain_merits.variable .∈ Ref(interesting_segments)
mean_brain_merits = mean_brain_merits[interesting_segments_idxes, :]

explore_importances_df = dropmissing(relative_brain_importances[:, vcat(["variable"], mean_brain_merits.variable) ])
rename!(explore_importances_df, :variable => :predictor)
explore_importances_df = stack(explore_importances_df, 2:ncol(explore_importances_df))
subset!(explore_importances_df, :predictor => (x -> x .!= "ageMonths"))
sort!(explore_importances_df, :value; rev = true)
@show unique(explore_importances_df.predictor)[1:20]

interesting_taxa_idxes = relative_brain_importances.variable .∈ Ref(interesting_taxa)
relative_brain_importances = dropmissing(relative_brain_importances[interesting_taxa_idxes, vcat(["variable"], mean_brain_merits.variable) ])

## Transpose the matrix to add information about the symmetric segments
transposedImportances = permutedims(relative_brain_importances, :variable)
insertcols!(transposedImportances, 1, :symmetricSegment => symmetric_segment_strings)
## Combine left and right of symmetric segments
combinedTransposedImportances = DataFrames.combine(
    groupby(transposedImportances, :symmetricSegment),
    Symbol.(names(transposedImportances))[3:end] .=> mean .=> Symbol.(names(transposedImportances))[3:end]
)
## Perform the actual HCA and store the order
dist_taxa = pairwise(Euclidean(), Matrix(combinedTransposedImportances[:, 2:end]); dims=2)
dist_segments = pairwise(Euclidean(), Matrix(combinedTransposedImportances[:, 2:end]); dims=1)
hcl_taxa = hclust(dist_taxa; linkage=:average, branchorder=:optimal)
hcl_segments = hclust(dist_segments; linkage=:average, branchorder=:optimal)
hclust_taxa_order = hcl_taxa.order
hclust_symmetric_segment_order = hcl_segments.order

reorder_segments_df = innerjoin(
    DataFrame(
        :original_segment => interesting_segments,
        :left_or_unique => repeat([1, 0], 17),
        :symmetric_segment => symmetric_segment_strings),
    DataFrame(
        :symmetric_segment => combinedTransposedImportances.symmetricSegment[hclust_symmetric_segment_order],
        :symmetric_order => hclust_symmetric_segment_order,
        :original_order => collect(1:length(combinedTransposedImportances.symmetricSegment[hclust_symmetric_segment_order]))),
    on = :symmetric_segment
)
insertcols!(reorder_segments_df, 1, :plot_order => collect(1:nrow(reorder_segments_df)))
sort!(reorder_segments_df, [:original_order, :left_or_unique])

plot_segments_order = reorder_segments_df.plot_order

mean_brain_merits = mean_brain_merits[plot_segments_order, :]
relative_brain_importances = relative_brain_importances[hclust_taxa_order, vcat( [ 1 ], (plot_segments_order .+ 1))]
#relative_brain_importances = relative_brain_importances[sortperm(map(mean, eachrow(relative_brain_importances[:, 2:end])); rev=true), vcat( [ 1 ], (plot_segments_order .+ 1))]
weighted_brain_importances = reduce(
    (x, y) -> outerjoin(x, y, on = :variable, makeunique=true),
        [ rename!(weighted_hpimportances(j; normalize_importances=true), 
                :weightedImportance => Symbol(split(j.name, '_')[2])) for (i, j) in brain_models ]
)

weighted_brain_importances = dropmissing(weighted_brain_importances[:, vcat(["variable"], interesting_segments) ])
weighted_noage_importances = weighted_brain_importances[2:end, :]

## -

AB_Subfig = Figure(; resolution=(1100, 700))
axA = Axis(
    AB_Subfig[1, 1];
    xlabel = "Correlation",
    yticks = (collect(1:4), [rich("taxa"; color=:green),rich("genes"; color=:green), rich("taxa"; color=:blue), rich("genes"; color=:blue)]),
    ylabel = "Model input composition",
    title = "Average correlation of model",
    yticklabelsize=16,
    alignmode=Outside(),
)

barplot!(axA,
    panelA_plot_df.xs,
    panelA_plot_df.Test_Cor_mean,
    dodge = panelA_plot_df.grp,
    color = panelA_plot_df.color,
    direction = :x
)

highlight_bugs = [
    "Anaerostipes_hadrus",
    "Bacteroides_vulgatus",
    "Blautia_wexlerae",
    "Eubacterium_eligens",
    "Coprococcus_comes",
    "Ruminococcus_torques",
    "Bifidobacterium_pseudocatenulatum",
    "Asaccharobacter_celatus",
    "Bacteroides_ovatus", 
    "Coprococcus_eutactus"
]

hlbugs_color = Dict(b=> c for (c,b) in zip(ColorSchemes.Set3_12, highlight_bugs))
@assert length(interesting_segments) % 2 == 0
@assert all(1:2:length(interesting_segments)) do i
    m1 = match(r"^(left)-(.+)$", interesting_segments[i])
    m2 = match(r"^(right)-(.+)$", interesting_segments[i+1])
    any(isnothing, (m1,m2)) && return false
    m1[2] == m2[2]
end

ableg = Legend(AB_Subfig[1,2],
    [
        [   MarkerElement(; marker=:rect, color=cm[1]),
            MarkerElement(; marker=:rect, color=cm[2]),
            MarkerElement(; marker=:rect, color=cm[3]),
            MarkerElement(; marker=:rect, color=cm[4])
        ],
        [   MarkerElement(; marker=:circle, color=:green),
            MarkerElement(; marker=:circle, color=:blue)
        ]
    ],
    [
        [ "Composite Score", "Expressive Language", "Gross Motor", "Visual Reception"],
        [ "-", "+" ]
    ],
    ["MSEL Subscale", "Demographics"];
    nbanks=2
)
save("/home/kevin/Downloads/resonance-rf-subscales.png", AB_Subfig)


#-
CD_Subfig = Figure(; resolution=(1100, 700), alignmode=Outside())

axC = Axis(
    CD_Subfig[1, 1];
    xlabel = "Correlation",
    yticks = (reverse(collect(1.5:2:length(interesting_segments)+0.5)),
              replace.(mean_brain_merits.variable[1:2:end], r"(right|left)-"=> "")),
    ylabel = "Target Variable",
    xticks = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
    title = "Mean RF correlations",
    yticklabelsize=16,
    alignmode=Inside(),
    #yticklabelrotation= -pi/2
)

tightlimits!(axC, Top())
tightlimits!(axC, Bottom())

barplot!(
    axC,
    reverse(collect(1:length(interesting_segments))),
    mean_brain_merits.Test_Cor,
    color = repeat( [ "blue", "red" ], 17)[plot_segments_order],
    direction=:x
)

cdleg = Legend(CD_Subfig[2,1],
    [MarkerElement(; marker=:rect, color=c) for c in ("blue", "red")],
    ["Left hemisphere", "Right hemisphere"];
    tellheight=false
)

nbugs_toplot = length(interesting_taxa) # the top 1/3 bugs (out of 129) From intersecting the top mean importances and top max importances

axD = Axis(
    CD_Subfig[1, 2];
    xlabel = "Predictor",
    xticks = (collect(1:nbugs_toplot),
            replace.(relative_brain_importances.variable[1:nbugs_toplot], "_"=>" ")),
    xticklabelsize=16,
    xticklabelrotation= pi/4,
    xticklabelfont="TeX Gyre Heros Makie Italic",
    yreversed=true,
    title = "Brain segmentation data variable importances"
)
hideydecorations!(axD)
tightlimits!(axD)

linkxaxes!(axD, axDticks)
hidexdecorations!(axDticks; ticklabels=false, label=false)
hideydecorations!(axDticks)
hidespines!(axDticks)
hm = CairoMakie.heatmap!(axD, Matrix(relative_brain_importances[1:nbugs_toplot, 2:end]), yflip=true)

Colorbar(CD_Subfig[1,3], hm; label= "Relative feature importance", ticks=0:0.01:0.04, minorticksvisible=true)
rowsize!(CD_Subfig.layout, 2, Relative(1/8))

save("/home/kevin/Downloads/resonance-brain-heatmap.png", CD_Subfig)

########
# MISC #
########

mdata = Resonance.load(Metadata())
seqs = subset(mdata, "sample"=> ByRow(!ismissing)) 

seqs.sample = [s for s in seqs.sample]
seqs.edfloat = map(x-> ismissing(x) ? missing : Float64(levelcode(x)), seqs.education)
taxa = Resonance.load(TaxonomicProfiles(); timepoint_metadata = seqs) # this can take a bit
species = filter(f-> taxrank(f) == :species, taxa)
relativeabundance!(species)
unirefs = Resonance.load(UnirefProfiles(); timepoint_metadata = seqs) # this can take a bit
#ecs = Resonance.load(ECProfiles(); timepoint_metadata = seqs)
kos = Resonance.load(KOProfiles(); timepoint_metadata = seqs)

ur_unmapped = abundances(relativeabundance(unirefs)[["UNMAPPED"], :])
ko_unmapped = abundances(relativeabundance(kos)[["UNMAPPED", "UNGROUPED"],:])
#-
let fig = Figure(; resolution=(1100, 700))
    ax1 = Axis(fig[1,1]; xticks=(1:2, ["unirefs", "kos"]), ylabel="N features")
    ax2 = Axis(fig[1,2]; xticks=(1:2, ["unirefs", "kos"]), ylabel="% classified")
    (nu, nk) = nfeatures(unirefs), nfeatures(kos)
    barplot!(ax1, [1,2], [nu, nk]; color=:gray)
    text!(ax1, Point.([1,2], [nu, nk]); text=string.([nu,nk]), offset=(0, 5), align=(:center,:bottom))
    barplot!(ax2, [1,2], 1 .- [mean(ur_unmapped), mean(sum(ko_unmapped; dims=1))]; color=:gray)
    ylims!(ax2, 0., 1.)
    tightlimits!(ax1, Bottom())
    save("/home/kevin/Downloads/resonance-gf-count.svg", current_figure())
end

hist(reshape(prevalence(unirefs), nfeatures(unirefs));
     axis=(; xlabel = "prevalence", ylabel="number of unirefs"), color=:gray
)
save("/home/kevin/Downloads/resonance-uniref-prevalence.svg", current_figure())

#-

lms_unirefs_18to120 = subset(CSV.read(tablefiles("figure2", "lms_unirefs_00to120.csv"), DataFrame),
                             "pvalue"=>ByRow(!isnan))

neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), lms_unirefs_18to120.feature))

#-

fig = Figure(; resolution=(1100, 700))
ax = Axis(fig[1,1]; xlabel="Z statistic", ylabel="count")
hist!(ax, lms_unirefs_18to120.stat; color=:gray)
vlines!(ax, [median(lms_unirefs_18to120[!, "stat"])]; color=:black, linestyle=:dash, linewidth=5)

vlines!(ax, lms_unirefs_18to120[neuroactive["Propionate synthesis"], "stat"]; color=:dodgerblue, linewidth=2)
save("/home/kevin/Downloads/resonance-fsea.png", current_figure())
save("/home/kevin/Downloads/resonance-fsea.svg", current_figure())


