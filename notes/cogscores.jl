include("startup.jl")

using GLM
using MixedModels
using AlgebraOfGraphics
using Statistics

lq, uq = quantile(skipmissing(get(species, :cogScore)), [0.25, 0.75])
for smp in samplenames(species)
    cs = get(species, smp, :cogScore)
    qrt = ismissing(cs) ? missing : 
          cs <= lq      ? "lower" :
          cs >= uq      ? "upper" : 
                          "middle"
                  
    set!(species, smp, :cogQuartile, qrt)
end

specdf = metadata(species) |> DataFrame

##

prevspecies = prevalence(species) |> vec .> 0.1
totest = map(x-> !ismissing(x) && x >=18, get(species, :ageMonths)) .&
         map(x-> !ismissing(x), get(species, :cogScore))

@info count(totest)

testcomm = species[prevspecies, totest]
lm_results = DataFrame()

##

for spc in featurenames(testcomm)
    @info spc
    ab = vec(abundances(testcomm[spc, :]))
    over0 = ab .> 0

    df = DataFrame(
        age = get(testcomm, :ageMonths)[over0],
        cogScore = get(testcomm, :cogScore)[over0],
        cogQuartile = categorical(get(testcomm, :cogQuartile)[over0], levels=["lower", "middle", "upper"]),
        bug = log.(collect(ab[over0]))
    )

    mod = lm(@formula(bug ~ age + cogScore), df)
    ct = DataFrame(coeftable(mod))
    ct.species .= spc
    ct.testvar .= "cogScore"
    append!(lm_results, ct)
    
    mod = lm(@formula(bug ~ age + levelcode(cogQuartile)), df)
    ct = DataFrame(coeftable(mod))
    ct.species .= spc
    ct.testvar .= "cogQuartile"
    append!(lm_results, ct)

end    

select!(lm_results, Cols(:species, :testvar, :))
rename!(lm_results, "Pr(>|t|)"=>"pvalue")
##
using MultipleTesting
using CategoricalArrays

@chain lm_results begin
    groupby(:testvar)
    transform!(:pvalue => (col-> adjust(collect(col), BenjaminiHochberg())) => :qvalue)
end


sort!(lm_results, :qvalue)

CSV.write("data/lm_results/species_cogscore.csv", lm_results)
##

fig = Figure(resolution=(1600,800))
cl = ["lower"=> :seagreen, "middle"=> :lightgray, "upper"=> :dodgerblue]

dt = data(specdf[completecases(specdf, [:ageMonths, :cogScore]), :]) * mapping(:ageMonths => "Age (months)", :cogScore=> "Cognitive function score", color = :cogQuartile)

draw!(fig[1,1], dt * visual(Scatter, strokewidth=0.5), palettes=(;color=cl))

##



bugs = ["s__Coprococcus_eutactus", "s__Eubacterium_eligens", "s__Ruminococcus_gnavus"]
bugdf = DataFrame()
for bug in bugs
    ab = vec(abundances(testcomm[spc, :]))
    over0 = ab .> 0
    df = DataFrame(
        age = get(testcomm, :ageMonths)[over0],
        cogScore = get(testcomm, :cogScore)[over0],
        cogQuartile = categorical(get(testcomm, :cogQuartile)[over0], levels=["lower", "middle", "upper"]),
        bug = log.(collect(ab[over0]))
    )

    df.bugname .= bug
    append!(bugdf, df)
end

pltr = data(bugdf) * mapping(:cogQuartile=>"Quartile", :bug=> "Log(relative abundance)", color=:cogQuartile, layout=:bugname)
draw!(subfig[1,1], pltr * (visual(BoxPlot, show_notch=true) + visual(Scatter)), palettes= (; color=cl))

fig

##

spc = "s__Coprococcus_eutactus"

ab = vec(abundances(testcomm[spc, :]))
over0 = ab .> 0

df = DataFrame(
    age = get(testcomm, :ageMonths)[over0],
    cogScore = get(testcomm, :cogScore)[over0],
    cogQuartile = categorical(get(testcomm, :cogQuartile)[over0], levels=["lower", "middle", "upper"]),
    bug = log.(collect(ab[over0]))
)


pltr = data(df) * mapping(:cogQuartile=>"Quartile", :bug=> "Log(relative abundance)", color=:cogQuartile) 
draw!(subfig[1,1], pltr * (visual(BoxPlot, show_notch=true) + visual(Scatter)), palettes= (; color=cl), axis = (;title=spc))

pltr = data(df) * mapping(:cogScore=>"Cognitive function score", :bug=> "Log(relative abundance)", color=:cogQuartile)
draw!(subfig[1,2], pltr * visual(Scatter), palettes= (; color=cl), axis = (;title=spc))

##

spc = "s__Eubacterium_eligens"

ab = vec(abundances(testcomm[spc, :]))
over0 = ab .> 0

df = DataFrame(
    age = get(testcomm, :ageMonths)[over0],
    cogScore = get(testcomm, :cogScore)[over0],
    cogQuartile = categorical(get(testcomm, :cogQuartile)[over0], levels=["lower", "middle", "upper"]),
    bug = log.(collect(ab[over0]))
)


pltr = data(df) * mapping(:cogQuartile=>"Quartile", :bug=> "Log(relative abundance)", color=:cogQuartile) 
draw!(subfig[2,1], pltr * (visual(BoxPlot, show_notch=true) + visual(Scatter)), palettes= (; color=cl), axis = (;title=spc))

pltr = data(df) * mapping(:cogScore=>"Cognitive function score", :bug=> "Log(relative abundance)", color=:cogQuartile)
draw!(subfig[2,2], pltr * visual(Scatter), palettes= (; color=cl), axis = (;title=spc))


##

spc = "s__Ruminococcus_gnavus"

ab = vec(abundances(testcomm[spc, :]))
over0 = ab .> 0

df = DataFrame(
    age = get(testcomm, :ageMonths)[over0],
    cogScore = get(testcomm, :cogScore)[over0],
    cogQuartile = categorical(get(testcomm, :cogQuartile)[over0], levels=["lower", "middle", "upper"]),
    bug = log.(collect(ab[over0]))
)


pltr = data(df) * mapping(:cogQuartile=>"Quartile", :bug=> "Log(relative abundance)", color=:cogQuartile) 
draw!(subfig[3,1], pltr * (visual(BoxPlot, show_notch=true) + visual(Scatter)), palettes= (; color=cl), axis = (;title=spc))

pltr = data(df) * mapping(:cogScore=>"Cognitive function score", :bug=> "Log(relative abundance)", color=:cogQuartile)
draw!(subfig[3,2], pltr * visual(Scatter), palettes= (; color=cl), axis = (;title=spc))

##

Legend(fig[2,1], [MarkerElement(color=:seagreen, marker=:circle, strokewidth=0.5),
                    MarkerElement(color=:lightgray, marker=:circle, strokewidth=0.5),
                    MarkerElement(color=:dodgerblue, marker=:circle, strokewidth=0.5)
                   ],
                   ["lower", "middle", "upper"];
                   orientation=:horizontal, tellwidth=false, tellheight=true,
                   title="Quartiles")
                   
save("figures/lms.png", fig)
fig
##

genemetab = DataFrame(
    gabasynth = map(sum, eachcol(abundances(unirefs[neuroactive["GABA synthesis"], [esmap[s] for s in samplenames(metaboverlap)]]))),
    gabadegr  = map(sum, eachcol(abundances(unirefs[neuroactive["GABA degradation"], [esmap[s] for s in samplenames(metaboverlap)]]))),
    gabagut   = log.(vec(abundances(metaboverlap[gaba, :]))),
    glutsynth = map(sum, eachcol(abundances(unirefs[neuroactive["Glutamate synthesis"], [esmap[s] for s in samplenames(metaboverlap)]]))),
    glutdegr  = map(sum, eachcol(abundances(unirefs[neuroactive["Glutamate degradation"], [esmap[s] for s in samplenames(metaboverlap)]]))),
    glutgut   = log.(vec(abundances(metaboverlap[glutamate, :]))),
    mc        = get(metaboverlap, :Mother_Child))

gabasynthlm = lm(@formula(gabagut ~ gabasynth + mc), genemetab)
gabadegrlm = lm(@formula(gabagut ~ gabadegr + mc), genemetab)
glutsynthlm = lm(@formula(glutgut ~ glutsynth + mc), genemetab)
glutdegrlm = lm(@formula(glutgut ~ glutdegr + mc), genemetab)

##

pred = DataFrame(gabasynth = repeat(range(extrema(genemetab.gabasynth)..., length=50), outer=2),
                 gabadegr  = repeat(range(extrema(genemetab.gabadegr)..., length=50), outer=2),
                 gabagut   = zeros(100),
                 glutsynth = repeat(range(extrema(genemetab.glutsynth)..., length=50), outer=2),
                 glutdegr  = repeat(range(extrema(genemetab.glutdegr)..., length=50), outer=2),
                 glutgut   = zeros(100),
                 mc        = repeat(["M", "C"], inner=50))

pred.gabasynth_pred = predict(gabasynthlm, pred)
pred.gabadegr_pred  = predict(gabadegrlm, pred)
pred.glutsynth_pred = predict(glutsynthlm, pred)
pred.glutdegr_pred  = predict(glutdegrlm, pred)

predgrp = groupby(pred, :mc)

##

fig = Figure(resolution=(800,800))
ax1 = Axis(fig[1,1], xlabel="GABA Synthesis (RPKM)", ylabel="GABA (log abundance)")
ax2 = Axis(fig[2,1], xlabel="GABA Degradation (RPKM)", ylabel="GABA (log abundance)")

scatter!(ax1, genemetab.gabasynth, genemetab.gabagut,
              color=[ismissing(x) ? :gray : x == "M" ? :dodgerblue : :orange for x in genemetab.mc])
              
for grp in groupby(pred, :mc)
    c = first(grp.mc) == "M" ? :dodgerblue : :orange
    lines!(ax1, grp.gabasynth, grp.gabasynth_pred, color=c)
end
            

scatter!(ax2, genemetab.gabadegr, genemetab.gabagut,
              color=[ismissing(x) ? :gray : x == "M" ? :dodgerblue : :orange for x in genemetab.mc])


for grp in groupby(pred, :mc)
    c = first(grp.mc) == "M" ? :dodgerblue : :orange
    lines!(ax2, grp.gabadegr, grp.gabadegr_pred, color=c)
end  

##

m1 = MarkerElement(color=:dodgerblue, marker=:circle)
l1 = LineElement(color=:dodgerblue)
m2 = MarkerElement(color=:orange, marker=:circle)
l2 = LineElement(color=:orange)
Legend(fig[1:2, 2], [[m1, l1], [m2, l2]], ["moms", "kids"])

save("figures/gaba_genes_metabolites.png", fig)
fig


