using Resonance
using DataFramesMeta
using CairoMakie

allmeta = CSV.read("data/wrangled.csv", DataFrame)
metabolites = CSV.read("data/metabolites.csv", DataFrame)
ms = [Resonance.Metabolite(row[:uid], row[:Metabolite], row[:MZ], row[:RT]) for row in eachrow(metabolites)]
metabolites = CommunityProfile(Matrix(metabolites[!, 9:end]), ms, MicrobiomeSample.(names(metabolites)[9:end]))
set!(metabolites, allmeta)

species = CSV.read("data/species.csv", DataFrame)
species = CommunityProfile(Matrix(species[!, 2:end]), Taxon.(species[!, 1]), MicrobiomeSample.(names(species)[2:end]))
set!(species, allmeta)


## 


met_pcoa = pcoa(metabolites)
spec_pcoa = pcoa(species)

## 
function getcolors(vals, clip=(0.1, 30); highclip=Makie.to_color(:white), lowclip=Makie.to_color(:black))
    kidscolor = Makie.to_colormap(:plasma)

    map(vals) do val
        if ismissing(val)
            return Makie.to_color(:gray)
        elseif val < clip[1]
            return lowclip
        elseif val > clip[2]
            return highclip
        else
            Makie.interpolated_getindex(kidscolor, val, clip)
        end
    end
end

##

fig = Figure(resolution=(800, 800))
ax1 = Axis(fig[1,1], xlabel="MDS1 ($(round(varexplained(met_pcoa)[1] * 100, digits=2))%)",
                     ylabel="MDS1 ($(round(varexplained(met_pcoa)[2] * 100, digits=2))%)")

scatter!(ax1, loadings(met_pcoa)[:,1], loadings(met_pcoa)[:,2],
        color=[ismissing(x) ? :grey : x == "M" ? :dodgerblue : :orange for x in get(metabolites, :Mother_Child)])

leg = Legend(fig[1,2], [MarkerElement(color=:orange, marker=:circle),
                        MarkerElement(color=:dodgerblue, marker=:circle)],
                        ["Kids", "Moms"])

ax2 = Axis(fig[2,1],
            xlabel="MDS1 ($(round(varexplained(met_pcoa)[1] * 100, digits=2))%)",
            ylabel="MDS1 ($(round(varexplained(met_pcoa)[2] * 100, digits=2))%)"
)

scatter!(ax2, loadings(met_pcoa)[:,1], loadings(met_pcoa)[:,2],
         color = getcolors(get(metabolites, :ageMonths), (0,25)))
cleg = Colorbar(fig[2, 2], limits=(0.1, 25), colormap=:plasma, highclip=:white, lowclip=:black, label="Age in Months")
save("figures/metabolites_pcoa.png", fig)

fig

##

fig = Figure(resolution=(800, 800))
ax1 = Axis(fig[1,1], xlabel="MDS1 ($(round(varexplained(spec_pcoa)[1] * 100, digits=2))%)",
                     ylabel="MDS1 ($(round(varexplained(spec_pcoa)[2] * 100, digits=2))%)")
scatter!(ax1, loadings(spec_pcoa)[:,1], loadings(spec_pcoa)[:,2],
        color=[ismissing(x) ? :grey : x == "M" ? :dodgerblue : :orange for x in get(species, :Mother_Child)])

leg = Legend(fig[1,2], [MarkerElement(color=:orange, marker=:circle),
                        MarkerElement(color=:dodgerblue, marker=:circle)],
                        ["Kids", "Moms"])

ax2 = Axis(fig[2,1],
            xlabel="MDS1 ($(round(varexplained(spec_pcoa)[1] * 100, digits=2))%)",
            ylabel="MDS1 ($(round(varexplained(spec_pcoa)[2] * 100, digits=2))%)"
)

scatter!(ax2, loadings(spec_pcoa)[:,1], loadings(spec_pcoa)[:,2],
         color = getcolors(get(species, :ageMonths) ./ 12, (0,10))
)
cleg = Colorbar(fig[2, 2], limits=(0.1, 10), colormap=:plasma, highclip=:white, lowclip=:black, label="Age in Years")
save("figures/species_pcoa.png", fig)

fig

##

gaba = findfirst(f-> commonname(f) === "gamma-Aminobutyric acid", features(metabolites))
glutamate = findfirst(f-> commonname(f) === "Glutamic acid", features(metabolites))

momsidx = string.(get(metabolites, :Mother_Child)) .=== "M"
kidsidx = string.(get(metabolites, :Mother_Child)) .=== "C"

##

fig = Figure(resolution=(1200,1200))
ax1 = Axis(fig[1:2,1:2], xlabel="GABA (log)", ylabel="Glutamate (log)")
ax2 = Axis(fig[0,1:2], height=200)
ax3 = Axis(fig[2:3, 3], width=200)


scmom = scatter!(ax1, log.(vec(abundances(metabolites[gaba, momsidx]))), log.(vec(abundances(metabolites[glutamate,momsidx]))))
hist!(ax2, log.(vec(metabolites[gaba, momsidx] |> abundances)))
hist!(ax3, log.(vec(metabolites[glutamate, momsidx] |> abundances)), direction=:x)

sckid = scatter!(ax1, log.(vec(metabolites[gaba,kidsidx] |> abundances)), log.(vec(metabolites[glutamate, kidsidx] |> abundances)))
hist!(ax2, log.(vec(metabolites[gaba, kidsidx] |> abundances)))
hist!(ax3, log.(vec(metabolites[glutamate, kidsidx] |> abundances)), direction=:x)


leg = Legend(fig[1,3], [scmom, sckid], ["Moms", "Kids"], tellwidth = false, tellheight = false)

save("figures/gaba-glutamate.png", fig)
fig


##

fig = Figure(resolution=(1200,1200))
ax1 = Axis(fig[1:2,1:2], xlabel="GABA (log)", ylabel="Glutamate (log)")
ax2 = Axis(fig[0,1:2], height=200)
ax3 = Axis(fig[2:3, 3], width=200)


scmom = scatter!(ax1, log.(vec(abundances(metabolites[gaba,momsidx]))), 
                      log.(vec(abundances(metabolites[glutamate, momsidx]))), color=:gray)
histmom = hist!(ax2, log.(vec(abundances(metabolites[gaba,momsidx]))), color=:gray)
hist!(ax3, log.(vec(abundances(metabolites[glutamate, momsidx]))), direction=:x, color=:gray)


sckid = scatter!(ax1, log.(vec(abundances(metabolites[gaba,kidsidx]))), 
                      log.(vec(abundances(metabolites[glutamate, kidsidx]))),
                      color = getcolors(get(metabolites[:, kidsidx], :ageMonths), (0,25)),
                      strokewidth=1)
histkid = hist!(ax2, log.(vec(abundances(metabolites[gaba,kidsidx]))))
hist!(ax3, log.(vec(abundances(metabolites[glutamate,kidsidx]))), direction=:x)



cleg = Colorbar(fig[2:3, 4], limits=(0.1, 25), colormap=:plasma, highclip=:white, lowclip=:black, label="Age in Months")
leg = Legend(fig[1,3], [histmom, histkid], ["Moms", "Kids"], tellwidth = false, tellheight = false)

save("figures/gaba-glutamate_age.png", fig)
fig

##


unirefs = Resonance.load_genefamilies()
set!(unirefs, allmeta)

neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(unirefs)))

metagrp = groupby(allmeta, [:subject, :timepoint])
esmap = Dict()
for grp in metagrp
    ss = unique(skipmissing(grp.sample))
    nrow(grp) > 1 || continue
    any(s-> contains(s, "FE"), ss) || continue
    for s in ss
        contains(s, "FE") || continue
        fgs = filter(s2-> contains(s2, "FG"), ss)
        isempty(fgs) && @info ss
        esmap[s] = isempty(fgs) ? missing : first(fgs)
    end
end

metagrp[(; subject=774, timepoint=2)]

for s in samples(metabolites)
    @show s.subject, s.timepoint
    break
end

metaboverlap = metabolites[:, findall(s-> haskey(esmap, s) && esmap[s] ∈ samplenames(unirefs), samplenames(metabolites))]

##

using GLM
using MixedModels
using AlgebraOfGraphics



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

##

##

fig = Figure(resolution=(850, 1100))
ax1 = Axis(fig[1,1], xlabel="Glutamate Synthesis (RPKM)", ylabel="Glutamate (log abundance)")
ax2 = Axis(fig[2,1], xlabel="Glutamate Degradation (RPKM)", ylabel="Glutamate (log abundance)")
ax3 = Axis(fig[3,1], xlabel="Glutamate Degradation (RPKM)", ylabel="Glutamate Glutamate Synthesis (RPKM)")


scatter!(ax1, genemetab.glutsynth, genemetab.glutgut,
              color=[ismissing(x) ? :gray : x == "M" ? :dodgerblue : :orange for x in genemetab.mc])
              
for grp in groupby(pred, :mc)
    c = first(grp.mc) == "M" ? :dodgerblue : :orange
    lines!(ax1, grp.glutsynth, grp.glutsynth_pred, color=c)
end
            

scatter!(ax2, genemetab.glutdegr, genemetab.glutgut,
              color=[ismissing(x) ? :gray : x == "M" ? :dodgerblue : :orange for x in genemetab.mc])


for grp in groupby(pred, :mc)
    c = first(grp.mc) == "M" ? :dodgerblue : :orange
    lines!(ax2, grp.glutdegr, grp.glutdegr_pred, color=c)
end  

scatter!(ax3, log.(1 .+ genemetab.glutdegr), log.(1 .+ genemetab.glutsynth),
              color=[ismissing(x) ? :gray : x == "M" ? :dodgerblue : :orange for x in genemetab.mc])
##

m1 = MarkerElement(color=:dodgerblue, marker=:circle)
l1 = LineElement(color=:dodgerblue)
m2 = MarkerElement(color=:orange, marker=:circle)
l2 = LineElement(color=:orange)
Legend(fig[4, 1], [[m1, l1], [m2, l2]], ["moms", "kids"], orientation=:horizontal, tellheight=true, tellwidth=false)

save("figures/glutamate_genes_metabolites.png", fig)
fig

##

scatter(log.(1 .+ genemetab.glutdegr), log.(1 .+ genemetab.glutsynth), genemetab.glutgut,
    axis=(; xlabel="degradataion", ylabel="synthesis", zlabel="concentration"))