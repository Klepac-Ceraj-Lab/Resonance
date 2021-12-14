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

gaba = findfirst(f-> commonname(f) == "gamma-Aminobutyric acid", features(metabolites))
glutamate = findfirst(f-> commonname(f) == "Glutamic acid", features(metabolites))

momsidx = get(metabolites, :Mother_Child) .=== "M"
kidsidx = get(metabolites, :Mother_Child) .=== "C"

##

fig = Figure(resolution=(1200,1200))
ax1 = Axis(fig[1:2,1:2], xlabel="GABA (log)", ylabel="Glutamate (log)")
ax2 = Axis(fig[0,1:2], height=200)
ax3 = Axis(fig[2:3, 3], width=200)


scmom = scatter!(ax1, vec(metabolites[gaba, momsidx]), vec(metabolites[glutamate,momsidx]))
hist!(ax2, vec(metabolites[gaba, momsidx]))
hist!(ax3, vec(metabolites[glutamate,momsidx]), direction=:x)

sckid = scatter!(ax1, vec(metabolites[gaba,kidsidx]), vec(metabolites[glutamate, kidsidx]))
hist!(ax2, gaba[kidsidx])
hist!(ax3, vec(metabolites[glutamate, kidsidx]), direction=:x)



save("figures/gaba-glutamate.png", fig)
fig


##

fig = Figure(resolution=(1200,1200))
ax1 = Axis(fig[1:2,1:2], xlabel="GABA (log)", ylabel="Glutamate (log)")
ax2 = Axis(fig[0,1:2], height=200)
ax3 = Axis(fig[2:3, 3], width=200)


scmom = scatter!(ax1, vec(metabolites[gaba,momsidx], vec(metabolites[glutamate, momsidx]), color=:gray))
histmom = hist!(ax2, vec(metabolites[gaba,momsidx], color=:gray))
hist!(ax3, vec(metabolites[glutamate, momsidx]), direction=:x, color=:gray)


sckid = scatter!(ax1, vec(metabolites[gaba,kidsidx]), vec(metabolites[glutamate, kidsidx]), color = getcolors(allmeta[kidsidx, :ageMonths], (0,25)),
                 strokewidth=1)
histkid = hist!(ax2, vec(metabolites[gaba,kidsidx]))
hist!(ax3, vec(metabolites[glutamate, kidsidx]), direction=:x)



cleg = Colorbar(fig[2:3, 4], limits=(0.1, 25), colormap=:plasma, highclip=:white, lowclip=:black, label="Age in Months")
leg = Legend(fig[1,3], [histmom, histkid], ["Moms", "Kids"], tellwidth = false, tellheight = false)

save("figures/gaba-glutamate.png", fig)
fig

##


unirefs = Resonance.load_genefamilies()
neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(unirefs)))

featurenames(unirefs)