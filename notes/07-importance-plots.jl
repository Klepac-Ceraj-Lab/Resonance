using Resonance
using CairoMakie
using CategoricalArrays
using Statistics
using HypothesisTests
using MultipleTesting
using KernelDensity
using MultivariateStats
using ThreadsX

#-

mdata = Resonance.load(Metadata())

#-

taxa = Resonance.load(TaxonomicProfiles(); timepoint_metadata = mdata)
species = filter(t-> taxrank(t) == :species, taxa)

unirefs = Resonance.load(UnirefProfiles(); timepoint_metadata = mdata) # this can take a bit
unirefs = filter(!hastaxon, unirefs) # don't use species stratification for summaries

metabolites = Resonance.load(MetabolicProfiles(); timepoint_metadata=mdata)
metabolites = metabolites[:, [!ismissing(a) && a < 14 for a in get(metabolites, :ageMonths)]]

brain = Resonance.load(Neuroimaging(); timepoint_metadata=mdata)

importance = CSV.read(datafiles("concatenated_importances_fromtaxa.csv"), DataFrame)
importance.model = categorical(importance.model)

igrp = groupby(importance, :Variable)
#-

x = vec(abundances(species[r"Gordonibacter_pamelaeae", :]))
y = get(species, :cogScore)
scatter(x, y; color = [x == 0 ? (:gray, 0.3) :
                                       age < 6 ? (:purple, 1.) :
                                       age < 12 ? (:dodgerblue, 1.) :
                                       age < 18 ? (:seagreen, 1.) : 
                                       age < 24 ? (:orchid, 1.) :
                                       (:gray, 0.3) for (x,y, age) in zip(x, y, get(species, :ageMonths))]
)

#

let fig = Figure()
    ax1 = Axis(fig[1,1])
    ax2 = Axis(fig[2,1])
    ax3 = Axis(fig[1,2])
    ax4 = Axis(fig[2,2])
    lv = levels(importance.model)

    for k in keys(igrp)
        sdf1 = subset(igrp[k], "model"=> ByRow(x-> contains(lv[levelcode(x)], "concurrent")))
        sdf2 = subset(igrp[k], "model"=> ByRow(x-> contains(lv[levelcode(x)], "future")))
        sdf3 = subset(igrp[k], "model"=> ByRow(x-> contains(lv[levelcode(x)], "regression")))
        sdf4 = subset(igrp[k], "model"=> ByRow(x-> contains(lv[levelcode(x)], "classification")))

        scatter!(ax1, levelcode.(sdf1.model), sdf1.Importance; color = (:gray, 0.3))
        lines!(ax1, levelcode.(sdf1.model), sdf1.Importance; color = (:gray, 0.3))
        scatter!(ax2, levelcode.(sdf2.model), sdf2.Importance; color = (:gray, 0.3))
        lines!(ax2, levelcode.(sdf2.model), sdf2.Importance; color = (:gray, 0.3))
        scatter!(ax3, levelcode.(sdf3.model), sdf3.Importance; color = (:gray, 0.3))
        lines!(ax3, levelcode.(sdf3.model), sdf3.Importance; color = (:gray, 0.3))
        scatter!(ax4, levelcode.(sdf4.model), sdf4.Importance; color = (:gray, 0.3))
        lines!(ax4, levelcode.(sdf4.model), sdf4.Importance; color = (:gray, 0.3))
    end
    fig
end