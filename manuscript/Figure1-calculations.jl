using Resonance
using FileIO
using CairoMakie # for plotting
using Distances
using MultivariateStats
using CategoricalArrays

#-

mdata = Resonance.load(Metadata())


species = Resonance.load(TaxonomicProfiles(); timepoint_metadata = mdata)
unirefs = Resonance.load(UnirefProfiles(); timepoint_metadata = mdata) # this can take a bit
# # unirefs = filter(!hastaxon, unirefs) # don't use species stratification for summaries

ecs = Resonance.load(ECProfiles(); timepoint_metadata = mdata)
ecs = ecs[:, findall(!ismissing, get(ecs, :ageMonths))]
ecs = filter(!hastaxon, ecs)

kos = Resonance.load(KOProfiles(); timepoint_metadata = mdata)
kos = kos[:, findall(!ismissing, get(kos, :ageMonths))]
kos = filter(!hastaxon, kos)

# metabolites = Resonance.load(MetabolicProfiles(); timepoint_metadata = mdata)
brain = Resonance.load(Neuroimaging(), timepoint_metadata = mdata, samplefield="omni")


# @assert all(samplenames(species) .== samplenames(unirefs))
@assert all(samplenames(species) .== samplenames(ecs) .== samplenames(kos))


#-

spedm = Microbiome.braycurtis(species)
CSV.write(scratchfiles("spedm.csv"), Tables.table(spedm))
# unidm = Microbiome.braycurtis(unirefs)
# CSV.write(scratchfiles("unidm.csv"), Tables.table(unidm))
ecsdm = Microbiome.braycurtis(ecs)
CSV.write(scratchfiles("ecsdm.csv"), Tables.table(ecsdm))
kosdm = Microbiome.braycurtis(kos)
CSV.write(scratchfiles("kosdm.csv"), Tables.table(kosdm))
# metdm = Microbiome.braycurtis(metabolites)
# CSV.write(scratchfiles("metdm.csv"), Tables.table(metdm))
brndm = Microbiome.braycurtis(brain)
CSV.write(scratchfiles("brndm.csv"), Tables.table(brndm))


#-

commlabels = ["taxa", "UniRef90s", "ECs", "KOs"]
mdlabels = ["Cog. score", "Age", "Sex", "Maternal Edu."]

p = permanovas([spedm[uidx, uidx], unidm[uidx, uidx], ecsdm[uidx, uidx], kosdm[uidx, uidx]], [
                        get(species, :cogScorePercentile)[uidx], 
                        get(species, :ageMonths)[uidx], 
                        get(species, :sex)[uidx], 
                        get(species, :education)[uidx]
            ]; commlabels, mdlabels
)
# p2 = permanovas(metdm, [
#                         get(metabolites, :cogScorePercentile), 
#                         get(metabolites, :ageMonths), 
#                         get(metabolites, :sex), 
#                         get(metabolites, :education)
#             ]; mdlabels
# )
# p2.label .= "metabolites"
# append!(p, p2)

p3 = permanovas(brndm[buidx, buidx], [
                        get(brain, :cogScorePercentile)[buidx], 
                        get(brain, :ageMonths)[buidx], 
                        get(brain, :sex)[buidx], 
                        get(brain, :education)[buidx]
            ]; mdlabels
)
p3.label .= "neuroimg"
append!(p, p3)
CSV.write(scratchfiles("permanovas_all.csv"), p)

#-


# idx = findall(<(6), get(species, :ageMonths))
# p = permanovas([spedm[idx, idx], unidm[idx, idx]], [
#                         get(species, :cogScorePercentile)[idx], 
#                         get(species, :ageMonths)[idx], 
#                         get(species, :sex)[idx], 
#                         get(species, :education)[idx]
#             ]; commlabels=["taxa", "UniRef90s"], mdlabels
# )

# idx = findall(<(6), get(metabolites, :ageMonths))
# p2 = permanovas(metdm[idx,idx], [
#                         get(metabolites, :cogScorePercentile)[idx], 
#                         get(metabolites, :ageMonths)[idx], 
#                         get(metabolites, :sex)[idx], 
#                         get(metabolites, :education)[idx]
#             ]; mdlabels
# )
# p2.label .= "metabolites"

# append!(p, p2)

idx = findall(<(6), get(brain, :ageMonths))
p3 = permanovas(brndm[idx, idx], [
                        get(brain, :cogScorePercentile)[idx], 
                        get(brain, :ageMonths)[idx], 
                        get(brain, :sex)[idx], 
                        get(brain, :education)[idx]
            ]; mdlabels
)
p3.label .= "neuroimg"
append!(p, p3)

CSV.write(scratchfiles("permanovas_u6mo.csv"), p)


#-

idx = findall(>(18), get(species, :ageMonths))
p = permanovas([spedm[idx, idx], unidm[idx, idx]], [
                        get(species, :cogScorePercentile)[idx], 
                        get(species, :ageMonths)[idx], 
                        get(species, :sex)[idx], 
                        get(species, :education)[idx]
            ]; commlabels=["taxa", "UniRef90s"], mdlabels
)
bidx = findall(>(18), get(brain, :ageMonths))
p3 = permanovas(brndm[bidx, bidx], [
                        get(brain, :cogScorePercentile)[bidx], 
                        get(brain, :ageMonths)[bidx], 
                        get(brain, :sex)[bidx], 
                        get(brain, :education)[bidx]
            ]; mdlabels
)
p3.label .= "neuroimg"
append!(p, p3)
CSV.write(scratchfiles("permanovas_o18mo.csv"), p)

#- Mantel -#

mdf = mantel([spedm, unidm, ecsdm, kosdm]; commlabels)

# (ol1, ol2) = stp_overlap(
#         collect(zip(get(species, :subject), get(species, :timepoint))),
#         collect(zip(get(metabolites, :subject), get(metabolites, :timepoint)))
# )
# m2 = DataFrame()
# for (i, dm1) in enumerate([spedm, unidm, ecsdm, kosdm])
#     m, p = mantel(dm1[ol1, ol1], metdm[ol2, ol2])
#     push!(m2, (; stat=m, pvalue=p, thing1=commlabels[i], thing2="metabolites"))
# end
# append!(mdf, m2)

(ol3, ol4) = Resonance.stp_overlap(
        collect(zip(get(species, :subject), get(species, :timepoint))),
        collect(zip(get(brain, :subject), get(brain, :timepoint)))
)
m3 = DataFrame()
for (i, dm1) in enumerate([spedm, unidm, ecsdm, kosdm])
    m, p = mantel(dm1[ol3, ol3], brndm[ol4, ol4])
    push!(m3, (; stat=m, pvalue=p, thing1=commlabels[i], thing2="neuroimg"))
end
append!(mdf, m3)

# (ol5, ol6) = stp_overlap(
#         collect(zip(get(metabolites, :subject), get(metabolites, :timepoint))),
#         collect(zip(get(brain, :subject), get(brain, :timepoint))),
# )

# m, p = mantel(metdm[ol5, ol5], brndm[ol6, ol6])
# push!(mdf, (; stat=m, pvalue=p, thing1="metabolites", thing2="neuroimg"))        

CSV.write(scratchfiles("mantel_all.csv"), mdf)

#- Under 6mo -#

speidx = get(species, :ageMonths) .< 6
# metidx = get(metabolites, :ageMonths) .< 6
brnidx = get(brain, :ageMonths) .< 6


speu6dm = spedm[speidx, speidx]
uniu6dm = unidm[speidx, speidx]
ecsu6dm = ecsdm[speidx, speidx]
kosu6dm = kosdm[speidx, speidx]
# metu6dm = metdm[metidx, metidx]
brnu6dm = brndm[brnidx, brnidx]


mdf = mantel([speu6dm, uniu6dm, ecsu6dm, kosu6dm]; commlabels)

# (ol1, ol2) = stp_overlap(
#         collect(zip(get(species, :subject)[speidx],
#                     get(species, :timepoint)[speidx])
#                 ),
#         collect(zip(get(metabolites, :subject)[metidx],
#                     get(metabolites, :timepoint)[metidx])
#                 )
# )



# m2 = DataFrame()
# for (i, dm1) in enumerate([speu6dm, uniu6dm, ecsu6dm, kosu6dm])
#     m, p = mantel(dm1[ol1, ol1], metu6dm[ol2, ol2])
#     push!(m2, (; stat=m, pvalue=p, thing1=commlabels[i], thing2="metabolites"))
# end
# append!(mdf, m2)


(ol3, ol4) = Resonance.stp_overlap(
        collect(zip(get(species, :subject)[speidx],
                    get(species, :timepoint)[speidx])
                ),
        collect(zip(get(brain, :subject)[brnidx],
                    get(brain, :timepoint)[brnidx])
                )
)
m3 = DataFrame()
for (i, dm1) in enumerate([speu6dm, uniu6dm, ecsu6dm, kosu6dm])
    m, p = mantel(dm1[ol3, ol3], brnu6dm[ol4, ol4])
    push!(m3, (; stat=m, pvalue=p, thing1=commlabels[i], thing2="neuroimg"))
end
append!(mdf, m3)

# (ol5, ol6) = stp_overlap(
#         collect(zip(get(metabolites, :subject)[metidx],
#                     get(metabolites, :timepoint)[metidx])
#                 ),
#         collect(zip(get(brain, :subject)[brnidx],
#                     get(brain, :timepoint)[brnidx])
#                 ),
# )

# m, p = mantel(metu6dm[ol5, ol5], brnu6dm[ol6, ol6])
# push!(mdf, (; stat=m, pvalue=p, thing1="metabolites", thing2="neuroimg"))        

CSV.write(scratchfiles("mantel_u6.csv"), mdf)

#- Over 18mo -#

speidx = get(species, :ageMonths) .> 18
brnidx = get(brain, :ageMonths) .> 18

speo18dm = spedm[speidx, speidx]
unio18dm = unidm[speidx, speidx]
ecso18dm = ecsdm[speidx, speidx]
koso18dm = kosdm[speidx, speidx]
brno18dm = brndm[brnidx, brnidx]


mdf = mantel([speo18dm, unio18dm, ecso18dm, koso18dm]; commlabels)


(ol3, ol4) = Resonance.stp_overlap(
        collect(zip(get(species, :subject)[speidx],
                    get(species, :timepoint)[speidx])
                ),
        collect(zip(get(brain, :subject)[brnidx],
                    get(brain, :timepoint)[brnidx])
                )
)
m3 = DataFrame()
for (i, dm1) in enumerate([speo18dm, unio18dm, ecso18dm, koso18dm])
    m, p = mantel(dm1[ol3, ol3], brno18dm[ol4, ol4])
    push!(m3, (; stat=m, pvalue=p, thing1=commlabels[i], thing2="neuroimg"))
end
append!(mdf, m3)   

CSV.write(scratchfiles("mantel_o18.csv"), mdf)
