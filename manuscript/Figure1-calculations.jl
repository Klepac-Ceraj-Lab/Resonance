using Resonance
using FileIO
using CairoMakie # for plotting
using Distances
using MultivariateStats
using CategoricalArrays

#-

mdata = Resonance.load(Metadata())
mdata.edfloat = map(x-> ismissing(x) ? missing : Float64(levelcode(x)), mdata.education)

species = Resonance.load(TaxonomicProfiles(); timepoint_metadata = mdata)
unirefs = Resonance.load(UnirefProfiles(); timepoint_metadata = mdata) # this can take a bit
ecs = Resonance.load(ECProfiles(); timepoint_metadata = mdata)
kos = Resonance.load(KOProfiles(); timepoint_metadata = mdata)

# metabolites = Resonance.load(MetabolicProfiles(); timepoint_metadata = mdata)
brain = Resonance.load(Neuroimaging(), timepoint_metadata = mdata, samplefield="sample")


# @assert all(samplenames(species) .== samplenames(unirefs))
@assert all(samplenames(species) .== samplenames(ecs) .== samplenames(kos) .== samplenames(unirefs))


#-

spedm = Microbiome.braycurtis(species)
CSV.write(scratchfiles("spedm.csv"), Tables.table(spedm))
unidm = Microbiome.braycurtis(unirefs)
CSV.write(scratchfiles("unidm.csv"), Tables.table(unidm))
ecsdm = Microbiome.braycurtis(ecs)
CSV.write(scratchfiles("ecsdm.csv"), Tables.table(ecsdm))
kosdm = Microbiome.braycurtis(kos)
CSV.write(scratchfiles("kosdm.csv"), Tables.table(kosdm))
# metdm = Microbiome.braycurtis(metabolites)
# CSV.write(scratchfiles("metdm.csv"), Tables.table(metdm))
brndm = pairwise(Euclidean(), abundances(brain))
CSV.write(scratchfiles("brndm.csv"), Tables.table(brndm))


#-

uidx = get(species, :filter_00to120)
buidx = get(brain, :filter_00to120)
commlabels = ["taxa", "UniRef90s", "ECs", "KOs"]
mdlabels = ["Cog. score", "Age", "Sex", "Maternal Edu.", "Race"]

p = permanovas([spedm[uidx, uidx], unidm[uidx, uidx], ecsdm[uidx, uidx], kosdm[uidx, uidx]], [
                        get(species, :cogScore)[uidx], 
                        get(species, :ageMonths)[uidx], 
                        get(species, :sex)[uidx], 
                        get(species, :edfloat)[uidx],
                        get(species, :race)[uidx]
            ]; commlabels, mdlabels
)
# p2 = permanovas(metdm, [
#                         get(metabolites, :cogScore), 
#                         get(metabolites, :ageMonths), 
#                         get(metabolites, :sex), 
#                         get(metabolites, :edfloat)
#             ]; mdlabels
# )
# p2.label .= "metabolites"
# append!(p, p2)

p3 = permanovas(brndm[buidx, buidx], [
                        get(brain, :cogScore)[buidx], 
                        get(brain, :ageMonths)[buidx], 
                        get(brain, :sex)[buidx], 
                        get(brain, :edfloat)[buidx],
                        get(brain, :race)[buidx]
            ]; mdlabels
)
p3.label .= "neuroimg"
append!(p, p3)
CSV.write(scratchfiles("permanovas_00t18to1200.csv"), p)

#-


idx = get(species, :filter_00to06)
p = permanovas([spedm[idx, idx], unidm[idx, idx]], [
                        get(species, :cogScore)[idx], 
                        get(species, :ageMonths)[idx], 
                        get(species, :sex)[idx], 
                        get(species, :edfloat)[idx],
                        get(species, :race)[idx],
            ]; commlabels=["taxa", "UniRef90s"], mdlabels
)

# idx = findall(<(6), get(metabolites, :ageMonths))
# p2 = permanovas(metdm[idx,idx], [
#                         get(metabolites, :cogScore)[idx], 
#                         get(metabolites, :ageMonths)[idx], 
#                         get(metabolites, :sex)[idx], 
#                         get(metabolites, :edfloat)[idx]
#             ]; mdlabels
# )
# p2.label .= "metabolites"

# append!(p, p2)

idx = get(brain, :filter_00to06)
p3 = permanovas(brndm[idx, idx], [
                        get(brain, :cogScore)[idx], 
                        get(brain, :ageMonths)[idx], 
                        get(brain, :sex)[idx], 
                        get(brain, :edfloat)[idx],
                        get(brain, :race)[idx]
            ]; mdlabels
)
p3.label .= "neuroimg"
append!(p, p3)

CSV.write(scratchfiles("permanovas_00to06.csv"), p)


#-

idx = get(species, :filter_18to120)
p = permanovas([spedm[idx, idx], unidm[idx, idx]], [
                        get(species, :cogScore)[idx], 
                        get(species, :ageMonths)[idx], 
                        get(species, :sex)[idx], 
                        get(species, :edfloat)[idx],
                        get(species, :race)[idx]
            ]; commlabels=["taxa", "UniRef90s"], mdlabels
)
bidx =  get(brain, :filter_18to120)
p3 = permanovas(brndm[bidx, bidx], [
                        get(brain, :cogScore)[bidx], 
                        get(brain, :ageMonths)[bidx], 
                        get(brain, :sex)[bidx], 
                        get(brain, :edfloat)[bidx],
                        get(brain, :race)[bidx]
            ]; mdlabels
)
p3.label .= "neuroimg"
append!(p, p3)
CSV.write(scratchfiles("permanovas_18to120.csv"), p)

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

CSV.write(scratchfiles("mantel_00to120.csv"), mdf)

#- Under 6mo -#

speidx = get(species, :filter_00to06)
# metidx = get(metabolites, :filter_00to06)
brnidx = get(brain, :filter_00to06)


spe00to06dm = spedm[speidx, speidx]
uni00to06dm = unidm[speidx, speidx]
ecs00to06dm = ecsdm[speidx, speidx]
kos00to06dm = kosdm[speidx, speidx]
# met00to06dm = metdm[metidx, metidx]
brn00to06dm = brndm[brnidx, brnidx]


mdf = mantel([spe00to06dm, uni00to06dm, ecs00to06dm, kos00to06dm]; commlabels)

# (ol1, ol2) = stp_overlap(
#         collect(zip(get(species, :subject)[speidx],
#                     get(species, :timepoint)[speidx])
#                 ),
#         collect(zip(get(metabolites, :subject)[metidx],
#                     get(metabolites, :timepoint)[metidx])
#                 )
# )



# m2 = DataFrame()
# for (i, dm1) in enumerate([spe00to06dm, uni00to06dm, ecs00to06dm, kos00to06dm])
#     m, p = mantel(dm1[ol1, ol1], met00to06dm[ol2, ol2])
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
for (i, dm1) in enumerate([spe00to06dm, uni00to06dm, ecs00to06dm, kos00to06dm])
    m, p = mantel(dm1[ol3, ol3], brn00to06dm[ol4, ol4])
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

# m, p = mantel(met00to06dm[ol5, ol5], brn00to06dm[ol6, ol6])
# push!(mdf, (; stat=m, pvalue=p, thing1="metabolites", thing2="neuroimg"))        

CSV.write(scratchfiles("mantel_00to06.csv"), mdf)

#- Over 18mo -#

speidx = get(species, :filter_18to120)
brnidx = get(brain, :filter_18to120)

spe18to120dm = spedm[speidx, speidx]
uni18to120dm = unidm[speidx, speidx]
ecs18to120dm = ecsdm[speidx, speidx]
kos18to120dm = kosdm[speidx, speidx]
brn18to120dm = brndm[brnidx, brnidx]


mdf = mantel([spe18to120dm, uni18to120dm, ecs18to120dm, kos18to120dm]; commlabels)


(ol3, ol4) = Resonance.stp_overlap(
        collect(zip(get(species, :subject)[speidx],
                    get(species, :timepoint)[speidx])
                ),
        collect(zip(get(brain, :subject)[brnidx],
                    get(brain, :timepoint)[brnidx])
                )
)
m3 = DataFrame()
for (i, dm1) in enumerate([spe18to120dm, uni18to120dm, ecs18to120dm, kos18to120dm])
    m, p = mantel(dm1[ol3, ol3], brn18to120dm[ol4, ol4])
    push!(m3, (; stat=m, pvalue=p, thing1=commlabels[i], thing2="neuroimg"))
end
append!(mdf, m3)   

CSV.write(scratchfiles("mantel_18to120.csv"), mdf)
