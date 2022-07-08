# # Taxonomic profiles
#
# This notebook just loads in all of the taxonomic profiles
# from the biobakery analysis folder (`$ANALYSIS_FILES`),
# subsets on the species,
# and then writes a new table to the data folder.
#
# This only needs to be run once each time new samples are added.

using Resonance

profiles = String[]
for (root, dirs, files) in walkdir(ENV["ANALYSIS_FILES"])
    contains(root, "metaphlan") || continue
    contains(root, "links") && continue
    filter!(f-> contains(f, r"^FG.+profile"), files)
    append!(profiles, joinpath.(Ref(root), files))
end

species = metaphlan_profiles(profiles, :species)
                                        
sns = map(samplenames(species)) do s
    replace(s, r"_S\d+_profile"=>"")
end

idx = unique([findfirst(==(x), sns) for x in sns])

species = CommunityProfile(abundances(species)[:, idx], features(species), map(samplenames(species)[idx]) do s
    MicrobiomeSample(replace(s, r"_S\d+_profile"=>""))
end)

CSV.write(datafiles("species.csv"), species)