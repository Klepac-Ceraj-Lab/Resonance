using Resonance

species = metaphlan_profiles(filter(f-> contains(f, "profile"), 
                                        readdir("/grace/echo/analysis/biobakery3/links/metaphlan", join=true)),
                                        :species)
                                        
sns = map(samplenames(species)) do s
    replace(s, r"_S\d+_profile"=>"")
end

idx = unique([findfirst(==(x), sns) for x in sns])

species = CommunityProfile(abundances(species)[:, idx], features(species), map(samplenames(species)[idx]) do s
    MicrobiomeSample(replace(s, r"_S\d+_profile"=>""))
end)

CSV.write("data/species.csv", species)