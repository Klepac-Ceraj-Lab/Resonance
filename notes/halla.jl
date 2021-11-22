using Resonance

metadata = CSV.read("data/wrangled.csv", DataFrame)
@rsubset!(metadata, !ismissing(:subject), 
                    !ismissing(:sample))

@rsubset!(metadata, !startswith(:sample, "C"),
                    !startswith(:sample, "z"))


volumes = metadata[:, r"^(left|right|sample)"i][:, 2:end]
@rsubset!(volumes, !ismissing(Symbol("Left-Lateral-Ventricle")))

met = metaphlan_profiles(filter(f-> contains(f, "profile"), readdir("/grace/echo/analysis/biobakery3/links/metaphlan", join=true)), :species)

sns = map(samplenames(met)) do s
    replace(s, r"_S\d+_profile"=>"")
end

idx = unique([findfirst(==(x), sns) for x in sns])

met = CommunityProfile(abundances(met)[:, idx], features(met), map(samplenames(met)[idx]) do s
    MicrobiomeSample(replace(s, r"_S\d+_profile"=>""))
end)

overlap = sort(samplenames(met) ∩ volumes.sample)

met = met[:, overlap]

volumes = @chain volumes begin
    @rsubset(:sample in overlap)
    @orderby(:sample)
end


CSV.write("data/halla_volumes.tsv", volumes, delim='\t')
CSV.write("data/halla_species.tsv", met, delim='\t')