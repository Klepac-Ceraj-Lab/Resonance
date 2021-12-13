using Resonance

metadata = CSV.read("data/wrangled.csv", DataFrame)
@rsubset!(metadata, !ismissing(:subject), 
                    !ismissing(:sample),
                    )
@rsubset!(metadata, !startswith(:sample, "C"),
                    !startswith(:sample, "z"),
                    )

volumes = metadata[:, r"^(left|right|sample)"i][:, 2:end]
dropmissing!(volumes)

spec = metaphlan_profiles(filter(f-> contains(f, "profile"), readdir("/grace/echo/analysis/biobakery3/links/metaphlan", join=true)), :species)

sns = map(samplenames(spec)) do s
    replace(s, r"_S\d+_profile"=>"")
end

idx = unique([findfirst(==(x), sns) for x in sns])

spec = CommunityProfile(abundances(spec)[:, idx], features(spec), map(samplenames(spec)[idx]) do s
    MicrobiomeSample(replace(s, r"_S\d+_profile"=>""))
end)

overlap = sort(samplenames(spec) ∩ volumes.sample)

brainspec = spec[:, overlap]

volumes = @chain volumes begin
    @rsubset(:sample in overlap)
    @orderby(:sample)
end

permutedims(volumes, "sample")

CSV.write("data/halla_volumes.tsv", permutedims(volumes, :sample), delim='\t')
CSV.write("data/halla_species.tsv", brainspec, delim='\t')

##

run(`halla -x data/halla_species.tsv -y data/halla_volumes.tsv --alla -o data/halla_brain/`)

##

metabs = CSV.read("data/metabolites.csv", DataFrame)

##

sample_meta = airtable_metadata()
@rsubset!(sample_meta, in(:sid_old, names(metabs)))

rename!(metabs, Dict(row.sid_old => row.sample for row in eachrow(sample_meta)))
rename!(sample_meta, "sample"=>"etoh_sample")

sample_meta = leftjoin(sample_meta[:, [:subject, :timepoint, :etoh_sample]], metadata, on=[:subject,:timepoint])

@rsubset!(sample_meta, !ismissing(:sample) && in(:sample, samplenames(spec)))
metabspec = let cmp = spec[:, string.(sample_meta.sample)]
    CommunityProfile(abundances(cmp), features(cmp), MicrobiomeSample.(string.(sample_meta.etoh_sample)))
end

CSV.write("data/halla_metabolites.tsv", select(metabs, ["Metabolite", sample_meta.etoh_sample...]), delim='\t')
CSV.write("data/halla_metab_species.tsv", metabspec, delim='\t')


##

run(`halla -x data/halla_metab_species.tsv -y data/halla_metabolites.tsv --alla -o data/halla_metab/`)

##

using CairoMakie

fig = Figure(resolution=(1200,1200))

ax, sc = scatter(fig[1,1], volumes."Right-Cerebral-Cortex", volumes."Left-Cerebral-Cortex")

fig

CSV.write("data/volumes_desc.csv", describe(volumes))