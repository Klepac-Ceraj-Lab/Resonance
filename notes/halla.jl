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

met = metaphlan_profiles(filter(f-> contains(f, "profile"), readdir("/grace/echo/analysis/biobakery3/links/metaphlan", join=true)), :species)

sns = map(samplenames(met)) do s
    replace(s, r"_S\d+_profile"=>"")
end

idx = unique([findfirst(==(x), sns) for x in sns])

met = CommunityProfile(abundances(met)[:, idx], features(met), map(samplenames(met)[idx]) do s
    MicrobiomeSample(replace(s, r"_S\d+_profile"=>""))
end)

overlap = sort(samplenames(met) ∩ volumes.sample)

brainmet = met[:, overlap]

volumes = @chain volumes begin
    @rsubset(:sample in overlap)
    @orderby(:sample)
end

permutedims(volumes, "sample")

CSV.write("data/halla_volumes.tsv", permutedims(volumes, :sample), delim='\t')
CSV.write("data/halla_species.tsv", brainmet, delim='\t')

##

# Note: had to delete empty "C1162_3E_1A" "M1162_3E_1A" columns from C8-pos (they were duplicated)
#
# Note: "M0855_1E_1A" missing from HILIC pos and neg


c18neg = load_metabolites("data/21_0924_VKC_C18-neg_results.xlsx")
c8pos = load_metabolites("data/21_0924_VKC_C8-pos_results.xlsx")
hneg = load_metabolites("data/21_0924_VKC_HILIC-neg_results.xlsx")
hpos = load_metabolites("data/21_0924_VKC_HILIC-pos_results.xlsx")

for df in [c18neg, c8pos, hneg, hpos]
    select!(df, Not(r"^PREF"))
end

c8pos = c8pos[:, names(c18neg)]
hneg = hneg[:, names(c18neg)]
hpos = hpos[:, names(c18neg)]

metabs = vcat(c18neg, c8pos, hneg, hpos)

sample_meta = airtable_metadata()
@rsubset!(sample_meta, in(:sid_old, names(metabs)))

rename!(metabs, [row.sid_old=> row.sample for row in eachrow(sample_meta)])
CSV.write("data/halla_metabolites.tsv", select(metabs, ["Metabolite", sample_meta.sample...]), delim='\t')
CSV.write("data/halla_metab_species.tsv", met[:, sample_meta.sample], delim='\t')
