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
metabs[!, :uid] = [join(p, "_") for p in zip(metabs.Method, metabs.Compound_ID)]


for row in eachrow(metabs)
    window = @view row[8:end-1]
    misscols = findall(ismissing, window)
    if !isempty(misscols)
        default = floor(Int, minimum(skipmissing(window)) / 2)
        window[misscols] .= default
    end
end
    
CSV.write("data/metabolites.csv", select(metabs, [:uid, Cols(:)]))

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