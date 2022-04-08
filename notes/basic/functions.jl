using Resonance
using SparseArrays

genefamilies = String[]
for (root, dirs, files) in walkdir("/grace/echo/analysis/biobakery3/")
    contains(root, "humann") || continue
    contains(root, "main") || continue
    contains(root, "links") && continue
    filter!(f-> contains(f, r"^FG.+genefamilies"), files)
    append!(genefamilies, joinpath.(Ref(root), files))
end

profiles = CommunityProfile[]
feature_set = Set(GeneFunction[])

for gff in genefamilies
    gf = humann_profile(gff; sample=replace(basename(gff), r"_S\d+_genefamilies.tsv"=>""))
    union!(feature_set, features(gf))
    push!(profiles, gf)
end

fidx = Dict(f=> i for (i,f) in enumerate(feature_set))

mat = spzeros(length(feature_set), length(profiles))

for (i, p) in enumerate(profiles)
    for (j, f) in enumerate(features(p))
    end
end

func = CommunityProfile(abundances(func)[:, idx], features(func), sns[idx])

CSV.write("data/genefamilies.csv", func)