# # Note - this execution takes a *very* long time. Probably need to re-think it.
# # Commenting out for now

using Resonance
using SparseArrays

genefamilies = String[]
for (root, dirs, files) in walkdir("/grace/echo/analysis/biobakery3/batch001")
    contains(root, "humann") || continue
    contains(root, "main") || continue
    contains(root, "links") && continue
    filter!(f-> contains(f, r"^FG.+genefamilies"), files)
    append!(genefamilies, joinpath.(Ref(root), files))
end

function builddf(files)
    df = DataFrame(feature=String[], value=Float64[], sample=String[])
    for f in files
        sdf = CSV.read(f, DataFrame; header=["feature", "value"], skipto=2)
        subset!(sdf, "feature"=> ByRow(f-> !contains(f, '|'))) # skip stratified features
        samplename = replace(splitext(basename(f))[1], r"_S\d+_genefamilies" => "")
        sdf.sample .= samplename
        append!(df, sdf)
    end
    return unique(df, [:sample, :feature])
end

function matrixize(longdf)
    fs = Set(longdf.feature)
    ss = Set(longdf.sample)
    features = Dict(f => i for (i, f) in enumerate(fs)) 
    samples = Dict(s => i for (i, s) in enumerate(ss))

    m = spzeros(length(fs), length(ss))
    for row in eachrow(longdf)
        rowidx = features[row.feature]
        colidx = samples[row.sample]
        m[rowidx, colidx] = row.value
    end
    return features, samples, m
end

function matrixize2(longdf)
    fs = Set(longdf.feature)
    ss = Set(longdf.sample)
    features = Dict(f => i for (i, f) in enumerate(fs)) 
    samples = Dict(s => i for (i, s) in enumerate(ss))

    longdf.fidx = [features[f] for f in longdf.feature]
    longdf.sidx = [samples[s] for s in longdf.sample]

    return features, samples, sparse(longdf.fidx, longdf.sidx, longdf.value)
end

# for gff in genefamilies
#     gf = humann_profile(gff; sample=replace(basename(gff), r"_S\d+_genefamilies.tsv"=>""))
#     union!(feature_set, features(gf))
#     push!(profiles, gf)
# end

# fidx = Dict(f=> i for (i,f) in enumerate(feature_set))

# mat = spzeros(length(feature_set), length(profiles))

# for (i, p) in enumerate(profiles)
#     for (j, f) in enumerate(features(p))
#     end
# end

# func = CommunityProfile(abundances(func)[:, idx], features(func), sns[idx])

# CSV.write("data/genefamilies.csv", func)