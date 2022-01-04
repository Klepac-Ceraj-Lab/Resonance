function load_genefamilies()

    allfiles = String[]

    @info "getting files"
    for (root, dirs, files) in walkdir("/grace/echo/analysis/biobakery3/")
        contains(root, "batch") || continue
        filter!(f-> contains(f, "genefamilies") && !contains(f, "relab") && !contains(f, "zymo") && !contains(f, r"^FE\d+"), files)
        append!(allfiles, joinpath.(Ref(root), files))
    end

    samples = unique(map(f-> first(split(basename(f), '_')) |> String, allfiles))
    unique!(allfiles) do f
        first(split(basename(f), '_')) |> String
    end

    features = Set{String}()
    @info "getting features"
    for f in allfiles
        @debug f
        union!(features, Iterators.filter(x-> !contains(x, "|"), Iterators.map(x-> x[1], CSV.File(f))))
    end

    abunds = zeros(length(features), length(samples))

    featuredict = Dict(f => i for (i, f) in enumerate((f for f in features)))
    sampledict = Dict(f => i for (i, f) in enumerate(samples))

    @info "Filling community profile"
    Threads.@threads for s in samples
        file = allfiles[findfirst(f-> contains(f,s), allfiles)]

        sidx = sampledict[first(split(basename(file), '_')) |> String]
        for row in CSV.File(file)
            contains(row[1], "|") && continue
            fidx = featuredict[row[1]]
            @inbounds abunds[fidx, sidx] = row[2]
        end
    end

    return CommunityProfile(abunds, GeneFunction.(features), MicrobiomeSample.(samples))
end

function get_neuroactive_kos(neuroactivepath="data/gbm.txt")
    neuroactive = Dictionary{String, Vector{String}}()
    desc = ""
    for line in eachline(neuroactivepath)
       line = split(line, r"[\t,]")
       if startswith(line[1], "MGB")
           (mgb, desc) = line
           desc = rstrip(replace(desc, r"\bI+\b.*$"=>""))
           desc = replace(desc, r" \([\w\s]+\)$"=>"")
           desc = replace(desc, r"^.+ \(([\w\-]+)\) (.+)$"=>s"\1 \2")
           @info "getting unirefs for $desc"
           !in(desc, keys(neuroactive)) && insert!(neuroactive, desc, String[])
       else
           filter!(l-> occursin(r"^K\d+$", l), line)
           append!(neuroactive[desc], String.(line))
       end
   end
   return neuroactive
end

function getneuroactive(features; neuroactivepath="data/gbm.txt", map_ko_uniref_path="data/map_ko_uniref90.txt.gz")
    neuroactivekos = get_neuroactive_kos(neuroactivepath)

    kos2uniref = Dictionary{String, Vector{String}}()
    for line in eachline(GzipDecompressorStream(open(map_ko_uniref_path)))
        line = split(line, '\t')
        insert!(kos2uniref, line[1], map(x-> String(match(r"UniRef90_(\w+)", x).captures[1]), line[2:end]))
    end

    neuroactive_index = Dictionary{String, Vector{Int}}()
    for na in keys(neuroactivekos)
        searchfor = Iterators.flatten([kos2uniref[ko] for ko in neuroactivekos[na] if ko in keys(kos2uniref)]) |> Set
        pos = findall(u-> u in searchfor, features)
        insert!(neuroactive_index, na, pos)
    end
    for k in keys(neuroactive_index)
        unique!(neuroactive_index[k])
    end
    return neuroactive_index
end

fsea(cors, pos) = (cors, pos, MannWhitneyUTest(cors[pos], cors[Not(pos)]))

function fsea(cors, allfeatures::AbstractVector, searchset::Set)
    pos = findall(x-> x in searchset, allfeatures)

    return fsea(cors, pos)
end

function fsea(occ::AbstractMatrix, metadatum::AbstractVector, pos::AbstractVector{<:Int})
    let notmissing = map(!ismissing, metadatum)
        occ = occ[:, notmissing]
        metadatum = metadatum[notmissing]
    end

    cors = cor(metadatum, occ, dims=2)'
    return fsea(cors, pos)
end