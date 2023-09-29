function write_gfs_arrow(; kind="genefamilies", names=false, stratified=false, sample_filter=nothing)
    root = analysisfiles("humann", names ? "rename" : 
                                    kind == "genefamilies" ? "main" : "regroup"
    )
    stripper = "_$kind" * (names ? "_rename.tsv" : ".tsv")
    filt = Regex(string(raw"SEQ\d+_S\d+_", kind))
    df = DataFrame(file = filter(f-> contains(f, filt), readdir(root, join=true)))
    df.sample = map(s-> replace(s, stripper => ""), basename.(df.file))
    df.seqid = map(s-> replace(s, r"_S\d+$"=> ""), df.sample)
    isnothing(sample_filter) || subset!(df, AsTable(["sample", "seqid"]) => ByRow(s-> s.sample ∈ sample_filter || s.seqid ∈ sample_filter))
    @debug "Found $(nrow(df)) files"

    knead = load(ReadCounts())
    leftjoin!(df, select(knead, "sample_uid"=>"sample", 
                                AsTable(["final pair1", "final pair2"])=> ByRow(row-> row[1]+row[2]) =>"read_depth");
                            on="sample"
    )
    
    @info "getting features"
    feats = mapreduce(union, eachrow(df)) do row
        fs = CSV.read(row.file, DataFrame; header=["feature", "value"], skipto=2, select=[1])[!,1]
        stratified || filter!(f->!contains(f, '|'), fs) # skip stratified feats
        Set(fs)
    end
    featuremap = Dict(f=> i for (i,f) in enumerate(feats))
    
    scratch = scratchfiles("genefunctions")
    isdir(scratch) || mkpath(scratch)

    @info "writing arrow file"
    open(joinpath(scratch, "$kind.arrow"), "w") do io
        tbls = Tables.partitioner(eachrow(df)) do row
            @debug "writing $(row.sample)"

            sdf = CSV.read(row.file, DataFrame; header=["feature", "value"], skipto=2)
            stratified || subset!(sdf, "feature"=> ByRow(f-> !contains(f, '|'))) # skip stratified feats
            sdf.sample .= row.sample
            sdf.sidx .= rownumber(row)
            sdf.fidx = ThreadsX.map(f-> featuremap[f], sdf.feature)

            sdf
        end

        Arrow.write(io, tbls; metadata=("features" => join(feats, '\n'), 
                                        "samples"  => join(df.sample, '\n'),
                                        "files"    => join(df.file, '\n'),
                                        "reads"    => join(df.read_depth, '\n')
                                        )
        )                                
    end


    return nothing
end

function read_gfs_arrow(; kind="genefamilies")
    @info "reading table"
    tbl = Arrow.Table(scratchfiles("genefunctions", "$kind.arrow"))
    @info "building sparse mat"
    mat = sparse(tbl.fidx, tbl.sidx, tbl.value)
    mdt =  Arrow.getmetadata(tbl)

    @info "getting features"
    fs = [genefunction(line) for line in eachline(IOBuffer(mdt["features"]))]
    @info "getting samples"
    mdt = DataFrame(
        sample = [MicrobiomeSample(line) for line in eachline(IOBuffer(mdt["samples"]))],
        read_depth = map(l-> l=="missing" ? missing : parse(Float64, l), eachline(IOBuffer(mdt["reads"]))),
        file = readlines(IOBuffer(mdt["files"]))
    )
    mdt.sample_base = replace.(name.(mdt.sample), r"_S\d+" => "")
    comm = CommunityProfile(mat, fs, mdt.sample)
    set!(comm, mdt)
    return comm
end


function get_neuroactive_kos(neuroactivepath=datafiles("gbm.txt"); consolidate=true)
    neuroactive = Dictionary{String, Vector{String}}()
    desc = ""
    for line in eachline(neuroactivepath)
       line = split(line, r"[\t,]")
       if startswith(line[1], "MGB")
            (mgb, desc) = line
            if consolidate
                desc = rstrip(replace(desc, r"\b[IV]+\b.*$"=>""))
                desc = replace(desc, r" \([\w\s\-]+\)"=>"")
                desc = replace(desc, r"^.+ \(([\w\-]+)\) (.+)$"=>s"\1 \2")
                desc = replace(desc, " (AA"=>"")
            end
            @info "getting unirefs for $desc"
            !in(desc, keys(neuroactive)) && insert!(neuroactive, desc, String[])
       else
           filter!(l-> occursin(r"^K\d+$", l), line)
           append!(neuroactive[desc], String.(line))
       end
   end
   return neuroactive
end

function getneuroactive(features; neuroactivepath=datafiles("gbm.txt"), map_ko_uniref_path=datafiles("map_ko_uniref90.txt.gz"), consolidate=true)
    neuroactivekos = get_neuroactive_kos(neuroactivepath; consolidate)

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
    cors = let notmissing = map(!ismissing, metadatum)
        occ = occ[:, notmissing]
        metadatum = metadatum[notmissing]
        cor(metadatum, occ, dims=2)'
    end

    return fsea(cors, pos)
end

function enrichment_score(cors, pos)
    srt = sortperm(cors; rev=true)
    ranks = invperm(srt)
    setranks = ranks[pos]    
    ng = length(pos)

    setscore =  1 / ng
    notscore = -1 / (length(cors) - ng)

    ys = cumsum(i ∈ setranks ? setscore : notscore for i in eachindex(ranks))
    lower, upper = extrema(ys)
    return abs(lower) > abs(upper) ? lower : upper
end

function fsea_permute(cors, pos; nperm=1000)
    es = enrichment_score(cors, pos)
    ng = length(pos)
    lt = ThreadsX.count(1:nperm) do _
        spos = rand(eachindex(cors), ng)
        ses = enrichment_score(cors, spos)
        es < 0 ? ses < es : ses > es
    end
    return (lt / nperm, es)
end
        
