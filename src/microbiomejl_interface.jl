function write_arrow(filename, cmp::CommunityProfile; unirefs = false)
    df = DataFrame(metadata(cmp))
    @debug "filtering profile"
    cmp = cmp[sort(ThreadsX.unique(first(findnz(abundances(cmp))))), :]
    @debug "getting features and samples"
    feats = features(cmp)
    featidx = Dict(string(f)=>i for (i, f) in enumerate(feats))
    smps = samples(cmp)
    smpidx = Dict(name(s)=>i for (i, s) in enumerate(smps))

    open(filename, "w") do io
        tbls = Tables.partitioner(smps) do smp
            sample = name(smp)
            sidx = smpidx[name(smp)]

            sdf = DataFrame([(; sample,
                                sidx,
                                feature = string(f),
                                fidx   = featidx[string(f)],
                                value = cmp[f, smp]
                     ) for f in feats if cmp[f, smp] > 0]
            )

            @debug "writing $sample"


            sdf
        end

        Arrow.write(io, tbls; metadata=("features" => unirefs ? "nothing" : join(string.(features(cmp)), '\n'), 
                                        "samples"  => join(name.(samples(cmp)), '\n'),
                                        "reads"    => join(df.read_depth, '\n')
                                        )
        )                                
    end
end

function read_arrow(filename; featurefunc = taxon, unirefs = false)
    @info "reading table"
    tbl = Arrow.Table(filename)
    @info "building sparse mat"
    mat = sparse(tbl.fidx, tbl.sidx, tbl.value)
    mdt =  Arrow.getmetadata(tbl)

    @info "getting features"
    fs = [featurefunc(line) for line in eachline(IOBuffer(mdt["features"]))]
    @info "getting samples"
    mdt = DataFrame(
        sample = [MicrobiomeSample(line) for line in eachline(IOBuffer(mdt["samples"]))],
        read_depth = map(l-> l=="missing" ? missing : parse(Float64, l), eachline(IOBuffer(mdt["reads"]))),
    )
    mdt.sample_base = replace.(name.(mdt.sample), r"_S\d+" => "")
    CommunityProfile(mat, fs, mdt.sample)
end

function comm2wide(cmp::CommunityProfile; feature_filter=identity, sample_filter=identity, feature_funct=name)
    cmp = cmp[feature_filter(features(cmp)), sample_filter(samples(cmp))]
    md = DataFrame(get(cmp))
    return hcat(md, DataFrame(collect(abundances(cmp)'), feature_funct.(features(cmp)), makeunique=true))
end