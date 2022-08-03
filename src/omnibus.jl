function permanovas(comm, metadatums; n = 1000)
    permdf = DataFrame()
    
    for md in metadatums
        com_md = get(comm, md)
        hasmd = findall(!ismissing, com_md)
        df = DataFrame(test = com_md[hasmd])
        disallowmissing!(df)
        size(df, 1) < 20 && @warn "Very small number of not missing $md"

        p = permanova(df, abundances(comm)[:, hasmd]', BrayCurtis, @formula(1 ~ test); n_perm=n)

        push!(permdf, (; metadatum = md, varexpl=varexpl(p), pvalue=pvalue(p)))
    end

    return permdf
end

function permanovas(comms::AbstractArray{<:CommunityProfile}, metadatums; n = 1000)
    permdf = DataFrame()
    
    for (i, c) in eumerate(comms)
        @info "Permanovas for community profile $i"
        append!(permdf, permanovas(c, metadatums); n)
    end
    
    return permdf
end

function mantel(mat1, mat2; n = 1000)
    size(mat1) == size(mat2) || throw(ArgumentError("Matrices must be the same size"))
    all(issymmetric, (mat1, mat2)) || error("Matries are not symetric")
    
    ord = size(mat1, 1)
    
    v1 = view(mat1, [CartesianIndex(i,j) for i in 1:ord for j in 1:ord if i < j])
    v2 = view(mat2, [CartesianIndex(i,j) for i in 1:ord for j in 1:ord if i < j])

    c_real = cor(v1, v2)
    cs = Vector{Float64}(undef, n)

    @inbounds for i in range(1, n)
        cs[i] = cor(view(v1, randperm(ord)), view(v2, randperm(ord)))
    end

    return (c_real, sum(>(c_real), cs) / n)

end

function mantel(comms::AbstractArray{<:CommunityProfile}, mdlabels; n = 1000)
    manteldf = DataFrame()

    for ((c1, l1), (c2, l2)) in combinations(collect(zip(comms, mdlabels)), 2)
        @info "Mantel for $l1 and $l2"
        c1, c2 = comm_overlap(c1, c2)
        c1 = c1[featurenames(c1) .!= "UNGROUPED", :]
        c2 = c2[featurenames(c2) .!= "UNGROUPED", :]

        dm1 = braycurtis(c1)
        dm2 = braycurtis(c2)

        m, p = mantel(dm1, dm2; n)
        push!(manteldf, (; stat = m, pvalue = p, thing1 = l1, thing2 = l2))      
    end

    return manteldf
end