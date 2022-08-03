function permanovas(comm, metadatums; n = 1000, mdlabels = String.(metadatums))
    permdf = DataFrame()
    
    for (md, lab) in zip(metadatums, mdlabels)
        com_md = get(comm, md)
        hasmd = findall(!ismissing, com_md)
        df = DataFrame(test = com_md[hasmd])
        disallowmissing!(df)
        size(df, 1) < 20 && @warn "Very small number of not missing $md"

        p = permanova(df, abundances(comm)[:, hasmd]', BrayCurtis, @formula(1 ~ test), n)

        push!(permdf, (; metadatum = lab, varexpl=varexpl(p), pvalue=pvalue(p)))
    end

    return permdf
end

function permanovas(comms::AbstractArray{<:CommunityProfile}, metadatums; n = 1000, commlabels = [], mdlabels = String.(metadatums))
    permdf = DataFrame()
    isempty(commlabels) && (commlabels = ["comm$i" for i in eachindex(comms)])
    
    for (i, c) in enumerate(comms)
        @info "Permanovas for $(commlabels[i])"
        df = permanovas(c, metadatums; n, mdlabels)
        df.label .= commlabels[i]
        append!(permdf, df)
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

function mantel(comms::AbstractArray{<:CommunityProfile}, commlabels = []; n = 1000)
    manteldf = DataFrame()
    isempty(commlabels) && (commlabels = ["comm$i" for i in eachindex(comms)])

    for ((c1, l1), (c2, l2)) in combinations(collect(zip(comms, commlabels)), 2)
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

@testset "Omnibus tests" begin
    c1 = CommunityProfile(hcat(rand(20,3), rand(20,2) .* 2), [Taxon("sp$i") for i in 1:20], [MicrobiomeSample("s$i") for i in 1:5])
    set!(c1, DataFrame(subject = 1:5, timepoint = ones(5),
                       sample = ["s$i" for i in 1:5], age = 1:5,
                       category=["young", "young", "young", "old", "old"]))
    c2 = CommunityProfile(hcat(rand(20,3), rand(20,2) .* 2), [Taxon("sp$i") for i in 1:20], [MicrobiomeSample("s$i") for i in 1:5])
    set!(c2, DataFrame(subject = 1:5, timepoint = ones(5),
                       sample = ["s$i" for i in 1:5], age = 1:5,
                       category=["young", "young", "young", "old", "old"]))

    @testset "Permanovas" begin
        p = permanovas(c1, [:age, :category])
        @test p isa DataFrame
        @test size(p) == (2,3)
        @test names(p) == ["metadatum", "varexpl", "pvalue"]
        
        p2 = permanovas([c1, c2], [:age, :category])
        @test p2 isa DataFrame
        @test size(p2) == (4,4)
        @test names(p2) == ["metadatum", "varexpl", "pvalue", "label"]
    end

    @testset "Mantel tests" begin
        m1 = mantel(braycurtis(c1), braycurtis(c2))
        @test m1 isa Tuple{Float64, Float64}
        @test all(<(1), m1)
        
        m2 = mantel([c1, c2])
        @test m2 isa DataFrame
        @test size(m2) == (1,4)
        @test names(m2) == ["stat", "pvalue", "thing1", "thing2"]
    end
end
