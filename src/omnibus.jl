function permanovas(dm::Matrix, metadatums::Vector{<:AbstractVector}; n = 1000, mdlabels = ["md$i" for i in eachindex(metadatums)])
    permdf = mapreduce(vcat, zip(metadatums, mdlabels)) do (md, lab)    
        hasmd = findall(!ismissing, md)
        df = DataFrame(test = md[hasmd])
        disallowmissing!(df)
        size(df, 1) < 20 && @warn "Very small number of not missing $lab"

        p = permanova(df, dm[hasmd, hasmd], @formula(1 ~ test), n)

        return (; metadatum = lab, varexpl=varexpl(p), pvalue=pvalue(p))
    end

    return DataFrame(permdf)
end

function permanovas(dms::AbstractArray{<:AbstractMatrix}, metadatums; n = 1000, commlabels = [], mdlabels = ["md$i" for i in eachindex(metadatums)])
    isempty(commlabels) && (commlabels = ["comm$i" for i in eachindex(dms)])
    
    mapreduce(vcat, enumerate(dms)) do (i, dm)
        @info "Permanovas for $(commlabels[i])"
        df = permanovas(dm, metadatums; n, mdlabels)
        df.label .= commlabels[i]
        return df
    end
end


function permanovas(comm::CommunityProfile, metadatums; n = 1000, mdlabels = String.(metadatums))
    permdf = mapreduce(vcat, zip(metadatums, mdlabels)) do (md, lab)    
        com_md = get(comm, md)
        hasmd = findall(!ismissing, com_md)
        df = DataFrame(test = com_md[hasmd])
        disallowmissing!(df)
        size(df, 1) < 20 && @warn "Very small number of not missing $lab"

        p = permanova(df, abundances(comm)[:, hasmd]', BrayCurtis, @formula(1 ~ test), n)

        return (; metadatum = lab, varexpl=varexpl(p), pvalue=pvalue(p))
    end

    return DataFrame(permdf)
end

function permanovas(comms::Vector{<:CommunityProfile}, metadatums; n = 1000, commlabels = [], mdlabels = String.(metadatums))
    isempty(commlabels) && (commlabels = ["comm$i" for i in eachindex(comms)])
    
    mapreduce(vcat, enumerate(comms)) do (i, c)
        @info "Permanovas for $(commlabels[i])"
        df = permanovas(c, metadatums; n, mdlabels)
        df.label .= commlabels[i]
        return df
    end
end

function mantel(mat1, mat2; n = 1000)
    size(mat1) == size(mat2) || throw(ArgumentError("Matrices must be the same size"))
    all(issymmetric, (mat1, mat2)) || error("Matries are not symetric")
    
    ord = size(mat1, 1)
    
    v1 = view(mat1, [CartesianIndex(i,j) for i in 1:ord for j in 1:ord if i < j])
    v2 = view(mat2, [CartesianIndex(i,j) for i in 1:ord for j in 1:ord if i < j])

    c_real = cor(v1, v2)
    cs = Vector{Float64}(undef, n)

    @inbounds ThreadsX.foreach(range(1, n)) do i
        cs[i] = cor(view(v1, randperm(ord)), view(v2, randperm(ord)))
    end

    return (c_real, sum(>(c_real), cs) / n)

end

function mantel(dms::Vector{<:Matrix}; commlabels = [], n = 1000)
    manteldf = DataFrame()
    isempty(commlabels) && (commlabels = ["comm$i" for i in eachindex(comms)])

    for (i,j) in combinations(eachindex(dms), 2)
        dm1, dm2 = dms[[i, j]]
        l1, l2 = commlabels[[i, j]]
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
