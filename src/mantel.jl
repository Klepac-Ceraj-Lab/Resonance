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