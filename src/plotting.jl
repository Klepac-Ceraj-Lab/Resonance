function count_set(df, cols, setcols)
    count(eachrow(df)) do row
        any(ismissing, row[cols[setcols]])        && return false # all setcols must be non-missing
        !all(row[cols[setcols]])                  && return false # all set colls must be true
        any(skipmissing(row[cols[Not(setcols)]])) && return false # no non-set cols can be true
        return true
    end
end

function upset_dots!(ax, colsets, nsets=maximum(Iterators.flatten(colsets)))
    xlims!(ax, (0, length(colsets) + 1))
    hidexdecorations!(ax)
    hideydecorations!(ax, ticklabels=false)
    hidespines!(ax)

    for (x, p) in enumerate(colsets)
        scatter!(ax, fill(x, nsets - length(p)), collect(1:nsets)[Not(p)], markersize=20, color=:lightgray)
        scatter!(ax, fill(x, length(p)), p, markersize=20, color=:black)
        length(p) > 1 && lines!(ax, [x,x], [extrema(p)...], color=:black, linewidth=3)
    end
end

# https://github.com/JuliaStats/MultivariateStats.jl/pull/162
function loadings(M::MultivariateStats.MDS)
    ev = eigvals(M)
    return ev' .* projection(M)[:, 1:length(ev)]
end

# https://github.com/JuliaStats/MultivariateStats.jl/pull/162
function loadings(M::MultivariateStats.MDS, dim)
    l = loadings(M)
    return l[:, dim]
end

varexplained(M::MultivariateStats.MDS) = eigvals(M) ./ sum(eigvals(M))

mdsaxis(M::MultivariateStats.MDS, dim::Int) = "MDS$dim ($(round(varexplained(M)[dim] * 100, digits=2))%)"

varexpl(p::PERMANOVA.PSummary) = p.results[1, 3] * 100
pvalue(p::PERMANOVA.PSummary) = p.results[1, 5]

function permanovas(comm, metadatums)
    permdf = DataFrame()
    for md in metadatums
        com_md = get(comm, md)
        hasmd = findall(!ismissing, com_md)
        df = DataFrame(test = com_md[hasmd])
        disallowmissing!(df)
        size(df, 1) < 20 && @warn "Very small number of not missing $md"

        p = permanova(df, abundances(comm)[:, hasmd]', BrayCurtis, @formula(1 ~ test))
        push!(permdf, (; metadatum = md, varexpl=varexpl(p), pvalue=pvalue(p)))
    end
    return permdf

end

function plot_permanovas(comms, metadatums; commlabels=[], mdlabels=[], colormap=:blues, colorrange=(0,10))
    fig = Figure()
    ax = Axis(fig[1,1])
    hm = plot_permanovas!(ax, comms, metadatums; commlabels, mdlabels, colormap, colorrange)

    return fig, ax, hm
end

function plot_permanovas!(ax, comms, metadatums; commlabels=[], mdlabels=[], colormap=:blues, colorrange=(0,10))
    ps = [permanovas(comm, metadatums) for comm in comms]
    vmat = mapreduce(df-> df.varexpl, hcat, ps)
    pmat = mapreduce(df-> df.pvalue, hcat, ps)
    
    hm = heatmap!(ax, vmat'; colormap, colorrange)

    for ci in CartesianIndices(pmat)
        text!(string(round(vmat[ci], digits=2), "%"); position=(ci[2],ci[1]), align=(:center, :center))
        p = pmat[ci]
        stars = p < 0.001 ? "***" : p < 0.01 ? "**" : p < 0.05 ? "*" : ""
        text!(stars; position=(ci[2],ci[1]), align=(:center, :bottom))
    end

    ax.xticks = (1:length(comms), isempty(commlabels) ? ["comm$i" for i in 1:length(comms)] : commlabels)
    ax.yticks = (1:length(metadatums), isempty(commlabels) ? string.(metadatums) : mdlabels)

    return hm
end