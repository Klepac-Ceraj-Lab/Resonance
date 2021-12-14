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
    ev' .* projection(M)[:, 1:length(ev)]
end

varexplained(M::MultivariateStats.MDS) = eigvals(M) ./ sum(eigvals(M))