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

function plot_upset(df, ycols, ylabels, intersects; figure=(;))
    fig = Figure(; figure...)

    intersection_ax = Axis(fig[1,1:2]; ylabel="intersection size", yautolimitmargin = (0, 0.15))
    dot_ax = Axis(fig[2,1:2], yticklabelpad = 10, yticks = (1:length(ylabels), ylabels))
    set_ax = Axis(fig[2,3], xlabel="set size", xticklabelrotation=π/4,
                    xautolimitmargin = (0, 0.25),  xgridvisible = false)


    barplot!(intersection_ax, 1:length(intersects), [count_set(df, ycols, i) for i in intersects],
                bar_labels=:y, color = :gray20,
                label_size = 14, label_formatter = x -> string(Int(x)))

    barplot!(set_ax, 1:length(ylabels), [count(x-> !ismissing(x) && x, df[!, col]) for col in ycols], 
                direction=:x, bar_labels=:y, label_size = 14, color = :gray20,
                label_formatter = x -> string(Int(x)))

                
    for i in 1:2:length(ycols)
        poly!(dot_ax,
        BBox(0, length(intersects) + 1, i-0.5, i+0.5),
        color = :gray95
        )
    end

    upset_dots!(dot_ax, intersects)

    hidexdecorations!(intersection_ax)
    hideydecorations!(set_ax)

    rowgap!(fig.layout, 0)
    linkyaxes!(dot_ax, set_ax)
    linkxaxes!(dot_ax, intersection_ax)
    hidespines!(intersection_ax, :t, :r, :b)
    hidespines!(set_ax, :t, :r, :l)
   
    return fig
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
        c = vmat[ci] < (colorrange[2] - colorrange[1]) / 2 ? :black : :lightgray
        text!(string(round(vmat[ci], digits=2), "%"); position=(ci[2],ci[1]), align=(:center, :center), color=c)
        p = pmat[ci]
        stars = p < 0.001 ? "***" : p < 0.01 ? "**" : p < 0.05 ? "*" : ""
        text!(stars; position=(ci[2],ci[1]), align=(:center, :bottom), color=c)
    end

    ax.xticks = (1:length(comms), isempty(commlabels) ? ["comm$i" for i in 1:length(comms)] : commlabels)
    ax.yticks = (1:length(metadatums), isempty(commlabels) ? string.(metadatums) : mdlabels)

    return hm
end

function plot_fsea(setcors, notcors; label="", bandres=5000)
    fullcors = [setcors; notcors]
    ncors = length(fullcors)

    srt = sortperm(fullcors; rev=true)
    ranks = invperm(srt)
    setranks = Set(ranks[1:length(setcors)])
    
    setscore =  1 / length(setcors)
    notscore = -1 / length(notcors)
    
    xs = 1:ncors
    ys = cumsum(i ∈ setranks ? setscore : notscore for i in eachindex(ranks))
    
    t = "Enrichment score - $(round(max(abs.(extrema(ys))...), digits=3))"
    !isempty(label) && (t = string(label, ": ", t))

    fig = Figure()
    ax1 = Axis(fig[1,1]; title=t, ylabel="enrichment score")
    hidexdecorations!(ax1)

    ax2 = Axis(fig[2,1])
    hidedecorations!(ax2)

    ax3 = Axis(fig[3,1]; ylabel="correlation", xlabel="rank")
    
    lines!(ax1, xs, ys)
    vlines!(ax2, ranks[1:length(setcors)]; color=:black)


    rn = ncors > bandres ? round.(Int, range(1, ncors; length=bandres)) : range(1, ncors)
    blow, bup = extrema(fullcors)
    band!(ax3, rn, fill(blow, length(rn)), fill(bup, length(rn)); color=fullcors[srt[rn]])
    
    lower = [x < 0 ? x : 0.0 for x in fullcors[srt]]
    upper = [x > 0 ? x : 0.0 for x in fullcors[srt]]
    band!(ax3, xs[rn], lower[rn], upper[rn]; color=:lightgray)

    
    rowsize!(fig.layout, 2, Relative(1/8))
    rowsize!(fig.layout, 3, Relative(1/4))
    
    linkxaxes!(ax1, ax2, ax3)
    tightlimits!.((ax1, ax2, ax3))
    fig
end


