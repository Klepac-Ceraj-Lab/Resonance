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

function plot_pcoa!(ax, M::MultivariateStats.MDS; dims=(1,2), kwargs...)
    ax.xlabel = get(kwargs, :xlabel, mdsaxis(M, dims[1]))
    ax.ylabel = get(kwargs, :ylabel, mdsaxis(M, dims[2]))
    
    scatter!(ax, loadings(M, dims[1]), loadings(M, dims[2]); kwargs...)
end


varexpl(p::PERMANOVA.PSummary) = p.results[1, 3] * 100
pvalue(p::PERMANOVA.PSummary) = p.results[1, 5]

function plot_permanovas(pdf; commlabels=[], mdlabels=[], colormap=:blues, colorrange=(0,10))
    fig = Figure()
    ax = Axis(fig[1,1])
    hm = plot_permanovas!(ax, pdf; commlabels, mdlabels, colormap, colorrange)

    return fig, ax, hm
end

function plot_permanovas!(ax, pdf; commlabels=unique(pdf.label), mdlabels=unique(pdf.metadatum), colormap=:blues, colorrange=(0,10))
    pgrp = groupby(pdf, :label)
    
    vmat = mapreduce(df-> df.varexpl, hcat, pgrp)
    pmat = mapreduce(df-> df.pvalue, hcat, pgrp)
    
    hm = heatmap!(ax, vmat'; colormap, colorrange)

    for ci in CartesianIndices(pmat)
        c = vmat[ci] < (colorrange[2] - colorrange[1]) / 2 ? :black : :lightgray
        text!(ax, string(round(vmat[ci], digits=2), "%"); position=(ci[2],ci[1]), align=(:center, :center), color=c)
        p = pmat[ci]
        stars = p < 0.001 ? "***" : p < 0.01 ? "**" : p < 0.05 ? "*" : ""
        text!(ax, stars; position=(ci[2],ci[1]), align=(:center, :bottom), color=c)
    end

    ax.xticks = (1:length(commlabels), commlabels)
    ax.yticks = (1:length(mdlabels), mdlabels)

    return hm
end

function plot_mantel(manteldf; commlabels=unique([manteldf.thing1; manteldf.thing2]), colormap=:deep, colorrange=(0,100))
    fig = Figure()
    ax = Axis(fig[1,1], title="Mantel tests")
    hm = plot_mantel!(ax, manteldf; commlabels, mdlabels, colormap, colorrange)

    return fig, ax, hm
end

function plot_mantel!(ax, manteldf; commlabels=unique([manteldf.thing1; manteldf.thing2]), colormap=:deep, colorrange=(0,100))
    n = length(commlabels)
    labidx = Dict(l=> i for (i, l) in enumerate(commlabels))

    vmat = zeros(n, n)
    pmat = ones(n, n)

    manteldf.idx = [CartesianIndex(i, j) for (i,j) in zip(map(t-> labidx[t], manteldf.thing1), map(t-> labidx[t], manteldf.thing2))]
    for row in eachrow(manteldf)
        vmat[row.idx] = row.stat
        pmat[row.idx] = row.pvalue
    end

    #-

    hm = heatmap!(ax, vmat[1:n-1, 2:n]'; colormap, colorrange = (0.01, 1), lowclip=:lightgray)

    for ci in manteldf.idx
        c = vmat[ci] < 0.5 ? :black : :lightgray 
        text!(ax, string(round(vmat[ci], digits=4)); position=(ci[2]-1,ci[1]), align=(:center, :center), color=c)
        p = pmat[ci]
        stars = p < 0.001 ? "***" : p < 0.01 ? "**" : p < 0.05 ? "*" : ""
        text!(ax, stars; position=(ci[2]-1,ci[1]), align=(:center, :bottom), color=c)
    end

    ax.xticks = (1:n-1, isempty(commlabels) ? ["comm$i" for i in 1:n-1] : commlabels[2:end])
    ax.yticks = (1:n-1, isempty(commlabels) ? ["comm$i" for i in 1:n] : commlabels[1:end-1])

    return hm
end


function format_species_label(lab)
    parts = split(String(lab), "_")
    length(parts) == 1 && return parts[1]
    parts[2] == "sp" && (parts[2] = "sp.")
    lab = rich(join(parts[1:2], " "); font=:italic)
    length(parts) > 2 ? rich(lab, " ", join(parts[3:end], " "); font=:regular) : lab
end

format_species_labels(labs) = format_species_label.(labs)

function plot_fsea(setcors, notcors; label="", ylabel="enrichment", kwargs...)
    fig = Figure()
    grid, ax1, ax2 = plot_fsea!(fig.layout, setcors, notcors; label, ylabel, kwargs...)
    fig, ax1, ax2
end

function plot_fsea!(grid, setcors, notcors; title="", ylabel="", linecolor=:gray, kwargs...)
    fullcors = [setcors; notcors]
    ncors = length(fullcors)

    srt = sortperm(fullcors)
    ranks = invperm(srt)
    setranks = Set(ranks[1:length(setcors)])
    
    setscore =  1 / length(setcors)
    notscore = -1 / length(notcors)
    
    xs = 1:ncors
    ys = cumsum(i ∈ setranks ? setscore : notscore for i in eachindex(ranks)) .* -1
    
    t = "ES: $(round(enrichment_score(fullcors, 1:length(setcors)), digits=3))"

    ax1 = Axis(grid[1,1]; title, ylabel=isempty(ylabel) ? t : ylabel, kwargs...)
    hidexdecorations!(ax1)

    ax2 = Axis(grid[2,1])
    hidedecorations!(ax2)
    
    lines!(ax1, xs, ys; color=linecolor)
    vlines!(ax2, ranks[1:length(setcors)]; color=:black)
    
    rowsize!(grid, 2, Relative(1/8))
    rowgap!(grid, 1, Fixed(0))

    linkxaxes!(ax1, ax2)
    tightlimits!.((ax1, ax2))
    grid, ax1, ax2
end

function plot_corrband!(ax, cors; bandres=5000)
    ax.xlabel="rank"
    ax.ylabel="correlation"
    srt = sortperm(cors)
    ncors = length(cors)
    rn = ncors > bandres ? round.(Int, range(1, ncors; length=bandres)) : range(1, ncors)
    blow, bup = extrema(cors)

    xs = 1:ncors
    band!(ax, xs[rn], fill(blow, length(rn)), fill(bup, length(rn)); color=cors[srt[rn]])
    
    lower = [x < 0 ? x : 0.0 for x in cors[srt]]
    upper = [x > 0 ? x : 0.0 for x in cors[srt]]
    band!(ax, rn, lower[rn], upper[rn]; color=:lightgray)
    tightlimits!(ax)
end

function treepositions(hc; useheight = true, orientation = :vertical)
    order = Dict(a=> i for (i,a) in enumerate(hc.order))
    nodepos = Dict(-i => (float(order[i]), 0.0) for i in hc.order)
    xs = []
    ys = []
    for i in 1:size(hc.merges, 1)
        x1, y1 = nodepos[hc.merges[i, 1]]
        x2, y2 = nodepos[hc.merges[i, 2]]
        xpos = (x1 + x2) / 2
        ypos = useheight ?  hc.heights[i] : (max(y1, y2) + 1)
        nodepos[i] = (xpos, ypos)
        push!(xs, [x1, x1, x2, x2])
        push!(ys, [y1, ypos, ypos, y2])
    end
    if orientation == :horizontal
        return ys, xs
    else
        return xs, ys
    end
end

function dendrogram(h; color = :black, kwargs...)
    ax = CairoMakie.Axis(Figure()[1, 1])
    dendrogram!(ax, h; color, kwargs...)
    current_figure()
end

function dendrogram!(ax, h; color = :black, kwargs...)
    for (x, y) in zip(treepositions(h; kwargs...)...)
        lines!(ax, x, y; color)
    end
end