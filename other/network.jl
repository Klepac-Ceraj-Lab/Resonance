using Resonance
using FileIO
using Graphs
using ThreadsX
using Dictionaries
using MetaGraphsNext
using NetworkLayout
using CairoMakie
using GraphMakie

#-

mdata = Resonance.load(Metadata())

mdata.sample

#-

taxa_files = filter(readdir("/grace/sequencing/processed/mgx/metaphlan/", join=true)) do f
    fname = basename(f)
    contains(fname, "profile") || return false
    m = match(r"^(FG\d+_S\d+)_", fname)
    isnothing(m) && return false
    return m[1] in mdata.sample
end

#-

g = MetaGraph(
    Graph();
    label_type       = Symbol,
    vertex_data_type = String,
    edge_data_type   = Float64,
    graph_data       = "Microbial community graph",
    weight_function  = identity,
    default_weight   = 0,
)

g[:root] = "root"

function add_taxstring!(graph, taxstring)
    spltax = Symbol.(split(taxstring, '|'))
    node = pop!(spltax)
    !haskey(graph, node) && (graph[node] = string(node))
    if isempty(spltax)
        g[:root, node] = 1.
    else
        parent = pop!(spltax)
        haskey(graph, parent) || (graph[parent] = string(parent))
        graph[parent, node] = 1.
    end
    return node
end


foreach(taxa_files) do f
    sample = Symbol(match(r"FG\d+_S\d+", basename(f)).match)
    df = CSV.read(f, DataFrame; header=["taxon", "ncbi", "abundance", "additional"], skipto=5)

    g[sample] = String(sample)
    for row in eachrow(df)
        node = add_taxstring!(g, row.taxon)
        startswith(String(node), "s__") && (g[sample, node] = row.abundance)
    end
end

#-

nv(g)

pathways_files = filter(readdir("/grace/sequencing/processed/mgx/humann/main", join=true)) do f
    fname = basename(f)
    contains(fname, "pathabundance") || return false
    m = match(r"^(FG\d+_S\d+)_", fname)
    isnothing(m) && return false
    return m[1] in mdata.sample
end

foreach(pathways_files) do f
    sample = Symbol(match(r"FG\d+_S\d+", basename(f)).match)
    df = CSV.read(f, DataFrame; header=["function", "abundance"], skipto=2)

    for row in eachrow(df)
        m = match(r"^(.+?)\|g__\w+\.(s__\w+)$", row.function)
        isnothing(m) && continue
        m[1] in ("UNINTEGRATED", "UNMAPPED") && continue
        (node, sp) = Symbol.(m.captures)
        !haskey(g, sp) && continue
        !haskey(g, node) && (g[node] = string(node))
        g[node, sp] = row.abundance / 1e6
        g[node, sample] = 1
    end

end

nv(g)
#-



function node_color(node)
    node = String(node)
    node == "root" && return Makie.to_color(:black)
    if contains(node, r"[kpcofgs]__") 
        return Makie.to_colormap(:viridis)[
            Dict('k'=> 1,'p'=> 41,'c'=> 81,'o'=> 121,'f'=> 161,'g'=> 201,'s'=> 256)[
                first(node)
            ]]
    elseif startswith(node, r"FG\d+")
        return Makie.to_color(:dodgerblue)
    else
        return Makie.to_color(:orchid)
    end
end

layout = Spring(; Ptype=Float32)
(f, ax, p) = graphplot(g; layout, node_color=[node_color(node) for node in labels(g)])


#-