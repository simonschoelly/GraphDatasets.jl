
using SimpleValueGraphs: AbstractValGraph, vertices, ne, vertexvals_type,
    edgevals_type, graphvals_type
using SimpleValueGraphs.AbstractTuples: AbstractTuple, typetuple


# ======================================
#   ValGraphCollection
# ======================================

"""
    ValGraphCollection{V <: Integer, V_VALS, E_VALS, G_VALS}

An immutable collection of multiple undirected value graphs.

### See also
[`ng`](@ref), [`ValGraphCollectionView`](@ref)
"""
struct ValGraphCollection{V <: Integer, V_VALS <: AbstractTuple, E_VALS <: AbstractTuple, G_VALS <: AbstractTuple}

    # This data structure is very similar to a sparse matrix in CSR or CSC format,
    # with an additional index array that denotes where a graph starts and ends.
    # One could imagine it like a block matrix.

    # TODO come up with better names
    graph_ids::Vector{Int}  # vertices for graph k are in graph_ids[k]:(graph_ids[k+1]-1)
    vertex_ids::Vector{Int} # outgoing edges of vertex j are in vertex_ids[k]:(vertex_ids[j+1]-1)
    edge_ids::Vector{V}     # edge_ids[i] contains the destination vertex  of edge i

    # TODO maybe change to tuple of arrays
    vertexvals::Vector{V_VALS}
    edgevals::Vector{E_VALS}
    graphvals::Vector{G_VALS}

    function ValGraphCollection(graph_ids, vertex_ids, edge_ids, vertexvals, edgevals, graphvals)

        V = eltype(edge_ids)
        V_VALS = eltype(vertexvals)
        E_VALS = eltype(edgevals)
        G_VALS = eltype(graphvals)

        result = new{V, V_VALS, E_VALS, G_VALS}(graph_ids, vertex_ids, edge_ids, vertexvals, edgevals, graphvals)
        _verify(result)
        return result
    end
end

function _verify(coll::ValGraphCollection{V, V_VALS, E_VALS, G_VALS}) where {V, V_VALS, E_VALS, G_VALS}

    graph_ids = coll.graph_ids
    vertex_ids = coll.vertex_ids
    edge_ids = coll.edge_ids
    vertexvals = coll.vertexvals
    edgevals = coll.edgevals
    graphvals = coll.graphvals

    # when there is no edge
    if length(coll.graph_ids) <= 1
        @assert length(vertex_ids) == 0
        @assert length(edge_ids) == 0
        @assert length(vertevals) == 0
        @assert length(edgevals) == 0
        @assert length(graphvals) == 0

        return nothing
    end

    @assert length(graph_ids) == length(graphvals) + 1
    @assert length(vertex_ids) == length(vertexvals) + 1
    @assert length(edge_ids) == length(edgevals)

    @assert length(graph_ids) >= 1
    @assert graph_ids[begin] >= 1
    @assert graph_ids[end] == length(vertex_ids)
    @assert issorted(graph_ids)

    @assert issorted(vertex_ids)

    N = length(graph_ids) - 1
    for k in 1:N
        v1 = graph_ids[k]
        v2 = graph_ids[k + 1]

        nvg = v2 - v1

        edges = Vector{@NamedTuple{src::Int, dst::Int, edge_index::Int}}()

        for v in v1:v2-1

            for e_id in vertex_ids[v]:vertex_ids[v+1]-1
                @assert edge_ids[e_id] ∈ 1:nvg
            end

            # TODO we should also verify here that the graph contains reverse edges
            # and that the weights of the reverse edges are the same
        end

    end

    return nothing
end


## ----------------------------------------------
##   ng
## ----------------------------------------------

"""
    ng(coll::ValGraphCollection)

Return the number of graphs in `coll`.

### See also
[`ValGraphCollection`](@ref), [`nv`](@ref), [`ne`](@ref)
"""
ng(coll::ValGraphCollection) = max(0, length(coll.graph_ids) - 1)


## ----------------------------------------------
##   iterable & indexable
## ----------------------------------------------

length(coll::ValGraphCollection) = ng(coll)

function getindex(coll::ValGraphCollection, i)

        return ValGraphCollectionView(coll, i)
end

firstindex(coll::ValGraphCollection) = 1
lastindex(coll::ValGraphCollection) = length(coll)


eltype(::Type{<:ValGraphCollection{V, V_VALS, E_VALS, G_VALS}}) where {V, V_VALS, E_VALS, G_VALS} =
    ValGraphCollectionView{V, V_VALS, E_VALS, G_VALS}

function iterate(coll::ValGraphCollection, state = 1)

    state > ng(coll) && return nothing
    return coll[state], (state + 1)
end

## ----------------------------------------------
##   show
## ----------------------------------------------

function show(io::IO, ::MIME"text/plain", coll::ValGraphCollection{V, V_VALS, E_VALS, G_VALS}) where {V, V_VALS, E_VALS, G_VALS}

    println("$(ng(coll))-element ValGraphCollection of graphs with")
    println(io, "              eltype: $V")
    println(io, "  vertex value types: $(typetuple(V_VALS))")
    println(io, "    edge value types: $(typetuple(E_VALS))")
    println(io, "   graph value types: $(typetuple(G_VALS))")
end

# ======================================
#   ValGraphCollectionView
# ======================================

"""
    ValGraphCollectionView{V, V_VALS, E_VALS, G_VALS} <: SimpleValueGraphs.AbstractValGraph

An immutable view on a single graph in a `ValGraphCollection`. Implements the graph
interface of an `AbstractValGraph`.

### See also
[`ValGraphCollection`](@ref)
"""
struct ValGraphCollectionView{V, V_VALS, E_VALS, G_VALS} <: AbstractValGraph{V, V_VALS, E_VALS, G_VALS}

    collection::ValGraphCollection{V, V_VALS, E_VALS, G_VALS}
    graph_num::Int

    function ValGraphCollectionView(coll::ValGraphCollection{V, V_VALS, E_VALS, G_VALS}, graph_num::Integer) where {V, V_VALS, E_VALS, G_VALS}

        @boundscheck checkbounds(1:ng(coll), graph_num)
        return new{V, V_VALS, E_VALS, G_VALS}(coll, graph_num)
    end

end

## ----------------------------------------------
##   implementation SimpleValueGraphs interface
## ----------------------------------------------

is_directed(::Type{<:ValGraphCollectionView}) = false

function nv(g::ValGraphCollectionView)

    graph_ids = g.collection.graph_ids
    @inbounds graph_id1 = graph_ids[g.graph_num]
    @inbounds graph_id2 = graph_ids[g.graph_num + 1]

    return eltype(g)(graph_id2 - graph_id1)
end

function has_edge(g::ValGraphCollectionView, s::Integer, d::Integer)

    if s ∉ vertices(g) || d ∉ vertices(g)
        return false
    end

    @inbounds graph_id = g.collection.graph_ids[g.graph_num]

    @inbounds v_id1 = g.collection.vertex_ids[graph_id + (s - 1)]
    @inbounds v_id2 = g.collection.vertex_ids[graph_id + s]

    # TODO maybe binary search
    return d ∈ @view g.collection.edge_ids[v_id1:v_id2-1]
end

function get_vertexval(g::ValGraphCollectionView, v::Integer, key::Integer)

    graph_id = g.collection.graph_ids[g.graph_num]
    return g.collection.vertexvals[graph_id + (v - 1)][key]
end

function get_edgeval(g::ValGraphCollectionView, s::Integer, d::Integer, key::Integer)

    graph_id = g.collection.graph_ids[g.graph_num]

    v_id1 = g.collection.vertex_ids[graph_id + (s - 1)]
    v_id2 = g.collection.vertex_ids[graph_id + s]

    # TODO error handling if no such edge
    # TODO maybe binary search
    for i ∈ v_id1:v_id2-1
        if g.collection.edge_ids[i] == d
            return g.collection.edgevals[i][key]
        end
    end
end

function get_graphval(g::ValGraphCollectionView, key::Integer)

    return g.collection.graphvals[g.graph_num][key]
end

## ----------------------------------------------
##   additional implementations for better performance
## ----------------------------------------------

function outneighbors(g::ValGraphCollectionView, u::Integer)

    graph_id = g.collection.graph_ids[g.graph_num]

    v_id1 = g.collection.vertex_ids[graph_id + (u - 1)]
    v_id2 = g.collection.vertex_ids[graph_id + u]

    # TODO maybe return an immutable view
    return @view g.collection.edge_ids[v_id1:v_id2-1]
end

function outedgevals(g::ValGraphCollectionView, u::Integer, ::Colon)

    graph_id = g.collection.graph_ids[g.graph_num]

    v_id1 = g.collection.vertex_ids[graph_id + (u - 1)]
    v_id2 = g.collection.vertex_ids[graph_id + u]

    # TODO maybe return an immutable view
    return @view g.collection.edgevals[v_id1:v_id2-1]
end

function outedgevals(g::ValGraphCollectionView, u::Integer, key::Integer)

    graph_id = g.collection.graph_ids[g.graph_num]

    v_id1 = g.collection.vertex_ids[graph_id + (u - 1)]
    v_id2 = g.collection.vertex_ids[graph_id + u]

    # TODO maybe there is some kind of iterator instead of an array
    # could also return a  generator but that looks ugly in the REPL
    return [tup[key] for tup ∈ @view g.collection.edgevals[v_id1:v_id2-1]]
end

## --------------------------------------
##      Converting to other graphs
## --------------------------------------

SimpleGraph(gv::ValGraphCollectionView) = SimpleGraph{eltype(gv)}(gv)

function SimpleGraph{T}(gv::ValGraphCollectionView) where {T}

    fadjlist = Vector{Vector{T}}(undef, nv(gv))
    for v ∈ vertices(gv)
        @inbounds fadjlist[v] = collect(T, outneighbors(gv, v))
    end

    return SimpleGraph(ne(gv), fadjlist)
end

ValGraph(gv::ValGraphCollectionView) = ValGraph{eltype(gv)}(gv)
# TODO this ugly constructor should be moved to SimpleValueGraphs.jl
function ValGraph{V}(gv::ValGraphCollectionView) where {V}

    nvg = Int(nv(gv))
    V_VALS = vertexvals_type(gv)
    E_VALS = edgevals_type(gv)
    G_VALS = graphvals_type(gv)

    fadjlist = Vector{Vector{V}}(undef, nvg)
    vertexvals = Tuple(Vector{VV}(undef, nvg) for VV in V_VALS.types)
    if V_VALS <: NamedTuple
        vertexvals = NamedTuple{Tuple(V_VALS.names)}(vertexvals)
    end
    edgevals = Tuple(Vector{Vector{EE}}(undef, nvg) for EE in E_VALS.types)
    if E_VALS <: NamedTuple
        edgevals = NamedTuple{Tuple(E_VALS.names)}(edgevals)
    end

    for v ∈ vertices(gv)
        @inbounds fadjlist[v] = collect(V, outneighbors(gv, v))
        vvalues = get_vertexval(gv, v, :)
        for (i, value) ∈ enumerate(vvalues)
            @inbounds vertexvals[i][v] = value
        end
        for (i, EE) in enumerate(E_VALS.types)
            @inbounds edgevals[i][v] = collect(EE, outedgevals(gv, v, i))
        end
    end

    V_VALS_C = typeof(vertexvals)
    E_VALS_C = typeof(edgevals)

    return ValGraph{V, V_VALS, E_VALS, G_VALS, V_VALS_C, E_VALS_C}(ne(gv), fadjlist, vertexvals, edgevals, get_graphval(gv, :))
end
