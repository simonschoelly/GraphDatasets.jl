module GraphDatasets

import Base: getindex, length, eltype, iterate, show, firstindex, lastindex

import Graphs: loadgraphs, SimpleGraph

import SimpleValueGraphs:
    nv, has_edge, is_directed,
    get_vertexval, get_edgeval, get_graphval,
    outedgevals, outneighbors,
    ValGraph

export
    list_datasets,

    loadgraphs,
    loadreadme,

    ValGraphCollection,
    ValGraphCollectionView,
    ng,

    TUDatasets,

    # reexport from Base
    getindex,
    length,
    eltype,
    iterate,
    show,
    firstindex,
    lastindex,


    # reexport from Graphs.jl & SimpleValueGraphs
    nv,
    has_edge,

    get_vertexval,
    get_edgeval,
    get_graphval,

    outedgevals,
    outneighbors,

    SimpleGraph,
    ValGraph

include("utils.jl")
include("graphdataset.jl")
include("valgraphcollection.jl")
include("TUDatasets.jl")

end
