
module TUDatasets

using DataDeps: @datadep_str
using CSV

using ..GraphDatasets
using ..GraphDatasets: GraphDataset, register_dataset, cat_tuple_types
import GraphDatasets: loadgraphs, loadreadme, dataset_root_directory_name, dataset_hash, dataset_url, dataset_message

# ======================================
#   TUDataset
# ======================================

"""
    TUDataset <: GraphDataset

Abstract graph dataset type for the graph collections from graphlearning.io.
"""
abstract type TUDataset <: GraphDataset end

## ----------------------------------------
##    methods to implement
## ----------------------------------------

"""
    dataset_name(::TUDataset)

The name of the dataset.

The URL and the directory where such dataset is stored is then
automatically derived from  that name.

This method should be implemented for subtypes of `TUDataset`.

### See also
[`TUDataset`](@ref)
"""
dataset_name(::TUDataset)

## ----------------------------------------
##    optional methods to implement
## ----------------------------------------

graph_eltype(::TUDataset) = Int8

node_labels_type(::TUDataset) = Tuple{}
node_labels_map(ds::TUDataset, i) = i
node_attributes_type(::TUDataset) = Tuple{}

edge_labels_type(::TUDataset) = Tuple{}
edge_labels_map(ds::TUDataset, i) = i
edge_attributes_type(::TUDataset) = Tuple{}

graph_labels_type(::TUDataset) = Tuple{}
graph_labels_map(ds::TUDataset, i) = i
graph_attributes_type(::TUDataset) = Tuple{}

readme_name(::TUDataset) = missing

dataset_references(::TUDataset) = Int[]

## ----------------------------------------
##    GraphDataset implementation
## ----------------------------------------

function dataset_message(ds::TUDataset)

    if isempty(dataset_references(ds))
        @warn "No references specified for dataset $ds"
    end

    graphlearning_reference = "Morris, Christopher, Nils M. Kriege, Franka Bause, Kristian Kersting, Petra Mutzel, and Marion Neumann. \"Tudataset: A collection of benchmark datasets for learning with graphs.\" arXiv preprint arXiv:2007.08663 (2020)."

    bibs = [graphlearning_reference; bibliography[dataset_references(ds)]]
    reference_list = join(map(x -> "[$(x[1])]: $(x[2])", enumerate(bibs)), "\n\n")

    return """
    The $(dataset_name(ds)) dataset is part of TUDatasets [1] that can be downloaded from www.graphlearning.io.

    == References ==
    $reference_list
    """
end

dataset_path(ds::TUDataset) = joinpath(@datadep_str(dataset_root_directory_name(ds)), dataset_name(ds))

dataset_root_directory_name(ds::TUDataset) = "tudatasets-$(dataset_name(ds))"

dataset_url(ds::TUDataset) = "https://www.chrsmrrs.com/graphkerneldatasets/$(dataset_name(ds)).zip"

# ======================================
#   loading
# ======================================


function loadreadme(ds::TUDataset)

    path = joinpath(dataset_path(ds), readme_name(ds))
    return Text(read(path, String))
end

full_vertexval_types(ds::TUDataset) = cat_tuple_types(node_labels_type(ds), node_attributes_type(ds))
full_edgeval_types(ds::TUDataset) = cat_tuple_types(edge_labels_type(ds), edge_attributes_type(ds))
full_graphval_types(ds::TUDataset) = cat_tuple_types(graph_labels_type(ds), graph_attributes_type(ds))

prefix(ds::TUDataset) = dataset_name(ds) * '_'

function load_edgelist(ds::TUDataset)

    # TODO find a way to specify correct integer type
    edgelist_path = joinpath(dataset_path(ds), prefix(ds) * "A.txt")
    return CSV.File(edgelist_path, header=false, strict=true, types=[Int32, Int32])
end

function load_graphindicator(ds::TUDataset)

    path = joinpath(dataset_path(ds), prefix(ds) * "graph_indicator.txt")
    return CSV.File(path, header=false, strict=true, type=Int32)
end

function load_node_labels(ds::TUDataset)

    path = joinpath(dataset_path(ds), prefix(ds) * "node_labels.txt")
    isfile(path) || return nothing
    return CSV.File(path, header=false, strict=true, type=Int32)
end

function load_node_attributes(ds::TUDataset)

    node_attributes_path = joinpath(dataset_path(ds), prefix(ds) * "node_attributes.txt")
    isfile(node_attributes_path) || return nothing
    return CSV.File(node_attributes_path, header=false, strict=true, types=[node_attributes_type(ds).types...])
end

function load_full_vertexvals(ds::TUDataset, n)

    node_labels = load_node_labels(ds::TUDataset)
    node_attributes = load_node_attributes(ds::TUDataset)

    V_VALS = full_vertexval_types(ds)

    vertexvals = Vector{V_VALS}(undef, n)

    for i ∈ 1:n
        label_i = node_labels == nothing ? () : tuple(node_labels_map(ds, node_labels[i][1]))
        attr_i = node_attributes == nothing ? () : (node_attributes[i]...,)
        vertexvals[i] = V_VALS((label_i..., attr_i...))
    end

    return vertexvals
end

function load_edge_labels(ds::TUDataset)

    path = joinpath(dataset_path(ds), prefix(ds) * "edge_labels.txt")
    isfile(path) || return nothing
    return CSV.File(path, header=false, strict=true, type=Int8)
end

function load_edge_attributes(ds::TUDataset)

    path = joinpath(dataset_path(ds), prefix(ds) * "edge_attributes.txt")
    isfile(path) || return nothing
    return CSV.File(path, header=false, strict=true, types=[edge_attributes_type(ds).types...])
end

function load_full_edgevals(ds::TUDataset, m)

    edge_labels = load_edge_labels(ds::TUDataset)
    edge_attributes = load_edge_attributes(ds::TUDataset)

    E_VALS = full_edgeval_types(ds)

    edgevals = Vector{E_VALS}(undef, m)

    for i ∈ 1:m
        label_i = edge_labels == nothing ? () : tuple(edge_labels_map(ds, edge_labels[i][1]))
        attr_i = edge_attributes == nothing ? () : (edge_attributes[i]...,)
        edgevals[i] = E_VALS((label_i..., attr_i...))
    end

    return edgevals
end

function load_graph_labels(ds::TUDataset)

    T = graph_labels_type(ds)
    path = joinpath(dataset_path(ds), prefix(ds) * "graph_labels.txt")
    if T == Tuple{}
        if isfile(path)
            @warn "TODO warn message"
        end
        return nothing
    end
    # TODO throw exception if there is no file but labels type defined
    return CSV.File(path, header=false, strict=true, type=Int8)
end

function load_graph_attributes(ds::TUDataset)

    path = joinpath(dataset_path(ds), prefix(ds) * "graph_attributes.txt")
    isfile(path) || return nothing
    return CSV.File(path, header=false, strict=true, delim=',', types=[graph_attributes_type(ds).types...])
end

function load_full_graphvals(ds::TUDataset, N)

    graph_labels = load_graph_labels(ds::TUDataset)
    graph_attributes = load_graph_attributes(ds::TUDataset)

    G_VALS = full_graphval_types(ds)

    graphvals = Vector{G_VALS}(undef, N)

    for i ∈ 1:N
        label_i = graph_labels == nothing ? () : tuple(graph_labels_map(ds, graph_labels[i][1]))
        attr_i = graph_attributes == nothing ? () : (graph_attributes[i]...,)
        graphvals[i] = G_VALS((label_i..., attr_i...))
    end

    return graphvals
end


# ======================================
#    datasets
# ======================================

function __init__()

    foreach(ds -> register_dataset(ds, [string(nameof(@__MODULE__))]), [
        AIDSDataset(),
        AspirinDataset(),
        BenzeneDataset(),
        BZRDataset(),
        BZR_MDDataset(),
        MutagenicityDataset(),
        MUTAGDataset(),
        DDDataset(),
        ENZYMESDataset(),
        COIL_DELDataset(),
        COIL_RAGDataset(),
        COLLABDataset(),
        DBLP_v1Dataset(),
        REDDIT_BINARYDataset(),
        COLORS_3Dataset(),
        QM9Dataset(),
        SYNTHETICDataset(),
        SYNTHETICnewDataset(),
        SynthieDataset(),
        TRIANGLESDataset()
       ])

    return nothing
end

## --------------------------------------
##    AIDS
## --------------------------------------

struct AIDSDataset <: TUDataset end

dataset_name(::AIDSDataset) = "AIDS"

dataset_hash(::AIDSDataset) = "ef65a8095846588ffd6e17f95c1968a247d8ada7295a61209c401bda23d19ab9"

dataset_references(::AIDSDataset) = [16, 17]

node_labels_type(::AIDSDataset) = NamedTuple{(:symbol,), Tuple{String}}
node_labels_map(::AIDSDataset, i) = (
               "C",  "O",  "N",  "Cl", "F",  "S",  "Se", "P", "Na", "I",  "Co", "Br", "Li",
                "Si", "Mg", "Cu", "As", "B",  "Pt", "Ru", "K",  "Pd", "Au", "Te", "W",  "Rh",
                "Zn", "Bi", "Pb", "Ge", "Sb", "Sn", "Ga", "Hg", "Ho", "Tl", "Ni", "Tb")[i + 1]
node_attributes_type(::AIDSDataset) = NamedTuple{(:chem, :charge, :x, :y), NTuple{4, Float64}}

edge_labels_type(::AIDSDataset) = NamedTuple{tuple(:valence), Tuple{Int8}}
edge_labels_map(::AIDSDataset, i) = (1, 2, 3)[i + 1]

graph_labels_type(::AIDSDataset) = NamedTuple{(:class,), Tuple{String}}
graph_labels_map(::AIDSDataset, i) = ("active", "inactive")[i + 1]

## --------------------------------------
##    BZR
## --------------------------------------

struct BZRDataset <: TUDataset end

dataset_name(::BZRDataset) = "BZR"

dataset_hash(::BZRDataset) = "584cb76f63cf5d0459cabfa254603b4e89823de0c03c96d60ab338f2d71a8b55"

readme_name(::BZRDataset) = "README.txt"

dataset_references(::BZRDataset) = [7]

graph_labels_type(::BZRDataset) = Tuple{Int8}

node_labels_type(::BZRDataset) = Tuple{Int8}
node_attributes_type(::BZRDataset) = NTuple{3, Float64}

## --------------------------------------
##    BZR_MD
## --------------------------------------

struct BZR_MDDataset <: TUDataset end

dataset_name(::BZR_MDDataset) = "BZR_MD"

dataset_hash(::BZR_MDDataset) = "f2ae267e19de998358bddeb20517768380ffb0594c594becbbd11aed538014d6"

readme_name(::BZR_MDDataset) = "README.txt"

dataset_references(::BZR_MDDataset) = [7, 23]

graph_labels_type(::BZR_MDDataset) = Tuple{Int8}

edge_labels_type(::BZR_MDDataset) = NamedTuple{(:bond_type, ), Tuple{String}}
edge_labels_map(::BZR_MDDataset, i) = ("aromatic", "no chemical bound", "single", "double", "triple")[i + 1]
edge_attributes_type(::BZR_MDDataset) = NamedTuple{(:distance,), Tuple{Float64}}

node_labels_type(::BZR_MDDataset) = NamedTuple{(:atom_type,), Tuple{String}}
node_labels_map(::BZR_MDDataset, i) = ("C", "N", "O", "F", "Cl", "S", "P", "BR")[i + 1]

## --------------------------------------
##    aspirin
## --------------------------------------

struct AspirinDataset <: TUDataset end

dataset_name(::AspirinDataset) = "aspirin"

dataset_hash(::AspirinDataset) = "b01b6b73841768670a3a79d04e04a79bb6027ee87a335b4cd3e073fbf6523fb4"

dataset_references(::AspirinDataset) = [36]

node_labels_type(::AspirinDataset) = NamedTuple{(:symbol,), Tuple{String}}
node_labels_map(::AspirinDataset, i) = ("C", "O", "H")[i + 1]

node_attributes_type(::AspirinDataset) = NamedTuple{(:x_coordinate, :y_coordinate, :z_coordinate, :atom_force_x, :atom_force_y, :atom_force_z), NTuple{6, Float64}}

graph_attributes_type(::AspirinDataset) = NamedTuple{(:total_energy,), Tuple{Float64}}

readme_name(::AspirinDataset) = "readme.txt"

## --------------------------------------
## benzene
## --------------------------------------

struct BenzeneDataset <: TUDataset end

dataset_name(::BenzeneDataset) = "benzene"

dataset_hash(::BenzeneDataset) = "d82d921cb20453b12c3e288d7b25e46d5a6c60052e8e276558078fb31e3503f1"

dataset_references(::BenzeneDataset) = [36]

readme_name(::BenzeneDataset) = "readme.txt"

node_attributes_type(::BenzeneDataset) = NamedTuple{(:x_coordinate, :y_coordinate, :z_coordinate, :atom_force_x, :atom_force_y, :atom_force_z), NTuple{6, Float64}}

node_labels_type(::BenzeneDataset) = NamedTuple{(:symbol,), Tuple{String}}
node_labels_map(::BenzeneDataset, i) = ("C", "O", "H")[i + 1]

graph_attributes_type(::BenzeneDataset) = NamedTuple{(:total_energy,), Tuple{Float64}}

## --------------------------------------
## Mutagenicity
## --------------------------------------

struct MutagenicityDataset <: TUDataset end

dataset_name(::MutagenicityDataset) = "Mutagenicity"

dataset_hash(::MutagenicityDataset) = "6230f94ba246b76834fb51ffa138370477b7bf8a784ade92c5e0586780d2ae0e"

dataset_references(::MutagenicityDataset) = [16, 20]

readme_name(::MutagenicityDataset) = "Mutagenicity_label_readme.txt"

node_labels_type(::MutagenicityDataset) = NamedTuple{(:chem,), Tuple{String}}
node_labels_map(::MutagenicityDataset, i) = ("C", "O", "Cl", "H", "N", "F", "Br", "S", "P", "I", "Na", "K", "Li", "Ca")[i + 1]

edge_labels_type(::MutagenicityDataset) = NamedTuple{(:valence,), Tuple{Int8}}
edge_labels_map(::MutagenicityDataset, i) = (1, 2, 3)[i + 1]

graph_labels_type(::MutagenicityDataset) = NamedTuple{(:class,), Tuple{String}}
graph_labels_map(::MutagenicityDataset, i) = ("mutagen", "nonmutagen")[i + 1]

## --------------------------------------
## MUTAG
## --------------------------------------

struct MUTAGDataset <: TUDataset end

dataset_name(::MUTAGDataset) = "MUTAG"

dataset_hash(::MUTAGDataset) = "c419bdc853c367d2d83da4973c45100954ae15e10f5ae2cddde6ca431f8207f6"

dataset_references(::MUTAGDataset) = [1, 23]

readme_name(::MUTAGDataset) = "README.txt"

node_labels_type(::MUTAGDataset) = NamedTuple{(:chem,), Tuple{String}}
node_labels_map(::MUTAGDataset, i) = ("C", "N", "O", "F", "I", "Cl", "Br")[i + 1]

edge_labels_type(::MUTAGDataset) = NamedTuple{(:bond_type,), Tuple{String}}
edge_labels_map(::MUTAGDataset, i) = ("aromatic", "single", "double", "triple")[i + 1]

graph_labels_type(::MUTAGDataset) = Tuple{Int8}


## --------------------------------------
##    QM9
## --------------------------------------

struct QM9Dataset <: TUDataset end

dataset_name(::QM9Dataset) = "QM9"

dataset_hash(::QM9Dataset) = "fd1e631ba58bc35cb11f67e0184e85237423e34c53511a5720fbf189e074a251"

readme_name(::QM9Dataset) = "README.txt"

dataset_references(::QM9Dataset) = [33, 34, 35]

# TODO  what these mean can probably found in https://arxiv.org/pdf/1704.01212.pdf
node_attributes_type(::QM9Dataset) = cat_tuple_types(
    NTuple{5, Bool},
    Tuple{Int8},
    NTuple{6, Bool},
    Tuple{Int},
    NTuple{3, Float64})

# TODO it is not really clear what these attributes signify
edge_attributes_type(::QM9Dataset) = NTuple{4, Bool}

graph_attributes_type(::QM9Dataset) = NamedTuple{(:μ, :α, :ϵ_HOMO, :ϵ_LUMO, :Δϵ, :electronic_spatial_energy,
        :ZPVE, :U_0, :U, :H, :G, :c_v, :UATOM_0, :UTAM, :HATOM, :GATOM, :A, :B, :C), NTuple{19, Float64}}

## --------------------------------------
##    DD
## --------------------------------------

struct DDDataset <: TUDataset end

dataset_name(::DDDataset) = "DD"

dataset_hash(::DDDataset) = "d033d7aeb1a48c4b2b47cf0390e7fc9671de70d98c8df5b11e458d2ec5515cd2"

dataset_references(::DDDataset) = [6, 22]

readme_name(::DDDataset) = "README.txt"

graph_eltype(::DDDataset) = Int16

node_labels_type(::DDDataset) = Tuple{Int8}
graph_labels_type(::DDDataset) = Tuple{Int8}

## --------------------------------------
##    ENZYMES
## --------------------------------------

struct ENZYMESDataset <: TUDataset end

dataset_name(::ENZYMESDataset) = "ENZYMES"

dataset_hash(::ENZYMESDataset) = "13d832eb6ffa084192daf6e5750250028a18437ee692c38d29a10cd60e18aaf4"

readme_name(::ENZYMESDataset) = "README.txt"

dataset_references(::ENZYMESDataset) = [4, 5]

graph_labels_type(::ENZYMESDataset) = Tuple{Int8}

node_labels_type(::ENZYMESDataset) = Tuple{Int8}
node_attributes_type(::ENZYMESDataset) = NTuple{18, Float64}

## --------------------------------------
##    COIL-DEL
## --------------------------------------

struct COIL_DELDataset <: TUDataset end

dataset_name(::COIL_DELDataset) = "COIL-DEL"

dataset_hash(::COIL_DELDataset) = "4dbaf79fc1d1b6fc90911b15faa086af0b5f2ce4600ecd93403af95a34e6e3e8"

dataset_references(::COIL_DELDataset) = [16, 18]

readme_name(::COIL_DELDataset) = "COIL-DEL_label_readme.txt"

node_attributes_type(::COIL_DELDataset) = NamedTuple{(:x, :y), Tuple{Float32, Float32}}

edge_labels_type(::COIL_DELDataset) = NamedTuple{(:valence,), Tuple{Int8}}
edge_labels_map(::COIL_DELDataset, i) = (2, 1)[i + 1] # 0 => 2, 1 => 1

graph_labels_type(::COIL_DELDataset) = Tuple{Int8} # TODO not sure what the labels mean

## --------------------------------------
##    COIL-RAG
## --------------------------------------

# TODO COIL-RAG contains multiple graphs without edges - one should analyze if that
# was done on purpose

struct COIL_RAGDataset <: TUDataset end

dataset_name(::COIL_RAGDataset) = "COIL-RAG"

dataset_hash(::COIL_RAGDataset) = "47f8be52a845c414299a4c9e6346acc074c5438fc12fe90b7d8da5c4613f667f"

dataset_references(::COIL_RAGDataset) = [16, 18]

readme_name(::COIL_RAGDataset) = "COIL-RAG_label_readme.txt"

node_attributes_type(::COIL_RAGDataset) = NTuple{64, Float64}

edge_attributes_type(::COIL_RAGDataset) = NamedTuple{(:boundary,), Tuple{Float32}}

graph_labels_type(::COIL_RAGDataset) = Tuple{Int8} # TODO not sure what the labels mean

## --------------------------------------
##    COLLAB
## --------------------------------------

struct COLLABDataset <: TUDataset end

dataset_name(::COLLABDataset) = "COLLAB"

dataset_hash(::COLLABDataset) = "cb5b772bc4f3a4690d0601378a3fc993b092aeafd63e2b0f2330137dda1cdcc6"

dataset_references(::COLLABDataset) = [14]

graph_eltype(::COLLABDataset) = UInt128

graph_labels_type(::COLLABDataset) = Tuple{Int8} # TODO not sure what the labels mean

## --------------------------------------
##    DBLP_v1
## --------------------------------------

struct DBLP_v1Dataset <: TUDataset end

dataset_name(::DBLP_v1Dataset) = "DBLP_v1"

dataset_hash(::DBLP_v1Dataset) = "67d8a383e8920e9f9e6d8afd55df4104619252f2d08772ee08ad8a76e57b547a"

dataset_references(::DBLP_v1Dataset) = [26]

dataset_readme(::DBLP_v1Dataset) = "readme.txt"

graph_labels_type(::DBLP_v1Dataset) = Tuple{Float32}
graph_labels_map(::DBLP_v1Dataset, i) = (1.0f0, -1.0f0)[i + 1]

edge_labels_type(::DBLP_v1Dataset) = Tuple{String}
edge_labels_map(::DBLP_v1Dataset, i) = ("P2P", "P2W", "W2W")[i + 1]

# TODO there is actually a node map with over 40000 entries defined in readme.txt
node_labels_type(::DBLP_v1Dataset) = Tuple{UInt16}

## --------------------------------------
##    REDDIT-BINARY
## --------------------------------------

struct REDDIT_BINARYDataset <: TUDataset end

dataset_name(::REDDIT_BINARYDataset) = "REDDIT-BINARY"

dataset_hash(::REDDIT_BINARYDataset) = "982dc64ddade42a6365ba9c93f5080d0b9d315df0df492b57fd9fbb7e40dbe16"

dataset_references(::REDDIT_BINARYDataset) = [14]

graph_eltype(::REDDIT_BINARYDataset) = Int16

graph_labels_type(::REDDIT_BINARYDataset) = Tuple{Int8}

## --------------------------------------
##    COLORS-3
## --------------------------------------

struct COLORS_3Dataset <: TUDataset end

dataset_name(::COLORS_3Dataset) = "COLORS-3"

dataset_hash(::COLORS_3Dataset) = "380d8f9e03c73455a2b280f859b136c53c80544eabd83c801a9b5c97d7830bfc"

dataset_references(::COLORS_3Dataset) = [27]

graph_eltype(::COLORS_3Dataset) = UInt128

graph_attributes_type(::COLORS_3Dataset) = Tuple{Int8} # TODO not sure what the labels mean

node_attributes_type(::COLORS_3Dataset) = NTuple{5, Bool} # TODO not sure what the attributes mean

## --------------------------------------
##    SYNTHETIC
## --------------------------------------

struct SYNTHETICDataset <: TUDataset end

dataset_name(::SYNTHETICDataset) = "SYNTHETIC"

readme_name(::SYNTHETICDataset) = "README.txt"

dataset_references(::SYNTHETICDataset) = [3]

dataset_hash(::SYNTHETICDataset) = "3c0344e0cd6518d8b3f52bf45d152b1a9a007a523f5e91fc2bda929cc353c84d"

graph_labels_type(::SYNTHETICDataset) = Tuple{Bool}
graph_labels_map(::SYNTHETICDataset, i) = Bool(i)

node_labels_type(::SYNTHETICDataset) = Tuple{Int8}
node_attributes_type(::SYNTHETICDataset) = Tuple{Float64}

## --------------------------------------
##    SYNTHETICnew
## --------------------------------------

struct SYNTHETICnewDataset <: TUDataset end

dataset_name(::SYNTHETICnewDataset) = "SYNTHETICnew"

dataset_hash(::SYNTHETICnewDataset) = "07e27d6ff1c25d036df5bf3593b1bf4676f51e656eb982e903dea4718634ed5e"

dataset_references(::SYNTHETICnewDataset) = [3, 10]

graph_labels_type(::SYNTHETICnewDataset) = Tuple{Int8}

node_attributes_type(::SYNTHETICnewDataset) = Tuple{Float64}

## --------------------------------------
##    Synthie
## --------------------------------------

struct SynthieDataset <: TUDataset end

dataset_name(::SynthieDataset) = "Synthie"

readme_name(::SynthieDataset) = "README.txt"

dataset_references(::SynthieDataset) = [21]

dataset_hash(::SynthieDataset) = "c5bada5ffe42b4a901d50e75f10c3f969fb962f94acaab0265f46097496154d5"

graph_labels_type(::SynthieDataset) = Tuple{Int8}

node_attributes_type(::SynthieDataset) = NTuple{15, Float64}

## --------------------------------------
##    TRIANGLES
## --------------------------------------

struct TRIANGLESDataset <: TUDataset end

dataset_name(::TRIANGLESDataset) = "TRIANGLES"

dataset_references(::TRIANGLESDataset) = [27]

readme_name(::TRIANGLESDataset) = "README.txt"

dataset_hash(::TRIANGLESDataset) = "d14094eecf75fd60cf08b9d18d33a7e8c7657ff07e554956f0230fba9bf63b60"

graph_labels_type(::TRIANGLESDataset) = Tuple{Int8}
graph_attributes_type(::TRIANGLESDataset) = Tuple{Int8}

node_attributes_type(::TRIANGLESDataset) = Tuple{Int8}

# ======================================
#   loadgraphs
# ======================================

#= The format for of the graphs here is describes in
 #  https://chrsmrrs.github.io/datasets/docs/format/, although some of the datasets
 #  deviate from that format. On the other hand certain, several assumptions are taken
 #  here that are not implied by that documentation. These assumptions should
 #  be documented and there should be some code that verifies that the datasets
 #  actually correspond to these assumptions
=#

# TODO this is quite ugly
# TODO fix type instabilities
# TODO check data for inconsistencies
function loadgraphs(ds::TUDataset)

    edgelist        = load_edgelist(ds)
    graph_indicator = load_graphindicator(ds)

    n = length(graph_indicator) # number of vertices
    m = length(edgelist) # number of edges
    N = maximum(row -> row[1], graph_indicator) # number of graphs

    vertex_ids = zeros(Int, n + 1)
    edge_ids = zeros(graph_eltype(ds), m)
    graph_ids = zeros(Int, N + 1)

    vertexvals = load_full_vertexvals(ds, n)
    edgevals = load_full_edgevals(ds, m)
    graphvals = load_full_graphvals(ds, N)

    @assert length(vertexvals) == n
    @assert length(edgevals) == m
    @assert length(graphvals) == N

    k = 0
    for (v_id, row) in enumerate(graph_indicator)

        k2 = row[1]
        @assert k2 >= max(k, 1)

        while k < k2
            k += 1
            graph_ids[k] = v_id
        end
    end
    graph_ids[end] = n + 1

    edges = Tuple{Int, Int}[] # src, dst, position in edge list
    i = 1
    j = 1
    k = 1
    delta = 0
    perm = Int[]
    for edge_num in 1:m+1

        if edge_num <= m
            src, dst = edgelist[edge_num][1], edgelist[edge_num][2]
            k_next = graph_indicator[src][1]
        end

        # processing k-th graph
        if edge_num == m+1 || k_next > k

            # at this point, edges contains exactly the edges of the k-th graph
            # but they are not in lexicographical ascending order and
            # there might be duplicates. We need to sort the edges and the edgevalues
            # in the same way.
            resize!(perm, length(edges))
            sortperm!(perm, edges)
            permute!(edges, perm)

            e_start = edge_num - length(edges)
            e_end = edge_num - 1 # edgevals[edge_num] is already part of the next graph

            @assert length(e_start:e_end) == length(edges) "k == $k : $(e_start:e_end)  != $(length(edges))"
            permute!(view(edgevals, e_start:e_end), perm)

            v_start = graph_ids[k]

            s_prev = 0
            d_prev = 0

            j = graph_ids[k] - 1
            for (s, d) in edges

                while s > j
                    j += 1
                    vertex_ids[j] = i - delta
                end

                # processing vertex j

                if (s == s_prev && d == d_prev)
                    delta += 1
                end
                s_prev = s
                d_prev = d
                edge_ids[i - delta] = d - (v_start - 1)
                edgevals[i - delta] = edgevals[i]
                i += 1
            end
            # in case the highest vertices of a graph are isolated
            while j < graph_ids[k + 1] - 1
                j += 1
                vertex_ids[j] = i - delta
            end

            if edge_num <= m
                k += 1
                empty!(edges)
            end
        end

        if edge_num <= m
            push!(edges, (src, dst))
        end
    end
    m -= delta
    vertex_ids[end] = m + 1

    if delta > 0
        @warn "Dataset $ds contains duplicate edges -> removing duplicates"

        resize!(edge_ids, m)
        resize!(edgevals, m)
    end

    graphs = ValGraphCollection(graph_ids, vertex_ids, edge_ids, vertexvals, edgevals, graphvals)

    return graphs

end

# ===============================
#     bibliography_entries
# ===============================

# TODO if we use markdown for entries, we could also embed URLs, or we could have
# a tuple where the second entry is an URL
# TODO use a dict with string identifiers would probably help to avoid mistakes if
# the order in the published biography ever changes
bibliography = [
    #= 1  =# "Debnath, A.K., Lopez de Compadre, R.L., Debnath, G., Shusterman, A.J., and Hansch, C. Structure-activity relationship of mutagenic aromatic and heteroaromatic nitro compounds. J. Med. Chem. 34(2):786-797 (1991).",
    #= 2  =# "Helma, C., King, R. D., Kramer, S., and Srinivasan, A. The Predictive Toxicology Challenge 2000–2001. Bioinformatics, 2001, 17, 107-108. The Predictive Toxicology Challenge",
    #= 3  =# "Feragen, A., Kasenburg, N., Petersen, J., de Bruijne, M., Borgwardt, K.M.: Scalable kernels for graphs with continuous attributes. In: C.J.C. Burges, L. Bottou, Z. Ghahramani, K.Q. Weinberger (eds.) NIPS, pp. 216-224 (2013).",
    #= 4  =# "K. M. Borgwardt, C. S. Ong, S. Schoenauer, S. V. N. Vishwanathan, A. J. Smola, and H. P. Kriegel. Protein function prediction via graph kernels. Bioinformatics, 21(Suppl 1):i47–i56, Jun 2005.",
    #= 5  =# "I. Schomburg, A. Chang, C. Ebeling, M. Gremse, C. Heldt, G. Huhn, and D. Schomburg. Brenda, the enzyme database: updates and major new developments. Nucleic Acids Research, 32D:431–433, 2004.",
    #= 6  =# "P. D. Dobson and A. J. Doig. Distinguishing enzyme structures from non-enzymes without alignments. J. Mol. Biol., 330(4):771–783, Jul 2003.",
    #= 7  =# "Sutherland, J. J.; O’Brien, L. A. & Weaver, D. F. Spline-fitting with a genetic algorithm: a method for developing classification structure-activity relationships. J. Chem. Inf. Comput. Sci., 2003, 43, 1906-1915.",
    #= 8  =# "N. Wale and G. Karypis. Comparison of descriptor spaces for chemical compound retrieval and classification. In Proc. of ICDM, pages 678–689, Hong Kong, 2006.",
    #= 9  =# "Pubchem",
    #= 10 =# "http://image.diku.dk/aasa/papers/graphkernels_nips_erratum.pdf",
    #= 11 =# "M. Neumann, P. Moreno, L. Antanas, R. Garnett, K. Kersting. Graph Kernels for Object Category Prediction in Task-Dependent Robot Grasping. Eleventh Workshop on Mining and Learning with Graphs (MLG-13), Chicago, Illinois, USA, 2013.",
    #= 12 =# "http://www.first-mm.eu/data.html",
    #= 13 =# "M. Neumann, R. Garnett, C. Bauckhage, and K. Kersting. Propagation kernels: efficient graph kernels from propagated information. Machine Learning, 102(2):209–245, 2016",
    #= 14 =# "Pinar Yanardag and S.V.N. Vishwanathan. 2015. Deep Graph Kernels. In Proceedings of the 21th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining, ACM, New York, NY, USA, 1365-1374.",
    #= 15 =# "Francesco Orsini, Paolo Frasconi, and Luc De Raedt. 2015 Graph invariant kernels. In Proceedings of the 24th International Conference on Artificial Intelligence (IJCAI’15), Qiang Yang and Michael Wooldridge (Eds.). AAAI Press 3756-3762.",
    #= 16 =# "Riesen, K. and Bunke, H.: IAM Graph Database Repository for Graph Based Pattern Recognition and Machine Learning. In: da Vitora Lobo, N. et al. (Eds.), SSPR&SPR 2008, LNCS, vol. 5342, pp. 287-297, 2008.",
    #= 17 =# "AIDS Antiviral Screen Data (2004)",
    #= 18 =# "S. A. Nene, S. K. Nayar and H. Murase. Columbia Object Image Library, Technical Report, Department of Computer Science, Columbia University CUCS-006-96, Feb. 1996.",
    #= 19 =# "NIST Special Database 4",
    #= 20 =# "Jeroen Kazius, Ross McGuire and, and Roberta Bursi. Derivation and Validation of Toxicophores for Mutagenicity Prediction, Journal of Medicinal Chemistry 2005 48 (1), 312-320",
    #= 21 =# "Christopher Morris, Nils M. Kriege, Kristian Kersting, Petra Mutzel. Faster Kernels for Graphs with Continuous Attributes via Hashing, IEEE International Conference on Data Mining (ICDM) 2016",
    #= 22 =# "Nino Shervashidze, Pascal Schweitzer, Erik Jan van Leeuwen, Kurt Mehlhorn, and Karsten M. Borgwardt. 2011. Weisfeiler-Lehman Graph Kernels. J. Mach. Learn. Res. 12 (November 2011), 2539-2561.",
    #= 23 =# "Nils Kriege, Petra Mutzel. 2012. Subgraph Matching Kernels for Attributed Graphs. International Conference on Machine Learning 2012.",
    #= 24 =# "Tox21 Data Challenge 2014",
    #= 25 =# "Nils M. Kriege, Matthias Fey, Denis Fisseler, Petra Mutzel, Frank Weichert. Recognizing Cuneiform Signs Using Graph Based Methods. International Workshop on Cost-Sensitive Learning (COST), SIAM International Conference on Data Mining (SDM) 2018, 31-44.",
    #= 26 =# "A Repository of Benchmark Graph Datasets for Graph Classification",
    #= 27 =# "Boris Knyazev, Graham W. Taylor, Mohamed R. Amer. Understanding Attention and Generalization in Graph Neural Networks. Neural Information Processing Systems (NeurIPS) 2019, 4204-4214.",
    #= 28 =# "Xifeng Yan, Hong Cheng, Jiawei Han, Philip S. Yu. Mining Significant Graph Patterns by Leap Search. ACM SIGMOD International Conference on Management of Data 2008, 433–444, Chemical Datasets.",
    #= 29 =# "Alchemy: A Quantum Chemistry Dataset for Benchmarking AI Models",
    #= 30 =# "An API Oriented Open-source Python Framework for Unsupervised Learning on Graphs",
    #= 31 =# "Xavier Bresson, Thomas Laurent. A Two-Step Graph Convolutional Decoder for Molecule Generation. Workshop on Machine Learning and the Physical Sciences. 2019.",
    #= 32 =# "Lutz Oettershagen, Nils Kriege, Christopher Morris, Petra Mutzel. Temporal Graph Kernels for Classifying Dissemination Processes. SIAM International Conference on Data Mining (SDM) 2020.",
    #= 33 =# "PyTorch Geometric datasets",
    #= 34 =# "Zhenqin Wu, Bharath Ramsundar, Evan Feinberg, Joseph Gomes, Caleb Geniesse, Aneesh Pappu, Karl Leswing, Vijay Pande. MoleculeNet: A Benchmark for Molecular Machine Learning. Chemical Science. 9. 2017.",
    #= 35 =# "Raghunathan Ramakrishnan, Pavlo Dral, Matthias Rupp, O. Anatole von Lilienfeld. Quantum Chemistry Structures and Properties of 134 kilo Molecules, Scientific Data 1: 140022, 2014.",
    #= 36 =# "Stefan Chmiela, Alexandre Tkatchenko, Huziel E. Sauceda, Igor Poltavsky, Kristof T. Schütt, Klaus-Robert Müller. Machine Learning of Accurate Energy-Conserving Molecular Force Fields. Science Advances 3(5). 2017.",
   ]

end # module