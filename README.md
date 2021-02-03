# GraphDatasets

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://simonschoelly.github.io/GraphDatasets.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://simonschoelly.github.io/GraphDatasets.jl/dev)
[![Build Status](https://github.com/simonschoelly/GraphDatasets.jl/workflows/CI/badge.svg)](https://github.com/simonschoelly/GraphDatasets.jl/actions)
[![Coverage](https://codecov.io/gh/simonschoelly/GraphDatasets.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/simonschoelly/GraphDatasets.jl)

*GraphDatasets.jl* is a package for downloading and working with graph datasets. It currently provides the
ability to download some of the datasets from [graphlearning.io](https://www.graphlearning.io).

## Quick example

```julia
# Load GraphDatasets
julia> using GraphDatasets

# List all available datasets
julia> list_datasets()
TUDatasets.AIDSDataset
TUDatasets.AspirinDataset
TUDatasets.BZRDataset
TUDatasets.BZR_MDDataset
TUDatasets.BenzeneDataset
TUDatasets.COIL_DELDataset
TUDatasets.COIL_RAGDataset
TUDatasets.COLLABDataset
TUDatasets.COLORS_3Dataset
TUDatasets.DBLP_v1Dataset
TUDatasets.DDDataset
TUDatasets.ENZYMESDataset
TUDatasets.QM9Dataset
TUDatasets.REDDIT_BINARYDataset
TUDatasets.SYNTHETICDataset
TUDatasets.SYNTHETICnewDataset
TUDatasets.SynthieDataset
TUDatasets.TRIANGLESDataset

# Load QM9 from TUDatasets. This dataset contains 129433 molecules represented as graphs.
# The resulting ValGraphCollection is an immutable collection of graphs.
julia> qm9 = loadgraphs(TUDatasets.QM9Dataset())
129433-element ValGraphCollection of graphs with
              eltype: Int8
  vertex value types: (Bool, Bool, Bool, Bool, Bool, Int8, Bool, Bool, Bool, Bool, Bool, Bool, Int64, Float64, Float64, Float64)
    edge value types: (Bool, Bool, Bool, Bool)
   graph value types: (μ = Float64, α = Float64, ϵ_HOMO = Float64, ϵ_LUMO = Float64, Δϵ = Float64, electronic_spatial_energy = Float64, ZPVE = Float64, U_0 = Float64, U = Float64, H = Float64, G = Float64, c_v = Float64, UATOM_0 = Float64, UTAM = Float64, HATOM = Float64, GATOM = Float64, A = Float64, B = Float64, C = Float64)

# We can have a look at the readme for this dataset
julia> loadreadme(TUDatasets.QM9Dataset())
README for dataset QM9
[..]

# A ValGraphCollection can be indexed to get a ValGraphCollectionView of a single graph.
julia> g = qm9[1234]
{19, 18} undirected ValGraphCollectionView with
              eltype: Int8
  vertex value types: (Bool, Bool, Bool, Bool, Bool, Int8, Bool, Bool, Bool, Bool, Bool, Bool, Int64, Float64, Float64, Float64)
    edge value types: (Bool, Bool, Bool, Bool)
   graph value types: (μ = Float64, α = Float64, ϵ_HOMO = Float64, ϵ_LUMO = Float64, Δϵ = Float64, electronig_spatial_energy = Float64, ZPVE = Float64, U_0 = Float64, U = Float64, H = Float64, G = Float64, c_v = Float64, UATOM_0 = Float64, UTAM = Float64, HATOM = Float64, GATOM = Float64, A = Float64, B = Float64, C = Float64)

# ValGraphCollectionView inherits from LightGraphs.AbstractGraph and SimpleValueGraphs.AbstractValGraph
# and can therefore be used like other graph types
julia> using LightGraphs: diameter

julia> diameter(g)
7

# We can also convert it to a SimpleGraph without keeping the metadata
julia> SimpleGraph(g)
{19, 18} undirected simple Int8 graph

# or to a ValGraph
julia> ValGraph(g)
{19, 18} undirected ValGraph with
              eltype: Int8
  vertex value types: (Bool, Bool, Bool, Bool, Bool, Int8, Bool, Bool, Bool, Bool, Bool, Bool, Int64, Float64, Float64, Float64)
    edge value types: (Bool, Bool, Bool, Bool)
   graph value types: (μ = Float64, α = Float64, ϵ_HOMO = Float64, ϵ_LUMO = Float64, Δϵ = Float64, electronig_spatial_energy = Float64, ZPVE = Float64, U_0 = Float64, U = Float64, H = Float64, G = Float64, c_v = Float64, UATOM_0 = Float64, UTAM = Float64, HATOM = Float64, GATOM = Float64, A = Float64, B = Float64, C = Float64)
```

## Alternatives

- [GraphMLDatasets.jl](https://github.com/yuehhua/GraphMLDatasets.jl)
- [SNAPDatasets.jl](https://github.com/JuliaGraphs/SNAPDatasets.jl)
- [LightGraphsExtras.j](https://github.com/JuliaGraphs/LightGraphsExtras.jl)
- [MatrixMarket.jl](https://github.com/JuliaSparse/MatrixMarket.jl)
- [MatrixDepot.jl](https://github.com/JuliaMatrices/MatrixDepot.jl)

## References

[1]: [Morris, Christopher, Nils M. Kriege, Franka Bause, Kristian Kersting, Petra Mutzel, and Marion Neumann. "Tudataset: A collection of benchmark datasets for learning with graphs." arXiv preprint arXiv:2007.08663 (2020).](https://arxiv.org/pdf/2007.08663.pdf)
