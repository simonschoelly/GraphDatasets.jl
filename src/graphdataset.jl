
using DataDeps: DataDep, register, unpack, @datadep_str

"""
   GraphDataset

Abstract graph dataset type.

## mandatory methods to implement
* [`loadgraph`](@ref) and/or [`loadgraphs`](@ref)
* [`dataset_root_directory_name`](@ref)
* [`dataset_url`](@ref)
* [`dataset_hash`](@ref)
* [`dataset_message`](@ref)

## optional methods to implement
* [`loadreadme`](@ref)
* [`dataset_path`](@ref)
"""
abstract type GraphDataset end

## ----------------------------------------
##    methods to implement
## ----------------------------------------

function dataset_root_directory_name end

function dataset_url end

dataset_hash(ds::GraphDataset) = nothing

function dataset_message(ds::GraphDataset)

    @warn "No message specified for dataset $ds"
    return ""
end

"""
    loadgraphs(ds::GraphDataset)

Loads multiple graphs from a dataset `ds`.

Either this method or `loadgraph` must be implement for new subytpes of `GraphDataset`.

### See also
[`loadgraph`](@Ref), `GraphDataset`](@Ref)
"""
loadgraphs(ds::GraphDataset)

## ----------------------------------------
##    optional methods to implement
## ----------------------------------------

function loadreadme end

dataset_path(ds::GraphDataset) = @datadep_str(dataset_root_directory_name(ds))

## ----------------------------------------
##    various methods
## ----------------------------------------

_registered_datasets = IdDict{Type{<:GraphDataset}, String}()

"""
    list_datasets()

List available graph datasets
"""
function list_datasets()

    foreach(println, sort(collect(values(_registered_datasets))))
    return nothing
end


function register_dataset(ds::GraphDataset, prefix_modules = String[])

    register(DataDep(
        dataset_root_directory_name(ds),
        dataset_message(ds::GraphDataset),
        dataset_url(ds),
        dataset_hash(ds),
        post_fetch_method = unpack
    ))

    name_in_list = join(map(s -> "$s.", prefix_modules)) * string(nameof(typeof(ds)))
    _registered_datasets[typeof(ds)] = name_in_list
end

