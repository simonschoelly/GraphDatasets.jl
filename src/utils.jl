


cat_tuple_types(T1::Type{<:Tuple}, T2::Type{<:Tuple}) = Tuple{T1.types..., T2.types...}
cat_tuple_types(T1::Type{<:NamedTuple}, ::Type{Tuple{}}) = T1
cat_tuple_types(::Type{Tuple{}}, T2::Type{<:NamedTuple}) = T2
cat_tuple_types(T1::Type{<:NamedTuple}, T2::Type{<:NamedTuple}) =
    NamedTuple{(T1.names..., T2.names...), Tuple{T1.types..., T2.types...}}
cat_tuple_types(T1, T2, T3, T...) = cat_tuple_types(cat_tuple_types(T1, T2), T3, T...)
