using GraphDatasets: cat_tuple_types

@testset "cat_tuple_types" begin

    @test cat_tuple_types(Tuple{}, Tuple{}) == Tuple{}
    @test cat_tuple_types(Tuple{}, Tuple{Int, String}) == Tuple{Int, String}
    @test cat_tuple_types(Tuple{Int, String}, Tuple{}) == Tuple{Int, String}

    @test cat_tuple_types(Tuple{}, @NamedTuple{a::String, b::Int}) == @NamedTuple{a::String, b::Int}
    @test cat_tuple_types(@NamedTuple{a::String, b::Int}, Tuple{}) == @NamedTuple{a::String, b::Int}

    @test cat_tuple_types(Tuple{Int, String}, Tuple{Char, Vector, Tuple}) == Tuple{Int, String, Char, Vector, Tuple}

    @test cat_tuple_types(@NamedTuple{a::Int, b::String}, @NamedTuple{c::Char, d::Vector, e::Tuple}) ==
        @NamedTuple{a::Int, b::String, c::Char, d::Vector, e::Tuple}

    @test cat_tuple_types(Tuple{Int}, Tuple{}, Tuple{Char, Int}, Tuple{String}) ==
            Tuple{Int, Char, Int, String}

    @test cat_tuple_types(@NamedTuple{a::String}, Tuple{}, @NamedTuple{b::Int}, @NamedTuple{c::String, e::Char}) ==
            @NamedTuple{a::String, b::Int, c::String, e::Char}
end
