using GraphDatasets
using Test

include("utils.jl")

if parse(Bool, get(ENV, "CI", "false"))
    @info "CI detected: skipping some test"
else
    include("TUDatasets.jl")
end

