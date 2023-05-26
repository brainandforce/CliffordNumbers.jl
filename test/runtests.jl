using CliffordNumbers
using Test

@testset "CliffordNumbers.jl" begin
    include("internals.jl")
    include("conversion.jl")
    include("operations.jl")
end
