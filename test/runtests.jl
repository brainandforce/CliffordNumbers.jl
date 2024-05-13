using CliffordNumbers
using Aqua, Test

Aqua.test_all(CliffordNumbers; unbound_args = false)

@testset "CliffordNumbers.jl" begin
    include("internals.jl")
    include("metrics.jl")
    include("construction.jl")
    include("indexing.jl")
    include("conversion.jl")
    include("operations.jl")
end
