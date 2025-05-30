using CliffordNumbers
using Aqua, Test
using LinearAlgebra, Quaternions, StaticArrays, Unitful

Aqua.test_all(CliffordNumbers; unbound_args = false)

# Define basis vectors for important algebras
@basis_vars(VGA(3), Int, 'σ')
@basis_vars(PGA(3), Int)
@basis_vars(STA)

# A subtype of Number that does not subtype Real or Complex
struct MockNumber <: Number
end

@testset "CliffordNumbers.jl" begin
    include("internals.jl")
    include("metrics.jl")
    include("construction.jl")
    include("indexing.jl")
    include("conversion.jl")
    include("operations.jl")
    include("ext/LinearAlgebra.jl")
    include("ext/Quaternions.jl")
    include("ext/StaticArraysCore.jl")
    include("ext/Unitful.jl")
end
