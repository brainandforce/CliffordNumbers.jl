using CliffordNumbers
using Aqua, Test
using LinearAlgebra, Quaternions, StaticArrays, StructArrays, Unitful
using Adapt
using Random
using ChainRulesCore

Aqua.test_all(CliffordNumbers; unbound_args = false)

# Define basis vectors for important algebras
@basis_vars(VGA(3), Int, 'σ')
@basis_vars(PGA(3), Int)
@basis_vars(STA)

# A subtype of Number that does not subtype Real or Complex
struct MockNumber <: Number
end

# The zero-allocation gates assert what the current stable compiler constant-folds. The 1.10 LTS
# misses some of those folds, so the gates are informative only there: results stay correct, the
# heap traffic is just not fully eliminated. Correctness and inference tests run on every version.
const ALLOCATION_GATES = VERSION >= v"1.12"

@testset "CliffordNumbers.jl" begin
    include("internals.jl")
    include("metrics.jl")
    include("construction.jl")
    include("indexing.jl")
    include("conversion.jl")
    include("bitindices.jl")
    include("operations.jl")
    include("determinant.jl")
    include("factorization.jl")
    include("pin_versors.jl")
    include("theorems.jl")
    include("custom_signature.jl")
    include("ext/Adapt.jl")
    include("ext/ChainRulesCore.jl")
    include("ext/LinearAlgebra.jl")
    include("ext/Quaternions.jl")
    include("ext/StaticArraysCore.jl")
    include("ext/StructArrays.jl")
    include("ext/Unitful.jl")
end
