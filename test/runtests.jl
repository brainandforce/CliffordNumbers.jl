using CliffordNumbers
using Aqua, Test, Unitful

Aqua.test_all(CliffordNumbers; unbound_args = false)

# Define basis vectors for important algebras
@basis_vars(VGA(3), Int, 'σ')
@basis_vars(PGA(3), Int)
@basis_vars(STA)

@testset "CliffordNumbers.jl" begin
    include("internals.jl")
    include("metrics.jl")
    include("construction.jl")
    include("indexing.jl")
    include("conversion.jl")
    include("operations.jl")
    if VERSION >= v"1.9.0"
        include("ext/Unitful.jl")
    end
end
