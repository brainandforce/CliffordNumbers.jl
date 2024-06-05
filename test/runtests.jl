using CliffordNumbers
using Aqua, Test

Aqua.test_all(CliffordNumbers; unbound_args = false)

# Define basis vectors for important algebras
const σ1 = KVector{1,VGA(3)}(0, 1, 0)
const σ2 = KVector{1,VGA(3)}(0, 0, 1)
const σ3 = KVector{1,VGA(3)}(0, 0, 0)

const e0 = KVector{1,PGA(3)}(1, 0, 0, 0)
const e1 = KVector{1,PGA(3)}(0, 1, 0, 0)
const e2 = KVector{1,PGA(3)}(0, 0, 1, 0)
const e3 = KVector{1,PGA(3)}(0, 0, 0, 1)

const γ0 = KVector{1,STA}(1, 0, 0, 0)
const γ1 = KVector{1,STA}(0, 1, 0, 0)
const γ2 = KVector{1,STA}(0, 0, 1, 0)
const γ3 = KVector{1,STA}(0, 0, 0, 1)

@testset "CliffordNumbers.jl" begin
    include("internals.jl")
    include("metrics.jl")
    include("construction.jl")
    include("indexing.jl")
    include("conversion.jl")
    include("operations.jl")
end
