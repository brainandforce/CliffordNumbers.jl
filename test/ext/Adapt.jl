@testset "Adapt extension" begin
    # Adapting to a typed array destination retypes the scalar coefficients, mirroring how array
    # storage is retyped on adaptation. The result is still an `isbits` Clifford number.
    k = KVector{2,PGA(3),Float64}(1, 2, 3, 4, 5, 6)
    k32 = adapt(Array{Float32}, k)
    @test k32 isa KVector{2,PGA(3),Float32,6}
    @test Tuple(k32) === (1.0f0, 2.0f0, 3.0f0, 4.0f0, 5.0f0, 6.0f0)
    @test isbits(k32)

    # A destination without a concrete scalar element type leaves the scalar type untouched.
    @test adapt(Array, k) === k
    @test adapt(Array{Float64}, k) === k

    # The behavior holds across the concrete subtypes.
    c = CliffordNumber{VGA(3),Int}(ntuple(identity, Val(8)))
    @test adapt(Array{Float32}, c) === CliffordNumber{VGA(3),Float32}(ntuple(identity, Val(8)))
    e = EvenCliffordNumber{VGA(3)}(1, 2, 3, 4)
    @test adapt(Array{Float32}, e) isa EvenCliffordNumber{VGA(3),Float32,4}
end
