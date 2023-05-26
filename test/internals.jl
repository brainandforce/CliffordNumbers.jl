@testset "Internals" begin
    @test CliffordNumbers.intlog2(3//4) === round(Int, log2(3//4))
    @test CliffordNumbers.intlog2(0.75) === round(Int, log2(0.75))
    @test CliffordNumbers.intlog2(1) === 0
    @test CliffordNumbers.intlog2(1//1) === 0
    @test CliffordNumbers.intlog2(1.0) === 0
    @test CliffordNumbers.intlog2(3//2) === round(Int, log2(3//2))
    @test CliffordNumbers.intlog2(1.5) === round(Int, log2(1.5))
    # Should always pass unless we specialize on integers
    @test all(CliffordNumbers.intlog2(n) === CliffordNumbers.intlog2(Float64(n)) for n in 0:65535)
end
