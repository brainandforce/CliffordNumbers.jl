@testset "StaticArraysCore extension" begin
    k = KVector{1,VGA(3)}(4, 2, 0)
    T = typeof(k)
    @test CliffordNumbers.similar_type(k, Float32) === StaticArrays.similar_type(k, Float32)
    @test CliffordNumbers.similar_type(T, Float32) === StaticArrays.similar_type(T, Float32)
end
