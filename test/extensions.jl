using LinearAlgebra, StaticArrays

@testset "LinearAlgebra extensions" begin
    k = KVector{1,VGA(3)}(4, 2, 0)
    l = KVector{2,VGA(3)}(0, 6, 9)
    @test CliffordNumbers.dot(k, l) = LinearAlgebra.dot(k, l)
end

@testset "StaticArraysCore extensions" begin
    k = KVector{1,VGA(3)}(4, 2, 0)
    @test CliffordNumbers.similar_type(k, Float32) = StaticArrays.similar_type(k, Float32)
end
