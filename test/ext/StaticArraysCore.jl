@testset "StaticArraysCore extension" begin
    k = KVector{1,VGA(3)}(4, 2, 0)
    @test CliffordNumbers.similar_type(k, Float32) === StaticArrays.similar_type(k, Float32)
end
