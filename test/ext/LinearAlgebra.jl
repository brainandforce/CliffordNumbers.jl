@testset "LinearAlgebra extensions" begin
    k = KVector{1,VGA(3)}(4, 2, 0)
    l = KVector{2,VGA(3)}(0, 6, 9)
    @test CliffordNumbers.dot(k, l) === LinearAlgebra.dot(k, l)
    @test CliffordNumbers.normalize(k) === LinearAlgebra.normalize(k)
end
