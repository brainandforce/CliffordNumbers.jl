@testset "Quadratic forms" begin
    @test dimension(QuadraticForm{3,2,1}) === 6
    @test CliffordNumbers.isdegenerate(QuadraticForm{3,2,1})
    @test !CliffordNumbers.isdegenerate(QuadraticForm{3,2,0})
    @test !CliffordNumbers.iseuclidean(QuadraticForm{3,2,0})
    @test CliffordNumbers.iseuclidean(QuadraticForm{3,0,0})
end