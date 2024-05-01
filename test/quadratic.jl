@testset "Quadratic forms" begin
    @test dimension(QuadraticForm{3,2,1}) === 6
    @test CliffordNumbers.is_degenerate(QuadraticForm{3,2,1})
    @test !CliffordNumbers.is_degenerate(QuadraticForm{3,2,0})
    @test !CliffordNumbers.is_positive_definite(QuadraticForm{3,2,0})
    @test CliffordNumbers.is_positive_definite(QuadraticForm{3,0,0})
end