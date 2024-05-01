@testset "Metrics" begin
    # Vanilla GAs
    @test Metrics.dimension(Metrics.VGA(3)) === 3
    @test Metrics.dimension(Metrics.VGA(3)) === length(Metrics.VGA(3))
    @test size(Metrics.VGA(3)) === tuple(3)
    @test axes(Metrics.VGA(3)) === tuple(1:3,)
    @test all(isone, Metrics.VGA(3))
    @test Metrics.basis_blades(Metrics.VGA(3)) == 8
    @test !Metrics.is_degenerate(Metrics.VGA(3))
    @test Metrics.is_positive_definite(Metrics.VGA(3))
    # Projective GAs
    @test Metrics.dimension(Metrics.PGA(3)) === 4
    @test Metrics.dimension(Metrics.PGA(3)) === length(Metrics.PGA(3))
    @test size(Metrics.PGA(3)) === tuple(4)
    @test axes(Metrics.PGA(3)) === tuple(0:3)
    @test iszero(first(Metrics.PGA(3)))
    @test all(isone, Metrics.PGA(3)[1:3])
    @test Metrics.is_degenerate(Metrics.PGA(3))
    @test !Metrics.is_positive_definite(Metrics.PGA(3))
    # Conformal GAs
    @test Metrics.dimension(Metrics.CGA(3)) === 5
    @test Metrics.dimension(Metrics.CGA(3)) === length(Metrics.CGA(3))
    @test size(Metrics.CGA(3)) === tuple(5)
    @test axes(Metrics.CGA(3)) === tuple(-1:3)
    @test eachindex(Metrics.CGA(3)) === -1:3
    @test Metrics.CGA(3)[-1] === Int8(-1)
    @test all(isone, Metrics.CGA(3)[0:3])
    @test !Metrics.is_degenerate(Metrics.CGA(3))
    @test !Metrics.is_positive_definite(Metrics.CGA(3))
    # Lorentzian GAs
    @test Metrics.dimension(Metrics.STAEast) === 4
    @test Metrics.dimension(Metrics.STAWest) === 4
    @test size(Metrics.STAEast) === tuple(4)
    @test size(Metrics.STAWest) === tuple(4)
    @test axes(Metrics.STAEast) === tuple(0:3)
    @test axes(Metrics.STAWest) === tuple(0:3)
    @test eachindex(Metrics.STAEast) === 0:3
    @test eachindex(Metrics.STAWest) === 0:3
    @test Metrics.STAEast[0] === Int8(-1)
    @test Metrics.STAWest[0] === Int8(+1)
    @test all(isequal(Int8(+1)), Metrics.STAEast[1:3])
    @test all(isequal(Int8(-1)), Metrics.STAWest[1:3])
    @test !Metrics.is_degenerate(Metrics.STA)
    @test !Metrics.is_positive_definite(Metrics.STA)
    @test Metrics.blade_symbol(Metrics.STA) === 'γ'
    # Exterior algebras
    @test Metrics.Exterior(6, 1) === Metrics.Exterior(6)
    @test Metrics.dimension(Metrics.Exterior(6, 1)) === 6
    @test Metrics.dimension(Metrics.Exterior(6, 0)) === 6
    @test Metrics.dimension(Metrics.Exterior(6, 1)) === length(Metrics.Exterior(6,1))
    @test Metrics.dimension(Metrics.Exterior(6, 0)) === length(Metrics.Exterior(6,0))
    @test size(Metrics.Exterior(6, 1)) == tuple(6)
    @test size(Metrics.Exterior(6, 0)) == tuple(6)
    @test axes(Metrics.Exterior(6, 1)) == tuple(1:6)
    @test axes(Metrics.Exterior(6, 0)) == tuple(0:5)
    @test all(iszero, Metrics.Exterior(6, 1))
    @test all(iszero, Metrics.Exterior(6, 0))
    @test Metrics.is_degenerate(Metrics.Exterior(6, 1))
    @test Metrics.is_degenerate(Metrics.Exterior(6, 0))
    @test !Metrics.is_positive_definite(Metrics.Exterior(6, 1))
    @test !Metrics.is_positive_definite(Metrics.Exterior(6, 0))
end
