@testset "Unitful extension" begin
    x = KVector{1,VGA(3)}(1, 0, 0)
    y = KVector{1,VGA(3)}(0, 1, 0)
    xu = x*u"m"
    yu = y*u"m"
    @test xu ∧ yu === KVector{2,VGA(3)}(1, 0, 0) * u"m^2"
    @test_throws Unitful.AffineError 2 ∧ 2u"°F"
    @test_throws Unitful.AffineError 2u"°F" ∧ 2
end
