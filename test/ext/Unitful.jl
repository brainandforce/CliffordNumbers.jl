@testset "Unitful extension" begin
    x = KVector{1,VGA(3)}(1, 0, 0)
    y = KVector{1,VGA(3)}(0, 1, 0)
    @test_throws ArgumentError KVector{1,VGA(3)}(1u"m", 0u"m", 0u"m")
    @test_throws ArgumentError KVector{1,VGA(3)}(0u"m", 1u"m/s", 0u"K")
    xu = x*u"m"
    yu = y*u"m"
    @test xu ∧ yu === KVector{2,VGA(3)}(1, 0, 0) * u"m^2"
    @test 2 ∧ xu === KVector{1,VGA(3)}(2, 0, 0) * u"m"
    @test yu ∧ 2 === KVector{1,VGA(3)}(0, 2, 0) * u"m"
    @test_throws Unitful.AffineError 2 ∧ 2u"°F"
    @test_throws Unitful.AffineError 2u"°F" ∧ 2
    # Printing
    @test Unitful.BracketStyle(x) === Unitful.RoundBrackets()
    @test string(xu) == "(1e₁) m"
    @test string(xu ∧ yu) == "(1e₁e₂) m^2"
    unitless = Quantity{KVector{1,VGA(3),Int,3}, NoDims, Unitful.FreeUnits{(), NoDims, nothing}}(y)
    @test string(unitless) == "1e₂ (no units)"
    @test repr("text/plain", unitless) == string(unitless)
end
