@testset "Conversion" begin
    # Conversion of scalar CliffordNumbers to Real subtypes
    @test convert(Int, CliffordNumber{APS,Float64}(1)) === Int(1)
    @test_throws InexactError convert(Int, CliffordNumber{APS}(1.5))
end
