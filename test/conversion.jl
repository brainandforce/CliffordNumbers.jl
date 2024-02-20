@testset "Conversion" begin
    # Conversion of scalar CliffordNumbers to Real subtypes
    @test convert(Int, CliffordNumber{APS,Float64}(1)) === Int(1)
    @test_throws InexactError convert(Int, CliffordNumber{APS}(1.5))
end

@testset "Promotion" begin
    @test promote_type(KVector{1,APS,Int}, KVector{1,APS,Float64}) === KVector{1,APS,Float64,3}
    @test promote_type(KVector{1,APS,Int}, KVector{2,APS,Int}) === CliffordNumber{APS,Int,8}
    @test promote_type(KVector{0,APS,Int}, Int) === KVector{0,APS,Int,1}
    @test promote_type(KVector{0,APS,Int}, Float64) === KVector{0,APS,Float64,1}
    @test promote_type(KVector{1,APS,Int}, CliffordNumber{APS,Int}) === CliffordNumber{APS,Int,8}
    @test promote_type(KVector{1,APS,Int}, CliffordNumber{APS,Float64}) ===
        CliffordNumber{APS,Float64,8}
    @test promote_type(KVector{1,APS,Float64}, CliffordNumber{APS,Int}) ===
        CliffordNumber{APS,Float64,8}
end
