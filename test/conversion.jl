@testset "Conversion" begin
    import CliffordNumbers: QFComplex
    # Conversion of scalar CliffordNumbers to Real subtypes
    @test convert(Int, CliffordNumber{APS,Float64}(1)) === Int(1)
    @test_throws InexactError convert(Int, CliffordNumber{APS}(1.5))
    # K-vectors of grade zero are scalars
    @test convert(Float64, KVector{0,APS}(3)) === Float64(3)
    @test convert(Complex{Float64}, CliffordNumber{QFComplex}(1, 2)) === 1.0 + 2.0im
    @test convert(CliffordNumber{QFComplex,Float64}, 1+2im) === CliffordNumber{QFComplex}(1.0, 2.0)
end

@testset "Promotion" begin
    @test promote_type(Int, CliffordNumber{APS}) === CliffordNumber{APS}
    @test promote_type(Int, CliffordNumber{APS,Float64}) === CliffordNumber{APS,Float64,8}
    @test promote_type(Float64, CliffordNumber{APS,Int}) === CliffordNumber{APS,Float64,8}
    @test promote_type(CliffordNumber{APS,Int}, CliffordNumber{APS,Float64}) ===
        CliffordNumber{APS,Float64,8}
    @test promote_type(KVector{1,APS,Int}, KVector{1,APS,Float64}) === KVector{1,APS,Float64,3}
    @test promote_type(KVector{0,APS,Int}, KVector{2,APS,Int}) === EvenCliffordNumber{APS,Int,4}
    @test promote_type(KVector{1,APS,Int}, KVector{3,APS,Int}) === OddCliffordNumber{APS,Int,4}
    @test promote_type(KVector{1,APS,Int}, KVector{2,APS,Int}) === CliffordNumber{APS,Int,8}
    @test promote_type(KVector{1,APS,Int}, EvenCliffordNumber{APS,Int}) === 
        CliffordNumber{APS,Int,8}
    @test promote_type(KVector{2,APS,Int}, EvenCliffordNumber{APS,Int}) === 
        EvenCliffordNumber{APS,Int,4}
    @test promote_type(KVector{2,APS,Int}, OddCliffordNumber{APS,Int}) === 
        CliffordNumber{APS,Int,8}
    @test promote_type(KVector{1,APS,Int}, OddCliffordNumber{APS,Int}) === 
        OddCliffordNumber{APS,Int,4}
    @test promote_type(KVector{0,APS,Int}, Int) === KVector{0,APS,Int,1}
    @test promote_type(KVector{0,APS,Int}, Float64) === KVector{0,APS,Float64,1}
    @test promote_type(KVector{1,APS,Int}, Int) === CliffordNumber{APS,Int,8}
    @test promote_type(KVector{1,APS,Int}, Float64) === CliffordNumber{APS,Float64,8}
    @test promote_type(KVector{1,APS,Int}, CliffordNumber{APS,Int}) === CliffordNumber{APS,Int,8}
    @test promote_type(KVector{1,APS,Int}, CliffordNumber{APS,Float64}) ===
        CliffordNumber{APS,Float64,8}
    @test promote_type(KVector{1,APS,Float64}, CliffordNumber{APS,Int}) ===
        CliffordNumber{APS,Float64,8}
end
