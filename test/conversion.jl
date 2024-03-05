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

@testset "Widening" begin
    @test widen(KVector{2,APS,Int32,3}) === KVector{2,APS,Int,3}
    @test widen(KVector{2,APS,Float16}(1, 2, 3)) === KVector{2,APS,Float32}(1, 2, 3)
    @test widen_grade(KVector) === CliffordNumber
    @test widen_grade(KVector{2}) === EvenCliffordNumber
    @test widen_grade(KVector{2,APS}) === EvenCliffordNumber{APS}
    @test widen_grade(KVector{2,APS,Int32}) === EvenCliffordNumber{APS,Int32}
    @test widen_grade(KVector{2,APS,Int32,3}) === EvenCliffordNumber{APS,Int32,4}
    @test widen_grade(KVector{1}) === OddCliffordNumber
    @test widen_grade(KVector{1,APS}) === OddCliffordNumber{APS}
    @test widen_grade(KVector{1,APS,Int32}) === OddCliffordNumber{APS,Int32}
    @test widen_grade(KVector{1,APS,Int32,3}) === OddCliffordNumber{APS,Int32,4}
    @test widen_grade(CliffordNumbers.Z2CliffordNumber) === CliffordNumber
    @test widen_grade(EvenCliffordNumber) === CliffordNumber
    @test widen_grade(OddCliffordNumber) === CliffordNumber
    @test widen_grade(EvenCliffordNumber{APS}) === CliffordNumber{APS}
    @test widen_grade(OddCliffordNumber{APS}) === CliffordNumber{APS}
    @test widen_grade(EvenCliffordNumber{APS,Int32}) === CliffordNumber{APS,Int32}
    @test widen_grade(OddCliffordNumber{APS,Int32}) === CliffordNumber{APS,Int32}
    @test widen_grade(EvenCliffordNumber{APS,Int32,4}) === CliffordNumber{APS,Int32,8}
    @test widen_grade(OddCliffordNumber{APS,Int32,4}) === CliffordNumber{APS,Int32,8}
    @test widen_grade(KVector{2,APS,Int32}(1, 2, 3)) === EvenCliffordNumber{APS,Int32}(0, 1, 2, 3)
    @test widen_grade(KVector{1,APS,Int}(4, 5, 6)) === OddCliffordNumber{APS,Int}(4, 5, 6, 0)
end
