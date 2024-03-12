@testset "Conversion" begin
    import CliffordNumbers: QFComplex
    k1 = KVector{1,APS}(1, 2, 3)
    k2 = KVector{2,APS}(1, 2, 3)
    # Conversion of scalar CliffordNumbers to Real subtypes
    @test convert(Int, CliffordNumber{APS,Float64}(1)) === Int(1)
    @test_throws InexactError convert(Int, CliffordNumber{APS}(1.5))
    # K-vectors of grade zero are scalars
    @test convert(Float64, KVector{0,APS}(3)) === Float64(3)
    @test convert(Complex{Float64}, CliffordNumber{QFComplex}(1, 2)) === 1.0 + 2.0im
    @test convert(CliffordNumber{QFComplex,Float64}, 1+2im) === CliffordNumber{QFComplex}(1.0, 2.0)
    @test convert(OddCliffordNumber, k1) === OddCliffordNumber{APS}(1, 2, 3, 0)
    @test convert(OddCliffordNumber{APS}, k1) === OddCliffordNumber{APS}(1, 2, 3, 0)
    @test convert(OddCliffordNumber{APS,Float64}, k1) === OddCliffordNumber{APS,Float64}(1, 2, 3, 0)
    @test convert(EvenCliffordNumber, k2) === EvenCliffordNumber{APS}(0, 1, 2, 3)
    @test_throws InexactError convert(EvenCliffordNumber, k1)
    @test_throws InexactError convert(OddCliffordNumber, k2)
    # Scalar conversion
    @test scalar_convert(Int, k1) === k1
    @test scalar_convert(Float32, k1) === KVector{1,APS,Float32}(1, 2, 3)
    @test scalar_convert(Int, 2) === 2
    @test scalar_convert(Float32, 2) === Float32(2)
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
    # Scalar promotion
    k = KVector{1,APS,Int16}(4, 2, 0)
    kk = KVector{1,APS,Float32}(4, 2, 0)
    l = KVector{2,APS}(0, 6, 9)
    ll = KVector{2,APS,Float32}(0, 6, 9)
    e = EvenCliffordNumber{APS}(1, 3, 3, 7)
    ee = EvenCliffordNumber{APS,Float32}(1, 3, 3, 7)
    @test scalar_promote() === tuple()
    @test scalar_promote(k) === tuple(k)
    @test scalar_promote(ee) === tuple(ee)
    @test scalar_promote(Float16(0)) === tuple(Float16(0))
    @test scalar_promote(k, l) === (KVector{1,APS}(4, 2, 0), l)
    @test scalar_promote(k, 2) === (KVector{1,APS}(4, 2, 0), 2)
    @test scalar_promote(k, l, ee) === (kk, ll, ee)
    @test scalar_promote(k, ll, e) === (kk, ll, ee)
    @test scalar_promote(kk, l, e) === (kk, ll, ee)
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
