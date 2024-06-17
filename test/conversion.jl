@testset "Conversion" begin
    k1 = KVector{1,VGA(3)}(1, 2, 3)
    k2 = KVector{2,VGA(3)}(1, 2, 3)
    # Conversion of scalar CliffordNumbers to Real subtypes
    @test convert(Int, CliffordNumber{VGA(3),Float64}(1)) === Int(1)
    @test_throws InexactError convert(Int, CliffordNumber{VGA(3)}(1.5))
    # K-vectors of grade zero are scalars
    @test convert(Float64, KVector{0,VGA(3)}(3)) === Float64(3)
    #=
    @test convert(Complex{Float64}, CliffordNumber{QFComplex}(1, 2)) === 1.0 + 2.0im
    @test convert(CliffordNumber{QFComplex,Float64}, 1+2im) === CliffordNumber{QFComplex}(1.0, 2.0)
    =#
    @test convert(OddCliffordNumber, k1) === OddCliffordNumber{VGA(3)}(1, 2, 3, 0)
    @test convert(OddCliffordNumber{VGA(3)}, k1) === OddCliffordNumber{VGA(3)}(1, 2, 3, 0)
    @test convert(OddCliffordNumber{VGA(3),Float64}, k1) === 
        OddCliffordNumber{VGA(3),Float64}(1, 2, 3, 0)
    @test convert(EvenCliffordNumber, k2) === EvenCliffordNumber{VGA(3)}(0, 1, 2, 3)
    @test_throws InexactError convert(EvenCliffordNumber, k1)
    @test_throws InexactError convert(OddCliffordNumber, k2)
    # Scalar conversion
    @test scalar_convert(Int, k1) === k1
    @test scalar_convert(Float32, k1) === KVector{1,VGA(3),Float32}(1, 2, 3)
    @test scalar_convert(Int, 2) === 2
    @test scalar_convert(Float32, 2) === Float32(2)
end

@testset "Promotion" begin
    @test promote_type(Int, CliffordNumber{VGA(3)}) === CliffordNumber{VGA(3)}
    @test promote_type(Int, CliffordNumber{VGA(3),Float64}) === CliffordNumber{VGA(3),Float64,8}
    @test promote_type(Float64, CliffordNumber{VGA(3),Int}) === CliffordNumber{VGA(3),Float64,8}
    @test promote_type(CliffordNumber{VGA(3),Int}, CliffordNumber{VGA(3),Float64}) ===
        CliffordNumber{VGA(3),Float64,8}
    @test promote_type(KVector{1,VGA(3),Int}, KVector{1,VGA(3),Float64}) === KVector{1,VGA(3),Float64,3}
    @test promote_type(KVector{0,VGA(3),Int}, KVector{2,VGA(3),Int}) === EvenCliffordNumber{VGA(3),Int,4}
    @test promote_type(KVector{1,VGA(3),Int}, KVector{3,VGA(3),Int}) === OddCliffordNumber{VGA(3),Int,4}
    @test promote_type(KVector{1,VGA(3),Int}, KVector{2,VGA(3),Int}) === CliffordNumber{VGA(3),Int,8}
    @test promote_type(KVector{1,VGA(3),Int}, EvenCliffordNumber{VGA(3),Int}) === 
        CliffordNumber{VGA(3),Int,8}
    @test promote_type(KVector{2,VGA(3),Int}, EvenCliffordNumber{VGA(3),Int}) === 
        EvenCliffordNumber{VGA(3),Int,4}
    @test promote_type(KVector{2,VGA(3),Int}, OddCliffordNumber{VGA(3),Int}) === 
        CliffordNumber{VGA(3),Int,8}
    @test promote_type(KVector{1,VGA(3),Int}, OddCliffordNumber{VGA(3),Int}) === 
        OddCliffordNumber{VGA(3),Int,4}
    @test promote_type(KVector{0,VGA(3),Int}, Int) === KVector{0,VGA(3),Int,1}
    @test promote_type(KVector{0,VGA(3),Int}, Float64) === KVector{0,VGA(3),Float64,1}
    @test promote_type(KVector{1,VGA(3),Int}, Int) === CliffordNumber{VGA(3),Int,8}
    @test promote_type(KVector{1,VGA(3),Int}, Float64) === CliffordNumber{VGA(3),Float64,8}
    @test promote_type(KVector{1,VGA(3),Int}, CliffordNumber{VGA(3),Int}) === CliffordNumber{VGA(3),Int,8}
    @test promote_type(KVector{1,VGA(3),Int}, CliffordNumber{VGA(3),Float64}) ===
        CliffordNumber{VGA(3),Float64,8}
    @test promote_type(KVector{1,VGA(3),Float64}, CliffordNumber{VGA(3),Int}) ===
        CliffordNumber{VGA(3),Float64,8}
    # Scalar promotion
    k = KVector{1,VGA(3),Int16}(4, 2, 0)
    kk = KVector{1,VGA(3),Float32}(4, 2, 0)
    l = KVector{2,VGA(3)}(0, 6, 9)
    ll = KVector{2,VGA(3),Float32}(0, 6, 9)
    e = EvenCliffordNumber{VGA(3)}(1, 3, 3, 7)
    ee = EvenCliffordNumber{VGA(3),Float32}(1, 3, 3, 7)
    @test scalar_promote() === tuple()
    @test scalar_promote(k) === tuple(k)
    @test scalar_promote(ee) === tuple(ee)
    @test scalar_promote(Float16(0)) === tuple(Float16(0))
    @test scalar_promote(k, l) === (KVector{1,VGA(3)}(4, 2, 0), l)
    @test scalar_promote(k, 2) === (KVector{1,VGA(3)}(4, 2, 0), 2)
    @test scalar_promote(k, l, ee) === (kk, ll, ee)
    @test scalar_promote(k, ll, e) === (kk, ll, ee)
    @test scalar_promote(kk, l, e) === (kk, ll, ee)
end

@testset "Widening" begin
    @test widen(KVector{2,VGA(3),Int32,3}) === KVector{2,VGA(3),Int,3}
    @test widen(KVector{2,VGA(3),Float16}(1, 2, 3)) === KVector{2,VGA(3),Float32}(1, 2, 3)
    @test widen_grade(KVector) === CliffordNumber
    @test widen_grade(KVector{2}) === EvenCliffordNumber
    @test widen_grade(KVector{2,VGA(3)}) === EvenCliffordNumber{VGA(3)}
    @test widen_grade(KVector{2,VGA(3),Int32}) === EvenCliffordNumber{VGA(3),Int32}
    @test widen_grade(KVector{2,VGA(3),Int32,3}) === EvenCliffordNumber{VGA(3),Int32,4}
    @test widen_grade(KVector{1}) === OddCliffordNumber
    @test widen_grade(KVector{1,VGA(3)}) === OddCliffordNumber{VGA(3)}
    @test widen_grade(KVector{1,VGA(3),Int32}) === OddCliffordNumber{VGA(3),Int32}
    @test widen_grade(KVector{1,VGA(3),Int32,3}) === OddCliffordNumber{VGA(3),Int32,4}
    @test widen_grade(CliffordNumbers.Z2CliffordNumber) === CliffordNumber
    @test widen_grade(EvenCliffordNumber) === CliffordNumber
    @test widen_grade(OddCliffordNumber) === CliffordNumber
    @test widen_grade(EvenCliffordNumber{VGA(3)}) === CliffordNumber{VGA(3)}
    @test widen_grade(OddCliffordNumber{VGA(3)}) === CliffordNumber{VGA(3)}
    @test widen_grade(EvenCliffordNumber{VGA(3),Int32}) === CliffordNumber{VGA(3),Int32}
    @test widen_grade(OddCliffordNumber{VGA(3),Int32}) === CliffordNumber{VGA(3),Int32}
    @test widen_grade(EvenCliffordNumber{VGA(3),Int32,4}) === CliffordNumber{VGA(3),Int32,8}
    @test widen_grade(OddCliffordNumber{VGA(3),Int32,4}) === CliffordNumber{VGA(3),Int32,8}
    @test widen_grade(KVector{2,VGA(3),Int32}(1, 2, 3)) === EvenCliffordNumber{VGA(3),Int32}(0, 1, 2, 3)
    @test widen_grade(KVector{1,VGA(3),Int}(4, 5, 6)) === OddCliffordNumber{VGA(3),Int}(4, 5, 6, 0)
end

@testset "Complex algebras" begin
    @test complex(CliffordNumber{VGA(1)}(1, 1)) === CliffordNumber{VGA(1),Complex{Int}}(1,1)
end

@testset "float() and big()" begin
    @test float(CliffordNumber{STA,Int}) <: CliffordNumber{STA,Float64}
    @test float(CliffordNumber{STA,Complex{Bool}}) <: CliffordNumber{STA,ComplexF64}
    @test big(CliffordNumber{STA,Int}) <: CliffordNumber{STA,BigInt}
    @test big(CliffordNumber{STA,Float64}) <: CliffordNumber{STA,BigFloat}
    @test float(KVector{1,VGA(3)}(1, 2, 3)) === KVector{1,VGA(3),Float64}(1, 2, 3)
    @test float(KVector{1,VGA(3),Complex{Int}}(1, 2, 3)) === KVector{1,VGA(3),ComplexF64}(1, 2, 3)
    @test big(KVector{1,VGA(3)}(1, 2, 3)) == KVector{1,VGA(3),BigInt}(1, 2, 3)
    @test big(KVector{1,VGA(3),Float64}(1, 2, 3)) == KVector{1,VGA(3),BigFloat}(1, 2, 3)
end
