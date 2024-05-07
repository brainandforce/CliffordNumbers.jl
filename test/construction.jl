@testset "Zero elements" begin
    @test zero(KVector{0,VGA(3)}) == KVector{0,VGA(3),Bool}(0)
    @test zero(KVector{1,VGA(3)}) == KVector{1,VGA(3),Bool}(0, 0, 0)
    @test zero(KVector{2,VGA(3)}) == KVector{2,VGA(3),Bool}(0, 0, 0)
    @test zero(KVector{3,VGA(3)}) == KVector{3,VGA(3),Bool}(0)
    @test zero(EvenCliffordNumber{VGA(3)}) == EvenCliffordNumber{VGA(3),Bool}(0, 0, 0, 0)
    @test zero(OddCliffordNumber{VGA(3)}) == OddCliffordNumber{VGA(3),Bool}(0, 0, 0, 0)
    @test zero(CliffordNumber{VGA(3)}) == CliffordNumber{VGA(3),Bool}(0, 0, 0, 0, 0, 0, 0, 0)
    @test zero(EvenCliffordNumber{VGA(3),Int}) === EvenCliffordNumber{VGA(3),Int}(0, 0, 0, 0)
    @test zero(OddCliffordNumber{VGA(3),Float32}) === OddCliffordNumber{VGA(3),Float32}(0, 0, 0, 0)
    @test zero(CliffordNumber{VGA(3),BigInt}) == CliffordNumber{VGA(3),BigInt}(0, 0, 0, 0, 0, 0, 0, 0)
end

@testset "Units" begin
    # one
    @test one(KVector{0,VGA(3)}) === KVector{0,VGA(3),Bool}(1)
    @test one(KVector{0,VGA(3),Float32}) === KVector{0,VGA(3),Float32}(1)
    @test one(KVector{0,VGA(3)}) == 1
    @test one(KVector{1,VGA(3)}) == KVector{0,VGA(3),Bool}(1)
    @test one(KVector{1,VGA(3),Float32}) === KVector{0,VGA(3),Float32}(1)
    @test one(KVector{1,VGA(3)}) == 1
    @test one(KVector{2,VGA(3)}) == KVector{0,VGA(3),Bool}(1)
    @test one(KVector{2,VGA(3),Float32}) === KVector{0,VGA(3),Float32}(1)
    @test one(KVector{2,VGA(3)}) == 1
    @test one(KVector{3,VGA(3)}) == KVector{0,VGA(3),Bool}(1)
    @test one(KVector{3,VGA(3),Float32}) === KVector{0,VGA(3),Float32}(1)
    @test one(KVector{3,VGA(3)}) == 1
    @test one(EvenCliffordNumber{VGA(3)}) == EvenCliffordNumber{VGA(3),Bool}(1, 0, 0, 0)
    @test one(EvenCliffordNumber{VGA(3)}) == 1
    @test one(OddCliffordNumber{VGA(3)}) == KVector{0,VGA(3),Bool}(1)
    @test one(OddCliffordNumber{VGA(3)}) == 1
    @test one(CliffordNumber{VGA(3)}) == CliffordNumber{VGA(3),Bool}(1, 0, 0, 0, 0, 0, 0, 0)
    @test one(CliffordNumber{VGA(3)}) == 1
    # oneunit
    @test oneunit(KVector{0,VGA(3)}) == KVector{0,VGA(3),Bool}(1)
    @test oneunit(KVector{0,VGA(3)}) == 1
    @test oneunit(EvenCliffordNumber{VGA(3)}) == EvenCliffordNumber{VGA(3),Bool}(1, 0, 0, 0)
    @test oneunit(EvenCliffordNumber{VGA(3)}) == 1
    @test oneunit(CliffordNumber{VGA(3)}) == CliffordNumber{VGA(3),Bool}(1, 0, 0, 0, 0, 0, 0, 0)
    @test oneunit(CliffordNumber{VGA(3)}) == 1
    # these should throw since they can't represent 1
    @test_throws InexactError oneunit(KVector{1,VGA(3)})
    @test_throws InexactError oneunit(KVector{2,VGA(3)})
    @test_throws InexactError oneunit(KVector{3,VGA(3)})
    @test_throws InexactError oneunit(OddCliffordNumber{VGA(3)})
end

@testset "Real and complex components" begin
    x = EvenCliffordNumber{VGA(3)}(0, 4, 2, 0)
    y = OddCliffordNumber{VGA(3)}(0, 6, 9, 0)
    @test real(x) === x
    @test complex(x) === EvenCliffordNumber{VGA(3),Complex{Int}}(0, 4, 2, 0)
    @test real(complex(y)) === y
    @test complex(y, y) === OddCliffordNumber{VGA(3)}(0, 6+6im, 9+9im, 0)
    @test complex(x, y) === CliffordNumber{VGA(3)}(0, 0, 6im, 4, 9im, 2, 0, 0)
end

@testset "Constructors" begin
    k1 = KVector{1,VGA(3)}(1, 2, 3)
    k2 = KVector{2,VGA(3)}(1, 2, 3)
    @test CliffordNumber{VGA(3)}(1337) === CliffordNumber{VGA(3)}(1337, 0, 0, 0, 0, 0, 0, 0)
    @test EvenCliffordNumber{VGA(3)}(1337) === EvenCliffordNumber{VGA(3)}(1337, 0, 0, 0)
    # Mixed types in input
    @test CliffordNumber{VGA(3)}(0.0, 0, 0, 0, 0, 0, 0, 0) === zero(CliffordNumber{VGA(3),Float64})
    @test CliffordNumber{VGA(3)}(0, 0, 0, 0, 0, 0, 0, 0.0) === zero(CliffordNumber{VGA(3),Float64})
    @test EvenCliffordNumber{VGA(3)}(0.0, 0, 0, 0) === zero(EvenCliffordNumber{VGA(3),Float64})
    @test EvenCliffordNumber{VGA(3)}(0, 0, 0, 0.0) === zero(EvenCliffordNumber{VGA(3),Float64})
    @test KVector{1,VGA(3)}(0, 0.0, 0) === zero(KVector{1,VGA(3),Float64})
    # Barest constructors
    @test CliffordNumber(k2) === CliffordNumber{VGA(3)}(0, 0, 0, 1, 0, 2, 3, 0)
    @test EvenCliffordNumber(k2) === EvenCliffordNumber{VGA(3)}(0, 1, 2, 3)
    @test OddCliffordNumber(k2) === OddCliffordNumber{VGA(3)}(0, 0, 0, 0)
    @test EvenCliffordNumber(k1) === EvenCliffordNumber{VGA(3)}(0, 0, 0, 0)
    @test OddCliffordNumber(k1) === OddCliffordNumber{VGA(3)}(1, 2, 3, 0)
    @test CliffordNumbers.Z2CliffordNumber(k1) === OddCliffordNumber{VGA(3)}(1, 2, 3, 0)
    @test CliffordNumbers.Z2CliffordNumber(k2) === EvenCliffordNumber{VGA(3)}(0, 1, 2, 3)
end

@testset "Similar types" begin
    import CliffordNumbers.similar_type
    @test similar(CliffordNumber{VGA(3)}, Int, Val(VGA(2))) isa CliffordNumber{VGA(2),Int}
    @test similar(EvenCliffordNumber{VGA(3)}, Int) isa EvenCliffordNumber{VGA(3),Int}
    @test similar(KVector, Int, Val(STA), Val(2)) isa KVector{2,STA,Int}
end

@testset "Scalar types" begin
    @test scalar_type(Int) === Int
    @test scalar_type(ComplexF16) === ComplexF16
    @test scalar_type(Int(420)) === Int
    @test scalar_type(ComplexF16(69)) === ComplexF16
    @test scalar_type(CliffordNumber{VGA(3),Float64,8}) === Float64
    @test scalar_type(zero(OddCliffordNumber{VGA(3),Int})) === Int
end

@testset "Metric signatures" begin
    @test signature(CliffordNumber{STAP}) === STAP
    @test signature(CliffordNumber{STAP,Float32}) === STAP
    @test signature(CliffordNumber{STAP,Float32,16}) === STAP
    @test signature(EvenCliffordNumber{STAP,Float32,8}) === STAP
    @test signature(OddCliffordNumber{STAP,Float32,8}) === STAP
    @test signature(KVector{2,STAP,Float32,6}) === STAP
    @test signature(zero(CliffordNumber{STAP})) === STAP
    @test signature(zero(EvenCliffordNumber{STAP})) === STAP
    @test signature(typeof(zero(CliffordNumber{STAP}))) === STAP
end
