@testset "Quaternions.jl extension" begin
    q = Quaternion(0, 1, 2, 3)
    k = KVector{1,VGA(3)}(4, 2, 0)
    l = KVector{2,VGA(3)}(0, 6, 9)
    # Conversion to AbstractCliffordNumber
    eq = EvenCliffordNumber{VGA(3)}(0, 1, 2, 3)
    @test AbstractCliffordNumber(q) === eq
    @test EvenCliffordNumber{VGA(3)}(q) === eq
    @test EvenCliffordNumber{VGA(3),Float64}(q) === scalar_convert(Float64, eq)
    @test KVector{1,VGA(3)}(q) === zero(KVector{1,VGA(3),Int})
    @test OddCliffordNumber{VGA(3)}(q) === zero(OddCliffordNumber{VGA(3),Int})
    @test convert(AbstractCliffordNumber, q) === eq
    @test convert(EvenCliffordNumber{VGA(3)}, q) === eq
    @test convert(EvenCliffordNumber{VGA(3),Float64}, q) === scalar_convert(Float64, eq)
    @test_throws InexactError convert(KVector{1,VGA(3)}, q)
    @test_throws InexactError convert(OddCliffordNumber{VGA(3)}, q)
    # Conversion to Quaternion
    @test Quaternion(l) === Quaternion(0, 0, 6, 9)
    @test Quaternion(float(l)) === Quaternion(0.0, 0.0, 6.0, 9.0)
    @test Quaternion{Float64}(l) === Quaternion(0.0, 0.0, 6.0, 9.0)
    @test QuaternionF32(l) === QuaternionF32(0, 0, 6, 9)
    @test convert(Quaternion, l) === Quaternion(0, 0, 6, 9)
    @test convert(Quaternion{Float64}, l) === Quaternion{Float64}(0, 0, 6, 9)
    @test Quaternion(k) === zero(Quaternion{scalar_type(k)})
    @test Quaternion(float(k)) === zero(Quaternion{Float64})
    @test Quaternion{Float64}(k) === zero(Quaternion{Float64})
    @test QuaternionF32(k) === zero(QuaternionF32)
    @test_throws InexactError convert(Quaternion, k)
    @test_throws InexactError convert(Quaternion{Float64}, k)
    @test q + k === CliffordNumber{VGA(3)}(0, 4, 2, 1, 0, 2, 3, 0)
    @test q + l === EvenCliffordNumber{VGA(3)}(0, 1, 8, 12)
    @test float(q) + k === CliffordNumber{VGA(3),Float64}(0, 4, 2, 1, 0, 2, 3, 0)
    @test float(q) + l === EvenCliffordNumber{VGA(3),Float64}(0, 1, 8, 12)
    @test q + float(k) === CliffordNumber{VGA(3),Float64}(0, 4, 2, 1, 0, 2, 3, 0)
    @test q + float(l) === EvenCliffordNumber{VGA(3),Float64}(0, 1, 8, 12)
    # Products
    @test q * k === convert(EvenCliffordNumber{VGA(3)}, q) * k
    @test q * l === convert(EvenCliffordNumber{VGA(3)}, q) * l
    @test k * q === k * convert(EvenCliffordNumber{VGA(3)}, q)
    @test l * q === l * convert(EvenCliffordNumber{VGA(3)}, q)
    # Spherical linear interpolation
    a = exp(1.09 * KVector{2,VGA(3)}(1, 2, 3))
    b = exp(1.11 * KVector{2,VGA(3)}(1, 2, 3))
    qa = Quaternion(a)
    qb = Quaternion(b)
    @test Quaternion(slerp(a, b, 0.5)) === slerp(qa, qb, 0.5)
    @test Quaternion(slerp(qa, b, 0.5)) === slerp(qa, qb, 0.5)
    @test Quaternion(slerp(a, qb, 0.5)) === slerp(qa, qb, 0.5)
end
