@testset "Equality" begin
    x = CliffordNumber{APS,Float64}(1, 2, 3, 4, 5, 6, 7, 8)
    @test convert(CliffordNumber{APS,Int}, x) isa CliffordNumber{APS,Int}
    @test x == convert(CliffordNumber{APS,Int}, x)
    @test x !== convert(CliffordNumber{APS,Int}, x)
    # Equality between disparate types
    @test KVector{2,APS}(4, 2, 0) == CliffordNumber{APS,Int}(0, 0, 0, 4, 0, 2, 0, 0)
    @test KVector{2,APS}(4, 2, 0) == CliffordNumber{APS,Float64}(0, 0, 0, 4, 0, 2, 0, 0)
end

@testset "Grade automorphisms" begin
    k0 = KVector{0,APS}(1)
    k1 = KVector{1,APS}(1, 2, 3)
    k2 = KVector{2,APS,Float64}(1, 2, 3)
    k3 = KVector{3,APS,Float64}(1)
    @test reverse(k0) === KVector{0,APS}(1)
    @test reverse(k1) === KVector{1,APS}(1, 2, 3)
    @test reverse(k2) === KVector{2,APS,Float64}(-1, -2, -3)
    @test reverse(k3) === KVector{3,APS,Float64}(-1)
    @test grade_involution(k0) === KVector{0,APS}(1)
    @test grade_involution(k1) === KVector{1,APS}(-1, -2, -3)
    @test grade_involution(k2) === KVector{2,APS,Float64}(1, 2, 3)
    @test grade_involution(k3) === KVector{3,APS,Float64}(-1)
    @test conj(k0) === KVector{0,APS}(1)
    @test conj(k1) === KVector{1,APS}(-1, -2, -3)
    @test conj(k2) === KVector{2,APS,Float64}(-1, -2, -3)
    @test conj(k3) === KVector{3,APS,Float64}(1)
end

@testset "Addition and subtraction" begin
    x = CliffordNumber{APS}(1, 2, 3, 4, 5, 6, 7, 8)
    y = CliffordNumber{APS,Float64}(9, -10, 11, -12, 13, 14, -15, 16)
    # Type promotion test
    @test x + y isa CliffordNumber{APS,Float64}
    # Equality test (mixed types)
    @test x + y == CliffordNumber{APS}(10, -8, 14, -8, 18, 20, -8, 24)
    @test -x === CliffordNumber{APS}(-1, -2, -3, -4, -5, -6, -7, -8)
    @test x - y === CliffordNumber{APS}(-8.0, 12.0, -8.0, 16.0, -8.0, -8.0, 22.0, -8.0)
end

@testset "muladd" begin
    e = EvenCliffordNumber{APS,Float32}(0, 4, 2, 0)
    ee = EvenCliffordNumber{APS,Float64}(0, 4, 2, 0)
    f = EvenCliffordNumber{APS,Float32}(0, 0, 6, 9)
    ff = EvenCliffordNumber{APS,Float64}(0, 0, 6, 9)
    @test muladd(2, ee, ff) === EvenCliffordNumber{APS,Float64}(0, 8, 10, 9)
    @test muladd(2, ff, ee) === EvenCliffordNumber{APS,Float64}(0, 4, 14, 18)
    @test muladd(2, e, f) === muladd(e, 2, f)
    @test muladd(2, e, ff) == muladd(2, ee, f)
end

@testset "Geometric product" begin
    x = CliffordNumber{APS}(0, 2, 0, 0, 0, 0, 0, 0)
    y = CliffordNumber{APS}(0, 3, 4, 0, 0, 0, 0, 0)
    k1 = KVector{1,APS}(4, 2, 0)
    k2 = KVector{2,APS}(4, 2, 0)
    l1 = KVector{1,APS}(0, 6, 9)
    l2 = KVector{2,APS}(0, 6, 9)
    # Signs included explicilty for clarity
    @test k1 * l1 === EvenCliffordNumber{APS}(+12, +24, +36, +18)
    @test l1 * k1 === EvenCliffordNumber{APS}(+12, -24, -36, -18)
    @test k1 * l1 === reverse(l1 * k1)
    @test k2 * l2 === EvenCliffordNumber{APS}(-12, -18, +36, -24)
    @test l2 * k2 === EvenCliffordNumber{APS}(-12, +18, -36, +24)
    @test k2 * l2 === reverse(l2 * k2)
    @test x * y == CliffordNumber{APS}(6, 0, 0, 8, 0, 0, 0, 0)
    # When multiplying vectors, the reverse should have a negative bivector
    @test y * x == CliffordNumber{APS}(6, 0, 0, -8, 0, 0, 0, 0)
end

@testset "Scalars" begin
    k = KVector{2,APS}(3, 4, 0)
    @test abs2(k) === 25
    @test abs(k) == 5
    @test normalize(k) == k / 5
end

@testset "Wedge product" begin
    x = KVector{1,APS}(1, 0, 0)
    y = KVector{1,APS}(0, 1, 0)
    z = KVector{1,APS}(0, 0, 1)
    k1 = KVector{1,APS}(3, 4, 0)
    k2 = KVector{2,APS}(4, 2, 0)
    five = 5 * one(CliffordNumber{APS})
    # Self wedges should be zero
    @test iszero(x ∧ x)
    @test iszero(y ∧ y)
    @test iszero(z ∧ z)
    # Reversing the order should flip the sign for vectors
    @test x ∧ y == -(y ∧ x)
    @test x ∧ z == -(z ∧ x)
    @test y ∧ z == -(z ∧ y)
    @test x ∧ y ∧ z == -(z ∧ y ∧ x)
    @test iszero(x ∧ y ∧ x)
    @test iszero(x ∧ y ∧ z ∧ x)
    @test KVector{1,APS}(x) ∧ KVector{1,APS}(y) == -(KVector{1,APS}(y) ∧ KVector{1,APS}(x))
    # Ensure the behavior of scalars (both CliffordNumber and normal) are correct
    @test 5 ∧ x == 5 * x
    @test five ∧ x == x ∧ five
    @test 5 ∧ 5 == 25
    # Turns out reversing the order shouldn't change anything
    @test k1 ∧ k2 == KVector{3,APS}(-8)
end

@testset "Exponentiation" begin
    e12 = CliffordNumber{APS}(0, 0, 0, 1, 0, 0, 0, 0)
    @test exppi(e12/2) == e12
    @test exptau(e12) == 1
    @test exppi(e12) ≈ exp(pi*e12)
    @test exptau(e12) ≈ exp(2*pi*e12)
    @test exp(KVector{1,APS}(1,0,0)) isa CliffordNumber
    @test exp(KVector{2,APS}(1,0,0)) isa EvenCliffordNumber
    k = KVector{2,APS}(1,0,0)
    @test CliffordNumbers.exp_taylor(pi/2 * k) ≈ exp(pi/2 * k)
    @test CliffordNumbers.exp_taylor(pi/2 * k) ≈ exppi(1//2 * k)
    @test CliffordNumbers.exp_taylor(pi/2 * k) ≈ exptau(1//4 * k)
    # Integer exponentiation of a multivector
    k1 = KVector{1,APS}(4, 2, 0)
    k2 = KVector{2,APS}(0, 6, 9)
    @test k1^2 isa CliffordNumber
    @test k2^2 isa EvenCliffordNumber
    @test k1^0 == one(k1)
    @test k1^1 == k1
    @test k1^2 == k1 * k1
    @test k2^0 == one(k2)
    @test k2^1 == k2
    @test k2^2 == k2 * k2
end
