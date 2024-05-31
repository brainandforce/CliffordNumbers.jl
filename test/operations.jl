@testset "(Approximate) equality" begin
    x = CliffordNumber{VGA(3),Float64}(1, 2, 3, 4, 5, 6, 7, 8)
    @test convert(CliffordNumber{VGA(3),Int}, x) isa CliffordNumber{VGA(3),Int}
    @test x == convert(CliffordNumber{VGA(3),Int}, x)
    @test x !== convert(CliffordNumber{VGA(3),Int}, x)
    # Equality between disparate types
    @test KVector{2,VGA(3)}(4, 2, 0) == CliffordNumber{VGA(3),Int}(0, 0, 0, 4, 0, 2, 0, 0)
    @test KVector{2,VGA(3)}(4, 2, 0) == CliffordNumber{VGA(3),Float64}(0, 0, 0, 4, 0, 2, 0, 0)
    @test one(EvenCliffordNumber{VGA(3)}) == 1
    @test 1 == one(EvenCliffordNumber{VGA(3)})
    @test one(EvenCliffordNumber{VGA(3)}) ≈ 1
    @test 1 ≈ one(EvenCliffordNumber{VGA(3)})
end

@testset "Grade automorphisms" begin
    k0 = KVector{0,VGA(3)}(1)
    k1 = KVector{1,VGA(3)}(1, 2, 3)
    k2 = KVector{2,VGA(3),Float64}(1, 2, 3)
    k3 = KVector{3,VGA(3),Float64}(1)
    @test reverse(k0) === KVector{0,VGA(3)}(1)
    @test reverse(k1) === KVector{1,VGA(3)}(1, 2, 3)
    @test reverse(k2) === KVector{2,VGA(3),Float64}(-1, -2, -3)
    @test reverse(k3) === KVector{3,VGA(3),Float64}(-1)
    @test grade_involution(k0) === KVector{0,VGA(3)}(1)
    @test grade_involution(k1) === KVector{1,VGA(3)}(-1, -2, -3)
    @test grade_involution(k2) === KVector{2,VGA(3),Float64}(1, 2, 3)
    @test grade_involution(k3) === KVector{3,VGA(3),Float64}(-1)
    @test conj(k0) === KVector{0,VGA(3)}(1)
    @test conj(k1) === KVector{1,VGA(3)}(-1, -2, -3)
    @test conj(k2) === KVector{2,VGA(3),Float64}(-1, -2, -3)
    @test conj(k3) === KVector{3,VGA(3),Float64}(1)
end

@testset "Addition and subtraction" begin
    x = CliffordNumber{VGA(3)}(1, 2, 3, 4, 5, 6, 7, 8)
    y = CliffordNumber{VGA(3),Float64}(9, -10, 11, -12, 13, 14, -15, 16)
    # Type promotion test
    @test x + y isa CliffordNumber{VGA(3),Float64}
    # Equality test (mixed types)
    @test x + y == CliffordNumber{VGA(3)}(10, -8, 14, -8, 18, 20, -8, 24)
    @test -x === CliffordNumber{VGA(3)}(-1, -2, -3, -4, -5, -6, -7, -8)
    @test x - y === CliffordNumber{VGA(3)}(-8.0, 12.0, -8.0, 16.0, -8.0, -8.0, 22.0, -8.0)
    @test 1 + KVector{0,VGA(3)}(1) === KVector{0,VGA(3)}(2)
    @test 1 - KVector{1,VGA(3)}(0, 0, 0) === CliffordNumber{VGA(3)}(1)
    @test OddCliffordNumber{VGA(3)}(0, 0, 0, 0) - 1 === CliffordNumber{VGA(3)}(-1)
    @test OddCliffordNumber{VGA(3)}(0, 0, 0, 0) - 1 === -CliffordNumber{VGA(3)}(1)
end

@testset "muladd" begin
    e = EvenCliffordNumber{VGA(3),Float32}(0, 4, 2, 0)
    ee = EvenCliffordNumber{VGA(3),Float64}(0, 4, 2, 0)
    f = EvenCliffordNumber{VGA(3),Float32}(0, 0, 6, 9)
    ff = EvenCliffordNumber{VGA(3),Float64}(0, 0, 6, 9)
    @test muladd(2, ee, ff) === EvenCliffordNumber{VGA(3),Float64}(0, 8, 10, 9)
    @test muladd(2, ff, ee) === EvenCliffordNumber{VGA(3),Float64}(0, 4, 14, 18)
    @test muladd(2, e, f) === muladd(e, 2, f)
    @test muladd(2, e, ff) == muladd(2, ee, f)
end

@testset "Geometric product" begin
    x = CliffordNumber{VGA(3)}(0, 2, 0, 0, 0, 0, 0, 0)
    y = CliffordNumber{VGA(3)}(0, 3, 4, 0, 0, 0, 0, 0)
    k1 = KVector{1,VGA(3)}(4, 2, 0)
    k2 = KVector{2,VGA(3)}(4, 2, 0)
    l1 = KVector{1,VGA(3)}(0, 6, 9)
    l2 = KVector{2,VGA(3)}(0, 6, 9)
    # Signs included explicilty for clarity
    @test k1 * l1 === EvenCliffordNumber{VGA(3)}(+12, +24, +36, +18)
    @test l1 * k1 === EvenCliffordNumber{VGA(3)}(+12, -24, -36, -18)
    @test k1 * l1 === reverse(l1 * k1)
    @test k2 * l2 === EvenCliffordNumber{VGA(3)}(-12, -18, +36, -24)
    @test l2 * k2 === EvenCliffordNumber{VGA(3)}(-12, +18, -36, +24)
    @test k2 * l2 === reverse(l2 * k2)
    @test x * y == CliffordNumber{VGA(3)}(6, 0, 0, 8, 0, 0, 0, 0)
    # When multiplying vectors, the reverse should have a negative bivector
    @test y * x == CliffordNumber{VGA(3)}(6, 0, 0, -8, 0, 0, 0, 0)
    @test x(y) === x * y
    @test (y)(x) === y * x
    @test_throws CliffordNumbers.AlgebraMismatch x * one(CliffordNumber{STA})
    # Check that scalar/pseudoscalar multiplications promote correctly
    @test k1 * KVector{0,VGA(3)}(2) === KVector{1,VGA(3)}(8, 4, 0)
    @test KVector{0,VGA(3)}(2) * k1 === KVector{1,VGA(3)}(8, 4, 0)
    # Test kernel directly for the scalar case (this is overridden)
    @test CliffordNumbers.mul(k1, KVector{0,VGA(3)}(2)) === KVector{1,VGA(3)}(8, 4, 0)
    @test CliffordNumbers.mul(KVector{0,VGA(3)}(2), k1) === KVector{1,VGA(3)}(8, 4, 0)
    @test k1 * KVector{3,VGA(3)}(2) === KVector{2,VGA(3)}(0, -4,  8)
    @test KVector{3,VGA(3)}(2) * k1 === KVector{2,VGA(3)}(0, -4,  8)
    # Test promotions of KVector with Z2CliffordNumber implicitly
    @test k1 * k2 * l2 isa OddCliffordNumber{VGA(3)}
    @test k1 * l1 * l2 isa EvenCliffordNumber{VGA(3)}
end

@testset "Scalars" begin
    k = KVector{2,VGA(3)}(3, 4, 0)
    @test abs2(k) === 25
    @test abs(k) == 5
    @test normalize(k) == k / 5
    # Testing with KVector{0}
    @test KVector{0,VGA(3)}(2) * k === 2 * k
    @test k * KVector{0,VGA(3)}(2) === k * 2
    @test KVector{0,VGA(3)}(2) * KVector{0,VGA(3)}(3) === KVector{0,VGA(3)}(6)
    @test isscalar(CliffordNumber{VGA(3)}(1))
    @test isscalar(EvenCliffordNumber{VGA(3)}(1, 0, 0, 0))
    @test !isscalar(EvenCliffordNumber{VGA(3)}(1, 2, 3, 4))
    @test !isscalar(OddCliffordNumber{VGA(3)}(1, 2, 3, 4))
    @test isscalar(OddCliffordNumber{VGA(3)}(0, 0, 0, 0))
    @test isscalar(KVector{1,STA}(0, 0, 0, 0))
    @test !isscalar(KVector{3,STA}(1, 3, 3, 7))
    @test isscalar(complex(420, 69))
end

@testset "Inverses and division" begin
    k = KVector{1,VGA(3)}(1, 2, 3)
    l = KVector{2,VGA(3)}(4, 5, 6)
    @test k * (1 / k) ≈ 1
    @test (1 / k) * k ≈ 1
    @test inv(k) ≈ 1 / k
    @test l * (1 / l) ≈ -1
    @test (1 / l) * l ≈ -1
    @test inv(l) ≈ -1 / l   # This one fails for exact equality, but I don't know why
    @test inv(KVector{0,VGA(2)}(2)) == KVector{0,VGA(2)}(1//2)
    @test inv(KVector{0,VGA(2)}(2)) === KVector{0,VGA(2)}(inv(2))
    @test_throws CliffordNumbers.InverseException inv(1 + KVector{1,VGA(2)}(1, 0))
end

@testset "Contractions and dot products" begin
    import CliffordNumbers: dot, hestenes_dot
    k = KVector{1,VGA(3)}(1, 2, 3)
    l = KVector{2,VGA(3)}(4, 5, 6)
    @test k ⨼ k === k ⨽ k
    @test l ⨼ l === l ⨽ l
    @test iszero(l ⨼ k)
    @test iszero(k ⨽ l)
    @test dot(k, l) === dot(l, k)
    @test dot(k, l) === k ⨼ l
    @test dot(l, k) === (k ⨼ l) * 1^(grade(l) * (grade(k) - grade(l)))
end

@testset "Wedge product" begin
    x = KVector{1,VGA(3)}(1, 0, 0)
    y = KVector{1,VGA(3)}(0, 1, 0)
    z = KVector{1,VGA(3)}(0, 0, 1)
    k1 = KVector{1,VGA(3)}(3, 4, 0)
    k2 = KVector{2,VGA(3)}(4, 2, 0)
    five = 5 * one(CliffordNumber{VGA(3)})
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
    @test KVector{1,VGA(3)}(x) ∧ KVector{1,VGA(3)}(y) == -(KVector{1,VGA(3)}(y) ∧ KVector{1,VGA(3)}(x))
    # Ensure the behavior of scalars (both CliffordNumber and normal) are correct
    @test 5 ∧ x == 5 * x
    @test five ∧ x == x ∧ five
    @test 5 ∧ 5 == 25
    # Turns out reversing the order shouldn't change anything
    @test k1 ∧ k2 == KVector{3,VGA(3)}(-8)
end

@testset "Exponentiation" begin
    e12 = CliffordNumber{VGA(3)}(0, 0, 0, 1, 0, 0, 0, 0)
    @test exppi(e12/2) == e12
    @test exptau(e12) == 1
    @test exppi(e12) ≈ exp(pi*e12)
    @test exptau(e12) ≈ exp(2*pi*e12)
    @test exp(KVector{1,VGA(3)}(1,0,0)) isa CliffordNumber
    @test exp(KVector{2,VGA(3)}(1,0,0)) isa EvenCliffordNumber
    k = KVector{2,VGA(3)}(1,0,0)
    @test CliffordNumbers.exp_taylor(pi/2 * k) ≈ exp(pi/2 * k)
    @test CliffordNumbers.exp_taylor(pi/2 * k) ≈ exppi(1//2 * k)
    @test CliffordNumbers.exp_taylor(pi/2 * k) ≈ exptau(1//4 * k)
    # Integer exponentiation of a multivector
    k = KVector{1,VGA(3)}(4, 2, 0)
    l = KVector{2,VGA(3)}(0, 6, 9)
    # Base.literal_pow tests
    @test k^-2 === inv(k) * inv(k)
    @test l^-2 === inv(l) * inv(l)
    @test k^-1 === inv(k)
    @test l^-1 === inv(l)
    @test k^0 === one(k)
    @test l^0 === one(l)
    @test k^1 === k
    @test l^1 === l
    @test k^2 isa EvenCliffordNumber
    @test l^2 isa EvenCliffordNumber
    @test k^3 isa OddCliffordNumber
    @test l^3 isa EvenCliffordNumber
end
