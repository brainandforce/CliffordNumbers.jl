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

@testset "Complements and duals" begin
    # We don't have support for the wedge product of BitIndex
    # But assuming an orthonormal basis, the geometric product suffices
    # VGA
    ps = σ1 ∧ σ2 ∧ σ3
    ps_index = BitIndex(ps, 1, 2, 3)
    @test left_complement(σ1) ∧ σ1  === ps
    @test σ1 ∧ right_complement(σ1) === ps
    @test left_complement(σ1 ∧ σ2) ∧ (σ1 ∧ σ2)  === ps
    @test (σ1 ∧ σ2) ∧ right_complement(σ1 ∧ σ2) === ps
    @test left_complement(BitIndex(σ1, 1)) * BitIndex(σ1, 1)  === ps_index
    @test BitIndex(σ1, 1) * right_complement(BitIndex(σ1, 1)) === ps_index
    @test left_complement(BitIndex(σ1 ∧ σ2, 1, 2)) * BitIndex(σ1 ∧ σ2, 1, 2)  === ps_index
    @test BitIndex(σ1 ∧ σ2, 1, 2) * right_complement(BitIndex(σ1 ∧ σ2, 1, 2)) === ps_index
    # STA
    ps = γ0 ∧ γ1 ∧ γ2 ∧ γ3
    ps_index = BitIndex(ps, 0, 1, 2, 3)
    @test left_complement(γ0) ∧ γ0  === ps
    @test γ0 ∧ right_complement(γ0) === ps
    @test left_complement(γ3) ∧ γ3  === ps
    @test γ3 ∧ right_complement(γ3) === ps
    @test left_complement(γ0 ∧ γ1) ∧ (γ0 ∧ γ1)  === ps
    @test (γ0 ∧ γ1) ∧ right_complement(γ0 ∧ γ1) === ps
    @test left_complement(γ2 ∧ γ3) ∧ (γ2 ∧ γ3)  === ps
    @test (γ2 ∧ γ3) ∧ right_complement(γ2 ∧ γ3) === ps
    @test left_complement(BitIndex(γ0, 0)) * BitIndex(γ0, 0)  === ps_index
    @test BitIndex(γ0, 0) * right_complement(BitIndex(γ0, 0)) === ps_index
    @test left_complement(BitIndex(γ1, 1)) * BitIndex(γ1, 1)  === ps_index
    @test BitIndex(γ1, 1) * right_complement(BitIndex(γ1, 1)) === ps_index
    @test left_complement(BitIndex(γ0 ∧ γ1, 0, 1)) * BitIndex(γ0 ∧ γ1, 0, 1)  === ps_index
    @test BitIndex(γ0 ∧ γ1, 0, 1) * right_complement(BitIndex(γ0 ∧ γ1, 0, 1)) === ps_index
    @test left_complement(BitIndex(γ2 ∧ γ3, 2, 3)) * BitIndex(γ2 ∧ γ3, 2, 3)  === ps_index
    @test BitIndex(γ2 ∧ γ3, 2, 3) * right_complement(BitIndex(γ2 ∧ γ3, 2, 3)) === ps_index
    # PGA
    ps = e0 ∧ e1 ∧ e2 ∧ e3
    ps_index = BitIndex(ps, 0, 1, 2, 3)
    @test left_complement(e0) ∧ e0  === ps
    @test e0 ∧ right_complement(e0) === ps
    @test left_complement(e3) ∧ e3  === ps
    @test e3 ∧ right_complement(e3) === ps
    @test left_complement(e0 ∧ e1) ∧ (e0 ∧ e1)  === ps
    @test (e0 ∧ e1) ∧ right_complement(e0 ∧ e1) === ps
    @test left_complement(e2 ∧ e3) ∧ (e2 ∧ e3)  === ps
    @test (e2 ∧ e3) ∧ right_complement(e2 ∧ e3) === ps
    @test left_complement(BitIndex(e0, 0)) * BitIndex(e0, 0)  === ps_index
    @test BitIndex(e0, 0) * right_complement(BitIndex(e0, 0)) === ps_index
    @test left_complement(BitIndex(e1, 1)) * BitIndex(e1, 1)  === ps_index
    @test BitIndex(e1, 1) * right_complement(BitIndex(e1, 1)) === ps_index
    @test left_complement(BitIndex(e0 ∧ e1, 0, 1)) * BitIndex(e0 ∧ e1, 0, 1)  === ps_index
    @test BitIndex(e0 ∧ e1, 0, 1) * right_complement(BitIndex(e0 ∧ e1, 0, 1)) === ps_index
    @test left_complement(BitIndex(e2 ∧ e3, 2, 3)) * BitIndex(e2 ∧ e3, 2, 3)  === ps_index
    @test BitIndex(e2 ∧ e3, 2, 3) * right_complement(BitIndex(e2 ∧ e3, 2, 3)) === ps_index
    # Test linear extension of the complement
    @test left_complement(6*σ1 + 9*σ2)  === 6*left_complement(σ1)  + 9*left_complement(σ2)
    @test right_complement(6*σ1 + 9*σ2) === 6*right_complement(σ1) + 9*right_complement(σ2)
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
    import CliffordNumbers.AlgebraMismatch as GAMismatch
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
    # Test non-positive-definite metrics
    @test e0 * e0 === EvenCliffordNumber{PGA(3)}(0)
    @test e1 * e1 === EvenCliffordNumber{PGA(3)}(1)
    @test e0 * e1 === -e1 * e0
    @test γ0 * γ0 === EvenCliffordNumber{STA}(1)
    @test γ1 * γ1 === EvenCliffordNumber{STA}(-1)
    @test γ0 * γ1 === -γ1 * γ0
    # Test where we attempt to multiply elements of different algebras
    @test_throws GAMismatch zero(EvenCliffordNumber{VGA(3),Int}) * zero(KVector{1,STA,Int})
end

@testset "Scalars and pseudoscalars" begin
    k = KVector{2,VGA(3)}(3, 4, 0)
    @test 2(k) === KVector{2,VGA(3)}(6, 8, 0)
    @test abs2(k) === 25
    @test abs(k) == 5
    @test normalize(k) == k / 5
    # Testing non-positive definite and degenerate metrics
    @test abs2(KVector{1,STA}(0, 2, 3, 4)) == -29
    @test abs2(KVector{1,PGA(3)}(1, 0, 0, 0)) == 0
    @test abs2(KVector{0,VGA(2)}(2)) === 4
    @test abs(KVector{0,VGA(2)}(-2)) === 2
    @test normalize(KVector{0,VGA(2),Float16}(3)) === one(KVector{0,VGA(2),Float16})
    # Testing with KVector{0}
    @test KVector{0,VGA(3)}(2) * k === 2 * k
    @test k * KVector{0,VGA(3)}(2) === k * 2
    @test KVector{0,VGA(3)}(2) * KVector{0,VGA(3)}(3) === KVector{0,VGA(3)}(6)
    # Extracting scalars
    @test scalar(EvenCliffordNumber{VGA(3)}(4, 3, 2, 1)) === 4
    @test scalar(KVector{0,VGA(3),Float64}(5)) === 5.0
    @test scalar(OddCliffordNumber{VGA(3),Float32}(6, 7, 8, 9)) === zero(Float32)
    # Testing if a multivector is a scalar
    @test isscalar(CliffordNumber{VGA(3)}(1))
    @test isscalar(EvenCliffordNumber{VGA(3)}(1, 0, 0, 0))
    @test !isscalar(EvenCliffordNumber{VGA(3)}(1, 2, 3, 4))
    @test !isscalar(OddCliffordNumber{VGA(3)}(1, 2, 3, 4))
    @test isscalar(OddCliffordNumber{VGA(3)}(0, 0, 0, 0))
    @test isscalar(KVector{1,STA}(0, 0, 0, 0))
    @test !isscalar(KVector{3,STA}(1, 3, 3, 7))
    @test isscalar(complex(420, 69))
    @test ispseudoscalar(KVector{4,STA}(69))
    @test ispseudoscalar(CliffordNumber{VGA(2)}(0, 0, 0, -69))
    @test !ispseudoscalar(OddCliffordNumber{VGA(3)}(1, 3, 3, 7))
    @test !ispseudoscalar(KVector{1,VGA(3)}(4, 2, 0))
end

@testset "Inverses and division" begin
    k = KVector{1,VGA(3)}(1, 2, 3)
    l = KVector{2,VGA(3)}(4, 5, 6)
    @test k * (1 / k) ≈ 1
    @test (1 / k) * k ≈ 1
    @test inv(k) ≈ 1 / k
    @test l * (1 / l) ≈ 1
    @test (1 / l) * l ≈ 1
    @test inv(l) ≈ 1 / l   # This one fails for exact equality, but I don't know why
    @test inv(KVector{0,VGA(2)}(2)) == KVector{0,VGA(2)}(1//2)
    @test inv(KVector{0,VGA(2)}(2)) === KVector{0,VGA(2)}(inv(2))
    @test 1 / γ0 == γ0
    @test 1 / γ1 == -γ1
    @test inv(γ0) == γ0
    @test inv(γ1) == -γ1
    @test inv(γ0 * γ1) == γ0 * γ1
    @test inv(γ1 * γ2) == -γ1 * γ2
    @test γ0 / γ1 == -γ0 * γ1
    @test γ0 \ γ1 == γ0 * γ1
    @test γ1 / γ0 == -γ0 * γ1
    @test γ1 \ γ0 == γ0 * γ1
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
    k1 = 3*σ1 + 4*σ2
    k2 = 4*(σ1 ∧ σ2) + 2*(σ1 ∧ σ3)
    five = CliffordNumber{VGA(3)}(5)
    # Self wedges should be zero
    @test iszero(σ1 ∧ σ1)
    @test iszero(σ2 ∧ σ2)
    @test iszero(σ3 ∧ σ3)
    # Reversing the order should flip the sign for vectors
    @test σ1 ∧ σ2 == -(σ2 ∧ σ1)
    @test σ1 ∧ σ3 == -(σ3 ∧ σ1)
    @test σ2 ∧ σ3 == -(σ3 ∧ σ2)
    @test σ1 ∧ σ2 ∧ σ3 == -(σ3 ∧ σ2 ∧ σ1)
    @test iszero(σ1 ∧ σ2 ∧ σ1)
    @test iszero(σ1 ∧ σ2 ∧ σ3 ∧ σ1)
    # Ensure the behavior of scalars (both CliffordNumber and normal) are correct
    @test 5 ∧ σ1 == 5 * σ1
    @test 5 ∧ σ1 == σ1 ∧ 5
    @test five ∧ σ1 == σ1 ∧ five
    @test 5 ∧ 5 == 25
    # Turns out reversing the order shouldn't change anything for odd k-vector results
    @test k1 ∧ k2 == KVector{3,VGA(3)}(-8)
    # Wedge products for non-positive-definite metrics
    @test iszero(γ0 ∧ γ0)
    @test iszero(γ1 ∧ γ1)
    @test γ0 ∧ γ1 === -γ1 ∧ γ0
    @test γ2 ∧ γ3 === -γ3 ∧ γ2
    @test γ0 ∧ γ1 ∧ γ2 ∧ γ3 === γ3 ∧ γ2 ∧ γ1 ∧ γ0
    @test iszero(e0 ∧ e0)
    @test iszero(e1 ∧ e1)
    @test e0 ∧ e1 === -e1 ∧ e0
    @test e2 ∧ e3 === -e3 ∧ e2
    @test e0 ∧ e1 ∧ e2 ∧ e3 === e3 ∧ e2 ∧ e1 ∧ e0
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
    # All must consist of float types to avoid issues with negative integer exponentiation
    k = KVector{1,VGA(3),Float64}(4, 2, 0)
    l = KVector{2,VGA(3),Float64}(0, 6, 9)
    m = KVector{0,VGA(3),Float64}(20)
    # Base.literal_pow tests
    if k^-2 === Base.literal_pow(^, k, Val(2)) && l^-2 === Base.literal_pow(^, l, Val(2))
        # These tests break on Julia 1.8: it seems that Base.literal_pow is not invoked
        @test k^-2 === inv(k) * inv(k)
        @test l^-2 === inv(l) * inv(l)
    else
        # As a fallback, check approximate equality
        @test k^-2 ≈ inv(k) * inv(k)
        @test l^-2 ≈ inv(l) * inv(l)
    end
    # These tests all use Base.literal_pow because the exponents are literal integers
    @test k^-1 === inv(k)
    @test l^-1 === inv(l)
    @test m^-1 === inv(m)
    @test k^0 === one(k)
    @test l^0 === one(l)
    @test m^0 === one(m)
    @test k^1 === k
    @test l^1 === l
    @test m^1 === m
    @test k^2 === KVector{0,VGA(3),Float64}(20)
    @test l^2 == -117
    @test m^2 === KVector{0,VGA(3),Float64}(400)
    # These examples are a little more complicated
    # 1-vectors raised to odd powers can always be represented as KVector{1}, so k^3 is a KVector{1}
    # But this is not true for arbitrary k-vectors, so l^3 is a Z2CliffordNumber
    # Also, multiplying by x as opposed to EvenCliffordNumber{Q}(x) can have unexpected effects on
    # the signs of the result
    # TODO: try to ensure signs are propagated correctly in these tests
    @test k^3 === k * 20 
    @test l^3 === EvenCliffordNumber(l) * EvenCliffordNumber{VGA(3)}(-117)
    @test m^3 === KVector{0,VGA(3),Float64}(8000)
    # Redo the same tests with an exponent variable
    # The types of some of the results should change due to type stability considerations
    let z
        z = -1
        @test_broken k^z === CliffordNumber(inv(k))
        @test_broken l^z === EvenCliffordNumber(inv(l))
        @test m^z === inv(m)
        z = 0
        @test k^z === one(CliffordNumber(k))
        @test l^z === one(EvenCliffordNumber(l))
        @test m^z === one(m)
        z = 1
        @test k^z === CliffordNumber(k)
        @test l^z === EvenCliffordNumber(l)
        @test m^z === m
        z = 2
        @test k^z === CliffordNumber{VGA(3),Float64}(20)
        @test l^z === EvenCliffordNumber{VGA(3),Float64}(-117)
        @test m^z === KVector{0,VGA(3),Float64}(400)
        z = 3
        @test k^z === CliffordNumber(k) * 20
        @test l^z === EvenCliffordNumber(l) * EvenCliffordNumber{VGA(3)}(-117)
        @test m^z === KVector{0,VGA(3)}(float(8000))
    end
    # Base.literal_pow only works for integers, apparently, but these are defined
    @test Base.literal_pow(^, k, Val(false)) === one(k)
    @test Base.literal_pow(^, l, Val(false)) === one(l)
    @test Base.literal_pow(^, k, Val(true)) === k
    @test Base.literal_pow(^, l, Val(true)) === l
end
