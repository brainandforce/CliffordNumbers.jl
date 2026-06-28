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
    # NaNs should be unequal with ==
    @test KVector{2,VGA(3)}(NaN, 6, 9) != KVector{2,VGA(3)}(NaN, 6, 9)
    # Signed zeros should be equal with ==
    @test KVector{2,VGA(3)}(0.0, 6, 9) == KVector{2,VGA(3)}(-0.0, 6, 9)
    # NaNs should be equal with isequal
    @test isequal(KVector{2,VGA(3)}(NaN, 6, 9), KVector{2,VGA(3)}(NaN, 6, 9))
    # Signed zeros should be unequal with ==
    @test !isequal(KVector{2,VGA(3)}(0.0, 6, 9), KVector{2,VGA(3)}(-0.0, 6, 9))
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

@testset "Sign automorphism parity across types" begin
    # The closed-form (KVector) and `@generated` (dense / Z2) implementations must agree with the
    # original generic indexing semantics `T(x[f.(BladeIndices(T))])` for every type and grade.
    # This also locks in `reverse(::KVector)` and the GPU-safe generated paths in src/math/duals.jl.
    refauto(f, x) = (T = typeof(x); T(x[f.(CliffordNumbers.BladeIndices(T))]))
    types = (
        CliffordNumber{VGA(3)}, EvenCliffordNumber{VGA(3)}, OddCliffordNumber{VGA(3)},
        KVector{0,VGA(3)}, KVector{1,VGA(3)}, KVector{2,VGA(3)}, KVector{3,VGA(3)},
        CliffordNumber{STA}, EvenCliffordNumber{PGA(3)}, KVector{2,PGA(3)},
    )
    for C in types
        x = C(ntuple(i -> Float64(i) - 2.5, nblades(C)))
        for f in (reverse, adjoint, conj, grade_involution)
            @test f(x) === refauto(f, x)
        end
    end
end

@testset "Complements and duals" begin
    # We don't have support for the wedge product of BladeIndex
    # But assuming an orthonormal basis, the geometric product suffices
    # VGA
    ps = σ1 ∧ σ2 ∧ σ3
    ps_index = BladeIndex(ps, 1, 2, 3)
    @test left_complement(σ1) ∧ σ1  === ps
    @test σ1 ∧ right_complement(σ1) === ps
    @test left_complement(σ1 ∧ σ2) ∧ (σ1 ∧ σ2)  === ps
    @test (σ1 ∧ σ2) ∧ right_complement(σ1 ∧ σ2) === ps
    @test left_complement(BladeIndex(σ1, 1)) * BladeIndex(σ1, 1)  === ps_index
    @test BladeIndex(σ1, 1) * right_complement(BladeIndex(σ1, 1)) === ps_index
    @test left_complement(BladeIndex(σ1 ∧ σ2, 1, 2)) * BladeIndex(σ1 ∧ σ2, 1, 2)  === ps_index
    @test BladeIndex(σ1 ∧ σ2, 1, 2) * right_complement(BladeIndex(σ1 ∧ σ2, 1, 2)) === ps_index
    # STA
    ps = γ0 ∧ γ1 ∧ γ2 ∧ γ3
    ps_index = BladeIndex(ps, 0, 1, 2, 3)
    @test left_complement(γ0) ∧ γ0  === ps
    @test γ0 ∧ right_complement(γ0) === ps
    @test left_complement(γ3) ∧ γ3  === ps
    @test γ3 ∧ right_complement(γ3) === ps
    @test left_complement(γ0 ∧ γ1) ∧ (γ0 ∧ γ1)  === ps
    @test (γ0 ∧ γ1) ∧ right_complement(γ0 ∧ γ1) === ps
    @test left_complement(γ2 ∧ γ3) ∧ (γ2 ∧ γ3)  === ps
    @test (γ2 ∧ γ3) ∧ right_complement(γ2 ∧ γ3) === ps
    @test left_complement(BladeIndex(γ0, 0)) * BladeIndex(γ0, 0)  === ps_index
    @test BladeIndex(γ0, 0) * right_complement(BladeIndex(γ0, 0)) === ps_index
    @test left_complement(BladeIndex(γ1, 1)) * BladeIndex(γ1, 1)  === ps_index
    @test BladeIndex(γ1, 1) * right_complement(BladeIndex(γ1, 1)) === ps_index
    @test left_complement(BladeIndex(γ0 ∧ γ1, 0, 1)) * BladeIndex(γ0 ∧ γ1, 0, 1)  === ps_index
    @test BladeIndex(γ0 ∧ γ1, 0, 1) * right_complement(BladeIndex(γ0 ∧ γ1, 0, 1)) === ps_index
    @test left_complement(BladeIndex(γ2 ∧ γ3, 2, 3)) * BladeIndex(γ2 ∧ γ3, 2, 3)  === ps_index
    @test BladeIndex(γ2 ∧ γ3, 2, 3) * right_complement(BladeIndex(γ2 ∧ γ3, 2, 3)) === ps_index
    # PGA
    ps = e0 ∧ e1 ∧ e2 ∧ e3
    ps_index = BladeIndex(ps, 0, 1, 2, 3)
    @test left_complement(e0) ∧ e0  === ps
    @test e0 ∧ right_complement(e0) === ps
    @test left_complement(e3) ∧ e3  === ps
    @test e3 ∧ right_complement(e3) === ps
    @test left_complement(e0 ∧ e1) ∧ (e0 ∧ e1)  === ps
    @test (e0 ∧ e1) ∧ right_complement(e0 ∧ e1) === ps
    @test left_complement(e2 ∧ e3) ∧ (e2 ∧ e3)  === ps
    @test (e2 ∧ e3) ∧ right_complement(e2 ∧ e3) === ps
    @test left_complement(BladeIndex(e0, 0)) * BladeIndex(e0, 0)  === ps_index
    @test BladeIndex(e0, 0) * right_complement(BladeIndex(e0, 0)) === ps_index
    @test left_complement(BladeIndex(e1, 1)) * BladeIndex(e1, 1)  === ps_index
    @test BladeIndex(e1, 1) * right_complement(BladeIndex(e1, 1)) === ps_index
    @test left_complement(BladeIndex(e0 ∧ e1, 0, 1)) * BladeIndex(e0 ∧ e1, 0, 1)  === ps_index
    @test BladeIndex(e0 ∧ e1, 0, 1) * right_complement(BladeIndex(e0 ∧ e1, 0, 1)) === ps_index
    @test left_complement(BladeIndex(e2 ∧ e3, 2, 3)) * BladeIndex(e2 ∧ e3, 2, 3)  === ps_index
    @test BladeIndex(e2 ∧ e3, 2, 3) * right_complement(BladeIndex(e2 ∧ e3, 2, 3)) === ps_index
    # Test linear extension of the complement
    @test left_complement(6*σ1 + 9*σ2)  === 6*left_complement(σ1)  + 9*left_complement(σ2)
    @test right_complement(6*σ1 + 9*σ2) === 6*right_complement(σ1) + 9*right_complement(σ2)
end

@testset "Regressive product and generated complements" begin
    # The generated complements must be bit-for-bit identical to the prior generic indexing path
    # `C(x[blade_complement.(BladeIndices(C))])` across every storage family and signature,
    # including mixed (CGA), degenerate (PGA), and Lorentzian (STA) metrics.
    function ref_left_complement(x)
        C = CliffordNumbers.complement_type(typeof(x))
        return C(x[right_complement.(CliffordNumbers.BladeIndices(C))])
    end
    function ref_right_complement(x)
        C = CliffordNumbers.complement_type(typeof(x))
        return C(x[left_complement.(CliffordNumbers.BladeIndices(C))])
    end
    for Q in (VGA(2), VGA(3), VGA(4), PGA(2), PGA(3), CGA(2), CGA(3), STA)
        D = CliffordNumbers.dimension(Q)
        types = Any[CliffordNumber{Q,Int}, EvenCliffordNumber{Q,Int}, OddCliffordNumber{Q,Int}]
        append!(types, [KVector{K,Q,Int} for K in 0:D])
        for T in types
            n = CliffordNumbers.nblades(T)
            x = T(ntuple(i -> i - 3, n))    # deterministic mixed-sign coefficients for === checks
            @test left_complement(x)  === ref_left_complement(x)
            @test right_complement(x) === ref_right_complement(x)
            # Round trip: right_complement ∘ left_complement and its mirror are the identity
            @test right_complement(left_complement(x)) === x
            @test left_complement(right_complement(x)) === x
        end
    end
    # Regressive product parity against the dense CliffordNumber path, over PGA(2)/PGA(3) grade pairs
    for Q in (PGA(2), PGA(3))
        D = CliffordNumbers.dimension(Q)
        for K1 in 0:D, K2 in 0:D
            n1 = binomial(D, K1); n2 = binomial(D, K2)
            a = KVector{K1,Q,Int}(ntuple(i -> i - 2, n1))
            b = KVector{K2,Q,Int}(ntuple(i -> 3 - i, n2))
            @test CliffordNumber(a ∨ b) === CliffordNumber(a) ∨ CliffordNumber(b)
        end
    end
    # The compact KVector ∨ KVector path allocates nothing
    join_op(a, b) = a ∨ b
    p1 = KVector{3,PGA(3)}(1.0, 2.0, 3.0, 4.0)
    p2 = KVector{3,PGA(3)}(4.0, 3.0, 2.0, 1.0)
    join_op(p1, p2)     # warm up
    @test @allocated(join_op(p1, p2)) == 0
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
    @test muladd(2,  ee, ff) === EvenCliffordNumber{VGA(3),Float64}(0, 8, 10, 9)
    @test muladd(2,  ff, ee) === EvenCliffordNumber{VGA(3),Float64}(0, 4, 14, 18)
    @test muladd(2,  e,  f) === muladd(e, 2, f)
    @test muladd(2,  e,  ff) == muladd(2, ee, f)
    @test muladd(ee, 2,  ff) === EvenCliffordNumber{VGA(3),Float64}(0, 8, 10, 9)
    @test muladd(ff, 2,  ee) === EvenCliffordNumber{VGA(3),Float64}(0, 4, 14, 18)
    @test muladd(e,  2,  f) === muladd(e, 2, f)
    @test muladd(e,  2,  ff) == muladd(2, ee, f)
end

@testset "Geometric product" begin
    import CliffordNumbers.AlgebraMismatch as GAMismatch
    x = CliffordNumber{VGA(3)}(0, 2, 0, 0, 0, 0, 0, 0)
    y = CliffordNumber{VGA(3)}(0, 3, 4, 0, 0, 0, 0, 0)
    @test KVector{0,VGA(3)}(2) * KVector{0,VGA(3)}(3) === KVector{0,VGA(3)}(6)
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

@testset "Rational kernels" begin
    import CliffordNumbers: mul, GradeFilter
    x = EvenCliffordNumber{VGA(3)}(1//2, 3//4, 5//6, 7//8)
    y = OddCliffordNumber{VGA(3)}(2//1, 4//3, 6//5, 8//7)
    z = EvenCliffordNumber{VGA(3)}(1, 3, 3, 7)
    @test x * y === CliffordNumbers.mul(x, y)
    @test x ∧ y === CliffordNumbers.mul(x, y, GradeFilter{:∧}())
    @test y * z === CliffordNumbers.mul(y, z)
    @test z * y === CliffordNumbers.mul(z, y)
end

@testset "Scalar and pseudoscalar components" begin
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

@testset "Scalar products and normalization" begin
    k = KVector{1,VGA(3)}(3, 4, 0)
    l = KVector{2,VGA(3)}(0, 5, 12)
    # In case we decide scalar_product should return a KVector{0}, the tests shouldn't break
    @test scalar_product(k, k) == 25
    @test scalar_product(OddCliffordNumber(k), CliffordNumber(k)) == 25
    @test scalar_product(l, l) == -169
    @test scalar_product(EvenCliffordNumber(l), EvenCliffordNumber(l)) == -169
    @test_throws CliffordNumbers.AlgebraMismatch scalar_product(one(KVector{0,STA}), k)
    @test iszero(scalar_product(k, l))  # no indices in commmon between k and l
    @test abs2(k) === 25
    @test abs(k) == 5
    @test abs2(l) === 169
    @test abs(l) == 13
    @test normalize(k) === k / 5
    @test normalize(l) === l / 13
    @test normalize(zero(k)) === zero(k)
    # Testing non-positive definite and degenerate metrics
    @test abs2(KVector{1,STA}(0, 2, 3, 4)) == -29
    @test abs2(KVector{1,PGA(3)}(1, 0, 0, 0)) == 0
    @test abs2(KVector{0,VGA(2)}(2)) === 4
    @test abs(KVector{0,VGA(2)}(-2)) === 2
    @test normalize(KVector{0,VGA(2),Float16}(3)) === one(KVector{0,VGA(2),Float16})
    @test normalize(KVector{1,STA}(3, 0, 0, 0)) === KVector{1,STA}(3, 0, 0, 0) / 3
    @test normalize(KVector{1,STA}(-3, 0, 0, 0)) === KVector{1,STA}(-3, 0, 0, 0) / 3
    @test normalize(KVector{1,STA}(0, 0, 0, 3)) === KVector{1,STA}(0, 0, 0, 3) / 3
    @test normalize(KVector{1,STA}(0, 0, 0, -3)) === KVector{1,STA}(0, 0, 0, -3) / 3
    @test normalize(KVector{1,STA}(2, 2, 0, 0)) === KVector{1,STA}(2, 2, 0, 0)
    @test normalize(KVector{0,STA,Float32}(420)) === one(KVector{0,STA,Float32})
    @test normalize(KVector{0,STA,Int8}(0)) === zero(KVector{0,STA,Int8})
end

@testset "Addition and multiplication with scalars" begin
    k = KVector{1,VGA(3)}(4, 2, 0)
    l = KVector{2,VGA(3)}(0, 6, 9)
    @test 1 + k === CliffordNumber{VGA(3)}(1, 4, 2, 0, 0, 0, 0, 0)
    @test k - 1 === CliffordNumber{VGA(3)}(-1, 4, 2, 0, 0, 0, 0, 0)
    @test 1 + l === EvenCliffordNumber{VGA(3)}(1, 0, 6, 9)
    @test l - 1 === EvenCliffordNumber{VGA(3)}(-1, 0, 6, 9)
    @test 2 * k === KVector{1,VGA(3)}(8, 4, 0)
    @test 2(k)  === KVector{1,VGA(3)}(8, 4, 0)
    @test k(2)  === KVector{1,VGA(3)}(8, 4, 0)
    @test k * 2 === k * KVector{0,VGA(3)}(2)
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

@testset "Even subalgebra fast paths" begin
    import CliffordNumbers: versor_inverse
    # Split (e₁₂² = +1) and dual (e₁₂² = 0) dimension-2 signatures, to cover every σ branch of the
    # closed-form spinor product alongside VGA(2) (e₁₂² = -1).
    Qsplit = CliffordNumbers.Metrics.Signature(2, 0b10, 0b00)   # e₂² = -1 ⟹ e₁₂² = +1
    Qdual = PGA(1)                                              # e₀² =  0 ⟹ e₀₁² =  0

    @testset "Closed-form spinor product / square" begin
        # Exact integer cases pin the (a·c + σ·b·d, a·d + b·c) formula per signature
        @test EvenCliffordNumber{VGA(2)}(2, 3) * EvenCliffordNumber{VGA(2)}(4, 5) ===
            EvenCliffordNumber{VGA(2)}(2*4 - 3*5, 2*5 + 3*4)        # σ = -1 → (-7, 22)
        @test EvenCliffordNumber{Qsplit}(2, 3) * EvenCliffordNumber{Qsplit}(4, 5) ===
            EvenCliffordNumber{Qsplit}(2*4 + 3*5, 2*5 + 3*4)        # σ = +1 → (23, 22)
        @test EvenCliffordNumber{Qdual}(2, 3) * EvenCliffordNumber{Qdual}(4, 5) ===
            EvenCliffordNumber{Qdual}(2*4, 2*5 + 3*4)               # σ =  0 → (8, 22)
        @test EvenCliffordNumber{VGA(2)}(2, 3)^2 === EvenCliffordNumber{VGA(2)}(2*2 - 3*3, 2*2*3)
        # The closed form must agree with the generic kernel (reached via CliffordNumber promotion)
        # and, for VGA(2), with ℂ. The wedge must still take the generic path, not the `*` fast path.
        for _ in 1:50
            (a, b, c, d) = randn(4)
            s1 = EvenCliffordNumber{VGA(2)}(a, b)
            s2 = EvenCliffordNumber{VGA(2)}(c, d)
            generic = CliffordNumber{VGA(2)}(s1) * CliffordNumber{VGA(2)}(s2)
            @test s1 * s2 ≈ EvenCliffordNumber{VGA(2)}(generic)
            @test all(Tuple(s1 * s2) .≈ reim(Complex(a, b) * Complex(c, d)))
            @test s1^2 ≈ EvenCliffordNumber{VGA(2)}(reim(Complex(a, b)^2)...)
        end
        @test EvenCliffordNumber{VGA(2)}(2, 3) ∧ EvenCliffordNumber{VGA(2)}(4, 5) ===
            EvenCliffordNumber{VGA(2)}(2*4, 2*5 + 3*4)              # wedge: e₁₂ ∧ e₁₂ = 0
    end

    @testset "Rotor (ℍ) geometric product / square" begin
        # Even VGA(3) × even VGA(3) (the rotor product) was not directly covered before. Exact
        # integer case pins the Hamilton-product structure.
        @test EvenCliffordNumber{VGA(3)}(2, 3, 5, 7) * EvenCliffordNumber{VGA(3)}(11, 13, 17, 19) ===
            EvenCliffordNumber{VGA(3)}(-235, 83, 55, 129)
        # Agreement with the 8-component product reached via CliffordNumber promotion.
        for _ in 1:50
            x = EvenCliffordNumber{VGA(3)}(ntuple(_ -> randn(), 4))
            y = EvenCliffordNumber{VGA(3)}(ntuple(_ -> randn(), 4))
            @test x * y ≈ EvenCliffordNumber{VGA(3)}(CliffordNumber{VGA(3)}(x) * CliffordNumber{VGA(3)}(y))
            @test x^2 ≈ EvenCliffordNumber{VGA(3)}(CliffordNumber{VGA(3)}(x)^2)
            @test x ∧ y ≈ EvenCliffordNumber{VGA(3)}(CliffordNumber{VGA(3)}(x) ∧ CliffordNumber{VGA(3)}(y))
        end
    end

    @testset "Closed-form abs2 (VGA(2) ≅ ℂ, VGA(3) ≅ ℍ)" begin
        @test abs2(EvenCliffordNumber{VGA(2)}(3, 4)) === 25
        @test abs2(EvenCliffordNumber{VGA(3)}(1, 2, 3, 4)) === 30
        for _ in 1:50
            s = EvenCliffordNumber{VGA(2)}(randn(), randn())
            r = EvenCliffordNumber{VGA(3)}(randn(), randn(), randn(), randn())
            @test abs2(s) ≈ scalar_product(s, s')
            @test abs2(r) ≈ scalar_product(r, r')
        end
    end

    @testset "Validation-free / closed-form inv and versor_inverse" begin
        # Generic versor formula, built from methods this PR does not specialize, as a reference
        generic_versor(x) = adjoint(x) / scalar_product(x, x')
        for _ in 1:50
            s = EvenCliffordNumber{VGA(2)}(randn(), randn())
            r = EvenCliffordNumber{VGA(3)}(randn(), randn(), randn(), randn())
            @test s * inv(s) ≈ one(s)
            @test inv(s) * s ≈ one(s)
            @test inv(s) === versor_inverse(s)              # inv delegates, no validation
            @test versor_inverse(s) ≈ generic_versor(s)     # closed form == x'/abs2(x)
            @test all(Tuple(inv(s)) .≈ reim(inv(Complex(Tuple(s)...))))
            @test r * inv(r) ≈ one(r)
            @test inv(r) === versor_inverse(r)
            @test versor_inverse(r) ≈ generic_versor(r)
        end
        # Integer input promotes to the same float type as the generic versor formula
        @test versor_inverse(EvenCliffordNumber{VGA(2)}(3, 4)) isa EvenCliffordNumber{VGA(2),Float64,2}
        @test typeof(versor_inverse(EvenCliffordNumber{VGA(3)}(1, 2, 3, 4))) ===
            typeof(generic_versor(EvenCliffordNumber{VGA(3)}(1, 2, 3, 4)))
    end
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
        @test k^z === inv(CliffordNumber(k))
        @test l^z === inv(EvenCliffordNumber(l))
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
    # These are not defined and should error
    @test_throws DomainError Base.literal_pow(^, k, Val(1//2))
    @test_throws DomainError Base.literal_pow(^, l, Val(1/3))
    @test_throws DomainError Base.literal_pow(^, CliffordNumber(m), Val(2*im))
    #= If KVector{0} passes through exponentiation, these *are* defined
    @test Base.literal_pow(^, m, Val(1//2)) === Base.literal_pow(^, scalar(m), Val(1//2))
    @test Base.literal_pow(^, m, Val(1/3)) === Base.literal_pow(^, scalar(m), Val(1/3))
    @test Base.literal_pow(^, m, Val(2*im)) === Base.literal_pow(^, scalar(m), Val(2*im))
    =#
    @test sqrt(KVector{0,VGA(3)}(1)) === KVector{0,VGA(3)}(sqrt(1))
    @test sqrt(KVector{0,VGA(3)}(2.0)) === KVector{0,VGA(3)}(sqrt(2.0))
    @test_throws DomainError sqrt(KVector{0,VGA(3)}(-1))
end

@testset "exp of mixed-signature bivectors" begin
    # In mixed-signature algebras (CGA, STA, ...) a bivector's `abs2` can be negative; when it lands
    # in (-1, 0) the power-of-two scaling exponent in `exp_taylor` went negative and `exp` threw a
    # DomainError. Construct such a bivector explicitly (abs2 ≈ -1/2) for each affected algebra.
    for Q in (CGA(3), STA)
        nb = binomial(CliffordNumbers.dimension(Q), 2)
        rng = MersenneTwister(0xCA75)
        local b
        # a bivector that squares to a negative scalar, rescaled so abs2 lands in the danger zone
        while true
            b = KVector{2,Q}(ntuple(_ -> randn(rng), nb))
            CliffordNumbers.abs2(b) < -0.1 && break
        end
        b *= sqrt(0.5 / -CliffordNumbers.abs2(b))
        @test -1 < CliffordNumbers.abs2(b) < 0          # confirm we are in the danger zone
        r = @test_nowarn exp(b)                         # formerly threw a DomainError
        @test all(isfinite, Tuple(r))
        @test isapprox(scalar(r * r'), 1; atol = 1e-8)  # exp(b) is a unit rotor
    end
end

@testset "Logarithms and square roots of rotors" begin
    rng = MersenneTwister(0x10AD)
    # `exp`/`log` and `exp`/`sqrt` round-trips across vanilla, projective, conformal, and
    # spacetime algebras, covering elliptic, hyperbolic, and degenerate (ideal) planes.
    for Q in (VGA(2), VGA(3), VGA(4), PGA(2), PGA(3), PGA(4), STA, CGA(3))
        d = CliffordNumbers.dimension(Q)
        nb = binomial(d, 2)
        for _ in 1:50
            B = KVector{2,Q}(ntuple(_ -> 0.5 * randn(rng), nb))
            R = exp(B)
            # `log` inverts `exp` (principal branch round-trips through `exp`)
            @test isapprox(exp(log(R)), R; atol = 1e-8)
            # `log` returns a bivector
            @test log(R) isa KVector{2,Q}
            # `sqrt` squares back to the rotor and halves the generator
            S = sqrt(R)
            @test isapprox(S * S, R; atol = 1e-8)
            @test S isa EvenCliffordNumber{Q}
        end
    end
    # Principal branch: `log(exp(B)) == B` for a bivector of small magnitude
    for Q in (VGA(3), PGA(3), STA, CGA(3))
        d = CliffordNumbers.dimension(Q)
        B = KVector{2,Q}(ntuple(i -> 0.1 / i, binomial(d, 2)))
        @test isapprox(log(exp(B)), B; atol = 1e-12)
    end
    # Isoclinic rotor (coincident invariant planes) is a documented edge case
    iso = KVector{2,VGA(4)}(0.7, 0, 0, 0, 0, 0.7)
    @test isapprox(exp(log(exp(iso))), exp(iso); atol = 1e-10)
    @test isapprox(log(exp(iso)), iso; atol = 1e-10)
    # Type stability and allocation-freedom on the simple (≅ ℍ) path
    R32 = exp(KVector{2,VGA(3),Float32}(0.3f0, 0.2f0, -0.1f0))
    @test log(R32) isa KVector{2,VGA(3),Float32}
    @test sqrt(R32) isa EvenCliffordNumber{VGA(3),Float32}
    @test (@allocated log(R32)) == 0
end

@testset "Miscellaneous properties" begin
    # Finiteness
    @test isfinite(EvenCliffordNumber{VGA(3)}(0, 1, 2, 3))
    @test !isfinite(EvenCliffordNumber{VGA(3)}(Inf, 1, 2, 3))
    @test isinf(EvenCliffordNumber{VGA(3)}(Inf, 1, 2, 3))
    @test !isinf(EvenCliffordNumber{VGA(3)}(0, 1, 2, 3))
    # NaN presence
    @test isnan(EvenCliffordNumber{VGA(3)}(NaN, 1, 2, 3))
    @test !isnan(EvenCliffordNumber{VGA(3)}(0, 1, 2, 3))
    # Equivalence to real values
    @test isreal(one(EvenCliffordNumber{VGA(3),Int}))
    @test isreal(one(EvenCliffordNumber{VGA(3),Float64}))
    @test isreal(one(EvenCliffordNumber{VGA(3),Complex{Float64}}))
    @test !isreal(EvenCliffordNumber{VGA(3)}(0, 1, 2, 3))
    @test !isreal(EvenCliffordNumber{VGA(3)}(im, 1, 2, 3))
    @test !isreal(EvenCliffordNumber{VGA(3)}(1 + 0im, 1, 2, 3))
    # Equivalence to integer values
    @test isinteger(zero(KVector{1,VGA(3)}))
    @test isinteger(one(EvenCliffordNumber{VGA(3),Int}))
    @test isinteger(one(EvenCliffordNumber{VGA(3),Float64}))
    @test isinteger(one(EvenCliffordNumber{VGA(3),Complex{Float64}}))
    @test !isinteger(π * one(EvenCliffordNumber{VGA(3),Float64}))
    @test !isinteger(EvenCliffordNumber{VGA(3)}(0, 1, 2, 3))
    @test !isinteger(KVector{0,VGA(3)}(1 + im))
    # Even and odd values
    @test iseven(KVector{0,VGA(3)}(0))
    @test !iseven(KVector{0,VGA(3)}(1))
    @test isodd(KVector{0,VGA(3)}(1))
    @test !isodd(KVector{0,VGA(3)}(0))
    @test !iseven(EvenCliffordNumber{VGA(3)}(1, 2, 3, 4))
    @test !isodd(EvenCliffordNumber{VGA(3)}(1, 2, 3, 4))
end
