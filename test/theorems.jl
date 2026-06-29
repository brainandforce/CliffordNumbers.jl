# Theorem-based tests drawn from the literature. These check that the package's products,
# involutions, contractions, and inverses satisfy published Clifford-algebra identities, using a
# rational-fixture convention (exact `Rational{BigInt}` coordinates, so structural grade
# content is checked with `==` rather than a tolerance; floats only where a transcendental `exp`
# enters). Sources, by testset:
#
#   • Hitzer & Sangwine, "Multivector and multivector matrix inverses in real Clifford algebras"
#     (2017): the scalar-of-conjugate-product structure (eqs 4.2, 5.3, 6.2, 7.2) and the universal
#     multivector inverse formulas for Cl(p,q), n ≤ 4 (eqs 4.3, 5.5, 6.10, 7.7). These invert
#     general (non-versor) multivectors that the package's own versor-only `inv` cannot.
#   • Shirokov, "On computing the determinant, other characteristic polynomial coefficients, and
#     inverse in Clifford algebras of arbitrary dimension" (2020): the cyclic scalar-part identities
#     ⟨UV⟩₀ = ⟨VU⟩₀ (any n) and ⟨UV⟩ₙ = ⟨VU⟩ₙ (odd n) of Lemma 1, and the involution
#     (anti)automorphism laws of eqs 8–9.
#   • Macfarlane, "Vector Analysis and Quaternions" (1906) and Lundholm & Svensson lecture notes —
#     the scalar triple product as a determinant, the BAC–CAB contraction identity, and reflections
#     as versor sandwiches (Cartan–Dieudonné). Theorem-oracle sources catalogued in the
#     ComputationalGeometricAlgebra.jl `wiki/bibliography/foundations.md`.

@testset "Theorem-based identities" begin
    using CliffordNumbers: dimension
    using LinearAlgebra: det

    # grade-k projection of any multivector, carried in the multivector's own scalar type
    gradepart(x, k) = KVector{k,signature(x),scalar_type(x)}(x)
    # Hitzer's grade-negation map m_{ks}: flip the sign of the listed grades, keep the rest. Built on
    # the dense representation so the grades the result spans are unconstrained.
    function negate_grades(x, ks)
        Q = signature(x); T = scalar_type(x)
        xd = CliffordNumber{Q,T}(x)
        for k in ks
            xd -= 2 * CliffordNumber{Q,T}(gradepart(xd, k))
        end
        return xd
    end

    randrat(rng, C) = C(ntuple(_ -> Rational{BigInt}(rand(rng, -3:3), rand(rng, 1:3)), nblades(C)))
    randf(rng, C) = C(ntuple(_ -> randn(rng), nblades(C)))

    # Non-degenerate signatures Cl(p,q) by dimension, for the dimension-specific Hitzer theorems.
    sig(p, q) = CliffordNumbers.Metrics.Signature(p + q, UInt(((1 << q) - 1) << p), UInt(0))
    dim1 = (VGA(1), sig(0, 1))                       # Cl(1,0), Cl(0,1)
    dim2 = (VGA(2), sig(1, 1), sig(0, 2))            # Cl(2,0), Cl(1,1), Cl(0,2)
    dim3 = (VGA(3), sig(2, 1), sig(1, 2), sig(0, 3)) # Cl(3,0), Cl(2,1), Cl(1,2), Cl(0,3)
    dim4 = (VGA(4), sig(3, 1), sig(2, 2), STA)       # Cl(4,0), Cl(3,1), Cl(2,2), Cl(1,3)

    @testset "Involutions are (anti)automorphisms" begin
        # Reversion and Clifford conjugation reverse products (anti-automorphisms); grade involution
        # preserves them (automorphism). Clifford conjugation is the composite of the other two, in
        # either order. (Hitzer eqs 2.4–2.6, Shirokov eqs 8–9.)
        rng = MersenneTwister(0x70F1)
        for Q in (VGA(3), STA, CGA(2), PGA(3)), _ in 1:4
            C = CliffordNumber{Q,Rational{BigInt}}
            x = randrat(rng, C); y = randrat(rng, C)
            @test grade_involution(x * y) == grade_involution(x) * grade_involution(y)
            @test reverse(x * y) == reverse(y) * reverse(x)
            @test conj(x * y) == conj(y) * conj(x)
            @test conj(x) == reverse(grade_involution(x))
            @test conj(x) == grade_involution(reverse(x))
            # Each is an involution: applying it twice is the identity.
            for f in (reverse, grade_involution, conj)
                @test f(f(x)) == x
            end
        end
    end

    @testset "Geometric product of vectors splits as a·b + a∧b" begin
        # The defining decomposition of the geometric product on grade-1 elements (Hestenes NFMP §1;
        # Dorst, Fontijne & Mann §6.1): the scalar part is the inner product, the bivector part is
        # the wedge, and the symmetric/antisymmetric halves recover each separately.
        rng = MersenneTwister(0x5D17)
        for _ in 1:25
            a = randrat(rng, KVector{1,VGA(3)}); b = randrat(rng, KVector{1,VGA(3)})
            ab = a * b
            @test scalar(ab) == scalar_product(a, b)
            @test gradepart(ab, 2) == a ∧ b
            @test (a * b + b * a) // 2 == scalar_product(a, b)   # symmetric half = inner product
            @test gradepart((a * b - b * a) // 2, 2) == a ∧ b    # antisymmetric half = wedge
        end
    end

    @testset "Macfarlane triple products" begin
        # Scalar triple product equals the 3×3 determinant of the component columns, read off the
        # σ₁₂₃ coefficient of a∧b∧c; and the BAC–CAB rule appears as the vector/bivector contraction
        # identity a⌋(b∧c) = (a·b)c − (a·c)b. (Macfarlane 1906, Ch. 7.)
        rng = MersenneTwister(0xBAC)
        for _ in 1:25
            a = randrat(rng, KVector{1,VGA(3)})
            b = randrat(rng, KVector{1,VGA(3)})
            c = randrat(rng, KVector{1,VGA(3)})
            triple = only(Tuple(gradepart(a ∧ b ∧ c, 3)))
            @test triple == det(hcat(collect.(Tuple.((a, b, c)))...))
            @test a ⨼ (b ∧ c) == scalar_product(a, b) * c - scalar_product(a, c) * b
        end
    end

    @testset "Reflections and rotor sandwiches" begin
        # A reflection of a vector v across the hyperplane ⊥ to a non-null vector u is the versor
        # sandwich −u v u⁻¹ (Cartan–Dieudonné; Lundholm & Svensson). It is again a vector, equals the
        # component formula v − 2(v·u)u/u², and is an isometry. A rotor sandwich (an even product of
        # two reflections) is likewise a length-preserving vector map.
        rng = MersenneTwister(0xEFEC)
        for _ in 1:25
            u = randrat(rng, KVector{1,VGA(3)}); v = randrat(rng, KVector{1,VGA(3)})
            iszero(abs2(u)) && continue
            vp = -u * v * inv(u)
            @test iszero(gradepart(vp, 3))                                  # stays grade 1
            @test gradepart(vp, 1) == v - 2 * scalar_product(v, u) * u / abs2(u)
            @test abs2(gradepart(vp, 1)) == abs2(v)                         # isometry
        end
        # Rotor sandwich R v R⁻¹ with R = exp(B): a unit transformation, so it preserves abs2 and
        # produces no trivector part. (Hathaway 1896; Dorst, Fontijne & Mann Ch. 7.)
        for _ in 1:25
            B = randf(rng, KVector{2,VGA(3)})
            R = exp(B)
            v = randf(rng, KVector{1,VGA(3)})
            vp = R * v * inv(R)
            @test isapprox(abs2(gradepart(vp, 1)), abs2(v); rtol = 1e-10)
            @test all(isapprox.(Tuple(gradepart(vp, 3)), 0; atol = 1e-10))
        end
    end

    @testset "Hitzer–Sangwine conjugate-product structure" begin
        # The scalar that discriminates divisors of zero from invertible elements is built by
        # multiplying x against grade-flipped copies of itself. Hitzer & Sangwine show which grades
        # survive at each dimension; everything else must cancel exactly.
        rng = MersenneTwister(0x4172)
        # n = 1: x·x̂ is a pure scalar (eq 4.2).
        for Q in dim1, _ in 1:4
            x = randrat(rng, CliffordNumber{Q,Rational{BigInt}})
            @test isscalar(x * grade_involution(x))
        end
        # n = 2: x·x̄ is a pure scalar and commutes, x·x̄ = x̄·x (eq 5.3).
        for Q in dim2, _ in 1:4
            x = randrat(rng, CliffordNumber{Q,Rational{BigInt}})
            @test isscalar(x * conj(x))
            @test x * conj(x) == conj(x) * x
        end
        # n = 3: x·x̄ is central (only grades 0 and 3), commutes, and (x·x̄)(x·x̄)~ is a pure scalar
        # (eqs 6.2, 6.4, 6.7).
        for Q in dim3, _ in 1:4
            x = randrat(rng, CliffordNumber{Q,Rational{BigInt}})
            p = x * conj(x)
            @test iszero(gradepart(p, 1)) && iszero(gradepart(p, 2))
            @test p == conj(x) * x
            @test isscalar(p * reverse(p))
        end
        # n = 4: x·x̄ has only grades 0, 3, 4 (eq 7.2).
        for Q in dim4, _ in 1:4
            x = randrat(rng, CliffordNumber{Q,Rational{BigInt}})
            p = x * conj(x)
            @test iszero(gradepart(p, 1)) && iszero(gradepart(p, 2))
        end
    end

    @testset "Hitzer–Sangwine universal multivector inverse" begin
        # Closed-form inverses for a general element of Cl(p,q), n ≤ 4 — including the non-versor
        # multivectors the package's own `inv` rejects. Each builds a scalar denominator from the
        # conjugate-product chain above; where it is non-zero, the formula is a true two-sided
        # inverse. Verified exactly over the rationals.
        rng = MersenneTwister(0x10F5)
        # n = 1: x⁻¹ = x̂ / scalar(x·x̂)   (eq 4.3)
        for Q in dim1, _ in 1:8
            x = randrat(rng, CliffordNumber{Q,Rational{BigInt}})
            xh = grade_involution(x); s = scalar(x * xh)
            iszero(s) && continue
            xinv = xh / s
            @test x * xinv == 1 && xinv * x == 1
        end
        # n = 2: x⁻¹ = x̄ / scalar(x·x̄)   (eq 5.5)
        for Q in dim2, _ in 1:8
            x = randrat(rng, CliffordNumber{Q,Rational{BigInt}})
            xb = conj(x); s = scalar(x * xb)
            iszero(s) && continue
            xinv = xb / s
            @test x * xinv == 1 && xinv * x == 1
        end
        # n = 3: x⁻¹ = x̄(x·x̄)~ / scalar(x·x̄·(x·x̄)~)   (Theorem 6.1, eq 6.10)
        for Q in dim3, _ in 1:8
            x = randrat(rng, CliffordNumber{Q,Rational{BigInt}})
            xb = conj(x); p = x * xb; pr = reverse(p); s = scalar(p * pr)
            iszero(s) && continue
            xinv = xb * pr / s
            @test x * xinv == 1 && xinv * x == 1
        end
        # n = 4: x⁻¹ = x̄·m₃₄(x·x̄) / scalar(x·x̄·m₃₄(x·x̄)), m₃₄ negating grades 3 and 4   (eq 7.7)
        for Q in dim4, _ in 1:8
            x = randrat(rng, CliffordNumber{Q,Rational{BigInt}})
            xb = conj(x); p = x * xb; m = negate_grades(p, (3, 4)); s = scalar(p * m)
            iszero(s) && continue
            xinv = xb * m / s
            @test x * xinv == 1 && xinv * x == 1
        end
    end

    @testset "Commutator and anticommutator products" begin
        using CliffordNumbers: commutator, anticommutator
        # The geometric product splits into its symmetric and antisymmetric halves,
        # xy = ⨰(x,y) + ×(x,y), with the anticommutator symmetric and the commutator antisymmetric;
        # the operators `⨰` and `×` are the named functions. The commutator of two bivectors is
        # again a bivector — the grade-2 subspace is closed under ½[·,·], the Lie algebra of the
        # spin group. (Doran, Hestenes, Sommen & Van Acker, "Lie Groups as Spin Groups", 1993;
        # Doran & Lasenby Ch. 11.)
        rng = MersenneTwister(0xC0AC)
        for Q in (VGA(3), STA, PGA(3)), _ in 1:5
            C = CliffordNumber{Q,Rational{BigInt}}
            x = randrat(rng, C); y = randrat(rng, C)
            @test x * y == anticommutator(x, y) + commutator(x, y)
            # `×` is qualified to avoid the LinearAlgebra.cross ambiguity in the test environment.
            @test commutator(x, y) == CliffordNumbers.:×(x, y)
            @test anticommutator(x, y) == x ⨰ y
            @test commutator(x, y) == -commutator(y, x)
            @test anticommutator(x, y) == anticommutator(y, x)
            B1 = randrat(rng, KVector{2,Q}); B2 = randrat(rng, KVector{2,Q})
            cb = commutator(B1, B2)
            @test gradepart(cb, 2) == cb        # bivector commutator stays a bivector
        end
    end

    @testset "Shirokov cyclic scalar/pseudoscalar parts" begin
        # The scalar part of a product is invariant under cyclic permutation, ⟨UV⟩₀ = ⟨VU⟩₀ for any
        # n; the pseudoscalar part shares the symmetry in odd dimension, ⟨UV⟩ₙ = ⟨VU⟩ₙ. (Shirokov
        # Lemma 1, eqs 3–4.) The scalar identity is the trace's cyclic invariance, since
        # tr(x) = 2ⁿ·scalar(x).
        rng = MersenneTwister(0x5417)
        for Q in (VGA(3), STA, PGA(3), CGA(2)), _ in 1:5
            C = CliffordNumber{Q,Rational{BigInt}}
            x = randrat(rng, C); y = randrat(rng, C)
            @test scalar(x * y) == scalar(y * x)
            n = dimension(Q)
            isodd(n) && @test gradepart(x * y, n) == gradepart(y * x, n)
        end
    end
end
