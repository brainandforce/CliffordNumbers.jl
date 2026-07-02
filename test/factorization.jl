@testset "Blade and versor structure" begin
    using CliffordNumbers: dimension

    function randrat(rng, ::Type{C}) where C
        return C(ntuple(_ -> Rational{Int}(rand(rng, -3:3), rand(rng, 1:3)), nblades(C)))
    end
    randf(rng, ::Type{C}) where C = C(ntuple(_ -> randn(rng), nblades(C)))
    # Non-degenerate signatures Cl(p,q) beyond the aliased ones
    sig(p, q) = CliffordNumbers.Metrics.Signature(p + q, UInt(((1 << q) - 1) << p), UInt(0))

    @testset "isblade" begin
        rng = MersenneTwister(0xB1AD)
        # Wedge products of random vectors are blades, exactly over the rationals
        for Q in (VGA(3), VGA(4), PGA(3), STA, CGA(2)), k in 1:min(dimension(Q), 4), _ in 1:4
            B = mapreduce(_ -> randrat(rng, KVector{1,Q,Rational{Int}}), ∧, 1:k)
            iszero(B) || @test isblade(B)
        end
        # ... and with floating-point noise from the products
        for Q in (VGA(4), PGA(3), STA), _ in 1:4
            B = randf(rng, KVector{1,Q,Float64}) ∧ randf(rng, KVector{1,Q,Float64})
            @test isblade(B)
        end
        # The canonical non-blade bivector, and its embeddings
        nonblade = KVector{2,VGA(4)}(1, 0, 0, 0, 0, 1)      # e₁₂ + e₃₄
        @test !isblade(nonblade)
        @test !isblade(EvenCliffordNumber{VGA(4)}(nonblade))
        @test !isblade(float(nonblade))
        @test !isblade(KVector{2,PGA(3)}(1, 0, 0, 0, 0, 1)) # degenerate metric, same relations
        # Scalars, vectors, pseudovectors, pseudoscalars, and zero are always blades
        @test isblade(KVector{0,VGA(3)}(3))
        @test isblade(randrat(rng, KVector{3,VGA(4),Rational{Int}}))
        @test isblade(pseudoscalar(CliffordNumber{VGA(4),Int}))
        @test isblade(zero(KVector{2,VGA(4),Int}))
        # Inhomogeneous multivectors are not blades; single-grade dense content is
        @test !isblade(CliffordNumber{VGA(3)}(1, 1, 0, 0, 0, 0, 0, 0))
        @test isblade(CliffordNumber{VGA(3)}(KVector{1,VGA(3)}(1, 2, 3)))
    end

    @testset "isversor" begin
        rng = MersenneTwister(0x7E45)
        # Products of invertible rational vectors are versors, exactly
        for Q in (VGA(3), VGA(4), STA, sig(2, 2)), k in 1:dimension(Q), _ in 1:4
            mirrors = KVector{1,Q,Rational{Int}}[]
            while length(mirrors) < k
                v = randrat(rng, KVector{1,Q,Rational{Int}})
                iszero(scalar_product(v, v)) || push!(mirrors, v)
            end
            @test isversor(prod(mirrors))
        end
        # PGA motors (products of non-null mirrors) qualify; a null mirror does not
        t = exp(KVector{2,PGA(3)}(0.0, 0.0, 0.0, 0.5, -0.25, 0.75))
        @test isversor(t)
        @test !isversor(KVector{1,PGA(3)}(1, 0, 0, 0))      # e₀ is null
        # Nonzero scalars are versors (empty products up to scale), zero is not
        @test isversor(KVector{0,VGA(3)}(3//2))
        @test !isversor(zero(EvenCliffordNumber{VGA(3),Int}))
        # Mixed parity, non-scalar x x̃, and isoclinic rotors
        @test !isversor(CliffordNumber{VGA(3)}(1, 1, 0, 0, 0, 0, 0, 0))
        @test !isversor(KVector{2,VGA(4)}(1, 0, 0, 0, 0, 1))
        @test isversor(exp(KVector{2,VGA(4)}(0.7, 0, 0, 0, 0, 0.7)))
    end

    @testset "factor_blade" begin
        rng = MersenneTwister(0xFB1A)
        # Rational reassembly is exact; float reassembly is approximate
        for Q in (VGA(3), VGA(4), PGA(3), STA), k in 2:min(dimension(Q), 4), _ in 1:4
            B = mapreduce(_ -> randrat(rng, KVector{1,Q,Rational{Int}}), ∧, 1:k)
            iszero(B) && continue
            fs = factor_blade(B)
            @test length(fs) == k
            @test all(f -> f isa KVector{1,Q}, fs)
            @test reduce(∧, fs) == B
        end
        for Q in (VGA(4), CGA(2)), _ in 1:4
            B = randf(rng, KVector{1,Q,Float64}) ∧ randf(rng, KVector{1,Q,Float64})
            @test reduce(∧, factor_blade(B)) ≈ B
        end
        # Degenerate directions factor exactly too
        B0 = KVector{1,PGA(3)}(1//1, 0, 0, 0) ∧ KVector{1,PGA(3)}(1//1, 2//1, 0, 3//1)
        @test reduce(∧, factor_blade(B0)) == B0
        # Error paths
        @test_throws DomainError factor_blade(KVector{0,VGA(3)}(2))
        @test_throws DomainError factor_blade(zero(KVector{2,VGA(3),Int}))
        @test_throws DomainError factor_blade(KVector{2,VGA(4)}(1, 0, 0, 0, 0, 1))
    end

    @testset "factor_versor" begin
        rng = MersenneTwister(0xFACE)
        # The peel is high-degree in the coefficients, so the exact corpus needs BigInt numerators
        for Q in (VGA(3), VGA(4), STA), k in 1:dimension(Q), _ in 1:4
            mirrors = KVector{1,Q,Rational{BigInt}}[]
            while length(mirrors) < k
                v = randrat(rng, KVector{1,Q,Rational{BigInt}})
                iszero(scalar_product(v, v)) || push!(mirrors, v)
            end
            V = prod(mirrors)
            fs = factor_versor(V)
            @test all(f -> f isa KVector{1,Q}, fs)
            @test length(fs) <= dimension(Q)
            @test prod(fs) == V                     # exact over the rationals
            @test iseven(length(fs)) == iseven(k)   # parity of the mirror count is preserved
        end
        # PGA translations peel as a parallel mirror pair; general motors reassemble
        for _ in 1:4
            M = exp(KVector{2,PGA(3)}(ntuple(_ -> 0.5 * randn(rng), 6)))
            fs = factor_versor(M)
            @test length(fs) <= dimension(PGA(3))
            @test prod(fs) ≈ EvenCliffordNumber{PGA(3),Float64}(M) atol = 1e-8
        end
        # STA boosts and rotations
        for _ in 1:4
            V = exp(KVector{2,STA}(ntuple(_ -> 0.4 * randn(rng), 6)))
            @test prod(factor_versor(V)) ≈ V atol = 1e-8
        end
        # Scalars factor as a parallel pair
        fs = factor_versor(KVector{0,VGA(3)}(3//2))
        @test length(fs) == 2 && prod(fs) == 3//2
        @test_throws DomainError factor_versor(CliffordNumber{VGA(3)}(1, 1, 0, 0, 0, 0, 0, 0))
    end

    @testset "bivector_decomposition" begin
        rng = MersenneTwister(0xB1FF)
        for Q in (VGA(4), PGA(3), STA, CGA(2), CGA(3)), _ in 1:25
            B = KVector{2,Q}(ntuple(_ -> randn(rng), binomial(dimension(Q), 2)))
            (b1, b2) = bivector_decomposition(B)
            @test b1 + b2 ≈ B
            @test b1 * b2 ≈ b2 * b1 atol = 1e-8 * sum(abs2, Tuple(B))
            @test isblade(b1)
            @test isblade(b2)
        end
        # Simple bivectors return themselves and zero; below four dimensions always
        (b1, b2) = bivector_decomposition(KVector{2,VGA(3)}(1.0, 2.0, 3.0))
        @test iszero(b2) && b1 ≈ KVector{2,VGA(3)}(1.0, 2.0, 3.0)
        simple = KVector{1,VGA(4)}(1, 1, 0, 0) ∧ KVector{1,VGA(4)}(0, 1, 1, 0)
        (b1, b2) = bivector_decomposition(simple)
        @test iszero(b2)
        # Isoclinic bivectors (repeated plane squares) still split into commuting simple planes
        iso = KVector{2,VGA(4)}(0.7, 0, 0, 0, 0, 0.7)
        (b1, b2) = bivector_decomposition(iso)
        @test b1 + b2 ≈ iso
        @test b1 * b2 ≈ b2 * b1
        @test isblade(b1) && isblade(b2)
        @test !iszero(b1) && !iszero(b2)
        # PGA motor generators split into the Euclidean plane and the ideal (null) plane
        (b1, b2) = bivector_decomposition(KVector{2,PGA(3)}(0.0, 0.0, 0.3, 0.0, 0.0, 1.2))
        @test b1 + b2 ≈ KVector{2,PGA(3)}(0.0, 0.0, 0.3, 0.0, 0.0, 1.2)
        @test isblade(b1) && isblade(b2)
        # Type stability and allocation-freedom
        B = KVector{2,VGA(4)}(ntuple(_ -> randn(rng), 6))
        bivector_decomposition(B)
        @test (@inferred bivector_decomposition(B)) isa NTuple{2,KVector{2,VGA(4),Float64,6}}
        ALLOCATION_GATES && @test (@allocated bivector_decomposition(B)) == 0
    end

    @testset "general log and sqrt" begin
        rng = MersenneTwister(0x106E)
        # Study numbers s + n with n² a scalar: hyperbolic (s > |n|), elliptic, and parabolic
        for _ in 1:10
            v = randf(rng, KVector{1,VGA(3),Float64})
            x = CliffordNumber{VGA(3)}(abs(v) + rand(rng) + 0.1) + v
            @test exp(log(x)) ≈ x rtol = 1e-6
            @test sqrt(x)^2 ≈ x rtol = 1e-6
        end
        for _ in 1:10
            v = randf(rng, KVector{1,sig(0, 3),Float64})        # v² < 0: elliptic, no scalar needed
            x = CliffordNumber{sig(0, 3)}(randn(rng)) + v
            @test exp(log(x)) ≈ x rtol = 1e-6
            @test sqrt(x)^2 ≈ x rtol = 1e-6
        end
        x = CliffordNumber{PGA(3)}(2.0) + KVector{1,PGA(3)}(0.5, 0, 0, 0)   # parabolic: e₀² = 0
        @test exp(log(x)) ≈ x
        @test sqrt(x)^2 ≈ x
        # Odd blades with negative square have elliptic logs and square roots
        I3 = KVector{3,VGA(3)}(2.0)
        @test sqrt(I3)^2 ≈ I3
        @test exp(log(I3)) ≈ CliffordNumber{VGA(3),Float64}(I3)
        # Scaled rotors and motors route through the rotor logarithm
        for Q in (VGA(3), VGA(4), PGA(3)), _ in 1:5
            R = exp(KVector{2,Q}(ntuple(_ -> 0.4 * randn(rng), binomial(dimension(Q), 2))))
            x = (1 + rand(rng)) * CliffordNumber{Q,Float64}(R)
            @test exp(log(x)) ≈ x rtol = 1e-6
            @test sqrt(x)^2 ≈ x rtol = 1e-6
        end
        # Results are returned as the exponential type; scalars keep KVector{0}
        @test log(KVector{1,sig(0, 3)}(1.0, 2.0, 0.0)) isa CliffordNumber{sig(0, 3),Float64}
        @test sqrt(KVector{3,VGA(3)}(2.0)) isa CliffordNumber{VGA(3),Float64}
        @test log(KVector{0,VGA(3)}(2.0)) === KVector{0,VGA(3)}(log(2.0))
        @test_throws DomainError log(KVector{0,VGA(3)}(-1.0))
        # No principal branch: hyperbolic without dominant scalar, negative parabolic, odd versors
        @test_throws DomainError log(KVector{1,VGA(3)}(1.0, 0.0, 0.0))
        @test_throws DomainError sqrt(CliffordNumber{VGA(2)}(-1.0, 0.0, 0.0, 0.0))
        @test_throws DomainError log(CliffordNumber{VGA(3)}(KVector{1,VGA(3)}(1.0, 2.0, 0.5)) +
                                     CliffordNumber{VGA(3)}(KVector{3,VGA(3)}(0.3)))
    end
end
