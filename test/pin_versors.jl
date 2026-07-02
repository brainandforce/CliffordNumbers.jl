@testset "Pin-group (odd versor) structure" begin
    using CliffordNumbers: dimension, versor_inverse

    function randrat(rng, ::Type{C}) where C
        return C(ntuple(_ -> Rational{BigInt}(rand(rng, -3:3), rand(rng, 1:3)), nblades(C)))
    end
    sig(p, q) = CliffordNumbers.Metrics.Signature(p + q, UInt(((1 << q) - 1) << p), UInt(0))

    # A rational versor of k mirrors, all with nonzero square
    function randversor(rng, Q, k)
        mirrors = KVector{1,Q,Rational{BigInt}}[]
        while length(mirrors) < k
            v = randrat(rng, KVector{1,Q,Rational{BigInt}})
            iszero(scalar_product(v, v)) || push!(mirrors, v)
        end
        return prod(mirrors)
    end

    @testset "parity-correct one" begin
        # The scalar identity of an odd carrier lives in the even subalgebra: odd·odd = even
        for Q in (VGA(3), VGA(4), STA, PGA(3))
            @test one(OddCliffordNumber{Q,Float64}) === one(EvenCliffordNumber{Q,Float64})
            @test one(OddCliffordNumber{Q}) isa EvenCliffordNumber{Q,Bool}
            x = OddCliffordNumber{Q}(ntuple(_ -> 1, div(CliffordNumbers.blade_count(Q), 2)))
            @test one(x) * x == x
            @test x * one(x) == x
            @test x^0 === one(x)
        end
        # KVector and even carriers keep their identities
        @test one(KVector{1,VGA(3),Int}) === KVector{0,VGA(3)}(1)
        @test one(EvenCliffordNumber{VGA(3),Int}) === EvenCliffordNumber{VGA(3)}(1, 0, 0, 0)
        # oneunit still throws: an odd carrier cannot hold the scalar identity
        @test_throws InexactError oneunit(OddCliffordNumber{VGA(3),Int})
    end

    @testset "mixed-parity products are type-stable" begin
        rng = MersenneTwister(0x0DD5)
        v = randversor(rng, VGA(3), 1)
        R = randversor(rng, VGA(3), 2)
        @test (@inferred prod((v, R))) isa OddCliffordNumber{VGA(3)}
        @test (@inferred prod((v, R, v))) isa EvenCliffordNumber{VGA(3)}
        @test (@inferred prod((R, v, R))) == R * v * R
    end

    @testset "mirror_peel" begin
        rng = MersenneTwister(0x9EE1)
        for Q in (VGA(2), VGA(3), VGA(4), STA, sig(2, 2)), k in (1, 3), _ in 1:8
            k > dimension(Q) && continue
            V = randversor(rng, Q, k)
            (v, R) = mirror_peel(V)
            @test v isa KVector{1,Q}
            @test R isa EvenCliffordNumber{Q}
            @test !iszero(scalar_product(v, v))
            @test v * R == V                            # exact over the rationals
            @test isversor(R)
        end
        # A PGA reflection composed with a motor peels to an invertible mirror
        p = KVector{1,PGA(3)}(0.0, 1.0, 0.0, 0.0)
        M = exp(KVector{2,PGA(3)}(0.0, 0.3, 0.0, 0.5, -0.25, 0.75))
        V = p * M
        (v, R) = mirror_peel(V)
        @test v * R ≈ OddCliffordNumber{PGA(3),Float64}(V)
        @test !isapprox(scalar_product(v, v), 0; atol = 1e-10)
        # A pure trivector has no grade-1 part; the peel falls back to a support basis vector
        V3 = KVector{3,VGA(4)}(2//1, 0, 0, 0)
        (v, R) = mirror_peel(V3)
        @test v * R == V3
        @test v == KVector{1,VGA(4)}(1//1, 0, 0, 0)
        # Even and mixed-parity inputs are rejected
        @test_throws DomainError mirror_peel(randversor(rng, VGA(3), 2))
        @test_throws DomainError mirror_peel(CliffordNumber{VGA(3)}(1, 1, 0, 0, 0, 0, 0, 0))
    end

    @testset "sandwich" begin
        rng = MersenneTwister(0x5A4D)
        # A single mirror acts as the proper reflection -u v u⁻¹, exactly over the rationals
        for Q in (VGA(3), STA, sig(2, 2)), _ in 1:8
            u = randversor(rng, Q, 1)
            v = randrat(rng, KVector{1,Q,Rational{BigInt}})
            @test sandwich(u, v) == -(u * v * versor_inverse(u))
            @test sandwich(u, sandwich(u, v)) == v      # reflections are involutions
        end
        # Odd versors preserve grade and the scalar product on every operand grade (exact)
        for Q in (VGA(3), VGA(4), STA), k in 1:3, _ in 1:4
            V = randversor(rng, Q, 3)
            x = randrat(rng, KVector{k,Q,Rational{BigInt}})
            y = randrat(rng, KVector{k,Q,Rational{BigInt}})
            xs = sandwich(V, x)
            @test CliffordNumber{Q}(KVector{k,Q,Rational{BigInt}}(xs)) == CliffordNumber{Q}(xs)
            @test scalar_product(xs, sandwich(V, y)) == scalar_product(x, y)
        end
        # The action is a homomorphism through the mirror-peel normal form
        for Q in (VGA(3), STA), _ in 1:4
            V = randversor(rng, Q, 3)
            (v, R) = mirror_peel(V)
            x = randrat(rng, KVector{2,Q,Rational{BigInt}})
            @test sandwich(V, x) == sandwich(v, sandwich(R, x))
        end
        # Runtime parity for dense carriers matches the typed dispatch
        V = randversor(rng, VGA(3), 3)
        x = randrat(rng, KVector{1,VGA(3),Rational{BigInt}})
        @test sandwich(CliffordNumber{VGA(3)}(V), x) == CliffordNumber{VGA(3)}(sandwich(V, x))
        # Type stability and allocation-freedom on the rotor path
        R = exp(KVector{2,VGA(3)}(0.3, 0.2, -0.1))
        w = KVector{1,VGA(3)}(1.0, 2.0, 3.0)
        sandwich(R, w)
        @test (@inferred sandwich(R, w)) isa OddCliffordNumber{VGA(3),Float64}
        ALLOCATION_GATES && @test (@allocated sandwich(R, w)) == 0
        # Mismatched algebras are rejected
        @test_throws CliffordNumbers.AlgebraMismatch sandwich(R, KVector{1,STA}(1.0, 0, 0, 0))
    end
end
