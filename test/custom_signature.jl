# Phase-space signatures (p, q, r) = (3, 3, 0) and (3, 3, 1): 3 position dimensions squaring to +1,
# 3 momentum dimensions squaring to −1, and r degenerate dimensions. (3,3,0) is the neutral
# phase-space metric; (3,3,1) is its projective extension. A downstream package defines such a
# metric as a `Signature` *value* — the supported, world-age-safe extension path (the family's own
# `STAP` is built the same way). The generated kernels read the metric from the value's `isbits`
# fields, so no downstream method dispatch happens inside a generator.
#
# `det` at n ≥ 6 reduces the explicit operator matrix and needs no geometric product, so it stays
# fast here. The *raw* geometric product at n = 7 compiles a very large generated kernel, so the one
# test that needs it is gated behind CLIFFORDNUMBERS_FULL_TESTS.

# A would-be custom AbstractSignature *subtype*. Its interface methods are correct, but they are
# defined here, downstream of CliffordNumbers — see the "known limitation" testset below.
struct PhaseSpace{R} <: CliffordNumbers.Metrics.AbstractSignature end
CliffordNumbers.Metrics.dimension(::PhaseSpace{R}) where R = 6 + R
Base.firstindex(::PhaseSpace) = 1
Base.@propagate_inbounds function Base.getindex(s::PhaseSpace, i::Int)
    @boundscheck checkbounds(s, i)
    return i <= 3 ? Int8(1) : (i <= 6 ? Int8(-1) : Int8(0))
end
CliffordNumbers.Metrics.is_degenerate(::PhaseSpace{R}) where R = R > 0
CliffordNumbers.Metrics.is_positive_definite(::PhaseSpace) = false

@testset "Downstream-defined signatures (phase space)" begin
    using CliffordNumbers: blade_signs, BladeIndices, blade_count, nblades, dimension

    full = get(ENV, "CLIFFORDNUMBERS_FULL_TESTS", "false") == "true"
    # negative mask 0x38 = bits for dims 4,5,6 (momenta); degenerate mask 0x40 = bit for dim 7
    PS330 = Signature(6, 0x38, 0x00)   # (3, 3, 0)
    PS331 = Signature(7, 0x38, 0x40)   # (3, 3, 1)

    refauto(f, w) = (T = typeof(w); T(w[f.(BladeIndices(T))]))

    @testset "signature interface" begin
        @test dimension(PS330) == 6 && blade_count(PS330) == 64
        @test dimension(PS331) == 7 && blade_count(PS331) == 128
        @test PS330[1] === Int8(1) && PS330[4] === Int8(-1)
        @test PS331[7] === Int8(0)
        @test !is_degenerate(PS330) && is_degenerate(PS331)
        @test count(==(Int8(1)), PS330) == 3 && count(==(Int8(-1)), PS330) == 3
    end

    # n = 7 automorphisms and blade_signs go through the (product-free) sign path, so this verifies
    # the (3,3,1) signature drives the generated blade machinery without paying the 7D mul compile.
    @testset "(3,3,1) signature machinery" begin
        C = CliffordNumber{PS331,Float64}
        x = C(ntuple(i -> Float64(i) - 64.5, nblades(C)))
        for f in (reverse, adjoint, conj, grade_involution)
            @test f(x) === refauto(f, x)
        end
        @test reverse(reverse(x)) === x
        @test blade_signs(C, Val(:reverse)) === map(b -> sign(reverse(b)), Tuple(BladeIndices(C)))
    end

    @testset "determinant (matrix path, product-free)" begin
        # det(λ) = λ^(2ⁿ), exact, on both phase spaces (n = 6 and n = 7)
        λ0 = Rational{BigInt}(3, 2)
        @test CliffordNumbers.det(CliffordNumber{PS330,Rational{BigInt}}(λ0)) == λ0^64
        λ1 = Rational{BigInt}(-2, 3)
        @test CliffordNumbers.det(CliffordNumber{PS331,Rational{BigInt}}(λ1)) == λ1^128
        # 1 + e₁ is a zero divisor (e₁² = +1), hence singular — built by addition, no product needed
        for Q in (PS330, PS331)
            D = CliffordNumber{Q,Rational{BigInt}}
            s = one(D) + KVector{1,Q,Rational{BigInt}}(ntuple(i -> Rational{BigInt}(i == 1), dimension(Q)))
            @test CliffordNumbers.det(s) == 0
        end
    end

    @testset "(3,3,0) products, automorphisms, inverse, det" begin
        C = CliffordNumber{PS330,Float64}
        rng = MersenneTwister(0xda7a)
        x = C(ntuple(_ -> randn(rng), nblades(C)))
        y = C(ntuple(_ -> randn(rng), nblades(C)))
        z = C(ntuple(_ -> randn(rng), nblades(C)))
        @test x * y isa CliffordNumber{PS330}
        @test (x * y) * z ≈ x * (y * z)        # associativity
        for f in (reverse, adjoint, conj, grade_involution)
            @test f(x) === refauto(f, x)
        end
        e1 = KVector{1,PS330,Float64}(ntuple(i -> Float64(i == 1), 6))
        @test inv(e1) * e1 ≈ one(C)            # e₁² = +1
        # Multiplicativity over ℚ (the product is the only n = 6 generated-kernel use here)
        D = CliffordNumber{PS330,Rational{BigInt}}
        rr = MersenneTwister(0x33a)
        xr = D(ntuple(_ -> Rational{BigInt}(rand(rr, -1:1)), nblades(D)))
        yr = D(ntuple(_ -> Rational{BigInt}(rand(rr, -1:1)), nblades(D)))
        @test CliffordNumbers.det(xr * yr) == CliffordNumbers.det(xr) * CliffordNumbers.det(yr)
    end

    if full
        # The only test that compiles the 7D geometric product kernel.
        @testset "(3,3,1) geometric product and determinant multiplicativity" begin
            D = CliffordNumber{PS331,Rational{BigInt}}
            rng = MersenneTwister(0xda7b)
            x = D(ntuple(_ -> Rational{BigInt}(rand(rng, -1:1)), nblades(D)))
            y = D(ntuple(_ -> Rational{BigInt}(rand(rng, -1:1)), nblades(D)))
            @test x * y isa CliffordNumber{PS331}
            @test CliffordNumbers.det(x * y) == CliffordNumbers.det(x) * CliffordNumbers.det(y)
        end
    end

    @testset "subtype signatures need a world-age-robust trait (known limitation)" begin
        # The subtype's interface is correct at normal world age...
        @test dimension(PhaseSpace{0}()) == 6
        @test PhaseSpace{0}()[4] === Int8(-1)
        # ...but it cannot yet drive the @generated kernels: their generators run in an older world
        # age and cannot see `dimension`/`getindex` defined downstream of CliffordNumbers, so a
        # geometric product falls back to the generic `signed(s.dimensions)` and throws. This is the
        # gap a world-age-robust signature trait interface would close.
        Q = PhaseSpace{0}()
        x = CliffordNumber{Q,Float64}(ntuple(_ -> 1.0, 64))
        @test_throws Exception x * x
    end
end
