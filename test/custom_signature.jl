# Phase-space signatures (p, q, r) = (D, D, R): D position dimensions squaring to +1, D momentum
# dimensions squaring to −1, and R degenerate dimensions. (D, D, 0) is the neutral phase-space metric;
# (D, D, 1) is its projective extension. A downstream package defines such a metric as a `Signature`
# *value* — the supported, world-age-safe extension path (the family's own `STAP` is built the same
# way). The generated kernels read the metric from the value's `isbits` fields, so no downstream
# method dispatch happens inside a generator. The tests below exercise both D = 2 (small enough to
# compile the geometric product) and D = 3 (where the 7D product is gated behind the full-test flag).

# A would-be custom AbstractSignature *subtype*, parameterized by the position/momentum count D and
# the degenerate count R. Its interface methods are correct, but they are defined here, downstream of
# CliffordNumbers — see the "world-age limitation" testset below.
struct PhaseSpace{D,R} <: CliffordNumbers.Metrics.AbstractSignature end
CliffordNumbers.Metrics.dimension(::PhaseSpace{D,R}) where {D,R} = 2D + R
Base.firstindex(::PhaseSpace) = 1
Base.@propagate_inbounds function Base.getindex(s::PhaseSpace{D}, i::Int) where D
    @boundscheck checkbounds(s, i)
    return i <= D ? Int8(1) : (i <= 2D ? Int8(-1) : Int8(0))
end
CliffordNumbers.Metrics.is_degenerate(::PhaseSpace{D,R}) where {D,R} = R > 0
CliffordNumbers.Metrics.is_positive_definite(::PhaseSpace) = false

@testset "Downstream-defined signatures (phase space)" begin
    using CliffordNumbers: BladeIndices, blade_count, nblades, dimension

    full = get(ENV, "CLIFFORDNUMBERS_FULL_TESTS", "false") == "true"
    # negative mask 0x38 = bits for dims 4,5,6 (momenta); degenerate mask 0x40 = bit for dim 7
    PS330 = Signature(6, 0x38, 0x00)   # (3, 3, 0)
    PS331 = Signature(7, 0x38, 0x40)   # (3, 3, 1)
    # 2D phase space: negative mask 0x0c = bits for dims 3,4 (momenta); degenerate mask 0x10 = dim 5
    PS220 = Signature(4, 0x0c, 0x00)   # (2, 2, 0)
    PS221 = Signature(5, 0x0c, 0x10)   # (2, 2, 1)

    refauto(f, w) = (T = typeof(w); T(w[f.(BladeIndices(T))]))

    @testset "signature interface" begin
        @test dimension(PS330) == 6 && blade_count(PS330) == 64
        @test dimension(PS331) == 7 && blade_count(PS331) == 128
        @test PS330[1] === Int8(1) && PS330[4] === Int8(-1)
        @test PS331[7] === Int8(0)
        @test !is_degenerate(PS330) && is_degenerate(PS331)
        @test count(==(Int8(1)), PS330) == 3 && count(==(Int8(-1)), PS330) == 3
    end

    # nblades of a KVector in a plain-Signature algebra must stay an Int, not a Float64: on Julia 1.10
    # a narrow `dimension` fed `binomial`'s Float64 gamma-function method and poisoned construction.
    @testset "nblades stays integer (Julia 1.10 binomial regression)" begin
        @test nblades(KVector{2,PS330,Float64}) === 15
        @test nblades(KVector{3,PS331,Float64}) === 35
        @test dimension(PS330) isa Int
    end

    # n = 7 automorphisms go through the (product-free) sign path, so this verifies the (3,3,1)
    # signature drives the generated blade machinery without paying the 7D mul compile.
    @testset "(3,3,1) signature machinery" begin
        C = CliffordNumber{PS331,Float64}
        x = C(ntuple(i -> Float64(i) - 64.5, nblades(C)))
        for f in (reverse, adjoint, conj, grade_involution)
            @test f(x) === refauto(f, x)
        end
        @test reverse(reverse(x)) === x
    end

    @testset "(3,3,0) products, automorphisms, inverse" begin
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
    end

    if full
        # The only test that compiles the 7D geometric product kernel.
        @testset "(3,3,1) geometric product" begin
            D = CliffordNumber{PS331,Rational{BigInt}}
            rng = MersenneTwister(0xda7b)
            x = D(ntuple(_ -> Rational{BigInt}(rand(rng, -1:1)), nblades(D)))
            y = D(ntuple(_ -> Rational{BigInt}(rand(rng, -1:1)), nblades(D)))
            @test x * y isa CliffordNumber{PS331}
        end
    end

    @testset "metric_tuple is the single source of metric data" begin
        # metric_tuple agrees with element-by-element indexing, and the masks derive from it.
        for Q in (PS220, PS221, PS330, PS331, STA, PGA(3), CGA(2), Signature(5, 0b11100, 0b00001, -1))
            m = metric_tuple(Q)
            @test m === ntuple(i -> Q[firstindex(Q) - 1 + i], dimension(Q))
            @test m isa NTuple{<:Any,Int8}
            neg = reduce(|, (UInt(m[i] === Int8(-1)) << (i - 1) for i in eachindex(m)); init = UInt(0))
            zer = reduce(|, (UInt(m[i] === Int8(0)) << (i - 1) for i in eachindex(m)); init = UInt(0))
            @test CliffordNumbers.negative_square_bits(Q) == neg
            @test CliffordNumbers.zero_square_bits(Q) == zer
        end
        # Downstream subtype: metric_tuple still reads its (run-time) interface correctly, for any D/R.
        @test metric_tuple(PhaseSpace{2,0}()) === (Int8(1), Int8(1), Int8(-1), Int8(-1))
        @test metric_tuple(PhaseSpace{3,0}()) === ntuple(i -> i <= 3 ? Int8(1) : Int8(-1), 6)
        @test metric_tuple(PhaseSpace{3,1}()) === (Int8(1), Int8(1), Int8(1), Int8(-1), Int8(-1), Int8(-1), Int8(0))
        # Allocation-free / constant-folded on a built-in.
        mt() = metric_tuple(STA)
        mt()
        ALLOCATION_GATES && @test @allocated(mt()) == 0
    end

    @testset "subtype signatures: world-age limitation and the Signature(s) bridge" begin
        # The subtype's interface is correct at normal world age, for any D/R...
        @test dimension(PhaseSpace{2,0}()) == 4
        @test dimension(PhaseSpace{3,0}()) == 6
        @test PhaseSpace{3,0}()[4] === Int8(-1)
        # ...but a bare subtype value cannot drive the @generated kernels: their generators run in an
        # older world age and cannot see `dimension`/`getindex` defined downstream of CliffordNumbers,
        # so a geometric product falls back to the generic `signed(s.dimensions)` and throws.
        Q = PhaseSpace{3,0}()
        x = CliffordNumber{Q,Float64}(ntuple(_ -> 1.0, 64))
        @test_throws Exception x * x

        # The supported bridge: reduce the subtype to a Signature *value* (computed at normal world
        # age, so it can read the subtype's interface), then parameterize Clifford numbers on it.
        @test Signature(PhaseSpace{2,0}()) === PS220
        @test Signature(PhaseSpace{2,1}()) === PS221
        @test Signature(PhaseSpace{3,0}()) === PS330
        @test Signature(PhaseSpace{3,1}()) === PS331
        @test Signature(PS330) === PS330                # idempotent on a Signature value

        # Over the bridged value the full kernel surface works: products, automorphisms, inverse.
        # 2D phase space (4 dims, 16 blades) is small enough to compile the real geometric product, so
        # this exercises the bridged subtype driving the @generated `mul` kernel end-to-end.
        @testset "bridged 2D phase space (2, 2, 0) drives the generated product" begin
            PS = Signature(PhaseSpace{2,0}())
            C = CliffordNumber{PS,Float64}
            rng = MersenneTwister(0x2d22)
            a = C(ntuple(_ -> randn(rng), nblades(C)))
            b = C(ntuple(_ -> randn(rng), nblades(C)))
            c = C(ntuple(_ -> randn(rng), nblades(C)))
            @test a * b isa CliffordNumber{PS}
            @test (a * b) * c ≈ a * (b * c)            # associativity through the real product
            for f in (reverse, adjoint, conj, grade_involution)
                @test f(a) === refauto(f, a)
            end
            e1 = KVector{1,PS,Float64}(ntuple(i -> Float64(i == 1), 4))
            @test inv(e1) * e1 ≈ one(C)                # e₁² = +1
        end

        # 3D phase space (64 blades): automorphisms go through the product-free sign path, so this
        # works on the bridge without paying the 6D product compile.
        @testset "bridged 3D phase space (3, 3, 0)" begin
            PS = Signature(PhaseSpace{3,0}())
            C = CliffordNumber{PS,Float64}
            rng = MersenneTwister(0x6a1d)
            a = C(ntuple(_ -> randn(rng), nblades(C)))
            for f in (reverse, adjoint, conj, grade_involution)
                @test f(a) === refauto(f, a)
            end
        end
    end
end
