@testset "StructArrays extension" begin
    @testset "Schema and column names" begin
        # Bit-pattern names for a bare grade-2 KVector in 3D
        kbiv = KVector{2,VGA(3)}(1.0, 2.0, 3.0)
        @test StructArrays.staticschema(typeof(kbiv)) ===
            NamedTuple{(:e12, :e13, :e23), NTuple{3,Float64}}
        # Grade-0 component is named :scalar; even subalgebra of VGA(3) ≅ ℍ
        eq = EvenCliffordNumber{VGA(3)}(1.0, 2.0, 3.0, 4.0)
        @test StructArrays.staticschema(typeof(eq)) ===
            NamedTuple{(:scalar, :e12, :e13, :e23), NTuple{4,Float64}}
        # Vectors get one column per dimension
        kvec = KVector{1,VGA(3)}(1.0, 2.0, 3.0)
        @test StructArrays.staticschema(typeof(kvec)) ===
            NamedTuple{(:e1, :e2, :e3), NTuple{3,Float64}}
        # CGA has a basis index of -1, so names fall back to positional :d1 … :dL
        kcga = KVector{1,CGA(2)}(1.0, 2.0, 3.0, 4.0)
        @test StructArrays.staticschema(typeof(kcga)) ===
            NamedTuple{(:d1, :d2, :d3, :d4), NTuple{4,Float64}}
    end

    @testset "SoA round-trip and columns" begin
        v = [KVector{1,VGA(3)}(Float64(i), 2i, 3i) for i in 1:5]
        sa = StructArray(v)
        # Each blade coefficient is its own column, not a single opaque NTuple column
        @test propertynames(sa) === (:e1, :e2, :e3)
        @test sa.e1 isa Vector{Float64}
        @test sa.e1 == Float64.(1:5)
        @test sa.e2 == Float64.(2 .* (1:5))
        @test sa.e3 == Float64.(3 .* (1:5))
        # Elementwise reconstruction is exact
        @test eltype(sa) === KVector{1,VGA(3),Float64,3}
        @test collect(sa) == v
        @test all(sa[i] === v[i] for i in eachindex(v))
    end

    @testset "Mutation through columns" begin
        v = [EvenCliffordNumber{VGA(3)}(Float64(i), 0.0, 0.0, 0.0) for i in 1:4]
        sa = StructArray(v)
        sa[2] = EvenCliffordNumber{VGA(3)}(10.0, 11.0, 12.0, 13.0)
        @test sa[2] === EvenCliffordNumber{VGA(3)}(10.0, 11.0, 12.0, 13.0)
        @test sa.scalar[2] == 10.0
        @test sa.e23[2] == 13.0
        # Column writes are reflected in reconstructed elements (e12 is the 2nd coefficient)
        sa.e12 .= -1.0
        @test all(Tuple(sa[i])[2] == -1.0 for i in eachindex(v))
    end

    @testset "Constructing a StructArray from columns" begin
        sa = StructArray{KVector{1,VGA(2),Float64,2}}((e1 = [1.0, 2.0], e2 = [3.0, 4.0]))
        @test sa[1] === KVector{1,VGA(2)}(1.0, 3.0)
        @test sa[2] === KVector{1,VGA(2)}(2.0, 4.0)
    end

    # Once each blade coefficient is its own column (Stage 1), the *default* StructArrays
    # broadcast machinery already operates on the SoA columns: geometric-operator broadcasts
    # over an SoA array stay SoA, agree elementwise with the plain-`Vector` path, and (because
    # these types are `isbits`) allocate only the output, never a per-element `KVector`. These
    # tests lock that behavior in.
    sandwich(rotor, p) = rotor * p * rotor'

    @testset "Columnar broadcast: SoA preserved and correct" begin
        R = EvenCliffordNumber{VGA(3)}(0.8, 0.1, 0.2, 0.3)
        pts = [KVector{1,VGA(3)}(float(i), 2.0i, 3.0i) for i in 1:6]
        sa = StructArray(pts)
        # A lone Clifford number broadcasts as a scalar; the result stays SoA. (A columnar
        # broadcast contracts FMAs in a different order than a scalar `*`, so products agree
        # with the per-element reference only up to floating-point rounding, hence `≈`.)
        @test R .* sa isa StructArray
        @test all(collect(R .* sa) .≈ [R * p for p in pts])
        # Fused product chains stay SoA and likewise agree up to rounding.
        @test (R .* sa .* R') isa StructArray
        @test all(collect(R .* sa .* R') .≈ [R * p * R' for p in pts])
        # Sandwich (rotation), rotor passed as a broadcast argument
        @test sandwich.(Ref(R), sa) isa StructArray
        @test all(collect(sandwich.(Ref(R), sa)) .≈ [sandwich(R, p) for p in pts])
        # Unary geometric automorphism
        @test collect(reverse.(sa)) == [reverse(p) for p in pts]
    end

    @testset "Columnar broadcast: no per-element allocation" begin
        # In-place broadcast cost must be O(1) in the number of elements, never O(N): if a
        # `KVector` were materialized per element, the cost would scale with N.
        function inplace_cost(n, R)
            sw(rotor, p) = rotor * p * rotor'
            sa = StructArray([KVector{1,VGA(3)}(float(i), 2.0i, 3.0i) for i in 1:n])
            out = sw.(Ref(R), sa)
            f() = (out .= sw.(Ref(R), sa))
            f(); f()    # warm up / compile
            return @allocated f()
        end
        R = EvenCliffordNumber{VGA(3)}(0.8, 0.1, 0.2, 0.3)
        small = inplace_cost(100, R)
        large = inplace_cost(5000, R)
        @test large <= small + 512
    end
end
