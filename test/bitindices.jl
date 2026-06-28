@testset "Compile-time blade data (S1)" begin
    using CliffordNumbers: bitindices, blade_signs, grade_negate, BladeIndices, nblades

    types = (
        CliffordNumber{VGA(3)}, EvenCliffordNumber{VGA(3)}, OddCliffordNumber{VGA(3)},
        KVector{0,VGA(3)}, KVector{1,VGA(3)}, KVector{2,VGA(3)}, KVector{3,VGA(3)},
        CliffordNumber{STA}, EvenCliffordNumber{PGA(3)}, KVector{2,PGA(3)}, KVector{2,CGA(3)},
    )

    @testset "bitindices matches blade enumeration" begin
        for C in types
            @test bitindices(C) === Tuple(BladeIndices(C))
            @test length(bitindices(C)) == nblades(C)
            # Every stored blade has a grade the type actually represents
            @test all(in(nonzero_grades(C)), grade.(bitindices(C)))
        end
        # Trait agrees on instances and types
        x = CliffordNumber{VGA(3),Float64}(ntuple(Float64, 8))
        @test bitindices(x) === bitindices(typeof(x))
    end

    @testset "blade_signs matches per-grade automorphism signs" begin
        # Independent closed forms for the sign each automorphism applies to a grade-g blade
        refsign(::Val{:reverse}, g) = (-1)^(g * (g - 1) ÷ 2)
        refsign(::Val{:adjoint}, g) = (-1)^(g * (g - 1) ÷ 2)
        refsign(::Val{:grade_involution}, g) = (-1)^g
        refsign(::Val{:conj}, g) = (-1)^(g * (g + 1) ÷ 2)
        for C in types, f in (:reverse, :adjoint, :grade_involution, :conj)
            expected = map(b -> Int8(refsign(Val(f), grade(b))), bitindices(C))
            @test blade_signs(C, Val(f)) === expected
        end
    end

    @testset "refactored automorphisms preserve semantics" begin
        # Locks the blade_signs rewrite to the generic indexing reference `T(x[f.(BladeIndices)])`
        refauto(f, x) = (T = typeof(x); T(x[f.(BladeIndices(T))]))
        for C in types
            x = C(ntuple(i -> Float64(i) - 2.5, nblades(C)))
            for f in (reverse, adjoint, conj, grade_involution)
                @test f(x) === refauto(f, x)
            end
        end
    end

    @testset "grade_negate" begin
        x = CliffordNumber{VGA(3),Int}(ntuple(i -> i, 8))
        # Negating no grades is the identity; negating every grade is plain negation
        @test grade_negate(x, Val(())) === x
        @test grade_negate(x, Val((0, 1, 2, 3))) === -x
        # Negating grades {2,3} flips exactly the grade-2 and grade-3 coefficients
        ref = CliffordNumber{VGA(3),Int}(map(b -> ifelse(grade(b) in (2, 3), -1, 1), bitindices(x)) .* Tuple(x))
        @test grade_negate(x, Val((2, 3))) === ref
    end

    @testset "allocation-free and inferred on Float64" begin
        x = CliffordNumber{VGA(3),Float64}(ntuple(Float64, 8))
        for f in (reverse, adjoint, conj, grade_involution)
            f(x)  # warm up
            @test @allocated(f(x)) == 0
            @test (@inferred f(x)) === f(x)
        end
        blade_signs(typeof(x), Val(:reverse))
        @test @allocated(blade_signs(typeof(x), Val(:reverse))) == 0
    end
end
