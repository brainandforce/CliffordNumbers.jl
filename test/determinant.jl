@testset "Determinant, trace, and characteristic polynomial" begin
    using CliffordNumbers: charpoly_coeffs, dimension
    using LinearAlgebra: det, tr

    # Dense 2ⁿ×2ⁿ matrix of the left-multiplication operator L_x : y ↦ x*y, the ground-truth oracle.
    function leftmul_matrix(x::AbstractCliffordNumber{Q,T}) where {Q,T}
        xd = CliffordNumber{Q,T}(x)
        N = CliffordNumbers.blade_count(Q)
        M = Matrix{T}(undef, N, N)
        for j in 1:N
            ej = CliffordNumber{Q,T}(ntuple(i -> T(i == j), N))
            M[:, j] = collect(Tuple(xd * ej))
        end
        return M
    end

    # Exact determinant by Gaussian elimination, independent of the Faddeev–LeVerrier path in `det`.
    function exact_det(A::AbstractMatrix{T}) where T
        M = Matrix{T}(A)
        n = size(M, 1)
        d = one(T)
        for k in 1:n
            p = findfirst(!iszero, @view M[k:n, k])
            p === nothing && return zero(T)
            p += k - 1
            if p != k
                M[[k, p], :] = M[[p, k], :]
                d = -d
            end
            d *= M[k, k]
            for i in (k + 1):n
                f = M[i, k] / M[k, k]
                M[i, :] .-= f .* M[k, :]
            end
        end
        return d
    end

    randmv(rng, ::Type{C}) where C = C(ntuple(_ -> Rational{BigInt}(rand(rng, -3:3), rand(rng, 1:3)), nblades(C)))
    randf(rng, ::Type{C}) where C = C(ntuple(_ -> randn(rng), nblades(C)))

    rng = MersenneTwister(0x5eed)
    small = (CliffordNumber{VGA(1)}, CliffordNumber{VGA(2)}, CliffordNumber{VGA(3)},
             CliffordNumber{PGA(2)}, CliffordNumber{STA}, CliffordNumber{CGA(2)})
    large = (CliffordNumber{VGA(4)}, CliffordNumber{PGA(3)}, CliffordNumber{CGA(3)})

    @testset "matches the dense left-multiplication matrix" begin
        # Exact rational oracle for the smaller algebras (n ≤ 4)
        for C in small, _ in 1:4
            x = randmv(rng, C{Rational{BigInt}})
            @test det(x) == exact_det(leftmul_matrix(x))
        end
        # Float oracle for the larger algebras (n = 4, 5)
        for C in large, _ in 1:3
            x = randf(rng, C{Float64})
            @test det(x) ≈ LinearAlgebra.det(leftmul_matrix(x)) rtol=1e-6
        end
    end

    @testset "det(λ) = λ^(2ⁿ) and multiplicativity" begin
        for C in (small..., large...)
            Q = signature(C)
            n2 = CliffordNumbers.blade_count(Q)
            λ = 3 // 2
            @test det(CliffordNumber{Q,Rational{BigInt}}(λ)) == λ^n2
        end
        for C in (CliffordNumber{VGA(2)}, CliffordNumber{VGA(3)}, CliffordNumber{STA}), _ in 1:3
            x = randmv(rng, C{Rational{BigInt}})
            y = randmv(rng, C{Rational{BigInt}})
            @test det(x * y) == det(x) * det(y)
        end
    end

    @testset "singular elements have zero determinant" begin
        # 1 + e₁ is a zero divisor wherever e₁² = +1: (1+e₁)(1−e₁) = 1 − e₁² = 0
        for Q in (VGA(2), VGA(3), VGA(4))
            x = one(CliffordNumber{Q,Rational{BigInt}}) +
                KVector{1,Q,Rational{BigInt}}(ntuple(i -> Rational{BigInt}(i == 1), dimension(Q)))
            @test det(x) == 0
        end
    end

    @testset "trace" begin
        for C in (CliffordNumber{VGA(3)}, EvenCliffordNumber{STA}, KVector{2,PGA(3)})
            x = randf(rng, C{Float64})
            @test tr(x) == CliffordNumbers.blade_count(signature(C)) * scalar(x)
        end
    end

    @testset "characteristic polynomial and Cayley–Hamilton" begin
        for C in (CliffordNumber{VGA(2)}, CliffordNumber{VGA(3)})
            Q = signature(C)
            N = CliffordNumbers.blade_count(Q)
            x = randmv(rng, C{Rational{BigInt}})
            c = charpoly_coeffs(x)
            @test length(c) == N
            @test -c[1] == tr(x)
            @test (isodd(N) ? -c[end] : c[end]) == det(x)
            # x satisfies its own characteristic polynomial: xᴺ + c₁xᴺ⁻¹ + ⋯ + c_N = 0, evaluated by
            # Horner from the monic leading term (h starts at the x⁰ coefficient 1)
            xd = CliffordNumber{Q,Rational{BigInt}}(x)
            h = one(xd)
            for k in 1:N
                h = muladd(c[k], one(xd), xd * h)
            end
            @test iszero(h)
        end
    end

    @testset "allocation-free and inferred determinant" begin
        for C in (KVector{2,VGA(3)}, EvenCliffordNumber{VGA(3)}, CliffordNumber{VGA(3)})
            x = randf(rng, C{Float64})
            det(x)  # warm up
            ALLOCATION_GATES && @test @allocated(CliffordNumbers.det(x)) == 0
            @test (@inferred CliffordNumbers.det(x)) === CliffordNumbers.det(x)
        end
    end

    @testset "matrix path agrees with the matrix-free path" begin
        # The n ≥ 6 dense-matrix fallback (forced here on small algebras) must match Faddeev–LeVerrier
        for C in (CliffordNumber{VGA(3)}, EvenCliffordNumber{STA}, KVector{2,CGA(2)}), _ in 1:3
            x = randmv(rng, C{Rational{BigInt}})
            @test CliffordNumbers._det(Val(false), x) == CliffordNumbers.det(x)
        end
    end

    @testset "LinearAlgebra forwarding" begin
        x = randf(rng, CliffordNumber{VGA(3),Float64})
        @test LinearAlgebra.det(x) === CliffordNumbers.det(x)
        @test LinearAlgebra.tr(x) === CliffordNumbers.tr(x)
    end
end
