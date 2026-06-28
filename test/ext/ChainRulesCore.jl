using Random
import CliffordNumbers: scalar_product, versor_inverse, nblades

# Build a Clifford number of concrete type `C` with random Float64 coefficients.
randc(::Type{C}) where C<:AbstractCliffordNumber = C(ntuple(_ -> randn(), Val(nblades(C))))

# Euclidean (coefficient-wise) linear functional ⟨W, Ω⟩, used to turn a Clifford-valued map into a
# scalar loss whose ChainRules cotangent at the output is exactly `W`.
eu(W::AbstractCliffordNumber, Ω::AbstractCliffordNumber) = sum(Tuple(W) .* Tuple(Ω))

# Central finite-difference gradient of a scalar loss `L` with respect to the coefficients of a
# Clifford number `x`, returned as a Clifford number of the same type.
function ngrad(L, x::AbstractCliffordNumber; h = 1.0e-5)
    t = Tuple(x)
    n = length(t)
    g = ntuple(n) do i
        xp = typeof(x)(ntuple(j -> j == i ? t[j] + h : t[j], n))
        xm = typeof(x)(ntuple(j -> j == i ? t[j] - h : t[j], n))
        (L(xp) - L(xm)) / 2h
    end
    return typeof(x)(g)
end

# Central finite-difference derivative of a scalar loss with respect to a scalar argument.
nderiv(L, a; h = 1.0e-5) = (L(a + h) - L(a - h)) / 2h

# Apply an rrule pullback and strip the leading `NoTangent()` (the seed cotangent is `one`).
pull(f, args...; seed = 1.0) = (y = rrule(f, args...); (y[1], y[2](seed)))

approx(a, b) = isapprox(a, b; atol = 1.0e-5, rtol = 1.0e-4)

@testset "ChainRulesCore.jl extension" begin
    Random.seed!(0x1234)

    # Representative type triples (X, Y) spanning Euclidean (VGA), indefinite non-degenerate (STA,
    # CGA) and degenerate (PGA) metrics, with k-vector, even/odd, and dense representations.
    product_cases = [
        (KVector{1,VGA(3),Float64,3},        KVector{1,VGA(3),Float64,3}),
        (EvenCliffordNumber{VGA(3),Float64,4}, KVector{1,VGA(3),Float64,3}),
        (EvenCliffordNumber{VGA(3),Float64,4}, EvenCliffordNumber{VGA(3),Float64,4}),
        (CliffordNumber{VGA(2),Float64,4},   CliffordNumber{VGA(2),Float64,4}),
        (KVector{1,STA,Float64,4},           KVector{1,STA,Float64,4}),
        (EvenCliffordNumber{STA,Float64,8},  KVector{1,STA,Float64,4}),
        (KVector{1,PGA(3),Float64,4},        KVector{1,PGA(3),Float64,4}),
        (EvenCliffordNumber{PGA(3),Float64,8}, EvenCliffordNumber{PGA(3),Float64,8}),
        (KVector{2,CGA(2),Float64,6},        KVector{1,CGA(2),Float64,4}),
    ]

    @testset "geometric product pullback ($(nameof(X)) * $(nameof(Y)))" for (X, Y) in product_cases
        x = randc(X)
        y = randc(Y)
        Ω, pb = rrule(*, x, y)
        @test Ω == x * y
        W = randc(typeof(Ω))
        (_, x̄, ȳ) = pb(W)
        # Pullback of ⟨W, x*y⟩ matches finite differences exactly, in every signature.
        @test approx(x̄, ngrad(xx -> eu(W, xx * y), x))
        @test approx(ȳ, ngrad(yy -> eu(W, x * yy), y))
        @test x̄ isa X
        @test ȳ isa Y
    end

    @testset "geometric product pullback reduces to reverse formula (positive-definite)" begin
        # In a Euclidean algebra the exact transpose equals the familiar (W*y', x'*W). Use dense
        # types so the outputs need no grade projection and the identity holds coefficient-for-
        # coefficient.
        x = randc(CliffordNumber{VGA(3),Float64,8})
        y = randc(CliffordNumber{VGA(3),Float64,8})
        Ω, pb = rrule(*, x, y)
        W = randc(typeof(Ω))
        (_, x̄, ȳ) = pb(W)
        @test approx(x̄, W * y')
        @test approx(ȳ, x' * W)
    end

    @testset "automorphism pullbacks ($f)" for f in (reverse, adjoint, conj, grade_involution)
        for C in (KVector{2,VGA(4),Float64,6}, EvenCliffordNumber{STA,Float64,8},
                  CliffordNumber{PGA(2),Float64,8})
            x = randc(C)
            y, pb = rrule(f, x)
            @test y == f(x)
            W = randc(C)
            (_, x̄) = pb(W)
            @test approx(x̄, ngrad(xx -> eu(W, f(xx)), x))
        end
    end

    @testset "scalar multiplication and division pullbacks" begin
        x = randc(EvenCliffordNumber{VGA(3),Float64,4})
        a = 1.7
        W = randc(typeof(x))
        # x * a
        (_, x̄, ā) = rrule(*, x, a)[2](W)
        @test approx(x̄, ngrad(xx -> eu(W, xx * a), x))
        @test ā ≈ nderiv(aa -> eu(W, x * aa), a)
        # a * x
        (_, ā2, x̄2) = rrule(*, a, x)[2](W)
        @test approx(x̄2, ngrad(xx -> eu(W, a * xx), x))
        @test ā2 ≈ nderiv(aa -> eu(W, aa * x), a)
        # x / a
        (_, x̄3, ā3) = rrule(/, x, a)[2](W)
        @test approx(x̄3, ngrad(xx -> eu(W, xx / a), x))
        @test ā3 ≈ nderiv(aa -> eu(W, x / aa), a)
    end

    @testset "scalar_product and abs2 pullbacks" begin
        for (X, Y) in [(KVector{1,STA,Float64,4}, KVector{1,STA,Float64,4}),
                       (EvenCliffordNumber{VGA(3),Float64,4}, EvenCliffordNumber{VGA(3),Float64,4})]
            x = randc(X)
            y = randc(Y)
            (_, x̄, ȳ) = rrule(scalar_product, x, y)[2](1.0)
            @test approx(x̄, ngrad(xx -> scalar_product(xx, y), x))
            @test approx(ȳ, ngrad(yy -> scalar_product(x, yy), y))
        end
        for C in (KVector{1,VGA(3),Float64,3}, EvenCliffordNumber{STA,Float64,8})
            x = randc(C)
            (_, x̄) = rrule(abs2, x)[2](1.0)
            @test approx(x̄, ngrad(abs2, x))
        end
    end

    @testset "additive pullbacks" begin
        x = randc(EvenCliffordNumber{VGA(3),Float64,4})
        y = randc(KVector{1,VGA(3),Float64,3})
        Ωp, pbp = rrule(+, x, y)
        @test Ωp == x + y
        W = randc(typeof(Ωp))
        (_, x̄, ȳ) = pbp(W)
        @test approx(x̄, ngrad(xx -> eu(W, xx + y), x))
        @test approx(ȳ, ngrad(yy -> eu(W, x + yy), y))
        (_, x̄m, ȳm) = rrule(-, x, y)[2](W)
        @test approx(x̄m, ngrad(xx -> eu(W, xx - y), x))
        @test approx(ȳm, ngrad(yy -> eu(W, x - yy), y))
        Wn = randc(typeof(x))
        (_, x̄n) = rrule(-, x)[2](Wn)
        @test approx(x̄n, ngrad(xx -> eu(Wn, -xx), x))
    end

    @testset "versor_inverse composes from reverse / abs2 / scalar division" begin
        # versor_inverse(x) = x' / abs2(x): differentiate by hand-composing the underlying rules and
        # check the result against finite differences of the whole expression.
        x = randc(EvenCliffordNumber{VGA(3),Float64,4})
        W = randc(typeof(x))
        f(z) = eu(W, versor_inverse(z))
        # Forward
        xr = x'
        a = abs2(x)
        # Reverse
        (_, x̄_div, ā) = rrule(/, xr, a)[2](W)
        (_, x̄_abs2) = rrule(abs2, x)[2](ā)
        (_, x̄_rev) = rrule(adjoint, x)[2](x̄_div)
        x̄ = x̄_rev + x̄_abs2
        @test approx(x̄, ngrad(f, x))
    end

    @testset "rotor sandwich gradient (sandwich = product + reverse)" begin
        # f(R) = ‖R p R' - q‖²  for a rotor R and vectors p, q in VGA(3). The gradient w.r.t. R is
        # obtained purely by composing the product, adjoint, subtraction, and abs2 rules — no
        # dedicated sandwich rule — and matches finite differences.
        R = randc(EvenCliffordNumber{VGA(3),Float64,4})
        p = randc(KVector{1,VGA(3),Float64,3})
        q = randc(KVector{1,VGA(3),Float64,3})

        f(r) = abs2(r * p * r' - q)

        # Forward pass, recording intermediates.
        A, pbA = rrule(*, R, p)        # A = R * p
        Radj, pbRadj = rrule(adjoint, R)
        B, pbB = rrule(*, A, Radj)     # B = A * R'
        D, pbD = rrule(-, B, q)        # D = B - q
        c, pbc = rrule(abs2, D)        # c = ‖D‖²
        @test c ≈ f(R)

        # Reverse pass.
        (_, D̄) = pbc(1.0)
        (_, B̄, q̄) = pbD(D̄)
        (_, Ā, R̄adj) = pbB(B̄)
        (_, R̄_from_adj) = pbRadj(R̄adj)
        (_, R̄_from_A, p̄) = pbA(Ā)
        R̄ = R̄_from_A + R̄_from_adj

        @test approx(R̄, ngrad(f, R))
    end

    @testset "ProjectTo maps cotangents onto the primal grade structure" begin
        x = randc(EvenCliffordNumber{VGA(3),Float64,4})
        p = ProjectTo(x)
        # A wider (dense) cotangent is projected back to the even subspace.
        dense = randc(CliffordNumber{VGA(3),Float64,8})
        @test p(dense) isa typeof(x)
        @test p(dense) == typeof(x)(dense)
        # AbstractZero passes through unchanged.
        @test p(ZeroTangent()) === ZeroTangent()
        @test p(NoTangent()) === NoTangent()
    end
end
