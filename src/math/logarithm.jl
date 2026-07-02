#---Logarithms and square roots of even multivectors-----------------------------------------------#
# These invert `exp` on the even subalgebra, mapping a unit rotor or motor
# (`x * x' ≈ 1`) back to its generating bivector. The algorithm is the "simple
# invariant decomposition" of Hadfield, Wieser, and Lasenby: below six dimensions
# a bivector splits into at most two commuting simple bivectors (`_bivector_split`
# in factorization.jl, exposed as `bivector_decomposition`), and the per-plane
# angles are recovered from the scalar, grade-2, and grade-4 parts of the rotor.
# This covers the elliptic, hyperbolic, and degenerate (ideal) cases uniformly.
# Non-unit dilators and even multivectors with an ideal part beyond the motor
# case are out of scope.

@inline _grade2(x::AbstractCliffordNumber{Q}) where Q = KVector{2,Q}(x)
@inline _grade4(x::AbstractCliffordNumber{Q}) where Q = KVector{4,Q}(x)

# Squared "length" of a bivector as a signed scalar: <b²>₀ = scalar_product(b, b). Its sign
# classifies the plane as elliptic (negative), hyperbolic (positive), or null (zero).
@inline _plane_square(b) = scalar_product(b, b)

"""
    CliffordNumbers._simple_log(s, b)

Logarithm of a simple rotor whose only nonzero parts are the scalar `s` and a *simple* bivector `b`
(one whose square is a scalar). Returns the generating bivector.
"""
function _simple_log(s, b)
    β = _plane_square(b)
    if β < zero(β)              # elliptic plane: s = cos θ, |b| = sin θ
        n = sqrt(-β)
        return (atan(n, s) / n) * b
    elseif β > zero(β)          # hyperbolic plane: s = cosh θ, |b| = sinh θ
        n = sqrt(β)
        return (asinh(n) / n) * b
    else                        # null/ideal plane (b² = 0), or the identity (b = 0)
        return iszero(s) ? b : b * inv(s)
    end
end

# Angle (or rapidity) of a plane from a 2-vector `(c, sdir)` aligned with `(cos, sin)` for elliptic
# planes, `(cosh, sinh)` for hyperbolic planes, and `(1, θ)` for null planes.
@inline function _plane_angle(c, sdir, ε)
    if ε < zero(ε)
        return atan(sdir, c)
    elseif ε > zero(ε)
        return atanh(clamp(sdir / c, -one(c), one(c)))
    else
        return sdir / c
    end
end

"""
    CliffordNumbers._general_log(R, s, b, b1, b2)

Logarithm of a rotor whose bivector part `b` splits into two distinct commuting simple bivectors
`b1` and `b2`. The per-plane angles are recovered from the rank-one matrix `[s n2; n1 ω]`, whose
columns are proportional to each plane's `(cos, sin)` pair.
"""
function _general_log(R, s, b, b1, b2)
    β1 = _plane_square(b1); β2 = _plane_square(b2)
    n1 = sqrt(abs(β1)); n2 = sqrt(abs(β2))
    # A plane is null (only possible in degenerate metrics) when its metric magnitude vanishes even
    # though the bivector itself does not. Its log contribution is simply b_null / s.
    tol = sqrt(eps(real(typeof(s)))) * (sum(abs2, Tuple(b)) + one(s))
    null1 = abs(β1) <= tol
    null2 = abs(β2) <= tol
    null1 && null2 && return b * inv(s)
    null1 && return b1 * inv(s) + _simple_log(s, b2)
    null2 && return b2 * inv(s) + _simple_log(s, b1)
    # signed coefficient of the unit 4-blade B̂1 ∧ B̂2 in the grade-4 part of R
    W  = _grade4(R)
    U4 = _grade4(b1 * b2)
    ω  = scalar_product(W, U4) / scalar_product(U4, U4) * (n1 * n2)
    col1 = (s, n1); col2 = (n2, ω)
    # both columns are proportional to (cos, sin) of plane 1; pick the larger for conditioning
    u = (col1[1]^2 + col1[2]^2) >= (col2[1]^2 + col2[2]^2) ? col1 : col2
    v = (s * u[1] + n1 * u[2], n2 * u[1] + ω * u[2])   # = Mᵀu, proportional to (cos, sin) of plane 2
    # `u` and `v` carry a shared sign gauge. Elliptic angles come from `atan` (sign-sensitive), so
    # fix the gauge from the hyperbolic plane, whose `cosh` must be positive (a hyperbolic plane
    # admits no π-shift). With no hyperbolic plane the shared sign cancels in the round trip.
    g = one(s)
    if β1 > zero(β1)
        g = ifelse(u[1] < zero(u[1]), -one(s), one(s))
    elseif β2 > zero(β2)
        g = ifelse(v[1] < zero(v[1]), -one(s), one(s))
    end
    u = (g * u[1], g * u[2])
    v = (g * v[1], g * v[2])
    θ1 = _plane_angle(u[1], u[2], sign(β1))
    θ2 = _plane_angle(v[1], v[2], sign(β2))
    return (θ1 / n1) * b1 + (θ2 / n2) * b2
end

@generated function _rotor_log(R::EvenCliffordNumber{Q,T}) where {Q,T}
    # In fewer than four dimensions every bivector is simple, so no decomposition is needed (and the
    # grade-4 machinery would not even type-check).
    dimension(Q) < 4 && return :(_simple_log(scalar(R), _grade2(R)))
    return quote
        s = scalar(R)
        b = _grade2(R)
        (b1, b2) = _bivector_split(b)
        iszero(b2) && return _simple_log(s, b)
        return _general_log(R, s, b, b1, b2)
    end
end

# Newton refinement polishes the approximate logarithm: if `exp(L) ≈ R`, then `R * exp(-L)` is a
# small rotor whose generator is approximately its grade-2 part over its scalar part. This restores
# full accuracy in the ill-conditioned neighborhood where the two invariant planes nearly coincide.
# Iteration continues only while the residual `‖R exp(-L) - 1‖` keeps shrinking, which both
# terminates at the achievable precision and guards against divergence for large boosts.
function _refine_log(R::EvenCliffordNumber, L)
    E = R * exp(-L)
    err = sum(abs2, Tuple(E - one(E)))
    for _ in 1:8
        δ = _grade2(E) * inv(scalar(E))
        Lnext = L + δ
        Enext = R * exp(-Lnext)
        errnext = sum(abs2, Tuple(Enext - one(Enext)))
        errnext >= err && break
        L = Lnext; E = Enext; err = errnext
    end
    return L
end

function _even_log(R::EvenCliffordNumber{Q}) where Q
    L = _rotor_log(R)
    # Below four dimensions every bivector is simple and the closed form is already exact.
    dimension(Q) < 4 && return L
    return _refine_log(R, L)
end

"""
    log(x::EvenCliffordNumber)

Returns the natural logarithm of a rotor or motor `x`, the bivector `B` such that `exp(B) ≈ x`.

The result is the principal logarithm: for a bivector `B` of sufficiently small magnitude,
`log(exp(B)) == B`. The implementation handles the elliptic (Euclidean rotation), hyperbolic
(Lorentzian boost), and degenerate (ideal/translation) planes that occur in vanilla, projective, and
conformal geometric algebras up to five dimensions.

`x` is assumed to be a unit versor of even grade (`x * x' ≈ 1`), as produced by [`exp`](@ref) of a
bivector.

See also: [`exp`](@ref), [`sqrt`](@ref).
"""
log(x::EvenCliffordNumber) = _even_log(float(x))

#---Square roots of even multivectors--------------------------------------------------------------#

@generated function _rotor_sqrt(R::EvenCliffordNumber{Q,T}) where {Q,T}
    # Positive-definite algebras with no grade-4 part (VGA(2) ≅ ℂ, VGA(3) ≅ ℍ) admit the exact
    # closed form sqrt(R) = (1 + R) / |1 + R|.
    if is_positive_definite(Q) && dimension(Q) < 4
        return :(normalize(one(R) + R))
    end
    # Otherwise fall back to the general identity sqrt(R) = exp(½ log(R)).
    return :(exp(_even_log(R) / 2))
end

"""
    sqrt(x::EvenCliffordNumber)

Returns the square root of a rotor or motor `x`: the versor `S` of even grade with `S * S ≈ x` and
the smaller rotation/boost.

For positive-definite algebras isomorphic to ℂ or ℍ this uses the closed form
`sqrt(R) = (1 + R) / |1 + R|`; otherwise it uses `sqrt(R) = exp(½ log(R))`, which is valid in the
mixed-signature and degenerate algebras (Minkowski, projective, conformal) as well.

`x` is assumed to be a unit versor of even grade (`x * x' ≈ 1`).

See also: [`exp`](@ref), [`log`](@ref).
"""
sqrt(x::EvenCliffordNumber) = _rotor_sqrt(float(x))

#---Logarithms and square roots of general multivectors--------------------------------------------#
# Beyond the even/versor case, closed forms exist whenever the non-scalar part `n` of `x = s + n`
# squares to a scalar: any vector, blade, or simple bivector plus a scalar. Such an `x` lives in a
# two-dimensional subalgebra isomorphic to ℂ, the split-complex numbers, or the dual numbers,
# according to the sign of `n²`, and the principal branch is the corresponding planar formula
# (the "study number" route of De Keninck & Roelfs, "Normalization, Square Roots, and the
# Exponential and Logarithmic Maps in Geometric Algebras of Less than 6D"). Even multivectors
# additionally reduce to the rotor logarithm after peeling their magnitude `√(x x̃)`.

"""
    CliffordNumbers._study_split(x::AbstractCliffordNumber) -> (s, n, β, flat)

Splits `x` into its scalar part `s` and non-scalar part `n`, with `β = ⟨n²⟩₀`. `flat` is `true`
when `n²` has no non-scalar residue, so that `x` lies in the planar (study number) subalgebra
where the closed-form logarithm and square root apply.
"""
function _study_split(x::AbstractCliffordNumber{Q,T}) where {Q,T}
    s = scalar(x)
    n = x - KVector{0,Q}(s)
    nn = n * n
    β = scalar(nn)
    tol = _structure_rtol(T) * sum(abs2, Tuple(x))
    flat = sum(abs2, Tuple(nn)) - β^2 <= tol^2
    return (s, n, β, flat)
end

"""
    CliffordNumbers._study_log(s, n, β)

Principal logarithm of the study number `s + n` with `β = n²` a scalar: the complex plane formula
for `β < 0`, the split-complex formula (requiring `s > |n|`) for `β > 0`, and the dual-number
formula (requiring `s > 0`) for `β = 0`. Throws a `DomainError` off the principal branch.
"""
function _study_log(s, n, β)
    if β < 0                    # elliptic plane: s + n ≅ the complex number s + |n| im
        m = sqrt(-β)
        return KVector{0,signature(n)}(log(hypot(s, m))) + (atan(m, s) / m) * n
    elseif β > 0                # hyperbolic plane: split-complex; the branch needs s > |n|
        m = sqrt(β)
        s > m || throw(DomainError(s, "no principal logarithm: s + n with n² > 0 needs s > |n|."))
        return KVector{0,signature(n)}(log(s^2 - β) / 2) + (atanh(m / s) / m) * n
    else                        # parabolic plane: dual number s (1 + n/s)
        s > zero(s) ||
            throw(DomainError(s, "no principal logarithm: s + n with n² = 0 needs s > 0."))
        return KVector{0,signature(n)}(log(s)) + n / s
    end
end

"""
    CliffordNumbers._study_sqrt(s, n, β)

Principal square root of the study number `s + n` with `β = n²` a scalar, on the same branches as
[`_study_log`](@ref CliffordNumbers._study_log). Throws a `DomainError` where no real square root
exists.
"""
function _study_sqrt(s, n, β)
    if β < 0
        m = sqrt(-β)
        ρ = sqrt(hypot(s, m))
        θ = atan(m, s) / 2
        return KVector{0,signature(n)}(ρ * cos(θ)) + (ρ * sin(θ) / m) * n
    elseif β > 0
        m = sqrt(β)
        s > m || throw(DomainError(s, "no real square root: s + n with n² > 0 needs s > |n|."))
        ρ = sqrt(sqrt(s^2 - β))
        θ = atanh(m / s) / 2
        return KVector{0,signature(n)}(ρ * cosh(θ)) + (ρ * sinh(θ) / m) * n
    else
        s > zero(s) || throw(DomainError(s, "no real square root: s + n with n² = 0 needs s > 0."))
        return KVector{0,signature(n)}(sqrt(s)) + n / (2 * sqrt(s))
    end
end

"""
    CliffordNumbers._even_content(y::AbstractCliffordNumber{Q,T})

Returns the even-grade content of `y` as an `EvenCliffordNumber{Q,T}`, or `nothing` if `y` has odd
coefficients. The odd energy is the difference between the total energy and the even projection's,
since the coefficients partition by parity.
"""
function _even_content(y::AbstractCliffordNumber{Q,T}) where {Q,T}
    e = EvenCliffordNumber{Q,T}(y)
    tol = _structure_rtol(T)^2 * sum(abs2, Tuple(y))
    return sum(abs2, Tuple(y)) - sum(abs2, Tuple(e)) <= tol ? e : nothing
end

"""
    CliffordNumbers._even_magnitude(e::EvenCliffordNumber)

Returns the squared magnitude `α = e ẽ` when it is a positive scalar — so that `e = √α (e/√α)`
peels into a magnitude and a unit rotor — or `nothing` otherwise.
"""
function _even_magnitude(e::EvenCliffordNumber{Q,T}) where {Q,T}
    r = e * e'
    α = scalar(r)
    tol = _structure_rtol(T) * sum(abs2, Tuple(e))
    (α > zero(α) && sum(abs2, Tuple(r)) - α^2 <= tol^2) || return nothing
    return α
end

"""
    log(x::AbstractCliffordNumber)

Principal natural logarithm of a general multivector, on the branches where a real logarithm
exists:
  * `x = s + n` with `n²` a scalar (a scalar plus a vector, blade, or simple bivector — in
    particular odd-grade blades): the planar closed form, elliptic (`n² < 0`, always defined),
    hyperbolic (`n² > 0`, defined for `s > |n|`), or parabolic (`n² = 0`, defined for `s > 0`).
  * even content with `x x̃` a positive scalar (a rotor or motor scaled by a positive magnitude):
    `log(x) = log(x x̃)/2 + log(x / sqrt(x x̃))` through the rotor logarithm.

Throws a `DomainError` outside these branches; in particular an odd *versor* has no logarithm (the
odd component of the Pin group is not connected to the identity). The result is returned as
`exponential_type(x)`.

See also: [`log(::EvenCliffordNumber)`](@ref), [`exp`](@ref), [`sqrt`](@ref).
"""
function log(x::AbstractCliffordNumber)
    y = float(x)
    C = exponential_type(y)
    (s, n, β, flat) = _study_split(y)
    flat && return convert(C, _study_log(s, n, β))
    e = _even_content(y)
    if !isnothing(e)
        α = _even_magnitude(e)
        if !isnothing(α)
            return convert(C, KVector{0,signature(y)}(log(α) / 2) + _even_log(e / sqrt(α)))
        end
    end
    throw(DomainError(x, "no principal logarithm exists for this multivector."))
end

log(x::KVector{0,Q}) where Q = KVector{0,Q}(log(scalar(x)))

"""
    sqrt(x::AbstractCliffordNumber)

Principal square root of a general multivector, on the branches where a real square root exists.
The branches mirror [`log(::AbstractCliffordNumber)`](@ref): the planar closed form for
`x = s + n` with `n²` a scalar, and `sqrt(x) = (x x̃)^(1/4) sqrt(x / sqrt(x x̃))` through the rotor
square root for even content with `x x̃` a positive scalar. Throws a `DomainError` outside these
branches. The result is returned as `exponential_type(x)`.

See also: [`sqrt(::EvenCliffordNumber)`](@ref), [`log(::AbstractCliffordNumber)`](@ref).
"""
function sqrt(x::AbstractCliffordNumber)
    y = float(x)
    C = exponential_type(y)
    (s, n, β, flat) = _study_split(y)
    flat && return convert(C, _study_sqrt(s, n, β))
    e = _even_content(y)
    if !isnothing(e)
        α = _even_magnitude(e)
        isnothing(α) || return convert(C, sqrt(sqrt(α)) * _rotor_sqrt(e / sqrt(α)))
    end
    throw(DomainError(x, "no real square root exists for this multivector."))
end
