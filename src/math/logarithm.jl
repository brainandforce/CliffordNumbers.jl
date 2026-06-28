#---Logarithms and square roots of even multivectors-----------------------------------------------#
# These operations invert `exp` on the even subalgebra: they map a rotor or motor (a versor of even
# grade, with `x * x' ≈ 1`) back to the bivector that generates it. The algorithms here follow the
# "simple invariant decomposition" of Hadfield, Wieser, and Lasenby: a bivector in fewer than six
# dimensions splits into at most two commuting simple bivectors, each of which exponentiates within
# its own plane. The decomposition is recovered from the scalar, grade-2, and grade-4 parts of the
# rotor, which lets `log` and `sqrt` cover the Euclidean (elliptic), Lorentzian (hyperbolic), and
# degenerate (parabolic/ideal) cases uniformly.
#
# Scope: these methods assume `x` is a unit rotor or motor (the image of `exp` on a bivector). They
# are not intended for general even multivectors that carry both a rotor and an independent ideal
# part beyond the motor case, nor for non-unit dilators.

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

"""
    CliffordNumbers._iso_log(s, b, P, G4)

Logarithm of an isoclinic rotor, where the two invariant planes rotate by equal angles so the
generating bivector is a scalar multiple `k * b` of the (non-simple) bivector part `b`. Here `P` is
`<b²>₀` and `G4` is `<b²>₄`.
"""
function _iso_log(s, b, P, G4)
    bg = _grade2(b * G4)
    # b * G4 is parallel to b for an isoclinic bivector: b*G4 = σ b.
    σ = scalar_product(bg, b) / _plane_square(b)
    μ = P + σ
    c = one(s) + (s - one(s)) * μ / P
    if μ < zero(μ)              # elliptic isoclinic
        ω = sqrt(-μ)
        return (acos(clamp(c, -one(c), one(c))) / ω) * b
    else                        # hyperbolic isoclinic
        ω = sqrt(μ)
        return (acosh(max(c, one(c))) / ω) * b
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
        s  = scalar(R)
        b  = _grade2(R)
        bb = b * b
        P  = scalar(bb)
        G4 = _grade4(bb)
        g  = scalar_product(G4, G4)
        disc = P * P - g
        scale = P * P + abs(g) + one(T)
        sq = sqrt(max(disc, zero(disc)))
        # Numerically stable quadratic roots (λ are the two plane squares, with λ1λ2 = g/4): take
        # the larger-magnitude root directly and the smaller from the product to avoid cancellation
        # when one plane nearly vanishes (a near-simple rotor).
        λ1 = (P + copysign(sq, P)) / 2
        λ2 = iszero(λ1) ? (P - sq) / 2 : (g / 4) / λ1
        if abs(λ1 - λ2) <= sqrt(eps(T)) * (abs(P) + one(T))
            # repeated eigenvalues: isoclinic if there is a grade-4 part, otherwise simple
            (abs(g) <= eps(T) * scale || iszero(P)) && return _simple_log(s, b)
            return _iso_log(s, b, P, G4)
        end
        bg = _grade2(b * G4)
        b1 = _grade2((λ1 * b - bg / 2) * inv(λ1 - λ2))
        b2 = _grade2((λ2 * b - bg / 2) * inv(λ2 - λ1))
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
    # Positive-definite algebras with no grade-4 part (VGA(2) ≅ ℂ, VGA(3) ≅ ℍ) admit the cheap,
    # exact closed form sqrt(R) = (1 + R) / |1 + R|.
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
