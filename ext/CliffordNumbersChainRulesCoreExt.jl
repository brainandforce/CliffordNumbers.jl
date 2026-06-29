module CliffordNumbersChainRulesCoreExt

using CliffordNumbers
using ChainRulesCore

# `rrule` and `ProjectTo` are extended here; the remaining names (`NoTangent`, `unthunk`,
# `AbstractZero`, `AbstractThunk`) are brought in by `using ChainRulesCore`.
import ChainRulesCore: rrule, ProjectTo

import CliffordNumbers: AbstractCliffordNumber, BaseNumber
import CliffordNumbers: versor_inverse, scalar_product
# Internal Cayley-product machinery, reused to build exact product pullbacks
import CliffordNumbers: BladeIndices, GradeFilter, bitindex_shuffle, mul_mask, mul_signs
import CliffordNumbers: to_index, nonzero_grades, nblades, scalar_type, similar_type, signature
import CliffordNumbers: sign_of_square

#---Cotangent space of a Clifford number-----------------------------------------------------------#
# The tangent/cotangent of an `AbstractCliffordNumber` is another Clifford number of the same type:
# these types are real (or complex) vector spaces with `NTuple` coefficient storage, and `+`, `-`,
# and scalar `*` are already defined. `ProjectTo` maps an arbitrary cotangent back onto the grade
# structure of the primal, dropping any blade coefficients the primal does not represent.

function ProjectTo(x::AbstractCliffordNumber)
    return ProjectTo{typeof(x)}()
end

(::ProjectTo{C})(dx::AbstractCliffordNumber) where C<:AbstractCliffordNumber = C(dx)
(::ProjectTo{C})(dx::BaseNumber) where C<:AbstractCliffordNumber = C(dx)
(::ProjectTo{C})(dx::AbstractZero) where C<:AbstractCliffordNumber = dx
(p::ProjectTo{C})(dx::AbstractThunk) where C<:AbstractCliffordNumber = p(unthunk(dx))

# Convert an incoming cotangent into the concrete Clifford type `C` (the primal output type).
_cotangent(::Type{C}, dΩ::AbstractCliffordNumber) where C<:AbstractCliffordNumber = C(dΩ)
_cotangent(::Type{C}, dΩ::BaseNumber) where C<:AbstractCliffordNumber = C(dΩ)

#---Exact transpose of the geometric product-------------------------------------------------------#
# For `Ω = x * y`, the geometric product is bilinear: `Ω[c] = Σ_a x[a] * ε(a,c) * y[idx(a,c)]`,
# where `ε` and `idx` are exactly the signs and shuffled indices the `mul` kernel generates. The
# Euclidean (coefficient-wise) pullbacks are the transposes of that bilinear map:
#
#     x̄[a] = Σ_c Ω̄[c] * ε(a,c) * y[idx(a,c)]
#     ȳ[b] = Σ_c Ω̄[c] * ε(b,c) * x[idx(b,c)]
#
# Building these directly from the Cayley table (rather than from a `reverse`-based identity) keeps
# them exact in every metric signature, including degenerate algebras (PGA, Exterior), where the
# reverse identity `⟨x*y, w⟩ = ⟨y, x̃*w⟩` does not account for null blades. The kernels mirror the
# two branches of `CliffordNumbers.mul`, so all index/sign resolution happens at compile time.

# Gradient with respect to the left operand `x` of `x * y`; mirrors the `mul` "else" branch.
@generated function _geometric_grad_left(
    Ω̄::AbstractCliffordNumber{Q},
    ::Type{X},
    y::AbstractCliffordNumber{Q},
) where {Q,X<:AbstractCliffordNumber{Q}}
    F = GradeFilter{:*}()
    Tz = promote_type(scalar_type(Ω̄), scalar_type(y))
    OutT = similar_type(X, Tz, Val(signature(X)))
    exprs = Any[:($(zero(Tz))) for _ in 1:nblades(X)]
    for a in BladeIndices(X)
        ia = to_index(X, a)
        inds = bitindex_shuffle(a, BladeIndices(Ω̄))
        x_mask = mul_mask(F, a, inds)
        y_mask = map(in, grade.(inds), ntuple(Returns(nonzero_grades(y)), Val(nblades(Ω̄))))
        mask = x_mask .& y_mask
        signs = mul_signs(F, a, inds)
        coeffs = signs .* mask
        terms = Any[]
        for c in 1:nblades(Ω̄)
            iszero(coeffs[c]) && continue
            iy = to_index(y, inds[c])
            push!(terms, :(Tuple(Ω̄)[$c] * $(coeffs[c]) * Tuple(y)[$iy]))
        end
        isempty(terms) || (exprs[ia] = length(terms) == 1 ? terms[1] : Expr(:call, :+, terms...))
    end
    return :($OutT(($(exprs...),)))
end

# Gradient with respect to the right operand `y` of `x * y`; mirrors the `mul` "if" branch.
@generated function _geometric_grad_right(
    x::AbstractCliffordNumber{Q},
    Ω̄::AbstractCliffordNumber{Q},
    ::Type{Y},
) where {Q,Y<:AbstractCliffordNumber{Q}}
    F = GradeFilter{:*}()
    Tz = promote_type(scalar_type(x), scalar_type(Ω̄))
    OutT = similar_type(Y, Tz, Val(signature(Y)))
    exprs = Any[:($(zero(Tz))) for _ in 1:nblades(Y)]
    for b in BladeIndices(Y)
        ib = to_index(Y, b)
        inds = bitindex_shuffle(BladeIndices(Ω̄), b)
        x_mask = map(in, grade.(inds), ntuple(Returns(nonzero_grades(x)), Val(nblades(Ω̄))))
        y_mask = mul_mask(F, inds, b)
        mask = x_mask .& y_mask
        signs = mul_signs(F, inds, b)
        coeffs = signs .* mask
        terms = Any[]
        for c in 1:nblades(Ω̄)
            iszero(coeffs[c]) && continue
            ix = to_index(x, inds[c])
            push!(terms, :(Tuple(Ω̄)[$c] * $(coeffs[c]) * Tuple(x)[$ix]))
        end
        isempty(terms) || (exprs[ib] = length(terms) == 1 ? terms[1] : Expr(:call, :+, terms...))
    end
    return :($OutT(($(exprs...),)))
end

#---Metric (coefficient) scaling by `sign_of_square`-----------------------------------------------#
# `_metric_scale(x)[i] = sign_of_square(blade i) * x[i]`. This is the diagonal map that turns the
# Euclidean coefficient inner product into the metric scalar product, and appears in the pullbacks
# of `scalar_product` and `abs2`, which intrinsically carry the metric.

@generated function _metric_scale(x::AbstractCliffordNumber)
    signs = map(sign_of_square, BladeIndices(x))
    return :(typeof(x)(Tuple(x) .* $signs))
end

# Plain (metric-free) coefficient dot product of two same-typed Clifford numbers.
_coeff_dot(a::C, b::C) where C<:AbstractCliffordNumber =
    mapreduce((u, v) -> u * conj(v), +, Tuple(a), Tuple(b); init = zero(promote_type(scalar_type(a), scalar_type(b))))

#---rrule: geometric product-----------------------------------------------------------------------#
"""
    rrule(*, x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Reverse-mode rule for the geometric product. The pullback is the exact transpose of the bilinear
Cayley product, so it is correct in every metric signature (including degenerate algebras). For a
positive-definite algebra it reduces to `Ω̄ -> (Ω̄ * y', x' * Ω̄)`.
"""
function rrule(
    ::typeof(*),
    x::AbstractCliffordNumber{Q},
    y::AbstractCliffordNumber{Q},
) where Q
    Ω = x * y
    C = typeof(Ω)
    X = typeof(x)
    Y = typeof(y)
    function times_pullback(dΩ)
        dΩ isa AbstractZero && return (NoTangent(), NoTangent(), NoTangent())
        Ω̄ = _cotangent(C, unthunk(dΩ))
        x̄ = _geometric_grad_left(Ω̄, X, y)
        ȳ = _geometric_grad_right(x, Ω̄, Y)
        return (NoTangent(), x̄, ȳ)
    end
    return Ω, times_pullback
end

#---rrule: scalar multiplication and division------------------------------------------------------#
function rrule(::typeof(*), x::AbstractCliffordNumber, a::BaseNumber)
    project_x = ProjectTo(x)
    project_a = ProjectTo(a)
    function scalar_mul_pullback(dΩ)
        dΩ isa AbstractZero && return (NoTangent(), NoTangent(), NoTangent())
        Ω̄ = unthunk(dΩ)
        x̄ = project_x(Ω̄ * conj(a))
        ā = project_a(_coeff_dot(Ω̄, x))
        return (NoTangent(), x̄, ā)
    end
    return x * a, scalar_mul_pullback
end

function rrule(::typeof(*), a::BaseNumber, x::AbstractCliffordNumber)
    project_x = ProjectTo(x)
    project_a = ProjectTo(a)
    function scalar_mul_pullback(dΩ)
        dΩ isa AbstractZero && return (NoTangent(), NoTangent(), NoTangent())
        Ω̄ = unthunk(dΩ)
        ā = project_a(_coeff_dot(Ω̄, x))
        x̄ = project_x(conj(a) * Ω̄)
        return (NoTangent(), ā, x̄)
    end
    return a * x, scalar_mul_pullback
end

function rrule(::typeof(/), x::AbstractCliffordNumber, a::BaseNumber)
    Ω = x / a
    project_x = ProjectTo(x)
    project_a = ProjectTo(a)
    function scalar_div_pullback(dΩ)
        dΩ isa AbstractZero && return (NoTangent(), NoTangent(), NoTangent())
        Ω̄ = unthunk(dΩ)
        x̄ = project_x(Ω̄ / conj(a))
        ā = project_a(-_coeff_dot(Ω̄, Ω) / conj(a))
        return (NoTangent(), x̄, ā)
    end
    return Ω, scalar_div_pullback
end

#---rrule: linear sign-changing automorphisms------------------------------------------------------#
# `reverse`, `adjoint`, `conj`, and `grade_involution` are diagonal sign flips in coefficient space,
# hence self-adjoint: the pullback applies the same automorphism. (For real scalar types `adjoint`
# coincides with `reverse`.)
for f in (:reverse, :adjoint, :conj, :grade_involution)
    @eval function rrule(::typeof($f), x::AbstractCliffordNumber)
        project_x = ProjectTo(x)
        function automorphism_pullback(dΩ)
            dΩ isa AbstractZero && return (NoTangent(), NoTangent())
            return (NoTangent(), project_x($f(unthunk(dΩ))))
        end
        return $f(x), automorphism_pullback
    end
end

#---rrule: scalar product and squared norm---------------------------------------------------------#
function rrule(
    ::typeof(scalar_product),
    x::AbstractCliffordNumber{Q},
    y::AbstractCliffordNumber{Q},
) where Q
    project_x = ProjectTo(x)
    project_y = ProjectTo(y)
    function scalar_product_pullback(ds)
        ds isa AbstractZero && return (NoTangent(), NoTangent(), NoTangent())
        s̄ = unthunk(ds)
        x̄ = project_x(s̄ * _metric_scale(y))
        ȳ = project_y(s̄ * _metric_scale(x))
        return (NoTangent(), x̄, ȳ)
    end
    return scalar_product(x, y), scalar_product_pullback
end

function rrule(::typeof(abs2), x::AbstractCliffordNumber)
    project_x = ProjectTo(x)
    function abs2_pullback(dc)
        dc isa AbstractZero && return (NoTangent(), NoTangent())
        c̄ = unthunk(dc)
        # d(abs2)/dx = 2 * sign_of_square ⊙ reverse(x)
        return (NoTangent(), project_x((2 * c̄) * _metric_scale(reverse(x))))
    end
    return abs2(x), abs2_pullback
end

#---rrule: additive structure----------------------------------------------------------------------#
function rrule(
    ::typeof(+),
    x::AbstractCliffordNumber{Q},
    y::AbstractCliffordNumber{Q},
) where Q
    project_x = ProjectTo(x)
    project_y = ProjectTo(y)
    function plus_pullback(dΩ)
        dΩ isa AbstractZero && return (NoTangent(), NoTangent(), NoTangent())
        Ω̄ = unthunk(dΩ)
        return (NoTangent(), project_x(Ω̄), project_y(Ω̄))
    end
    return x + y, plus_pullback
end

function rrule(
    ::typeof(-),
    x::AbstractCliffordNumber{Q},
    y::AbstractCliffordNumber{Q},
) where Q
    project_x = ProjectTo(x)
    project_y = ProjectTo(y)
    function minus_pullback(dΩ)
        dΩ isa AbstractZero && return (NoTangent(), NoTangent(), NoTangent())
        Ω̄ = unthunk(dΩ)
        return (NoTangent(), project_x(Ω̄), project_y(-Ω̄))
    end
    return x - y, minus_pullback
end

function rrule(::typeof(-), x::AbstractCliffordNumber)
    project_x = ProjectTo(x)
    function negate_pullback(dΩ)
        dΩ isa AbstractZero && return (NoTangent(), NoTangent())
        return (NoTangent(), project_x(-unthunk(dΩ)))
    end
    return -x, negate_pullback
end

# `versor_inverse(x) = x' / abs2(x)` differentiates by composition through the `reverse`, `abs2`,
# and scalar-`/` rules above, so it needs no dedicated rule. The sandwich product `x * p * x'` (and
# hence `rotate`) likewise composes from the product and `adjoint`/`reverse` rules.

end
