#---Random sampling of 0-blades passes through sampling of the scalar type-------------------------#

rand(rng::AbstractRNG, ::SamplerType{<:KVector{0,Q,T}}) where {Q,T} = KVector{0,Q}(rand(rng, T))
randn(rng::AbstractRNG, ::SamplerType{<:KVector{0,Q,T}}) where {Q,T} = KVector{0,Q}(randn(rng, T))

function randexp(rng::AbstractRNG, ::SamplerType{<:KVector{0,Q,T}}) where {Q,T}
    return KVector{0,Q}(randexp(rng, T))
end

# Default to Float64 scalars
rand(rng::AbstractRNG, ::SamplerType{KVector{0,Q}}) where Q = KVector{0,Q}(rand(rng, Float64))
randn(rng::AbstractRNG, ::SamplerType{KVector{0,Q}}) where Q = KVector{0,Q}(randn(rng, Float64))
randexp(rng::AbstractRNG, ::SamplerType{KVector{0,Q}}) where Q = KVector{0,Q}(randexp(rng, Float64))

#---Random sampling of unit 1-blades---------------------------------------------------------------#
"""
    rand(rng::AbstractRNG, K::Type{<:KVector{0,Q,[T = Float64]}}

Generates a uniformly distributed random 1-blade sampled from a uniform distribution on the unit
sphere centered at the origin of the modeled space.

In a positive-definite metric (VGA), the 1-blade coefficients are drawn from a standard normal
distribution and the resulting vector is normalized so that it squares to 1.

In a projective geometric algebra (PGA), the e₀ component is set equal to 1 and the remaining
coefficients are drawn from a standard normal distribution, then normalized so that the result
squares to 1.

In a conformal geometric algebra (CGA), the n₀ component is set equal to 1, the coordinate
coefficients are drawn from a standard normal distribution and normalized so they square to 1, and
the n∞ component is chosen so that the result squares to 0.
"""
rand(rng::AbstractRNG, K::Type{<:KVector{1,Q,T}}) where {Q,T} = rand_by_algebra(rng, Q, K)

# Default to Float64 scalars
rand(rng::AbstractRNG, ::Type{KVector{1,Q}}) where Q = rand(rng, KVector{1,Q,Float64})

"""
    randn([rng=default_rng()], ::Type{KVector{1,Q,[T = Float64]}})

Generates a random 1-blade sampled from a radially symmetric multivariate normal distribution.

In a positive-definite metric (VGA), the 1-blade coefficients are drawn from a standard normal 
distribution (μ = 0, σ = 1). Their norms are distributed according to the [chi distribution], and
their moduli distributed according to the [chi-squared distribution].

In a projective geometric algebra (PGA), the resulting 1-blades have their projective component
fixed to `one(T)`, and the remaining components are drawn from a standard normal distribution.

In a conformal geometric algebra (CGA), the n₀ coeffient is `one(T)`, the coordinate coefficients
are drawn from a standard normal distribution, and the n∞ component is chosen so that the result
squares to 0.

[chi distribution]:         https://en.wikipedia.org/wiki/Chi_distribution
[chi-squared distribution]: https://en.wikipedia.org/wiki/Chi-squared_distribution
"""
randn(rng::AbstractRNG, K::Type{<:KVector{1,Q,T}}) where {Q,T} =  randn_by_algebra(rng, Q, K)

# Default to Float64 scalars
randn(rng::AbstractRNG, ::Type{KVector{1,Q}}) where Q = randn(rng, KVector{1,Q,Float64})

#---Defining random number generation by algebra---------------------------------------------------#
"""
    CliffordNumbers.randn_coordinate(rng::AbstractRNG, ::Type{NTuple{N,T}}; normalize = false)

Generates a coordinate in `N` dimensions by drawing each of its coefficients from a standard
multivariate normal distribution.
"""
function randn_coordinate(
    rng::AbstractRNG,
    ::Type{NTuple{N,T}}
) where {N,T<:AbstractFloat}
    return ntuple(_ -> randn(rng, T), Val(N))
end

# NOTE: this is one of the reasons why I think the algebra module needs to be reimplemented.
# Managing dispatch is way too complicated here.
"""
    CliffordNumbers.randn_by_algebra(
        ::Metrics.Signature,
        rng::AbstractRNG,
        ::Type{<:KVector{1,Q,T}}
    ) -> KVector{1,Q,T}

Generates a random 1-blade sampled from a multivariate normal distribution. For more information 
about the public-facing implementation, see [`randn`](@ref).

This internal function is needed because we cannot directly specialize on the type of the algebra
parameter.
"""
function randn_by_algebra(rng::AbstractRNG, ::VGA, K::Type{<:KVector{1,Q,T}}) where {Q,T}
    Q isa VGA || error(LazyString("Algebra ", Q, " is not a VGA."))
    return KVector{1,Q,T}(randn_coordinate(rng, NTuple{dimension(Q), T})...)
end

function randn_by_algebra(rng::AbstractRNG, ::PGA, K::Type{<:KVector{1,Q,T}}) where {Q,T}
    Q isa PGA || error(LazyString("Algebra ", Q, " is not a PGA."))
    return KVector{1,Q,T}(one(T), randn_coordinate(rng, NTuple{dimension(Q) - 1, T})...)
end

function randn_by_algebra(rng::AbstractRNG, ::CGA, K::Type{<:KVector{1,Q,T}}) where {Q,T}
    Q isa CGA || error(LazyString("Algebra ", Q, " is not a CGA."))
    data = randn_coordinate(rng, NTuple{dimension(Q) - 2, T})
    sq = sum(abs2, data)
    return KVector{1,Q,T}(1//2 * (sq + 1), 1//2 * (sq - 1), data...)
end

function randn_by_algebra(::AbstractRNG, ::Any, ::Type{<:KVector{1,Q,T}}) where {Q,T}
    error("Custom algebras require their own implementation of CliffordNumbers.randn_by_algebra.")
end

function rand_coordinate(
    rng::AbstractRNG,
    ::Type{NTuple{N,T}}
) where {N,T<:AbstractFloat}
    result = randn_coordinate(rng, NTuple{N,T})
    return result ./ hypot(result...)
end

"""
    CliffordNumbers.rand_by_algebra(
        ::Metrics.Signature,
        rng::AbstractRNG,
        ::Type{<:KVector{1,Q,T}}
    ) -> KVector{1,Q,T}

Generates a random 1-blade with unit distance from the origin. For more information about the
public-facing implementation, see [`rand`](@ref).

This internal function is needed because we cannot directly specialize on the type of the algebra
parameter.
"""
function rand_by_algebra(rng::AbstractRNG, ::VGA, K::Type{<:KVector{1,Q,T}}) where {Q,T}
    return normalize(randn(rng, K))
end

function rand_by_algebra(rng::AbstractRNG, ::PGA, K::Type{<:KVector{1,Q,T}}) where {Q,T}
    Q isa PGA || error(LazyString("Algebra ", Q, " is not a PGA."))
    return KVector{1,Q,T}(one(T), rand_coordinate(rng, NTuple{dimension(Q) - 1, T})...)
end

function rand_by_algebra(rng::AbstractRNG, ::CGA, K::Type{<:KVector{1,Q,T}}) where {Q,T}
    Q isa CGA || error(LazyString("Algebra ", Q, " is not a CGA."))
    return KVector{1,Q,T}(one(T), zero(T), rand_coordinate(rng, NTuple{dimension(Q) - 2, T})...)
end

function rand_by_algebra(::AbstractRNG, ::Any, ::Type{<:KVector{1,Q,T}}) where {Q,T}
    error("Custom algebras require their own implementation of CliffordNumbers.rand_by_algebra.")
end
