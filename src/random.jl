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
    randn([rng=default_rng()], ::Type{KVector{1,Q,[T = Float64]}})

Generates a random 1-blade sampled from a radially symmetric multivariate normal distribution.

In a positive-definite metric (VGA), the 1-blade coefficients are drawn from a standard normal 
distribution (μ = 0, σ = 1). Their norms are distributed according to the [chi distribution], and
their moduli distributed according to the [chi-squared distribution].

In a projective geometric algebra (PGA), the resulting 1-blades have their projective component
fixed to `one(T)`, and the remaining components are drawn from a standard normal distribution.

[chi distribution]:         https://en.wikipedia.org/wiki/Chi_distribution
[chi-squared distribution]: https://en.wikipedia.org/wiki/Chi-squared_distribution
"""
randn(rng::AbstractRNG, K::Type{KVector{1,Q,T}}) where {Q,T} =  randn_by_algebra(rng, Q, K)

# Default to Float64 scalars
randn(rng::AbstractRNG, ::Type{KVector{1,Q}}) where Q = randn(rng, KVector{1,Q,Float64})

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
    return KVector{1,Q,T}(ntuple(_ -> randn(rng, T), Val(dimension(Q))))
end

function randn_by_algebra(rng::AbstractRNG, ::PGA, K::Type{<:KVector{1,Q,T}}) where {Q,T}
    Q isa PGA || error(LazyString("Algebra ", Q, " is not a PGA."))
    return KVector{1,Q,T}(one(T), ntuple(_ -> randn(rng, T), Val(dimension(Q) - 1))...)
end

function randn_by_algebra(rng::AbstractRNG, ::CGA, K::Type{<:KVector{1,Q,T}}) where {Q,T}
    Q isa CGA || error(LazyString("Algebra ", Q, " is not a CGA."))
    coords = ntuple(_ -> randn(rng, T), Val(dimension(Q) - 2))
    sq = sum(abs2, coords)
    return KVector{1,Q,T}(1//2 * (sq + 1), 1//2 * (sq - 1), coords...)
end

function randn_by_algebra(::AbstractRNG, ::Any, ::Type{<:KVector{1,Q,T}}) where {Q,T}
    error("Custom algebras require their own implementation of CliffordNumbers.randn_by_algebra.")
end

#---Functions for other distributions not in Base.Random-------------------------------------------#

#=
function rand_by_algebra(
    rand_f,
    algebra::Metrics.Signature
)
=#
