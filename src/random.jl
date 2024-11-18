#---Random sampling of 0-blades passes through sampling of the scalar type-------------------------#

rand(rng::AbstractRNG, ::SamplerType{<:KVector{0,Q,T}}) where {Q,T} = KVector{0,Q}(rand(rng, T))
randn(rng::AbstractRNG, ::SamplerType{<:KVector{0,Q,T}}) where {Q,T} = KVector{0,Q}(randn(rng, T))

function randexp(rng::AbstractRNG, ::SamplerType{<:KVector{0,Q,T}}) where {Q,T}
    return KVector{0,Q}(randexp(rng, T))
end

#---Random sampling of unit 1-blades---------------------------------------------------------------#
"""
    randn([rng=default_rng()], ::Type{KVector{1,Q,T}})

Generates a random 1-blade whose coefficients are sampled from an `n`-spherically symmetric
multivariate normal distribution, where `n = dimension(Q)`.

In a positive-definite metric, the resulting 1-blades have their norms distributed according to
the [chi distribution], and their moduli distributed according to the [chi-squared distribution].

[chi distribution]:         https://en.wikipedia.org/wiki/Chi_distribution
[chi-squared distribution]: https://en.wikipedia.org/wiki/Chi-squared_distribution
"""
function randn(rng::AbstractRNG, ::SamplerType{<:KVector{1,Q,T}}) where {Q,T}
    return KVector{K,Q}(ntuple(_ -> randn(rng, T), Val(dimension(Q))))
end
