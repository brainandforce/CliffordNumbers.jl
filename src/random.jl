#---Random sampling of 0-blades passes through sampling of the scalar type-------------------------#

rand(rng::AbstractRNG, ::SamplerType{<:KVector{0,Q,T}}) where {Q,T} = KVector{0,Q}(rand(rng, T))
randn(rng::AbstractRNG, ::SamplerType{<:KVector{0,Q,T}}) where {Q,T} = KVector{0,Q}(randn(rng, T))

function randexp(rng::AbstractRNG, ::SamplerType{<:KVector{0,Q,T}}) where {Q,T}
    return KVector{0,Q}(randexp(rng, T))
end
