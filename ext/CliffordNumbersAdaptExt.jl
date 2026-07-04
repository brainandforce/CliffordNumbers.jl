module CliffordNumbersAdaptExt

using CliffordNumbers
using CliffordNumbers: BaseNumber
using Adapt

"""
    adapt_scalar_type(to, ::Type{T}) where T<:CliffordNumbers.BaseNumber

Returns the scalar type implied by an `Adapt` destination `to` for a Clifford number whose current
scalar type is `T`. When `to` is an array type carrying a concrete scalar element type (e.g.
`Array{Float32}` or `CuArray{Float32}`), the coefficients are retyped to match, mirroring how array
storage is retyped on adaptation. Otherwise the existing scalar type is preserved.
"""
adapt_scalar_type(to, ::Type{T}) where T<:BaseNumber = T
adapt_scalar_type(::Type{<:AbstractArray{U}}, ::Type{<:BaseNumber}) where U<:BaseNumber = U

# `AbstractCliffordNumber` subtypes are `isbits` with `NTuple{L,T}` storage, so there is no device
# pointer to swap: the only meaningful adaptation is retyping the scalar coefficients to match the
# destination's element type. The number itself moves to the device inside an adapted array.
function Adapt.adapt_structure(to, x::AbstractCliffordNumber)
    return scalar_convert(adapt_scalar_type(to, scalar_type(x)), x)
end

end
