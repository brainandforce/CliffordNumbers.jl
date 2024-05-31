# Default conversion should check for exact representability
function convert(T::Type{<:AbstractCliffordNumber}, x::AbstractCliffordNumber)
    result = T(x)::T
    return (has_grades_of(x, result) ? result : throw(InexactError(:convert, T, x)))
end

function convert(::Type{T}, x::AbstractCliffordNumber) where T<:BaseNumber
    return isscalar(x) ? T(scalar(x)) : throw(InexactError(:convert, T, x))
end

#---Specialized conversion methods for certain representations and signatures----------------------#

#= TODO: fix this once metric signature interface stabilizes
function convert(::Type{T}, x::AbstractCliffordNumber{QFComplex,<:Real}) where T<:BaseNumber
    return convert(T, x[scalar_index(x)] + x[pseudoscalar_index(x)] * im)
end

function convert(::Type{T}, z::Complex) where T<:AbstractCliffordNumber{QFComplex,<:Real}
    return convert(T, CliffordNumber{QFComplex}(real(z), imag(z)))
end
=#

# k-vectors of grade 0 are scalars
convert(::Type{T}, k::KVector{0}) where T<:BaseNumber = convert(T, only(k.data))

#---Convert only the scalar portion of an AbstractCliffordNumber-----------------------------------#

float(::Type{C}) where C<:AbstractCliffordNumber = similar_type(C, float(scalar_type(C)))
float(x::AbstractCliffordNumber) = convert(float(typeof(x)), x)

big(::Type{C}) where C<:AbstractCliffordNumber = similar_type(C, big(scalar_type(C)))
big(x::AbstractCliffordNumber) = convert(big(typeof(x)), x)

"""
    scalar_convert(T::Type{<:Union{Real,Complex}}, x::AbstractCliffordNumber) -> T
    scalar_convert(T::Type{<:Union{Real,Complex}}, x::Union{Real,Complex}) -> T

If `x` is an `AbstractCliffordNumber`, converts the scalars of `x` to type `T`.

If `x` is a `Real` or `Complex`, converts `x` to `T`.

# Examples
```julia-repl
julia> scalar_convert(Float32, KVector{1,APS}(1, 2, 3))
3-element KVector{1, VGA(3), Float32}:
1.0σ₁ + 2.0σ₂ + 3.0σ

julia> scalar_convert(Float32, 2)
2.0f0
```
"""
scalar_convert(::Type{T}, x::AbstractCliffordNumber) where T<:BaseNumber = similar_type(x, T)(x)
scalar_convert(::Type{T}, x::AbstractCliffordNumber{<:Any,T}) where T<:BaseNumber = x

scalar_convert(::Type{T}, x::BaseNumber) where T<:BaseNumber = convert(T, x)
scalar_convert(::Type{T}, x::T) where T<:BaseNumber = x
