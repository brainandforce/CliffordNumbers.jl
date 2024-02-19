#---Abstract type for all Clifford numbers---------------------------------------------------------#
"""
    AbstractCliffordNumber{Q,T} <: Number

An element of a Clifford algebra, often referred to as a multivector, with quadratic form `Q` and
element type `T`.
"""
abstract type AbstractCliffordNumber{Q<:QuadraticForm,T<:BaseNumber} <: Number
end

#---Get type parameters----------------------------------------------------------------------------#

QuadraticForm(::Type{<:AbstractCliffordNumber{Q}}) where Q = Q
QuadraticForm(::AbstractCliffordNumber{Q}) where Q = Q

"""
    numbertype(::Type{<:AbstractCliffordNumber{Q,T}}) = T
    numbertype(::Type{<:AbstractCliffordNumber}) = Union{Real,Complex}

Returns the numeric type associated with an `AbstractCliffordNumber` instance. 

# Why not define `eltype`?

`AbstractCliffordNumber` instances behave like numbers, not arrays. If `collect()` is called on a
Clifford number of type `T`, it should not construct a vector of coefficients; instead it should
return an `Array{T,0}`. Similarly,

For subtypes `T` of `Number`, `eltype(T) === T`, and this is true for `AbstractCliffordNumber`.
"""
numbertype(::Type{<:AbstractCliffordNumber{Q,T}}) where {Q,T} = T
numbertype(::Type{<:AbstractCliffordNumber}) = Union{Real,Complex}
