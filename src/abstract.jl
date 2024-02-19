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
    numeric_type(::Type{<:AbstractCliffordNumber{Q,T}}) = T
    numeric_type(T::Type{<:Union{Real,Complex}}) = T
    numeric_type(x) = numeric_type(typeof(x))

Returns the numeric type associated with an `AbstractCliffordNumber` instance. For subtypes of
`Real` and `Complex`, or their instances, this simply returns the input type or instance type. For 
incompletely instantiated types lacking information about the backing tuple, `Base.Bottom` is
returned.

# Why not define `eltype`?

`AbstractCliffordNumber` instances behave like numbers, not arrays. If `collect()` is called on a
Clifford number of type `T`, it should not construct a vector of coefficients; instead it should
return an `Array{T,0}`. Similarly, a broadcasted multiplication should return the same result as
normal multiplication, as is the case with complex numbers.

For subtypes `T` of `Number`, `eltype(T) === T`, and this is true for `AbstractCliffordNumber`.
"""
numeric_type(::Type) = BaseNumber
numeric_type(::Type{<:AbstractCliffordNumber{Q,T}}) where {Q,T} = T
numeric_type(T::Type{<:BaseNumber}) = T
numeric_type(x) = numeric_type(typeof(x))

#---Construct similar types------------------------------------------------------------------------#
"""
    CliffordNumbers.similar_type(
        C::Type{<:AbstractCliffordNumber},
        [N::Type{<:BaseNumber} = numeric_type(C)],
        [Q::Type{<:QuadraticForm} = QuadraticForm(C)]
    ) -> Type{<:AbstractCliffordNumber{Q,N}}

Constructs a type similar to `T` but with numeric type `N` and quadratic form `Q`.

This function must be defined with all its arguments for each concrete type subtyping
`AbstractCliffordNumber`.
"""
function similar_type(x::AbstractCliffordNumber, T::Type{<:BaseNumber}, Q::Type{<:QuadraticForm})
    return similar_type(typeof(x), T, Q)
end

similar_type(x, T::Type{<:BaseNumber}) = similar_type(x, T, QuadraticForm(x))
similar_type(x, Q::Type{<:QuadraticForm}) = similar_type(x, numeric_type(x), Q)

function similar(x, T::Type{<:BaseNumber}, Q::Type{<:QuadraticForm})
    return zero(similar_type(x))
end

similar(C::Type{<:AbstractCliffordNumber}, args...) = zero(similar_type(C, args...))
similar(x::AbstractCliffordNumber, args...) = zero(similar_type(C, args...))

#---Default constructor for zeros------------------------------------------------------------------#
import Base: zero

zero(C::Type{<:AbstractCliffordNumber{Q,T}}) where {Q,T} = C(_ -> ntuple(zero(T), Val(length(C))))
zero(C::Type{<:AbstractCliffordNumber}) = C(ntuple(_ -> zero(Bool), Val(length(C))))
zero(x::AbstractCliffordNumber) = zero(typeof(x))
