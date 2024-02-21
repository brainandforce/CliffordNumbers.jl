#---Abstract type for all Clifford numbers---------------------------------------------------------#
"""
    AbstractCliffordNumber{Q,T} <: Number

An element of a Clifford algebra, often referred to as a multivector, with quadratic form `Q` and
element type `T`.

# Interface

## Required implementation

All subtypes `C` of `AbstractCliffordNumber{Q}` must implement the following functions:
  * `Base.length(x::C)` should return the number of nonzero basis elements represented by `x`.
  * `CliffordNumbers.similar_type(::Type{C}, ::Type{T}, ::Type{Q}) where {C,T,Q}` should construct a
new type similar to `C` which subtypes `AbstractCliffordNumber{Q,T}` that may serve as a
constructor.
  * `Base.getindex(x::C, b::BitIndex{Q})` should allow one to recover the coefficients associated
with each basis blade represented by `C`.

## Required implementation for static types
  * `Base.length(::Type{C})` should be defined, with `Base.length(x::C) = length(typeof(x))`.
  * `Base.Tuple(x::C)` should return the tuple used to construct `x`. The fallback is
`getfield(x, :data)::Tuple`, so any type declared with a `NTuple` field named `data` should have
this defined automatically.
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
`Real` and `Complex`, or their instances, this simply returns the input type or instance type.

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

#---Get underlying tuple---------------------------------------------------------------------------#

Base.Tuple(x::AbstractCliffordNumber) = getfield(x, :data)::Tuple

#---Zero multivectors------------------------------------------------------------------------------#
import Base: zero

zero(C::Type{<:AbstractCliffordNumber{Q,T}}) where {Q,T} = C(_ -> ntuple(zero(T), Val(length(C))))
zero(C::Type{<:AbstractCliffordNumber}) = C(ntuple(_ -> zero(Bool), Val(length(C))))
zero(x::AbstractCliffordNumber) = zero(typeof(x))

#---Construct similar types------------------------------------------------------------------------#
import Base.similar

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

similar(C::Type{<:AbstractCliffordNumber}, args...) = zero(similar_type(C, args...))
similar(x::AbstractCliffordNumber, args...) = zero(similar_type(x, args...))

#---Error checking---------------------------------------------------------------------------------#
"""
    CliffordNumbers.check_element_count(f, Q::Type{<:QuadraticForm}, data)

Ensures that the number of elements in `data` is the same as the result of `f(Q)`, where `f` is a
function that generates the expected number of elements for the type.

This function is used in the inner constructors of subtypes of `AbstractCliffordNumber{Q}` to ensure
that the input is the correct length.
"""
function check_element_count(f, Q::Type{QuadraticForm{X,Y,Z}}, data) where {X,Y,Z}
    @assert length(data) == f(Q) "Expected $(f(Q)) scalars from input, got $(length(data))"
    return nothing
end

"""
    CliffordNumbers.check_length_parameter(f, Q::Type{<:QuadraticForm}, L)

Ensures that the length type parameter `L` is a `Int` equal to `f(Q)`, where `f` where `f` is a
function that generates the expected number of elements for the type.

This function is used in the inner constructors of subtypes of `AbstractCliffordNumber{Q}` to ensure
that the input is the correct length.
"""
function check_length_parameter(f, Q::Type{QuadraticForm{X,Y,Z}}, L) where {X,Y,Z}
    @assert L isa Int || "Length type parameter must be an Int (got $(typeof(L)))."
    @assert L = f(Q) || "Length type parameter must equal $(f(Q)) (got $L)."
    return nothing
end
