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
abstract type AbstractCliffordNumber{Q,T<:BaseNumber} <: Number
end

#---Default varargs constructors for types---------------------------------------------------------#

(::Type{T})(x::Vararg{BaseNumber}) where {Q,T<:AbstractCliffordNumber{Q}} = T(x)

#---Get type parameters----------------------------------------------------------------------------#

QuadraticForm(::Type{<:AbstractCliffordNumber{Q}}) where Q = Q
QuadraticForm(::AbstractCliffordNumber{Q}) where Q = Q

"""
    signature(T::Type{<:AbstractCliffordNumber{Q}}) = Q
    signature(x::AbstractCliffordNumber{Q}) = Q

Returns the metric signature object associated with an `AbstractCliffordNumber` `x` or its type `T`.
"""
signature(::Type{<:AbstractCliffordNumber{Q}}) where Q = Q
signature(::AbstractCliffordNumber{Q}) where Q = Q

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

#---Additive and multiplicative identities---------------------------------------------------------#
"""
    CliffordNumbers.zero_tuple(::Type{T}, ::Val{L}) -> NTuple{L,T}

Generates a `Tuple` of length `L` with all elements being `zero(T)`.
"""
zero_tuple(::Type{T}, ::Val{L}) where {T,L} = ntuple(Returns(zero(T)), Val(L))

"""
    CliffordNumbers.zero_tuple(::Type{C<:AbstractCliffordNumber})
        -> NTuple{length(C),numeric_type(C)}

Generates a `Tuple` that can be used to construct `zero(C)`.
"""
zero_tuple(::Type{C}) where C<:AbstractCliffordNumber = zero_tuple(numeric_type(C), Val(length(C)))

zero(C::Type{<:AbstractCliffordNumber{Q,T}}) where {Q,T} = C(zero_tuple(C))
zero(C::Type{<:AbstractCliffordNumber}) = C(zero_tuple(Bool, Val(length(C))))
zero(x::AbstractCliffordNumber) = zero(typeof(x))

# The default defintion assumes oneunit(T) = T(one(x))
# But this doesn't work here, because T(one(x)) doesn't do any error checking
# Only the explicit conversion does the error checking
oneunit(::Union{T,Type{T}}) where T<:AbstractCliffordNumber = convert(T, one(T))

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

similar(C::Type{<:AbstractCliffordNumber}, args...) = zero(similar_type(C, args...))
similar(x::AbstractCliffordNumber, args...) = zero(similar_type(x, args...))

#---real() and complex()---------------------------------------------------------------------------#
"""
    real(x::AbstractCliffordNumber{Q,T})

Gets the real portion of each coefficient of `x`. For `T<:Real` this operation does nothing; for
`T<:Complex{S}` this an `AbstractCliffordNumber{Q,S}`.

Note that this does not return the scalar (grade 0) coefficient of `x`. Use `real(scalar(x))` to
obtain this result.
"""
function Base.real(x::AbstractCliffordNumber)
    T = similar_type(x, real(numeric_type(x)))
    return T(real.(Tuple(x)))
end

Base.real(x::AbstractCliffordNumber{<:Any,<:Real}) = x

"""
    complex(x::AbstractCliffordNumber, [y::AbstractCliffordNumber = zero(typeof(x))])

For a single argument `x`, converts the type of each coefficient to a suitable complex type.

For two arguments `x` and `y`, which are both real Clifford numbers, performs the sum `x + y*im`,
constructing a complex Clifford number.

Note that this operation does not isolate a scalar (grade 0) coefficient of `x` or `y`. Use
`complex(scalar(x), [scalar(y)])` to obtain this result.
"""
function Base.complex(x::AbstractCliffordNumber)
    C = similar_type(x, complex(numeric_type(x)))
    return C(complex.(Tuple(x)))
end

Base.complex(x::AbstractCliffordNumber{<:Any,<:Complex}) = x

function Base.complex(x::T, y::T) where T<:AbstractCliffordNumber{<:QuadraticForm,<:Real}
    C = similar_type(x, complex(numeric_type(x)))
    return C(complex.(Tuple(x), Tuple(y)))
end

Base.complex(x::AbstractCliffordNumber, y::AbstractCliffordNumber) = complex(promote(x, y)...)

#---Error checking---------------------------------------------------------------------------------#
"""
    CliffordNumbers.check_element_count(sz, Q::Type{<:QuadraticForm}, [L], data)

Ensures that the number of elements in `data` is the same as the result of `f(Q)`, where `f` is a
function that generates the expected number of elements for the type. This function is used in the
inner constructors of subtypes of `AbstractCliffordNumber{Q}` to ensure that the input has the
correct length.

If provided, the length type parameter `L` can be included as an argument, and it will be checked
for type (must be an `Int`) and value (must be equal to `sz`).

This function returns nothing, but throws an `AssertionError` for failed checks.
"""
@inline function check_element_count(sz, data)
    @assert length(data) == sz "Expected $sz scalars from input, got $(length(data))"
    return nothing
end

@inline function check_element_count(sz, L, data)
    @assert L isa Int "Length type parameter must be an Int (got $(typeof(L)))."
    @assert L == sz "Length type parameter must equal $sz (got $L)."
    check_element_count(sz, data)
end

#---Printed representations------------------------------------------------------------------------#
"""
    CliffordNumbers.short_typename(T::Type{<:AbstractCliffordNumber})
    CliffordNumbers.short_typename(x::AbstractCliffordNumber})
    
Returns a type with a shorter name than `T`, but still constructs `T`. This is achieved by removing
dependent type parameters; often this includes the length parameter.
"""
short_typename(x::AbstractCliffordNumber) = short_typename(typeof(x))

function show(io::IO, x::AbstractCliffordNumber)
    print(io, short_typename(x), (isone(length(x)) ? string('(', only(Tuple(x)), ')') : Tuple(x)))
end

function summary(io::IO, x::AbstractCliffordNumber)
    println(io, length(x), "-element ", short_typename(x), ":")
end
