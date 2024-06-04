#---Abstract type for all Clifford numbers---------------------------------------------------------#
"""
    AbstractCliffordNumber{Q,T} <: Number

An element of a Clifford algebra, often referred to as a multivector, with quadratic form `Q` and
element type `T`. These are statically size and therefore should be able to be stored inline in
arrays or other data structures.

# Interface

## Required implementation

All subtypes `C` of `AbstractCliffordNumber{Q}` must implement the following functions:
  * `CliffordNumbers.similar_type(::Type{C}, ::Type{T}, ::Type{Q}) where {C,T,Q}` should construct a
new type similar to `C` which subtypes `AbstractCliffordNumber{Q,T}` that may serve as a
constructor.
  * `Base.getindex(x::C, b::BitIndex{Q})` should allow one to recover the coefficients associated
with each basis blade represented by `C`.
  * `nblades(::Type{C})` should be defined to return the number of basis blades represented by the
type. By default, `nblades(x::AbstractCliffordNumber) = nblades(typeof(x))`.
  * `Base.Tuple(x::C)` should return the tuple used to construct `x`. The fallback is
`getfield(x, :data)::Tuple`, so any type declared with a `NTuple` field named `data` should have
this defined automatically.
"""
abstract type AbstractCliffordNumber{Q,T<:BaseNumber} <: Number
end

#---Default varargs constructors for types---------------------------------------------------------#

(::Type{T})(x::Vararg{BaseNumber}) where {Q,T<:AbstractCliffordNumber{Q}} = T(x)

#---Number of blades-------------------------------------------------------------------------------#
"""
    nblades(::Type{<:Number}) -> Int
    nblades(x::Number)

Returns the number of blades represented by a `Number` subtype or instance. For subtypes of `Number`
that are not `AbstractCliffordNumber`, this is always 1.

This function is separate from `Base.length` since `AbstractCliffordNumber` is a scalar type for
which `collect()` returns a zero-dimensional array. For consistency, `length(x)` should always equal
`length(collect(x))`.
"""
nblades(::Type{<:Number}) = 1
nblades(x::Number) = nblades(typeof(x))

#---Get type parameters----------------------------------------------------------------------------#
"""
    signature(T::Type{<:AbstractCliffordNumber{Q}}) = Q
    signature(x::AbstractCliffordNumber{Q}) = Q

Returns the metric signature object associated with an `AbstractCliffordNumber` `x` or its type `T`.
"""
signature(::Type{<:AbstractCliffordNumber{Q}}) where Q = Q
signature(::AbstractCliffordNumber{Q}) where Q = Q

"""
    scalar_type(::Type{<:AbstractCliffordNumber{Q,T}}) = T
    scalar_type(T::Type{<:Union{Real,Complex}}) = T
    scalar_type(x) = scalar_type(typeof(x))

Returns the numeric type associated with an `AbstractCliffordNumber` instance. For subtypes of
`Real` and `Complex`, or their instances, this simply returns the input type or instance type.

# Why not define `eltype`?

`AbstractCliffordNumber` instances behave like numbers, not arrays. If `collect()` is called on a
Clifford number of type `T`, it should not construct a vector of coefficients; instead it should
return an `Array{T,0}`. Similarly, a broadcasted multiplication should return the same result as
normal multiplication, as is the case with complex numbers.

For subtypes `T` of `Number`, `eltype(T) === T`, and this is true for `AbstractCliffordNumber`.
"""
scalar_type(::Type) = BaseNumber
scalar_type(::Type{<:AbstractCliffordNumber{Q,T}}) where {Q,T} = T
scalar_type(T::Type{<:BaseNumber}) = T
scalar_type(x) = scalar_type(typeof(x))

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
        -> NTuple{nblades(C),scalar_type(C)}

Generates a `Tuple` that can be used to construct `zero(C)`.
"""
zero_tuple(::Type{C}) where C<:AbstractCliffordNumber = zero_tuple(scalar_type(C), Val(nblades(C)))

zero(::Type{C}) where C<:AbstractCliffordNumber = C(zero_tuple(Bool, Val(nblades(C))))
zero(x::AbstractCliffordNumber) = zero(typeof(x))

# The default defintion assumes oneunit(T) = T(one(x))
# But this doesn't work here, because T(one(x)) doesn't do any error checking
# Only the explicit conversion does the error checking
oneunit(::Union{T,Type{T}}) where T<:AbstractCliffordNumber = convert(T, one(T))

#---Construct similar types------------------------------------------------------------------------#
"""
    CliffordNumbers.similar_type(
        C::Type{<:AbstractCliffordNumber},
        [N::Type{<:BaseNumber} = scalar_type(C)],
        [Q::Val = Val(signature(C))]
    ) -> Type{<:AbstractCliffordNumber{Q,N}}

Constructs a type similar to `T` but with numeric type `N` and quadratic form `Q`. The quadratic
form must be wrapped in a `Val` to preserve type stability.

This function must be defined with all its arguments for each concrete type subtyping
`AbstractCliffordNumber`.
"""
function similar_type(x::AbstractCliffordNumber, T::Type{<:BaseNumber}, Q::Val)
    return similar_type(typeof(x), T, Q)
end

similar_type(x, T::Type{<:BaseNumber}) = similar_type(x, T, Val(signature(x)))
similar_type(x, Q::Val) = similar_type(x, scalar_type(x), Q)

similar(C::Type{<:AbstractCliffordNumber}, args...) = zero(similar_type(C, args...))
similar(x::AbstractCliffordNumber, args...) = zero(similar_type(x, args...))

"""
    CliffordNumbers.complement_type(C::Type{<:AbstractCliffordNumber})
    CliffordNumbers.complement_type(x::AbstractCliffordNumber)

Constructs a type capable of storing the complementary grades of `x`. For the types provided by this
package:
  * The complement type of `KVector{K,Q,T}` is `KVector{dimension(Q)-K,Q,T}`.
  * The complement type of `EvenCliffordNumber{Q,T}` and `OddCliffordNumber{Q,T}` is the same type
    if `dimension(Q)` is even, or the type of opposite grade parity if `dimension(Q)` is odd.
  * All other subtypes of `AbstractCliffordNumber{Q,T}`, including `CliffordNumber{Q,T}`, have
    complement type `CliffordNumber{Q,T}` (unless defined otherwise).
"""
complement_type(x::AbstractCliffordNumber) = complement_type(typeof(x))

#---real() and complex()---------------------------------------------------------------------------#
"""
    real(x::AbstractCliffordNumber{Q,T})

Gets the real portion of each coefficient of `x`. For `T<:Real` this operation does nothing; for
`T<:Complex{S}` this an `AbstractCliffordNumber{Q,S}`.

Note that this does not return the scalar (grade 0) coefficient of `x`. Use `scalar(x)` to obtain
this result in general, or `real(scalar(x))` if only the real portion is desired.
"""
function Base.real(x::AbstractCliffordNumber)
    T = similar_type(x, real(scalar_type(x)))
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
    C = similar_type(x, complex(scalar_type(x)))
    return C(complex.(Tuple(x)))
end

Base.complex(x::AbstractCliffordNumber{<:Any,<:Complex}) = x

function Base.complex(x::T, y::T) where T<:AbstractCliffordNumber{<:Any,<:Real}
    C = similar_type(x, complex(scalar_type(x)))
    return C(complex.(Tuple(x), Tuple(y)))
end

Base.complex(x::AbstractCliffordNumber, y::AbstractCliffordNumber) = complex(promote(x, y)...)

#---Error checking---------------------------------------------------------------------------------#
"""
    CliffordNumbers.check_element_count(sz, [L], data)

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
    
Returns a type with a shorter name than `T`, but still constructs an instance of `T`. This is
achieved by removing dependent type parameters; often this includes the length parameter.
"""
short_typename(::Type{C}) where C<:AbstractCliffordNumber = C
short_typename(x::AbstractCliffordNumber) = short_typename(typeof(x))

function show(io::IO, x::AbstractCliffordNumber)
    print(io, short_typename(x), (isone(nblades(x)) ? string('(', only(Tuple(x)), ')') : Tuple(x)))
end

function summary(io::IO, x::AbstractCliffordNumber)
    println(io, nblades(x), "-element ", short_typename(x), ":")
end

#---Algebra mismatch errors------------------------------------------------------------------------#
"""
    CliffordNumbers.AlgebraMismatch(f, args::Tuple)

The arguments to function `f` expected all of it arguments `args` to be part of the same algebra.
"""
struct AlgebraMismatch <: Exception
    f::Any
    args::Tuple{Vararg{Any}}
end

function Base.showerror(io::IO, ex::AlgebraMismatch)
    print(io, "Arguments to ", ex.f, " were not Clifford numbers of the same algebra.")
    for (n,x) in enumerate(ex.args)
        if x isa AbstractCliffordNumber
            print(io, "\n  Argument ", n, " has signature parameter ", signature(x))
        end
    end
end
