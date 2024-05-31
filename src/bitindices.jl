#---Sets of `BitIndex` objects that index specific types-------------------------------------------#
"""
    AbstractBitIndices{Q,C<:AbstractCliffordNumber{Q}} <: AbstractVector{BitIndex{Q}}

Supertype for vectors containing all valid `BitIndex{Q}` objects for the basis elements represented
by `C`.
"""
abstract type AbstractBitIndices{Q,C<:AbstractCliffordNumber{Q}} <: AbstractVector{BitIndex{Q}}
end

size(::Type{<:AbstractBitIndices{Q,C}}) where {Q,C} = tuple(nblades(C))
size(b::AbstractBitIndices) = size(typeof(b))

length(::Type{<:AbstractBitIndices{Q,C}}) where {Q,C} = nblades(C)
length(::T) where T<:AbstractBitIndices = length(T)

# Conversion to tuple
# Base.Tuple(b::T) where T<:AbstractBitIndices = ntuple(i -> b[i], Val(length(T)))

"""
    CliffordNumbers.bitindices_type(C::Type{<:AbstractCliffordNumber{Q,T}})

Removes extraneous type parameters from `C`, converting it to the least parameterized type that
can be used to parameterize an `AbstractBitIndices{Q,C}` object. This is to avoid issues with the
proliferation of type parameters that would construct identical `BitIndices` objects otherwise:
for instance, `BitIndices{VGA(3),EvenCliffordNumber{VGA(3),Float32,4}}()` and
`BitIndices{VGA(3),EvenCliffordNumber{VGA(3),Int}}()` have identical elements, and are equal when
compared with `==`, but are not the same object.

For types defined in this package, this strips the scalar type parameter `T` and any length
parameters present.

# Examples
```julia-repl
julia> CliffordNumbers.bitindices_type(CliffordNumber{VGA(3),Float32,8})
CliffordNumber{VGA(3)}

julia> CliffordNumbers.bitindices_type(KVector{2,STA,Bool})
KVector{2,STA}
```
"""
bitindices_type(::Type{AbstractCliffordNumber{Q,<:Any}}) where Q = AbstractCliffordNumber{Q}
bitindices_type(x::AbstractCliffordNumber) = bitindices_type(typeof(x))

#---Clifford number iteration----------------------------------------------------------------------#
"""
    BitIndices{Q,C<:AbstractCliffordNumber{Q,<:Any}} <: AbstractVector{BitIndex{Q}}

Represents a range of valid `BitIndex` objects for the nonzero components of a given multivector of
algebra `Q`.

For a generic `AbstractCliffordNumber{Q}`, this returns `BitIndices{Q,CliffordNumber{Q}}`, which
contains all possible indices for a multivector associated with the algebra parameter `Q`. For
sparser representations, such as `KVector{K,Q}`, the object only contains the indices of the
nonzero elements of the multivector.

# Construction

`BitIndices` can be constructed by calling the type constructor with the Clifford number or its
type. This will automatically strip some type parameters so that identical `BitIndices` objects are
constructed regardless of the scalar type.

For this reason, you should not use `BitIndices{Q,C}()`; instead use `BitIndices(C)`.

# Indexing

`BitIndices` always uses one-based indexing like most Julia arrays. Although it is more natural in
the dense case to use zero-based indexing, as the basis blades are naturally encoded in the indices
for the dense representation of `CliffordNumber`, one-based indexing is used by the tuples which
contain the data associated with this package's implementations of Clifford numbers.

# Broadcasting

Because `BitIndices(x)` only lazily references the indices of `x`, we define a new type,
[`TransformedBitIndices`](@ref), which allows for a function `f` to be lazily associated with the
`BitIndices` object, and this type is constructed when a `f.(BitIndices(x))` is called.

# Interfaces for new subtypes of `AbstractCliffordNumber`

When defining the behavior of `BitIndices` for new subtypes `T` of `AbstractCliffordNumber`, 
`Base.getindex(::BitIndices{Q,T}, i::Integer)` should be defined so that all indices of T that are
not constrained to be zero are returned.

You should also define `CliffordNumbers.bitindices_type(::Type{T})` so that type parameters that do
not affect the construction of the `BitIndices` object are stripped.
"""
struct BitIndices{Q,C<:AbstractCliffordNumber{Q}} <: AbstractBitIndices{Q,C}
end

BitIndices{Q}(::Type{C}) where {Q,C<:AbstractCliffordNumber} = BitIndices{Q,bitindices_type(C)}()
BitIndices(::Type{C}) where C<:AbstractCliffordNumber = BitIndices{signature(C)}(C)
(::Type{B})(x::AbstractCliffordNumber) where B<:BitIndices = B(typeof(x))

# TODO: more efficient defintion of equality

function getindex(b::BitIndices{Q}, i::Integer) where Q
    @boundscheck checkbounds(b, i)
    return BitIndex{Q}(signbit(i-1), unsigned(i-1))
end

# Very efficient tuple generation
@generated function _Tuple(::B) where B<:BitIndices
    data = ntuple(i -> B()[i], Val(length(B)))
    return :($data)
end

# Avoid generating more methods if C contains extraneous parameters
Base.Tuple(::BitIndices{Q,C}) where {Q,C} = _Tuple(BitIndices{Q,bitindices_type(C)}())

Base.map(f, b::AbstractBitIndices) = map(f, Tuple(b))

#---Range of valid indices for CliffordNumber------------------------------------------------------#

Base.keys(x::AbstractCliffordNumber) = keys(typeof(x))  # only need to define on types
Base.keys(::Type{T}) where T<:AbstractCliffordNumber = BitIndices(T)

#---Transformed BitIndices-------------------------------------------------------------------------#
"""
    TransformedBitIndices{Q,C,F} <: AbstractBitIndices{Q,C}

Lazy representation of `BitIndices{Q,C}` with some function of type `f` applied to each element.
These objects can be used to perform common operations which act on basis blades or grades, such as
the reverse or grade involution.
"""
struct TransformedBitIndices{Q,C<:AbstractCliffordNumber{Q},F} <: AbstractBitIndices{Q,C}
    f::F
end

TransformedBitIndices{Q,C}(f) where {Q,C} = TransformedBitIndices{Q,C,typeof(f)}(f)

function TransformedBitIndices(f, ::BitIndices{Q,C}) where {Q,C}
    return TransformedBitIndices{Q,bitindices_type(C)}(f)
end

TransformedBitIndices(f, x) = TransformedBitIndices(f, BitIndices(x))

function getindex(b::TransformedBitIndices{Q,C}, i::Integer) where {Q,C}
    @boundscheck checkbounds(b, i)
    return b.f(BitIndices{Q,C}()[i])
end

Base.Tuple(b::TransformedBitIndices{Q,C}) where {Q,C} = map(b.f, BitIndices{Q,C}())

#---Broadcasting implementation--------------------------------------------------------------------#

# AbstractBitIndices types behave like Tuples when broadcast
Broadcast.BroadcastStyle(::Type{<:AbstractBitIndices}) = Broadcast.Style{Tuple}()

Broadcast.broadcasted(::Broadcast.BroadcastStyle, f, b::BitIndices) = TransformedBitIndices(f, b)

function Broadcast.broadcasted(
    ::Broadcast.BroadcastStyle,
    f,
    b::TransformedBitIndices{Q,C}
) where {Q,C}
    return TransformedBitIndices{Q,C}(x -> f(b.f(x)))
end

#---Aliases for commonly used cases----------------------------------------------------------------#

const ReversedBitIndices{Q,C} = TransformedBitIndices{Q,C,typeof(reverse)}
ReversedBitIndices(x) = TransformedBitIndices(reverse, x)

const GradeInvolutedBitIndices{Q,C} = TransformedBitIndices{Q,C,typeof(grade_involution)}
GradeInvolutedBitIndices(x) = TransformedBitIndices(grade_involution, x)

const ConjugatedBitIndices{Q,C} = TransformedBitIndices{Q,C,typeof(conj)}
ConjugatedBitIndices(x) = TransformedBitIndices(conj, x)

#---Indexing an AbstractCliffordNumber with BitIndices of a type-----------------------------------#

# Indexing with an NTuple returns an NTuple of the coefficients at the BitIndex members
@inline function getindex(x::AbstractCliffordNumber{Q}, B::NTuple{L,BitIndex{Q}}) where {L,Q}
    return map((@inline b -> x[b]), B)
end

"""
    CliffordNumbers.getindex_as_tuple(x::AbstractCliffordNumber{Q}, B::BitIndices{Q,C})

An efficient method for indexing the elements of `x` into a tuple that can be used to construct a
Clifford number of type `C`, or a similar type from `CliffordNumbers.similar_type`.
"""
@generated function getindex_as_tuple(x::AbstractCliffordNumber{Q}, ::BitIndices{Q,C}) where {Q,C}
    inds = to_index.(x, BitIndices(C))
    # The mask is needed for the KVector case.
    # BitIndex objects that don't map to a tuple index just get to turned to index 1
    mask = in.(BitIndices(C), tuple(BitIndices(x)))
    tup = :(getfield.(getindex.(tuple(Tuple(x)), $inds), :value))
    return :(map(*, $tup, $mask))
end

@inline function getindex(x::AbstractCliffordNumber{Q}, B::BitIndices{Q,C}) where {Q,C}
    T = similar_type(C, scalar_type(x))
    return T(getindex_as_tuple(x, B))
end

function getindex(x::AbstractCliffordNumber{Q}, B::AbstractBitIndices{Q,C}) where {Q,C}
    T = similar_type(C, scalar_type(x))
    return T(x[Tuple(B)])
end

# Constructors can follow similar logic
# The type bound is required here to get the expected dispatch behavior
function (C::Type{<:AbstractCliffordNumber{Q1}})(x::AbstractCliffordNumber{Q2}) where {Q1,Q2}
    @assert Q1 === Q2 string(
        "Cannot construct a Clifford number from another Clifford number of a different algebra."
    )
    return C(getindex_as_tuple(x, BitIndices(C)))
end

function (C::Type{<:AbstractCliffordNumber})(x::AbstractCliffordNumber)
    T = similar_type(C, scalar_type(x), Val(signature(x)))
    return T(getindex_as_tuple(x, BitIndices(T)))
end
