#---Sets of `BitIndex` objects that index specific types-------------------------------------------#
"""
    AbstractBitIndices{Q,C<:AbstractCliffordNumber{Q}} <: AbstractVector{BitIndex{Q}}

Supertype for vectors containing all valid `BitIndex{Q}` objects for the basis elements represented
by `C`.
"""
abstract type AbstractBitIndices{Q,C<:AbstractCliffordNumber{Q}} <: AbstractVector{BitIndex{Q}}
end

Base.size(::Type{<:AbstractBitIndices{Q,C}}) where {Q,C} = tuple(length(C))
Base.size(b::AbstractBitIndices) = size(typeof(b))

#---Broadcasting-----------------------------------------------------------------------------------#

#=
struct BitIndicesStyle{Q,C<:AbstractCliffordNumber{Q}} <: Base.Broadcast.AbstractArrayStyle{1}
end

Base.Broadcast.BroadcastStyle(::Type{<:AbstractBitIndices{Q,C}}) = BitIndicesStyle{Q,C}()
Base.Broadcast.BroadcastStyle(::BitIndicesStyle, s::Base.Broadcast.DefaultArrayStyle) = s
Base.Broadcast.BroadcastStyle(s::Base.Broadcast.DefaultArrayStyle, ::BitIndicesStyle) = s
=#

#---Clifford number iteration----------------------------------------------------------------------#
"""
    BitIndices{Q<:QuadraticForm,C<:AbstractCliffordNumber{Q,<:Any}} <: AbstractVector{BitIndex{Q}}

Represents a range of valid `BitIndex` objects for the nonzero components of a given multivector 
with quadratic form `Q`.

For a generic `AbstractCliffordNumber{Q}`, this returns `BitIndices{CliffordNumber{Q}}`, which
contains all possible indices for a multivector associated with the quadratic form `Q`. This may 
also be constructed with `BitIndices(Q)`.

For sparse representations, such as `KVector{K,Q}`, the object only contains the indices of the
nonzero elements of the multivector.

# Construction

`BitIndices` can be constructed by calling the type constructor with either the multivector or its
type.

# Indexing

`BitIndices` always uses one-based indexing like most Julia arrays. Although it is more natural in
the dense case to use zero-based indexing, as the basis blades are naturally encoded in the indices
for the dense representation of `CliffordNumber`, one-based indexing is used by the tuples which
contain the data associated with this package's implementations of Clifford numbers.

# Interfaces for new subtypes of `AbstractCliffordNumber`

When defining the behavior of `BitIndices` for new subtypes `T` of `AbstractCliffordNumber`, 
`Base.getindex(::BitIndices{Q,T}, i::Integer)` should be defined so that all indices of T that are
not constrained to be zero are returned.
"""
struct BitIndices{Q,C<:AbstractCliffordNumber{Q}} <: AbstractBitIndices{Q,C}
end

BitIndices{Q}(x::AbstractCliffordNumber) where Q = BitIndices{Q,typeof(x)}()
# Ensure that the second type parameter is a concrete type
# It doesn't actually need to be, but doing this should simplify things
BitIndices{Q}(C::Type{<:AbstractCliffordNumber}) where Q = BitIndices{Q}(zero(C))
BitIndices(x) = BitIndices{QuadraticForm(x)}(x)

# TODO: more efficient defintion of equality

function Base.getindex(b::BitIndices{Q}, i::Integer) where Q
    @boundscheck checkbounds(b, i)
    return BitIndex{Q}(signbit(i-1), unsigned(i-1))
end

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
struct TransformedBitIndices{Q,C,F} <: AbstractBitIndices{Q,C}
    f::F
end

TransformedBitIndices{Q,C}(f) where {Q,C} = TransformedBitIndices{Q,C,typeof(f)}(f)
TransformedBitIndices(f, ::BitIndices{Q,C}) where {Q,C} = TransformedBitIndices{Q,C}(f)
TransformedBitIndices(f, x) = TransformedBitIndices(BitIndices(x), f)

function Base.getindex(b::TransformedBitIndices{Q}, i::Integer) where Q
    @boundscheck checkbounds(b, i)
    return b.f(BitIndex{Q}(signbit(i-1), unsigned(i-1)))
end

#---Broadcasting implementation--------------------------------------------------------------------#

Base.Broadcast.broadcasted(f, b::BitIndices) = TransformedBitIndices(f, b)

function Base.Broadcast.broadcasted(f, b::TransformedBitIndices{Q,C}) where {Q,C}
    return TransformedBitIndices{Q,C}(b.f |> f)
end

#---Aliases for commonly used cases----------------------------------------------------------------#

const ReversedBitIndices{Q,C} = TransformedBitIndices{Q,C,typeof(reverse)}
ReversedBitIndices(x) = TransformedBitIndices(reverse, x)

const GradeInvolutedBitIndices{Q,C} = TransformedBitIndices{Q,C,typeof(grade_involution)}
GradeInvolutedBitIndices(x) = TransformedBitIndices(grade_involution, x)

const ConjugatedBitIndices{Q,C} = TransformedBitIndices{Q,C,typeof(conj)}
ConjugatedBitIndices(x) = TransformedBitIndices(conj, x)

#---Indexing an AbstractCliffordNumber with BitIndices of a type-----------------------------------#

function Base.getindex(x::AbstractCliffordNumber{Q}, b::AbstractBitIndices{Q,C}) where {Q,C}
    data = ntuple(i -> x[b[i]], Val(length(C)))
    return similar_type(C, numeric_type(x))(data)
end

# Constructors can follow similar logic, just don't use similar_type
function (T::Type{<:AbstractCliffordNumber{Q}})(x::AbstractCliffordNumber{Q}) where Q
    data = ntuple(i -> x[BitIndices(T)[i]], Val(length(T)))
    return T(data)
end
