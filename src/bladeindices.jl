"""
    CliffordNumbers.blade_indices_type(C::Type{<:AbstractCliffordNumber{Q,T}})

Removes extraneous type parameters from `C`, converting it to the least parameterized type that
can be used to parameterize an `AbstractBitIndices{Q,C}` object. This is to avoid issues with the
proliferation of type parameters that would construct identical `BladeIndices` objects otherwise:
for instance, `BladeIndices{VGA(3),EvenCliffordNumber{VGA(3),Float32,4}}()` and
`BladeIndices{VGA(3),EvenCliffordNumber{VGA(3),Int}}()` have identical elements, and are equal when
compared with `==`, but are not the same object.

Preventing proliferation of type parameters is critical when working with generated functions that
have no generic implementation (i.e. an `if @generated` block).

For types defined in this package, this strips the scalar type parameter `T` and any length
parameter present, but leaves the algebra type parameter `Q`.

# Examples
```julia-repl
julia> CliffordNumbers.bitindices_type(CliffordNumber{VGA(3),Float32,8})
CliffordNumber{VGA(3)}

julia> CliffordNumbers.bitindices_type(KVector{2,STA,Bool})
KVector{2,STA}
```
"""
blade_indices_type(::Type{AbstractCliffordNumber{Q,<:Any}}) where Q = AbstractCliffordNumber{Q}
blade_indices_type(x::AbstractCliffordNumber) = bitindices_type(typeof(x))

#---Type and aliases-------------------------------------------------------------------------------#
"""
    CliffordNumbers.CGNBladeIndices{Q,C<:AbstractCliffordNumber{Q},N,I,R}

Compactly represents all valid indices of `C` in a singleton type. The last three type parameters 
are Boolean values which represent the application of negation (`N`), grade involution (`I`), or
reversion (`R`) to the indices.

# Broadcasting

Custom broadcast implementations are used to change the values of the type parameters `N`, `I`, and
`R`, so that the types do not need to be managed directly with their constructors.

Printing of each variant type is done using broadcasted functions:
```julia_repl
julia> show(CliffordNumbers.CGNBladeIndices{VGA(3),KVector{2,VGA(3)},true,true,true}())
-grade_involution.(reverse.(BladeIndices(KVector{2, VGA(3)})))
```
# Aliases

For convenience, the alias [`BladeIndices`](@ref) is provided, which produces none of these
transformations.
"""
struct CGNBladeIndices{Q,C<:AbstractCliffordNumber{Q},N,I,R} <: AbstractVector{BladeIndex{Q}}
    function CGNBladeIndices{Q,C,N,I,R}() where {Q,C,N,I,R}
        @assert (N,I,R) isa NTuple{3,Bool} "The final type parameters must be of type Bool."
        return new()
    end
end

"""
    BladeIndices{Q,C<:AbstractCliffordNumber{Q}} <: AbstractVector{BladeIndex{Q}}

Represents all valid nonzero indices of an `AbstractCliffordNumber` in a lazily generated array of
static size.

New subtypes of `AbstractCliffordNumber{Q}` should define `Base.getindex(::BladeIndex{Q,C}, ::Int)`.
For performance, this definition should be annotated with the macros `@inline` and
`Base.@propagate_inbounds` (to preserve the `@inbounds` context).

# Construction

`BladeIndices` can be constructed by calling the type constructor with the Clifford number or its
type. This will automatically strip some type parameters so that identical `BladeIndices` objects
are constructed regardless of the scalar type.

For this reason, you should not use `BitIndices{Q,C}()`; instead use `BitIndices(C)`.

# Implementation

This type is implemented as an alias of 
`CliffordNumbers.CGNBladeIndices{Q,C<:AbstractCliffordNumber{Q},false,false,false}`. More 
information is available at the docstring for [`CliffordNumbers.CGNBladeIndices`](@ref).
"""
const BladeIndices{Q,C} = CGNBladeIndices{Q,C,false,false,false}

#---Constructors-----------------------------------------------------------------------------------#

function BladeIndices{Q}(::Type{C}) where {Q,C<:AbstractCliffordNumber}
    return BladeIndices{Q,blade_indices_type(C)}()
end

function BladeIndices(::Type{C}) where {Q,C<:AbstractCliffordNumber{Q}}
    return BladeIndices{Q,blade_indices_type(C)}()
end

BladeIndices{Q}(x::AbstractCliffordNumber) where Q = BladeIndices{Q}(typeof(x))
BladeIndices(x::AbstractCliffordNumber{Q}) where Q = BladeIndices{Q}(typeof(x))

#---Size and indexing------------------------------------------------------------------------------#

Base.size(::CGNBladeIndices{Q,C}) where {Q,C} = tuple(nblades(C))
Base.size(::Type{<:CGNBladeIndices{Q,C}}) where {Q,C} = tuple(nblades(C))

Base.IndexStyle(::Type{<:CGNBladeIndices}) = IndexLinear()

# In practice, new types should only need to define `Base.getindex(::BladeIndices{Q,C}, ::Int)`
# This definition transforms a `BladeIndices{Q,C}` object with a `CyclicGradeNegation{N,I,R}`
@inline function getindex(::CGNBladeIndices{Q,C,N,I,R}, i::Int) where {Q,C,N,I,R}
    return CyclicGradeNegation{N,I,R}()(BladeIndices{Q,C}()[i])
end

function getindex(::BladeIndices{Q,C}, i::Int) where {Q,C<:AbstractCliffordNumber{Q}}
    T = blade_indices_type(C)
    error(
        "The method Base.getindex(::BladeIndices{Q,$T}, ::Int) must be defined."
    )
end

#---Conversion to a Tuple--------------------------------------------------------------------------#

@generated function _Tuple(
    ::CGNBladeIndices{Q,C,N,I,R}
) where {Q,C<:AbstractCliffordNumber{Q},N,I,R}
    data = ntuple(i -> CGNBladeIndices{Q,C,N,I,R}()[i], Val(nblades(C)))
    return :($data)
end

function Base.Tuple(::CGNBladeIndices{Q,C,N,I,R}) where {Q,C,N,I,R}
    return _Tuple(CGNBladeIndices{Q,blade_indices_type(C),N,I,R}())
end

Base.map(f, inds::CGNBladeIndices{Q,C,N,I,R}) where {Q,C,N,I,R} = map(f, Tuple(inds))

#---Custom broadcasting behavior-------------------------------------------------------------------#

Broadcast.BroadcastStyle(::Type{<:CGNBladeIndices}) = Broadcast.Style{Tuple}()
Broadcast.BroadcastStyle(::Type{<:BladeIndices}) = Broadcast.Style{Tuple}()

-(::CGNBladeIndices{Q,C,N,I,R}) where {Q,C,N,I,R} = CGNBladeIndices{Q,C,!N,I,R}()
Broadcast.broadcasted(::typeof(-), inds::CGNBladeIndices) = -inds

function Broadcast.broadcasted(
    ::typeof(grade_involution),
    ::CGNBladeIndices{Q,C,N,I,R}
) where {Q,C,N,I,R}
    return CGNBladeIndices{Q,C,N,!I,R}()
end

function Broadcast.broadcasted(
    ::Union{typeof(adjoint), typeof(reverse)},
    ::CGNBladeIndices{Q,C,N,I,R}
) where {Q,C,N,I,R}
    return CGNBladeIndices{Q,C,N,I,!R}()
end

function Broadcast.broadcasted(::typeof(conj), ::CGNBladeIndices{Q,C,N,I,R}) where {Q,C,N,I,R}
    return CGNBladeIndices{Q,C,N,!I,!R}()
end

#---Indexing AbstractCliffordNumber instances------------------------------------------------------#
"""
    CliffordNumbers.getindex_as_tuple(x::AbstractCliffordNumber{Q}, B::BladeIndices{Q,C})

An efficient method for indexing the elements of `x` into a tuple that can be used to construct a
Clifford number of type `C`, or a similar type from `CliffordNumbers.similar_type`.
"""
@generated function getindex_as_tuple(x::AbstractCliffordNumber{Q}, ::BladeIndices{Q,C}) where {Q,C}
    inds = to_index.(x, BladeIndices(C))
    # The mask is needed for the KVector case.
    # BladeIndex objects that don't map to a tuple index just get to turned to index 1
    mask = in.(BladeIndices(C), tuple(BladeIndices(x)))
    return :(map(*, getindex.(tuple(Tuple(x)), $inds), $mask))
end

@inline function getindex(x::AbstractCliffordNumber{Q}, B::BladeIndices{Q,C}) where {Q,C}
    T = similar_type(C, scalar_type(x))
    return T(getindex_as_tuple(x, B))
end

#---Show methods to avoid printing too much type detail--------------------------------------------#

function Base.show(io::IO, ::BladeIndices{Q,C}) where {Q,C}
    if C === blade_indices_type(C)
        print(io, BladeIndices, '(', C, ')')
    else
        print(io, BladeIndices{Q,C}, "()")
    end
end

function Base.show(io::IO, ::CGNBladeIndices{Q,C,N,I,R}) where {Q,C,N,I,R}
    print(io, '-'^N)
    (I && !R) && print(io, grade_involution)
    (R && !I) && print(io, reverse)
    (I && R) && print(io, conj)
    (I || R) && print(io, ".(")
    show(io, BladeIndices{Q,C}())
    print(io, ')'^(I || R))
end
