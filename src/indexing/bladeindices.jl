"""
    CliffordNumbers.blade_indices_type(C::Type{<:AbstractCliffordNumber{Q,T}})

Removes extraneous type parameters from `C`, converting it to the least parameterized type that
can be used to parameterize an `AbstractBladeIndices{Q,C}` object. This is to avoid issues with the
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
julia> CliffordNumbers.blade_indices_type(CliffordNumber{VGA(3),Float32,8})
CliffordNumber{VGA(3)}

julia> CliffordNumbers.blade_indices_type(KVector{2,STA,Bool})
KVector{2,STA}
```
"""
blade_indices_type(::Type{AbstractCliffordNumber{Q,<:Any}}) where Q = AbstractCliffordNumber{Q}
blade_indices_type(x::AbstractCliffordNumber) = blade_indices_type(typeof(x))

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

For this reason, you should not use `BladeIndices{Q,C}()`; instead use `BladeIndices(C)`.

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

#---Fast equality checks---------------------------------------------------------------------------#

==(::B, ::B) where B<:CGNBladeIndices = true

function ==(::CGNBladeIndices{Q,S,N,I,R}, ::CGNBladeIndices{Q,T,N,I,R}) where {Q,S,T,N,I,R}
    return blade_indices_type(S) === blade_indices_type(T)
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

#---Range of valid indices for CliffordNumber------------------------------------------------------#

Base.keys(::Type{T}) where T<:AbstractCliffordNumber = BladeIndices(T)
Base.keys(x::AbstractCliffordNumber) = keys(typeof(x))  # only need to define on types

#---Custom broadcast and map implementations-------------------------------------------------------#

Broadcast.BroadcastStyle(::Type{<:CGNBladeIndices}) = Broadcast.Style{Tuple}()
Broadcast.BroadcastStyle(::Type{<:BladeIndices}) = Broadcast.Style{Tuple}()

+(inds::CGNBladeIndices) = inds
Broadcast.broadcasted(::Union{typeof(+),typeof(identity)}, inds::CGNBladeIndices) = inds

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

Base.map(::Union{typeof(+),typeof(identity)}, inds::CGNBladeIndices) = inds
Base.map(::typeof(-), inds::CGNBladeIndices) = -inds
Base.map(::typeof(grade_involution), inds::CGNBladeIndices) = grade_involution.(inds)
Base.map(::Union{typeof(adjoint), typeof(reverse)}, inds::CGNBladeIndices) = reverse.(inds)
Base.map(::typeof(conj), inds::CGNBladeIndices) = conj.(inds)

#---Indexing AbstractCliffordNumber instances------------------------------------------------------#

# Indexing with an NTuple returns an NTuple of the coefficients at the BladeIndex members
@inline function getindex(x::AbstractCliffordNumber{Q}, B::NTuple{L,BladeIndex{Q}}) where {L,Q}
    return map((@inline b -> x[b]), B)
end

"""
    CliffordNumbers.getindex_as_tuple(x::AbstractCliffordNumber{Q}, B::BladeIndices{Q,C})

An efficient method for indexing the elements of `x` into a tuple that can be used to construct a
Clifford number of type `C`, or a similar type from `CliffordNumbers.similar_type`.
"""
@generated function getindex_as_tuple(
    x::AbstractCliffordNumber{Q},
    ::B
) where {Q,C,B<:CGNBladeIndices{Q,C}}
    return :(getindex.(tuple(x), B()))
end

@inline function getindex(x::AbstractCliffordNumber{Q}, B::CGNBladeIndices{Q,C}) where {Q,C}
    T = similar_type(C, scalar_type(x))
    return T(getindex_as_tuple(x, B))
end

# Constructors can follow similar logic
# The type bound is required here to get the expected dispatch behavior
function (C::Type{<:AbstractCliffordNumber{Q1}})(x::AbstractCliffordNumber{Q2}) where {Q1,Q2}
    @assert Q1 === Q2 string(
        "Cannot construct a Clifford number from another Clifford number of a different algebra."
    )
    return C(getindex_as_tuple(x, BladeIndices(C)))
end

function (C::Type{<:AbstractCliffordNumber})(x::AbstractCliffordNumber)
    T = similar_type(C, scalar_type(x), Val(signature(x)))
    return T(getindex_as_tuple(x, BladeIndices(T)))
end

#---Show methods to avoid printing too much type detail--------------------------------------------#

Base.summary(io::IO, i::BladeIndices) = print(io, length(i), "-element ", typeof(i))

function Base.summary(io::IO, i::CGNBladeIndices{Q,C,N,I,R}) where {Q,C,N,I,R}
    automorphism = ifelse(I, ifelse(R, "Clifford conjugated","grade involuted"), "reversed"^R)
    opstring = "negated"^N * ", "^(N && !isempty(automorphism)) * automorphism
    print(
        io,
        length(i), "-element ", BladeIndices{Q,C}, " (", opstring, ')'
    )
end

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
