"""
    CliffordNumbers.Z2CliffordNumber{P,Q,T,L} <: AbstractCliffordNumber{Q,T}

A Clifford number whose only nonzero grades are even or odd. Clifford numbers of this form naturally
arise as versors, the geometric product of 1-vectors.

The type parameter `P` is constrained to be a `Bool`: `true` for odd grade Clifford numbers, and
`false` for even grade Clifford numbers, corresponding to the Boolean result of each grade modulo 2.

# Type aliases

This type is not exported, and usually you will want to refer to the following aliases:
    
    const EvenCliffordNumber{Q,T,L} = Z2CliffordNumber{false,Q,T,L}
    const OddCliffordNumber{Q,T,L} = Z2CliffordNumber{true,Q,T,L}
"""
struct Z2CliffordNumber{P,Q,T<:BaseNumber,L} <: AbstractCliffordNumber{Q,T}
    data::NTuple{L,T}
    function Z2CliffordNumber{P,Q,T,L}(x::Tuple) where {P,Q,T,L}
        @assert P isa Bool "The first type parameter must be a Bool (got $P)."
        check_element_count(div(blade_count(Q), 2), L, x)
        return new(x)
    end
end

"""
    EvenCliffordNumber{P,Q,T,L} (alias for CliffordNumbers.Z2CliffordNumber{false,Q,T,L})

A Clifford number whose only nonzero grades are even. These are the natural choice of representation
for rotors and motors (Euclidean isometries preserving orientation, or "proper" isometries), as well
as their composition with dilations.
"""
const EvenCliffordNumber{Q,T<:BaseNumber,L} = Z2CliffordNumber{false,Q,T,L}

"""
    OddCliffordNumber{P,Q,T,L} (alias for CliffordNumbers.Z2CliffordNumber{true,Q,T,L})

A Clifford number whose only nonzero grades are odd. These are the natural choice of representation
for reflections, as well as their compositions with rotors and motors (Euclidean isometries 
preserving orientation, or "proper" isometries), as well as their composition with dilations.
"""
const OddCliffordNumber{Q,T<:BaseNumber,L} = Z2CliffordNumber{true,Q,T,L}

#---Constructors-----------------------------------------------------------------------------------#

function Z2CliffordNumber{P,Q,T}(x::Tuple{Vararg{BaseNumber,L}}) where {P,Q,T,L}
    return Z2CliffordNumber{P,Q,T,L}(x)
end

Z2CliffordNumber{P,Q}(x::Tuple{Vararg{T}}) where {P,Q,T<:BaseNumber} = Z2CliffordNumber{P,Q,T}(x)

# Automatically convert arguments to a common type
function Z2CliffordNumber{P,Q}(x::Tuple{Vararg{BaseNumber}}) where {P,Q}
    return Z2CliffordNumber{P,Q}(promote(x...))
end

# Convert real/complex numbers to CliffordNumber
(::Type{T})(x::BaseNumber) where {T<:Z2CliffordNumber} = T(ntuple(i -> x*isone(i), Val(nblades(T))))

function (::Type{T})(x::BaseNumber) where {T<:OddCliffordNumber}
    isone(nblades(T)) && return T(ntuple(i -> x*isone(i), Val(nblades(T))))
    throw(
        DomainError(
            x,
            "An OddCliffordNumber cannot be constructed from a scalar unless the algebra " *
            "contains exactly one odd element."
        )
    )
end

#---Number of elements-----------------------------------------------------------------------------#

nblades(::Type{<:Z2CliffordNumber{P,Q}}) where {P,Q} = div(blade_count(Q), 2)
nonzero_grades(::Type{<:Z2CliffordNumber{P,Q}}) where {P,Q} = P:2:dimension(Q)

#---Indexing---------------------------------------------------------------------------------------#

bitindices_type(::Type{<:Z2CliffordNumber{P,Q}}) where {P,Q} = Z2CliffordNumber{P,Q}

@inline function getindex(b::BitIndices{Q,<:Z2CliffordNumber{P,Q}}, i::Integer) where {P,Q}
    @boundscheck checkbounds(b, i)
    n = number_of_parity(i, P)
    return BitIndex{Q}(signbit(n), unsigned(n))
end

@inline to_index(::Type{<:Z2CliffordNumber{P,Q}}, b::BitIndex{Q}) where {P,Q} = div(Int(b), 2) + 1
@inline to_index(x::Z2CliffordNumber{P,Q}, b::BitIndex{Q}) where {P,Q} = to_index(typeof(x), b)

@inline function getindex(x::Z2CliffordNumber{P,Q}, b::BitIndex{Q}) where {P,Q}
    return xor(iseven(grade(b)), P) * flipsign((@inbounds x.data[to_index(x, b)]), b)
end

#---Multiplicative identity------------------------------------------------------------------------#

one(C::Type{<:EvenCliffordNumber{Q}}) where Q = C(ntuple(isone, Val(nblades(C))))

#---Similar types----------------------------------------------------------------------------------#

function similar_type(::Type{<:Z2CliffordNumber{P}}, ::Type{T}, ::Val{Q}) where {P,Q,T<:BaseNumber}
    return Z2CliffordNumber{P,Q,T,div(blade_count(Q),2)}
end

function complement_type(::Type{Z2CliffordNumber{P,Q,T,L}}) where {P,Q,T,L}
    return Z2CliffordNumber{xor(P,isodd(dimension(Q))),Q,T,L}
end

function complement_type(::Type{Z2CliffordNumber{P,Q,T}}) where {P,Q,T}
    return Z2CliffordNumber{xor(P,isodd(dimension(Q))),Q,T}
end

function complement_type(::Type{Z2CliffordNumber{P,Q}}) where {P,Q}
    return Z2CliffordNumber{xor(P,isodd(dimension(Q))),Q}
end

#---Show methods-----------------------------------------------------------------------------------#

short_typename(::Type{<:Z2CliffordNumber{P,Q,T}}) where {P,Q,T} = Z2CliffordNumber{P,Q,T}
