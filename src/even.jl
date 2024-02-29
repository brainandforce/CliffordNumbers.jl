"""
    CliffordNumbers.Z2CliffordNumber{P,Q,T,L}

A Clifford number whose only nonzero grades are even or odd. Clifford numbers of this form naturally
arise as versors, the geometric product of 1-vectors.

The type parameter `P` is constrained to be a `Bool`: `true` for odd grade Clifford numbers, and
`false` for even grade Clifford numbers, corresponding to the Boolean result of each grade modulo 2.

# Type aliases

This type is not exported, and usually you will want to refer to the following aliases:
    
    const EvenCliffordNumber{Q,T,L} = Z2CliffordNumber{false,Q,T,L}
    const OddCliffordNumber{Q,T,L} = Z2CliffordNumber{true,Q,T,L}
"""
struct Z2CliffordNumber{P,Q<:QuadraticForm,T<:BaseNumber,L} <: AbstractCliffordNumber{Q,T}
    data::NTuple{L,T}
    function Z2CliffordNumber{P,Q,T,L}(x::Tuple) where {P,Q,T,L}
        @assert P isa Bool "The first type parameter must be a Bool (got $P)."
        check_element_count(q -> div(elements(q), 2), Q, L, x)
        return new(x)
    end
end

const EvenCliffordNumber{Q<:QuadraticForm,T<:BaseNumber,L} = Z2CliffordNumber{false,Q,T,L}
const OddCliffordNumber{Q<:QuadraticForm,T<:BaseNumber,L} = Z2CliffordNumber{true,Q,T,L}

#---Constructors-----------------------------------------------------------------------------------#

function Z2CliffordNumber{P,Q,T}(x::Tuple{Vararg{BaseNumber,L}}) where {P,Q,T,L}
    return Z2CliffordNumber{P,Q,T,L}(x)
end

Z2CliffordNumber{P,Q}(x::Tuple{Vararg{T}}) where {P,Q,T<:BaseNumber} = Z2CliffordNumber{P,Q,T}(x)

# Automatically convert arguments to a common type
function Z2CliffordNumber{P,Q}(x::Tuple{Vararg{BaseNumber}}) where {P,Q}
    return Z2CliffordNumber{P,Q}(promote(x...))
end

# Allow varargs arguments
(::Type{T})(x::Vararg{BaseNumber}) where {P,Q,T<:Z2CliffordNumber{P,Q}} = T(x)

# Convert real/complex numbers to CliffordNumber
(::Type{T})(x::BaseNumber) where {T<:Z2CliffordNumber} = T(ntuple(i -> x*isone(i), Val(length(T))))

#---Number of elements-----------------------------------------------------------------------------#

Base.length(::Type{<:Z2CliffordNumber{P,Q}}) where {P,Q} = div(elements(Q), 2)
Base.length(x::Z2CliffordNumber) = length(typeof(x))

nonzero_grades(::Type{<:Z2CliffordNumber{P,Q}}) where {P,Q} = P:2:dimension(Q)

#---Indexing---------------------------------------------------------------------------------------#

function Base.getindex(b::BitIndices{Q,<:Z2CliffordNumber{P,Q}}, i::Integer) where {P,Q}
    @boundscheck checkbounds(b, i)
    n = number_of_parity(i, P)
    return BitIndex{Q}(signbit(n), unsigned(n))
end

function Base.getindex(x::Z2CliffordNumber{P,Q}, b::BitIndex{Q}) where {P,Q}
    return sign(b) * (@inbounds x.data[div(b.blade, 2) + 1]) * xor(isevil(b.blade), P)
end

#---Multiplicative identity------------------------------------------------------------------------#

Base.one(C::Type{<:EvenCliffordNumber{Q}}) where Q = C(ntuple(isone, Val(length(C))))

#---Similar types----------------------------------------------------------------------------------#

function similar_type(
    ::Type{<:Z2CliffordNumber{P,<:QuadraticForm}},
    T::Type{<:BaseNumber}, 
    Q::Type{<:QuadraticForm}
) where P
    return Z2CliffordNumber{P,Q,T,div(elements(Q),2)}
end

#---Show methods-----------------------------------------------------------------------------------#

short_typename(::Type{<:Z2CliffordNumber{P,Q,T}}) where {P,Q,T} = Z2CliffordNumber{P,Q,T}
