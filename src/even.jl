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

Z2CliffordNumber{P,Q,T}(x::NTuple{L,<:BaseNumber}) where {P,Q,T,L} = Z2CliffordNumber{P,Q,T,L}(x)
Z2CliffordNumber{P,Q,T}(x::Vararg{BaseNumber,L}) where {P,Q,T,L} = Z2CliffordNumber{P,Q,T,L}(x)

function Z2CliffordNumber{P,Q}(x::NTuple{L,<:BaseNumber}) where {P,Q,L}
    data = promote(x...)
    return Z2CliffordNumber{P,Q,eltype(data),L}(data)
end

Z2CliffordNumber{P,Q}(x::Vararg{BaseNumber,L}) where {P,Q,L} = Z2CliffordNumber{P,Q}(x)

const EvenCliffordNumber{Q<:QuadraticForm,T<:BaseNumber,L} = Z2CliffordNumber{false,Q,T,L}
const OddCliffordNumber{Q<:QuadraticForm,T<:BaseNumber,L} = Z2CliffordNumber{true,Q,T,L}

#---Convert scalars to even Clifford numbers-------------------------------------------------------#

function EvenCliffordNumber{Q,T,L}(x::BaseNumber) where {Q,T,L}
    return EvenCliffordNumber{Q,T,L}(ntuple(i -> T(isone(i) * x), Val(L)))
end

EvenCliffordNumber{Q,T}(x::BaseNumber) where {Q,T} = EvenCliffordNumber{Q,T,div(elements(Q), 2)}(x)
EvenCliffordNumber{Q}(x::T) where {Q,T<:BaseNumber} = EvenCliffordNumber{Q,T,div(elements(Q), 2)}(x)

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
    return sign(b) * x.data[div(b.blade, 2) + 1] * xor(isevil(b.blade), P)
end

#---Multiplicative identity for EvenCliffordNumber-------------------------------------------------#

Base.oneunit(C::Type{<:EvenCliffordNumber{Q}}) where Q = C(ntuple(isone, Val(length(C))))
Base.one(C::Type{<:EvenCliffordNumber{Q}}) where Q = oneunit(C)

#---Similar types----------------------------------------------------------------------------------#

function similar_type(
    ::Type{<:Z2CliffordNumber{P,<:QuadraticForm}},
    T::Type{<:BaseNumber}, 
    Q::Type{<:QuadraticForm}
) where P
    return Z2CliffordNumber{P,Q,T,div(elements(Q),2)}
end

#---Show methods-----------------------------------------------------------------------------------#

function Base.show(io::IO, x::Z2CliffordNumber{P}) where P
    print(io, Z2CliffordNumber{P}, "{", QuadraticForm(x), ",", numeric_type(x), "}", x.data)
end

function Base.summary(io::IO, x::Z2CliffordNumber{P}) where P
    println(io, Z2CliffordNumber{P}, "{", QuadraticForm(x), ",", numeric_type(x), "}:")
end
