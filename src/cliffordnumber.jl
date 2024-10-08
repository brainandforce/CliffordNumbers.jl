#---Dense representation of Clifford numbers-------------------------------------------------------#
"""
    CliffordNumber{Q,T,L} <: AbstractCliffordNumber{Q,T}

A dense multivector (or Clifford number), with quadratic form `Q`, element type `T`, and length `L`
(which depends entirely on `Q`).

The coefficients are ordered by taking advantage of the natural binary structure of the basis. The
grade of an element is given by the Hamming weight of its index. For the algebra of physical space,
the order is: 1, e₁, e₂, e₁₂, e₃, e₁₃, e₂₃, e₁₂₃ = i. This order allows for more aggressive SIMD
optimization when calculating the geometric product.
"""
struct CliffordNumber{Q,T<:BaseNumber,L} <: AbstractCliffordNumber{Q,T}
    data::NTuple{L,T}
    function CliffordNumber{Q,T,L}(x::Tuple) where {Q,T,L}
        check_element_count(blade_count(Q), L, x)
        return new{Q,T,L}(x)
    end
end

#---Constructors-----------------------------------------------------------------------------------#

CliffordNumber{Q,T}(x::Tuple{Vararg{BaseNumber,L}}) where {Q,T,L} = CliffordNumber{Q,T,L}(x)
CliffordNumber{Q}(x::Tuple{Vararg{T}}) where {Q,T<:BaseNumber} = CliffordNumber{Q,T}(x)

# Automatically convert arguments to a common type
CliffordNumber{Q}(x::Tuple{Vararg{BaseNumber}}) where Q = CliffordNumber{Q}(promote(x...))

# Convert real/complex numbers to CliffordNumber
(::Type{T})(x::BaseNumber) where T<:CliffordNumber = T(ntuple(i -> x*isone(i), Val(nblades(T))))

#---Number of elements-----------------------------------------------------------------------------#

nblades(::Type{<:CliffordNumber{Q}}) where Q = blade_count(Q)
nonzero_grades(::Type{<:CliffordNumber{Q}}) where Q = 0:dimension(Q)

#---Default BitIndices construction should include all possible BitIndex objects-------------------#

bitindices_type(::Type{<:CliffordNumber{Q}}) where Q = CliffordNumber{Q}

BitIndices{Q}() where Q = BitIndices{Q,CliffordNumber{Q}}()
BitIndices(Q::Metrics.AbstractSignature) = BitIndices{Q}()

#---Clifford number indexing-----------------------------------------------------------------------#

@inline to_index(::Type{<:CliffordNumber{Q}}, i::BitIndex{Q}) where Q = (Int(i) % blade_count(Q) + 1)
@inline to_index(x::CliffordNumber{Q}, i::BitIndex{Q}) where Q = to_index(typeof(x), i)

@inline function getindex(x::CliffordNumber{Q}, b::BitIndex{Q}) where Q 
    return flipsign((@inbounds x.data[to_index(x, b)]), b)
end

#---Multiplicative identity------------------------------------------------------------------------#

one(C::Type{<:CliffordNumber{Q}}) where Q = C(ntuple(isone, Val(nblades(C))))

#---Similar types----------------------------------------------------------------------------------#

function similar_type(::Type{<:CliffordNumber}, ::Type{T}, ::Val{Q}) where {Q,T<:BaseNumber}
    return CliffordNumber{Q,T,blade_count(Q)}
end

complement_type(::Type{<:AbstractCliffordNumber{Q,T}}) where {Q,T<:BaseNumber} = CliffordNumber{Q,T}
complement_type(::Type{<:AbstractCliffordNumber{Q}}) where Q = CliffordNumber{Q}
complement_type(::Type{<:AbstractCliffordNumber}) = CliffordNumber

complement_type(::Type{CliffordNumber{Q,T,L}}) where {Q,T,L} = CliffordNumber{Q,T,L}

#---Show methods-----------------------------------------------------------------------------------#

short_typename(::Type{<:CliffordNumber{Q,T}}) where {Q,T} = CliffordNumber{Q,T}
