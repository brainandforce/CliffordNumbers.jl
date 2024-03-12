#---Efficient multiplication kernels---------------------------------------------------------------#

"""
    CliffordNumbers.bitindex_shuffle(a::BitIndex{Q}, B::NTuple{L,BitIndex{Q}})
    CliffordNumbers.bitindex_shuffle(a::BitIndex{Q}, B::BitIndices{Q})
    
    CliffordNumbers.bitindex_shuffle(B::NTuple{L,BitIndex{Q}}, a::BitIndex{Q})
    CliffordNumbers.bitindex_shuffle(B::BitIndices{Q}, a::BitIndex{Q})

Performs the multiplication `-a * b` for each element of `B` for the above ordering, or `-b * a` for
the below ordering, generating a reordered `NTuple` of `BitIndex{Q}` objects suitable for
implementing a geometric product.
"""
@inline bitindex_shuffle(a::BitIndex{Q}, B::NTuple{L,BitIndex{Q}}) where {L,Q} = map(b -> -a*b, B)
bitindex_shuffle(a::BitIndex{Q}, B::BitIndices{Q}) where Q = map(b -> -a*b, Tuple(B))

@inline bitindex_shuffle(B::NTuple{L,BitIndex{Q}}, a::BitIndex{Q}) where {L,Q} = map(b -> -b*a, B)
bitindex_shuffle( B::BitIndices{Q}, a::BitIndex{Q}) where Q = map(b -> -a*b, Tuple(B))

"""
    CliffordNumbers.nondegenerate_mask(a::BitIndex{Q}, B::NTuple{L,BitIndex{Q}})

Constructs a Boolean mask which is `false` for any multiplication that squares a degenerate blade;
`true` otherwise.
"""
function nondegenerate_mask(a::BitIndex{Q}, B::NTuple{L,BitIndex{Q}}) where {L,Q}
    return map(b -> nondegenerate_mult(a, b), B)
end

"""
    CliffordNumbers.widen_for_mul(x::AbstractCliffordNumber)

Widens `x` to an `EvenCliffordNumber`, `OddCliffordNumber`, or `CliffordNumber` as appropriate
for the fast multiplication kernel.
"""
widen_grade_for_mul(x::Union{CliffordNumber,Z2CliffordNumber}) = x
widen_grade_for_mul(k::KVector{K}) where K = Z2CliffordNumber{isodd(K)}(k)

# Generic fallback for future user-defined types
function widen_grade_for_mul(x::AbstractCliffordNumber)
    all(iseven, nonzero_grades(x)) && return EvenCliffordNumber(x)
    all(iodd, nonzero_grades(x)) && return OddCliffordNumber(x)
    return CliffordNumber(x)
end

#---Geometric product------------------------------------------------------------------------------#
"""
    CliffordNumbers.geometric_product_type(::Type{S}, ::Type{T})

Returns the type of the result of the geometric product of the input types.
"""
@generated function geometric_product_type(
    ::Type{C1},
    ::Type{C2}
) where {Q,C1<:AbstractCliffordNumber{Q},C2<:AbstractCliffordNumber{Q}}
    c1_odd = all(isodd, nonzero_grades(C1))
    c2_odd = all(isodd, nonzero_grades(C2))
    c1_even = all(iseven, nonzero_grades(C1))
    c2_even = all(iseven, nonzero_grades(C2))
    P = (c1_odd && c2_even) || (c1_even && c2_odd)
    T = promote_numeric_type(C1,C2)
    if (!c1_odd && !c1_even) || (!c2_odd && !c2_even)
        return :(CliffordNumber{Q,$T,elements(Q)})
    else
        return :(Z2CliffordNumber{$P,Q,$T,div(elements(Q), 2)})
    end
end

"""
    CliffordNumbers.mul(x::AbstractCliffordNumber{Q,T}, y::AbstractCliffordNumber{Q,T})

A fast geometric product implementation using generated functions for specific cases, and generic
methods which either convert the arguments or fall back to other methods.

The arguments to this function should all agree in scalar type `T`. The `*` function, which exposes
the fast geometric product implementation, promotes the scalar types of the arguments before
utilizing this kernel.

# Notes (for internal use)

Testing with `KVector` instances shows an extremely strong dependency on the multiplication order,
with `x::KVector` and `y::CliffordNumber` being 10x faster than the opposite order.

This could potentially be solved with kernels specific to those cases, but for sufficiently small
multiplications this may be best solved by simply converting KVector arguments using `widen_grade`.
"""
@generated function mul(
    x::Union{CliffordNumber{Q,T},Z2CliffordNumber{<:Any,Q,T}},
    y::Union{CliffordNumber{Q,T},Z2CliffordNumber{<:Any,Q,T}}
) where {Q,T}
    C = geometric_product_type(x, y)
    # Construct the 
    BC = Tuple(BitIndices(C))
    z = zero_tuple(C)
    ex = :($z)
    for a in BitIndices(x)
        inds = bitindex_shuffle(a, BC)
        mask = nondegenerate_mask(a, BC)
        ex = :(map(muladd, x[$a] .* $mask, y[$inds], $ex))
    end
    return :($C($ex))
end
