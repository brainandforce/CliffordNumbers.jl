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
@inline function bitindex_shuffle(a::BitIndex{Q}, B::NTuple{L,BitIndex{Q}}) where {L,Q}
    return map(b -> reverse(a) * b, B)
end

@inline function bitindex_shuffle(B::NTuple{L,BitIndex{Q}}, a::BitIndex{Q}) where {L,Q}
    return map(b -> b * reverse(a), B)
end

bitindex_shuffle(a::BitIndex{Q}, B::BitIndices{Q}) where Q = bitindex_shuffle(a, Tuple(B))
bitindex_shuffle(B::BitIndices{Q}, a::BitIndex{Q}) where Q = bitindex_shuffle(Tuple(B), a)

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

#---Grade filters----------------------------------------------------------------------------------#
"""
    CliffordNumbers.GradeFilter{S}

A type that can be used to filter certain products of blades in a geometric product multiplication.
The type parameter `S` must be a `Symbol`. The single instance of `GradeFilter{S}` is a callable
object which implements a function that takes two or more `BitIndex{Q}` objects `a` and `b` and
returns `false` if the product of the blades indexed is zero.

To implement a grade filter for a product function `f`, define the following method:
    (::GradeFilter{:f})(::BitIndex{Q}, ::BitIndex{Q})
    # Or if the definition allows for more arguments
    (::GradeFilter{:f})(::BitIndex{Q}...) where Q
"""
struct GradeFilter{S}
    GradeFilter{S}() where S = (@assert S isa Symbol "Type parameter must be a Symbol."; new())
end

(::GradeFilter{S})(args...) where S = error("This filter has not been implemented.")

(::GradeFilter{:*})(args::BitIndex{Q}...) where Q = true

(::GradeFilter{:∧})(args::BitIndex{Q}...) where Q = has_wedge(args...)

(::GradeFilter{:⨼})(a::BitIndex{Q}, b::BitIndex{Q}) where Q = (grade(b) - grade(a)) == grade(a*b)
(::GradeFilter{:⨽})(a::BitIndex{Q}, b::BitIndex{Q}) where Q = (grade(a) - grade(b)) == grade(a*b)

function (::GradeFilter{:dot})(a::BitIndex{Q}, b::BitIndex{Q}) where Q
    return abs(grade(a) - grade(b)) == grade(a*b)
end

const ContractionGradeFilters = Union{GradeFilter{:⨼},GradeFilter{:⨽},GradeFilter{:dot}}

"""
    CliffordNumbers.mul_mask(F::GradeFilter, a::BitIndex{Q}, B::NTuple{L,BitIndices{Q}})
    CliffordNumbers.mul_mask(F::GradeFilter, B::NTuple{L,BitIndices{Q}}, a::BitIndex{Q})

    CliffordNumbers.mul_mask(F::GradeFilter, a::BitIndex{Q}, B::BitIndices{Q})
    CliffordNumbers.mul_mask(F::GradeFilter, B::BitIndices{Q}, a::BitIndex{Q})

Generates a `NTuple{L,Bool}` which is `true` whenever the multiplication of the blade indexed by `a`
and blades indexed by `B` is nonzero. `false` is returned if the grades multiply to zero due to the
squaring of a degenerate component, or if they are filtered by `F`.
"""
function mul_mask(F::GradeFilter, a::BitIndex{Q}, B::NTuple{L,BitIndex{Q}}) where {L,Q}
    return map(b -> F(a,b) & nondegenerate_mult(a,b), B)
end

function mul_mask(F::GradeFilter, B::NTuple{L,BitIndex{Q}}, a::BitIndex{Q}) where {L,Q}
    return map(b -> F(b,a) & nondegenerate_mult(b,a), B)
end

mul_mask(F::GradeFilter, a::BitIndex{Q}, B::BitIndices{Q}) where Q = mul_mask(F, a, Tuple(B))
mul_mask(F::GradeFilter, B::BitIndices{Q}, a::BitIndex{Q}) where Q = mul_mask(F, Tuple(B), a)

#---Product return types---------------------------------------------------------------------------#
"""
    CliffordNumbers.product_return_type(::Type{X}, ::Type{Y}, [::GradeFilter{S}])

Returns a suitable type for representing the product of Clifford numbers of types `X` and `Y`. The
`GradeFilter{S}` argument allows for the return type to be changed depending on the type of product.
Without specialization on `S`, a type suitable for the geometric product is returned.
"""
@generated function product_return_type(
    ::Type{C1},
    ::Type{C2},
    ::GradeFilter{<:Any}
) where {Q,C1<:AbstractCliffordNumber{Q},C2<:AbstractCliffordNumber{Q}}
    c1_odd = all(isodd, nonzero_grades(C1))
    c2_odd = all(isodd, nonzero_grades(C2))
    c1_even = all(iseven, nonzero_grades(C1))
    c2_even = all(iseven, nonzero_grades(C2))
    # Parity: true for odd multivectors, false for even multivectors
    P = (c1_odd && c2_even) || (c1_even && c2_odd)
    T = promote_numeric_type(C1,C2)
    if (!c1_odd && !c1_even) || (!c2_odd && !c2_even)
        return :(CliffordNumber{Q,$T,elements(Q)})
    else
        return :(Z2CliffordNumber{$P,Q,$T,div(elements(Q), 2)})
    end
end

function product_return_type(
    X::Type{<:KVector{K1,Q}},
    Y::Type{<:KVector{K2,Q}},
    ::GradeFilter{:wedge}
) where {Q,K1,K2}
    K = min(K1 + K2, dimension(Q))
    return KVector{K, Q, promote_numeric_type(X, Y), binomial(dimension(Q), k)}
end

function product_return_type(
    X::Type{<:KVector{K1,Q}},
    Y::Type{<:KVector{K2,Q}},
    ::ContractionGradeFilters
) where {Q,K1,K2}
    K = abs(K1 - K2)
    return KVector{K, Q, promote_numeric_type(X, Y), binomial(dimension(Q), k)}
end

function product_return_type(
    x::AbstractCliffordNumber,
    y::AbstractCliffordNumber,
    F::GradeFilter = GradeFilter{:*}()
)
    return product_return_type(typeof(x), typeof(y), F)
end

#---Geometric product------------------------------------------------------------------------------#
"""
    CliffordNumbers.geometric_product_type(::Type{S}, ::Type{T})

Returns the type of the result of the geometric product of the input types.
"""
function geometric_product_type(
    ::Type{C1},
    ::Type{C2}
) where {Q,C1<:AbstractCliffordNumber{Q},C2<:AbstractCliffordNumber{Q}}
    return product_return_type(C1, C2, GradeFilter{:*}())
end

"""
    CliffordNumbers.mul(
        x::Union{CliffordNumber{Q,T},Z2CliffordNumber{<:Any,Q,T}},
        y::Union{CliffordNumber{Q,T},Z2CliffordNumber{<:Any,Q,T}},
        [F::GradeFilter = GradeFilter{:*}()]
    )

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
    y::Union{CliffordNumber{Q,T},Z2CliffordNumber{<:Any,Q,T}},
    F::GradeFilter = GradeFilter{:*}()
) where {Q,T}
    C = product_return_type(x, y, F())
    ex = :($(zero_tuple(C)))
    for a in BitIndices(x)
        inds = bitindex_shuffle(a, Tuple(BitIndices(C)))
        mask = mul_mask(F(), a, Tuple(BitIndices(C)))
        ex = :(map(muladd, x[$a] .* $mask, y[$inds], $ex))
    end
    return :($C($ex))
end
