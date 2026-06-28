#---Efficient multiplication kernels---------------------------------------------------------------#
"""
    CliffordNumbers.bitindex_shuffle(a::BladeIndex{Q}, B::NTuple{L,BladeIndex{Q}})
    CliffordNumbers.bitindex_shuffle(a::BladeIndex{Q}, B::BladeIndices{Q})
    
    CliffordNumbers.bitindex_shuffle(B::NTuple{L,BladeIndex{Q}}, a::BladeIndex{Q})
    CliffordNumbers.bitindex_shuffle(B::BladeIndices{Q}, a::BladeIndex{Q})

Performs the multiplication `-a * b` for each element of `B` for the above ordering, or `-b * a` for
the below ordering, generating a reordered `NTuple` of `BladeIndex{Q}` objects suitable for
implementing a geometric product.
"""
@inline function bitindex_shuffle(a::BladeIndex{Q}, B::NTuple{L,BladeIndex{Q}}) where {L,Q}
    return map(b -> _inv(a) * b, B)
end

@inline function bitindex_shuffle(B::NTuple{L,BladeIndex{Q}}, a::BladeIndex{Q}) where {L,Q}
    return map(b -> b * _inv(a), B)
end

bitindex_shuffle(a::BladeIndex{Q}, B::BladeIndices{Q}) where Q = bitindex_shuffle(a, Tuple(B))
bitindex_shuffle(B::BladeIndices{Q}, a::BladeIndex{Q}) where Q = bitindex_shuffle(Tuple(B), a)

"""
    CliffordNumbers.nondegenerate_mask(a::BladeIndex{Q}, B::NTuple{L,BladeIndex{Q}})

Constructs a Boolean mask which is `false` for any multiplication that squares a degenerate blade;
`true` otherwise.
"""
function nondegenerate_mask(a::BladeIndex{Q}, B::NTuple{L,BladeIndex{Q}}) where {L,Q}
    return map(b -> nondegenerate_mult(a, b), B)
end

#=
"""
    CliffordNumbers.widen_grade_for_mul(x::AbstractCliffordNumber)

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
=#

#---Grade filters----------------------------------------------------------------------------------#
"""
    CliffordNumbers.GradeFilter{S}

A type that can be used to filter certain products of blades in a geometric product multiplication.
The type parameter `S` must be a `Symbol`. The single instance of `GradeFilter{S}` is a callable
object which implements a function that takes two or more `BladeIndex{Q}` objects `a` and `b` and
returns `false` if the product of the blades indexed is zero.

To implement a grade filter for a product function `f`, define the following method:
    (::GradeFilter{:f})(::BladeIndex{Q}, ::BladeIndex{Q})
    # Or if the definition allows for more arguments
    (::GradeFilter{:f})(::BladeIndex{Q}...) where Q
"""
struct GradeFilter{S}
    GradeFilter{S}() where S = (@assert S isa Symbol "Type parameter must be a Symbol."; new())
end

(::GradeFilter{S})(args...) where S = error("This filter has not been implemented.")

(::GradeFilter{:*})(args::BladeIndex{Q}...) where Q = true

(::GradeFilter{:∧})(args::BladeIndex{Q}...) where Q = has_wedge(args...)

(::GradeFilter{:⨼})(a::T, b::T) where {Q,T<:BladeIndex{Q}} = (grade(b) - grade(a)) == grade(a*b)
(::GradeFilter{:⨽})(a::T, b::T) where {Q,T<:BladeIndex{Q}} = (grade(a) - grade(b)) == grade(a*b)

function (::GradeFilter{:dot})(a::T, b::T) where {Q,T<:BladeIndex{Q}}
    return abs(grade(a) - grade(b)) == grade(a*b)
end

const ContractionGradeFilters = Union{GradeFilter{:⨼},GradeFilter{:⨽},GradeFilter{:dot}}

"""
    CliffordNumbers.mul_mask(F::GradeFilter, a::BladeIndex{Q}, B::NTuple{L,BladeIndices{Q}})
    CliffordNumbers.mul_mask(F::GradeFilter, B::NTuple{L,BladeIndices{Q}}, a::BladeIndex{Q})

    CliffordNumbers.mul_mask(F::GradeFilter, a::BladeIndex{Q}, B::BladeIndices{Q})
    CliffordNumbers.mul_mask(F::GradeFilter, B::BladeIndices{Q}, a::BladeIndex{Q})

Generates a `NTuple{L,Bool}` which is `true` whenever the multiplication of the blade indexed by `a`
and blades indexed by `B` is nonzero. `false` is returned if the grades multiply to zero due to the
squaring of a degenerate component, or if they are filtered by `F`.
"""
function mul_mask(F::GradeFilter, a::BladeIndex{Q}, B::NTuple{L,BladeIndex{Q}}) where {L,Q}
    return map(b -> F(a,b) & nondegenerate_mult(a,b), B)
end

function mul_mask(F::GradeFilter, B::NTuple{L,BladeIndex{Q}}, a::BladeIndex{Q}) where {L,Q}
    return map(b -> F(b,a) & nondegenerate_mult(b,a), B)
end

mul_mask(F::GradeFilter, a::BladeIndex{Q}, B::BladeIndices{Q}) where Q = mul_mask(F, a, Tuple(B))
mul_mask(F::GradeFilter, B::BladeIndices{Q}, a::BladeIndex{Q}) where Q = mul_mask(F, Tuple(B), a)

"""
    CliffordNumbers.mul_signs(F::GradeFilter, a::BladeIndex{Q}, B::NTuple{L,BladeIndices{Q}})
    CliffordNumbers.mul_signs(F::GradeFilter, B::NTuple{L,BladeIndices{Q}}, a::BladeIndex{Q})

    CliffordNumbers.mul_signs(F::GradeFilter, a::BladeIndex{Q}, B::BladeIndices{Q})
    CliffordNumbers.mul_signs(F::GradeFilter, B::BladeIndices{Q}, a::BladeIndex{Q})

Generates an `NTuple{L,Int8}` which represents the sign associated with the multiplication needed to
calculate components of a multiplication result.

This is equivalent to `sign.(B)` unless `F === CliffordNumbers.GradeFilter{:dot}()`.
"""
mul_signs(::GradeFilter, ::BladeIndex{Q}, B::NTuple{L,BladeIndex{Q}}) where {L,Q} = sign.(B)
mul_signs(::GradeFilter, B::NTuple{L,BladeIndex{Q}}, ::BladeIndex{Q}) where {L,Q} = sign.(B)

function mul_signs(::GradeFilter{:dot}, a::BladeIndex{Q}, B::NTuple{L,BladeIndex{Q}}) where {L,Q}
    return sign.(B) .* Int8(-1).^(grade.(B) .* (grade(a) .- grade.(B)))
end

function mul_signs(::GradeFilter{:dot}, B::NTuple{L,BladeIndex{Q}}, a::BladeIndex{Q}) where {L,Q}
    return sign.(B) .* Int8(-1).^(grade(a) .* (grade.(B) .- grade(a)))
end

mul_signs(F::GradeFilter, a::BladeIndex{Q}, B::BladeIndices{Q}) where Q = mul_signs(F, a, Tuple(B))
mul_signs(F::GradeFilter, B::BladeIndices{Q}, a::BladeIndex{Q}) where Q = mul_signs(F, Tuple(B), a)

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
    T = promote_scalar_type(C1,C2)
    if (!c1_odd && !c1_even) || (!c2_odd && !c2_even)
        return :(CliffordNumber{Q})
    else
        return :(Z2CliffordNumber{$P,Q})
    end
end

# Extra implementation for k-vectors: special handling scalar and pseudoscalar arguments
# TODO: can we integrate this into the above function?
@generated function product_return_type(
    ::Type{C1},
    ::Type{C2},
    ::GradeFilter{<:Any}
) where {K1,K2,Q,C1<:KVector{K1,Q},C2<:KVector{K2,Q}}
    D = dimension(Q)
    # Handle the scalar and pseudoscalar cases
    if isone(nblades(C1))
        iszero(K1) && return :(KVector{K2,Q})
        K1 == D && return :(KVector{$D-K2,Q})
    elseif isone(nblades(C2))
        iszero(K2) && return :(KVector{K1,Q})
        K2 == D && return :(KVector{$D-K1,Q})
    end
    # Fall back to a Z2CliffordNumber with the right parity
    P = isodd(xor(K1, K2))
    return :(Z2CliffordNumber{$P,Q})
end

function product_return_type(
    ::Type{<:KVector{K1,Q}},
    ::Type{<:KVector{K2,Q}},
    ::GradeFilter{:∧}
) where {Q,K1,K2}
    # TODO: do we need this minimizing behavior?
    # Maybe we can let K exceed the expected value and return a zero-element multivector.
    K = min(K1 + K2, dimension(Q))
    return KVector{K,Q}
end

function product_return_type(
    ::Type{<:KVector{K1,Q}},
    ::Type{<:KVector{K2,Q}},
    ::ContractionGradeFilters
) where {Q,K1,K2}
    K = abs(K1 - K2)
    return KVector{K,Q}
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
    CliffordNumbers.mul(
        x::AbstractCliffordNumber{Q},
        y::AbstractCliffordNumber{Q},
        [F::GradeFilter = GradeFilter{:*}()]
    )

A fast geometric product implementation using generated functions for specific cases, and generic
methods which either convert the arguments or fall back to other methods.

The arguments to this function should all agree in scalar type `T`. The `*` function, which exposes
the fast geometric product implementation, promotes the scalar types of the arguments before
utilizing this kernel. The scalar multiplication operations are implemented using [`muladd`](@ref), 
allowing for hardware fma operations to be used when available.

The `GradeFilter` `F` allows for some blade multiplications to be excluded if they meet certain
criteria. This is useful for implementing products besides the geometric product, such as the wedge
product, which excludes multiplications between blades with shared vectors. Without a filter, this
kernel just returns the geometric product.
"""
@generated function mul(
    x::AbstractCliffordNumber{Q,T},
    y::AbstractCliffordNumber{Q,T},
    F::GradeFilter = GradeFilter{:*}()
) where {Q,T<:BaseNumber}
    C = product_return_type(x, y, F())
    ex = :($(zero_tuple(C)))
    # Generate the expression differently if the first argument is longer
    # This helps leverage SIMD and cuts down on the number of loop iterations
    # However, it seems that a smaller first argument is still *slightly* slower, why?
    # TODO: can we unify the cases and not have to repeat so much code?
    if nblades(x) > nblades(y) 
        for b in BladeIndices(y)
            # Permute the indices of x so that the result coefficients are correctly ordered
            inds = bitindex_shuffle(BladeIndices(C), b)
            # Filter out indexing operations that automatically go to zero
            # This must be done manually since we want to work directly with tuples
            x_mask = map(in, grade.(inds), ntuple(Returns(nonzero_grades(x)), Val(nblades(C))))
            # Filter out multiplications which necessarily go to zero
            y_mask = mul_mask(F(), inds, b)
            mask = x_mask .& y_mask
            if any(mask)
                # Resolve BladeIndex to an integer here to avoid calling Base.to_index at runtime
                # This function cannot be inlined or unrolled for KVector arguments
                # But all values are known at compile time, so interpolate them into expressions
                ib = to_index(y, b)
                tuple_inds = to_index.(x, inds)
                signs = mul_signs(F(), inds, b)
                # Construct the tuples that contribute to the product
                x_tuple_ex = :(getindex.(tuple(Tuple(x)), $tuple_inds))
                y_tuple_ex = :(Tuple(y)[$ib] .* $(signs .* mask))
                # Combine the tuples using muladd operations
                ex = :(map(muladd, $x_tuple_ex, $y_tuple_ex, $ex))
            end
        end
    else
        for a in BladeIndices(x)
            # Permute the indices of y so that the result coefficients are correctly ordered
            inds = bitindex_shuffle(a, BladeIndices(C))
            # Filter out multiplications which necessarily go to zero
            x_mask = mul_mask(F(), a, inds)
            # Filter out indexing operations that automatically go to zero
            # This must be done manually since we want to work directly with tuples
            y_mask = map(in, grade.(inds), ntuple(Returns(nonzero_grades(y)), Val(nblades(C))))
            mask = x_mask .& y_mask
            # Don't append operations that won't actually do anything
            if any(mask)
                # Resolve BladeIndex to an integer here to avoid calling Base.to_index at runtime
                # This function cannot be inlined or unrolled for KVector arguments
                # But all values are known at compile time, so interpolate them into expressions
                ia = to_index(x, a)
                tuple_inds = to_index.(y, inds)
                signs = mul_signs(F(), a, inds)
                # Construct the tuples that contribute to the product
                x_tuple_ex = :(Tuple(x)[$ia] .* $(signs .* mask))
                y_tuple_ex = :(getindex.(tuple(Tuple(y)), $tuple_inds))
                # Combine the tuples using muladd operations
                ex = :(map(muladd, $x_tuple_ex, $y_tuple_ex, $ex))
            end
        end
    end
    return :(($C)($ex))
end

"""
    CliffordNumbers.mul(x::EvenCliffordNumber{Q,T,2}, y::EvenCliffordNumber{Q,T,2}, ::GradeFilter{:*})

Closed-form geometric product for the even subalgebra of a dimension-2 algebra, whose elements
`a + b·e₁₂` form a 2-component spinor (≅ ℂ for `VGA(2)`, the split/dual analogues for other
signatures). Writing `σ = e₁₂²` (the sign of the pseudoscalar's square, `±1` or `0`),

    (a + b·e₁₂)(c + d·e₁₂) = (a·c + σ·b·d) + (a·d + b·c)·e₁₂

The generic kernel lowers this to a `map(muladd, …)` chain seeded with `zero_tuple`, whose leading
`muladd(0, 0, x)` cannot be folded away under IEEE signed-zero rules; emitting the closed form
directly removes that overhead. `x^2` benefits transitively, since `literal_pow(^, x, Val(2))` is
`x * x`. Other products (`∧`, contractions) keep the generic kernel by dispatching only on
`GradeFilter{:*}`.
"""
@generated function mul(
    x::EvenCliffordNumber{Q,T,2},
    y::EvenCliffordNumber{Q,T,2},
    ::GradeFilter{:*} = GradeFilter{:*}()
) where {Q,T<:BaseNumber}
    σ = sign_of_square(last(BladeIndices(EvenCliffordNumber{Q,T,2})))
    re = if σ > 0
        :(muladd(a, c, b * d))
    elseif σ < 0
        :(muladd(a, c, -(b * d)))
    else
        :(a * c)
    end
    return quote
        (a, b) = Tuple(x)
        (c, d) = Tuple(y)
        return EvenCliffordNumber{Q,T,2}(($re, muladd(a, d, b * c)))
    end
end

# NOTE: a closed-form Hamilton product for `EvenCliffordNumber{VGA(3),T,4}` (≅ ℍ) was tried and
# dropped. Unlike the 2-component spinor, the generic kernel below already SIMD-vectorizes the
# 4-component product (`map(muladd, …)` over 4-tuples), so removing its `zero_tuple` seed gave no
# benefit and a scalar/broadcast closed form measured *slower* in the compute-bound reduce case
# (chain-compose ≈ 5.3 ns/mul generic vs ≈ 7.4–7.9 ns/mul closed). See bench/results_spinor.md.

#= Update (2024-05-07)
    The performance bottlenecks for KVector arguments have been (mostly) resolved.
    This was done by resolving all BladeIndex objects to tuple indices at compile time, avoiding
    all calls to Base.to_index(::KVector, ::Int) at runtime (since this function cannot be unrolled
    or inlined at compile time).

    There might be some lingering performance issues with the order of differently sized arguments.
=#

function mul(
    x::AbstractCliffordNumber{Q},
    y::AbstractCliffordNumber{Q},
    F::GradeFilter = GradeFilter{:*}()
) where Q
    return @inline mul(scalar_promote(x, y)..., F)
end

#= Rational multiplication is slow without some optimization

function gcd_rescale(t::Tuple{Vararg{Rational{T}}}) where T
    g = lcm(denominator.(t)...)
    return (convert.(T, t .* g), g)
end

gcd_rescale(t::Tuple{Vararg{T}}) where T<:Integer = (t, one(T))

function gcd_rescale(x::AbstractCliffordNumber{Q,Rational{T}}) where {Q,T}
    C = similar_type(x, T)
    (t, g) = gcd_rescale(x)
    return (C(t), g)
end
=#

function mul(
    x::AbstractCliffordNumber{Q,Rational{S}},
    y::AbstractCliffordNumber{Q,Rational{T}},
    F::GradeFilter = GradeFilter{:*}()
) where {Q,S<:Integer,T<:Integer}
    gx = lcm(denominator.(Tuple(x))...)
    gy = lcm(denominator.(Tuple(y))...)
    sx = scalar_convert(promote_type(S, T), x * gx)
    sy = scalar_convert(promote_type(S, T), y * gy)
    return (@inline mul(sx, sy, F)) // (gx * gy)
end

function mul(
    x::AbstractCliffordNumber{Q,Rational{S}},
    y::AbstractCliffordNumber{Q,T},
    F::GradeFilter = GradeFilter{:*}()
) where {Q,S<:Integer,T<:Integer}
    gx = lcm(denominator.(Tuple(x))...)
    sx = scalar_convert(promote_type(S, T), x * gx)
    return (@inline mul(sx, y, F)) // gx
end

function mul(
    x::AbstractCliffordNumber{Q,S},
    y::AbstractCliffordNumber{Q,Rational{T}},
    F::GradeFilter = GradeFilter{:*}()
) where {Q,S<:Integer,T<:Integer}
    gy = lcm(denominator.(Tuple(y))...)
    sy = scalar_convert(promote_type(S, T), y * gy)
    return (@inline mul(x, sy, F)) // gy
end

# Throw an error if the algebras are different

function mul(
    x::AbstractCliffordNumber,
    y::AbstractCliffordNumber,
    ::GradeFilter{S} = GradeFilter{:mul}()
) where S
    throw(AlgebraMismatch(S, (x, y)))
end
