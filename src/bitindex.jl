"""
    CliffordNumbers.signmask([T::Type{<:Integer} = UInt], [signbit::Bool = true]) -> T

Generates a signmask, or a string of bits where the only 1 bit is the sign bit. If `signbit` is set
to false, this returns zero (or whatever value is represented by all bits being 0).
"""
signmask(T::Type{<:Integer}, signbit::Bool = true) = bitreverse(T(signbit))
signmask(x::Integer, signbit::Bool = true) = bitreverse(typeof(x)(signbit))
signmask(signbit::Bool = true) = signmask(UInt, signbit)

"""
    BitIndex{Q<:QuadraticForm}

A representation of an index corresponding to a basis blade of the geometric algebra with quadratic
form `Q`.
"""
struct BitIndex{Q<:QuadraticForm}
    i::UInt
    BitIndex{Q}(i::Unsigned) where Q = new(UInt(i % elements(Q)) | (UInt(i) & signmask(i)))
end

# Construct with the sign bit as a separate argument
@inline function BitIndex{Q}(signbit::Bool, blade::Unsigned) where Q
    return BitIndex{Q}(UInt(blade % elements(Q)) | signmask(UInt, signbit))
end

Base.UInt(i::BitIndex) = getfield(i, :i)
Base.Int(i::BitIndex) = Int(UInt(i) & ~signmask(Int))

const GenericBitIndex = BitIndex{QuadraticForm}

#---Treat as scalar for the sake of broadcasting---------------------------------------------------#

Broadcast.broadcastable(b::BitIndex) = tuple(b)

#---Convenience constructors-----------------------------------------------------------------------#
"""
    CliffordNumbers._sort_with_parity!(v::AbstractVector{<:Real}) -> Tuple{typeof(v),Bool}

Performs a parity-tracking insertion sort of `v`, which modifies `v` in place. The function returns
a tuple containing `v` and the parity, which is `true` for an odd permutation and `false` for an
even permutation. This is implemented with a modified insertion sort algorithm.
"""
function _sort_with_parity!(v::AbstractVector{<:Real})
    swaps = 0
    # Algorithm adapted from Julia Base.Sort._sort!
    for i in axes(v,1)[2:end]
        j = i
        x = v[i]
        # j will decrement until it equals the first index of v
        while j > firstindex(v)
            # Look behind element i
            y = v[j-1]
            # It's sorted when x is greater than or equal to y
            x < y || break
            # Insert y where x was
            v[j] = y
            # Move to the previous index for comparison
            j -= 1
            # Track number of swaps
            swaps += 1
        end
        # Complete the swap of x and y values
        v[j] = x
    end
    return (v, isodd(swaps))
end

function _sort_with_parity(t::NTuple{L}) where L
    (v, sign) = _sort_with_parity!(collect(t))
    return (ntuple(i -> v[i], Val(L)), sign)
end

function _bitindex(Q::Type{<:QuadraticForm}, t::NTuple)
    @assert all(in(1:elements(Q)), t) "1-vector indices are between 1 and $(elements(Q))."
    (t, parity) = _sort_with_parity(t)
    i = signmask(UInt, parity)
    for x in t
        # Pairs of identical elements will be removed with xor
        i = xor(i, 2^(x-1))
    end
    return BitIndex{Q}(i)
end

"""
    BitIndex(Q::Type{<:QuadraticForm}, i::Integer...)
    BitIndex(x, i::Integer...) = BitIndex(QuadraticForm(x), i...)

Constructs a `BitIndex{Q}` from a list of integers that represent the basis 1-vectors of the space.
`Q` can be determined from the `QuadraticForm` associated with `x`, whether it be a type or object.

This package uses a lexicographic convention for basis blades: in the algebra of physical space, the
basis bivectors are {e₁e₂, e₁e₃, e₂e₃}. The sign of the `BitIndex{Q}` is negative when the parity of
the basis vector permutation is odd.
"""
BitIndex(::Type{Q}, i::Integer...) where Q<:QuadraticForm = _bitindex(Q, promote(i...))
BitIndex(x, i::Integer...) = BitIndex(QuadraticForm(x), i...)

#---Show method------------------------------------------------------------------------------------#

function show(io::IO, b::BitIndex{Q}) where Q
    print(io, "-"^signbit(b), BitIndex, "(", Q, iszero(Int(b)) ? ")" : ", ")
    iszero(UInt(b)) && return nothing
    found_first_vector = false
    for a in 1:min(dimension(Q), 8*sizeof(UInt) - 1)
        if !iszero(UInt(b) & 2^(a-1))
            found_first_vector && print(io, ", ")
            print(io, a)
            found_first_vector = true
        end
    end
    print(io, ")")
end

#---Other useful functions-------------------------------------------------------------------------#

signbit(i::BitIndex) = !iszero(UInt(i) & signmask(UInt))
sign(i::BitIndex) = Int8(-1)^signbit(i)

-(i::BitIndex) = typeof(i)(xor(signmask(UInt), UInt(i)))
abs(i::BitIndex) = typeof(i)(UInt(i) & ~signmask(UInt))

grade(i::BitIndex) = count_ones(UInt(i) & ~signmask(UInt))

"""
    CliffordNumbers.is_same_blade(a::BitIndex{Q}, b::BitIndex{Q})

Checks if `a` and `b` perform identical indexing up to sign.
"""
is_same_blade(a::T, b::T) where T<:BitIndex = (UInt(a) << 1) == (UInt(b) << 1)
is_same_blade(a::BitIndex, b::BitIndex) = false

"""
    scalar_index(x::AbstractCliffordNumber{Q}) -> BitIndex{Q}()

Constructs the `BitIndex` used to obtain the scalar (grade zero) portion of `x`.
"""
scalar_index(::Type{Q}) where Q<:QuadraticForm = BitIndex{Q}(false, UInt(0))
scalar_index(x) = scalar_index(QuadraticForm(x))

"""
    pseudoscalar_index(x::AbstractCliffordNumber{Q}) -> BitIndex{Q}

Constructs the `BitIndex` used to obtain the pseudoscalar (highest grade) portion of `x`.
"""
pseudoscalar_index(::Type{Q}) where Q<:QuadraticForm = BitIndex{Q}(false, typemax(UInt))
pseudoscalar_index(x) = pseudoscalar_index(QuadraticForm(x))

#---Grade dependent sign inversion-----------------------------------------------------------------#
"""
    reverse(i::BitIndex) -> BitIndex
    reverse(x::AbstractCliffordNumber) -> typeof(x)

Performs the reverse operation on the basis blade indexed by `b` or the Clifford number `x`. The 
sign of the reverse depends on the grade, and is positive for `g % 4 in 0:1` and negative for
`g % 4 in 2:3`.
"""
reverse(i::BitIndex) = typeof(i)(xor(signbit(i), !iszero(grade(i) & 2)), UInt(i))

"""
    grade_involution(i::BitIndex) -> BitIndex
    grade_involution(x::AbstractCliffordNumber) -> typeof(x)

Calculates the grade involution of the basis blade indexed by `b` or the Clifford number `x`. This
effectively reflects all of the basis vectors of the space along their own mirror operation, which
makes elements of odd grade flip sign.
"""
grade_involution(i::BitIndex) = typeof(i)(xor(signbit(i), isodd(grade(i))), UInt(i))

"""
    conj(i::BitIndex) -> BitIndex
    conj(x::AbstractCliffordNumber) -> typeof(x)

Calculates the Clifford conjugate of the basis blade indexed by `b` or the Clifford number `x`. This
is equal to `grade_involution(reverse(x))`.
"""
conj(i::BitIndex) = typeof(i)(xor(signbit(i), !iszero(grade(i)+1 & 2)), UInt(i))

#---Multiplication tools---------------------------------------------------------------------------#
"""
    CliffordNumbers.signbit_of_mult(a::Integer, [b::Integer]) -> Bool
    CliffordNumbers.signbit_of_mult(a::BitIndex, [b::BitIndex]) -> Bool

Calculates the sign bit associated with multiplying basis elements indexed with bit indices supplied
as either integers or `BitIndex` instances. The sign bit flips when the order of `a` and `b` are
reversed, unless `a === b`. 

As with `Base.signbit()`, `true` represents a negative sign and `false` a positive sign. However,
in degenerate metrics (such as those of projective geometric algebras) the sign bit may be
irrelevant as the multiplication of those basis blades would result in zero.
"""
function signbit_of_mult(a::Unsigned, b::Unsigned)
    a = abs(a) >>> 1
    sum = 0
    while !iszero(a)
        sum += count_ones(a & abs(b))
        a = a >>> 1
    end
    return !iszero(sum & 1)
end

# Account for the sign bits of signed integers
signbit_of_mult(a::Integer, b::Integer) = xor(signbit_of_mult(unsigned.(a,b)...), signbit(xor(a,b)))
signbit_of_mult(a::GenericBitIndex, b::GenericBitIndex) = signbit_of_mult(UInt(a), UInt(b))

function signbit_of_mult(
    a::BitIndex{QuadraticForm{P,Q,R}},
    b::BitIndex{QuadraticForm{P,Q,R}}
) where {P,Q,R}
    base_signbit = xor(signbit_of_mult(UInt(a), UInt(b)), signbit(a), signbit(b))
    iszero(Q) && return base_signbit
    # Only perform this test for pseudo-Riemannian metrics
    q = sum(UInt(2)^(n-1) for n in P .+ (1:Q); init=0)
    return xor(base_signbit, !isevil(UInt(a) & UInt(b) & q))
end

signbit_of_mult(i) = signbit_of_mult(i,i)

"""
    CliffordNumbers.nondegenerate_mult(a::T, b::T) where T<:BitIndex{QuadraticForm{P,Q,R}} -> Bool

Returns `false` if the product of `a` and `b` is zero due to the squaring of a degenerate component,
`true` otherwise. This function always returns `true` if `R === 0`.
"""
function nondegenerate_mult(
    a::BitIndex{QuadraticForm{P,Q,R}},
    b::BitIndex{QuadraticForm{P,Q,R}}
) where {P,Q,R}
    # This mask filters out the nondegenerate components, which are the highest bits
    mask = -UInt(2)^(P+Q)
    return iszero(UInt(a) & UInt(b) & mask)
end

nondegenerate_mult(a::BitIndex{Q}, b::BitIndex{Q}) where Q<:QuadraticForm{<:Any,<:Any,0} = true

"""
    CliffordNumbers.sign_of_mult(a::T, b::T) where T<:BitIndex{QuadraticForm{P,Q,R}} -> Int8

Returns an `Int8` that carries the sign associated with the multiplication of two basis blades of
Clifford/geometric algebras of the same quadratic form.
"""
function sign_of_mult(
    a::BitIndex{QuadraticForm{P,Q,R}},
    b::BitIndex{QuadraticForm{P,Q,R}}
) where {P,Q,R}
    base_signbit = signbit_of_mult(UInt(a), UInt(b))
    # For Euclidean spaces no further processing is needed
    iszero(R) && return Int8(-1)^base_signbit
    # If any dimension squares to zero, just return zero
    r = sum(UInt(2)^(n-1) for n in (P + Q) .+ (1:R); init=0)
    return iszero(UInt(a) & UInt(b) & r) ? Int8(-1)^base_signbit : Int8(0)
end

sign_of_mult(a::GenericBitIndex, b::GenericBitIndex) = Int8(-1)^signbit_of_mult(a,b)
sign_of_mult(i) = sign_of_mult(i,i)

#---Multiplication and duals-----------------------------------------------------------------------#
"""
    *(a::BitIndex{Q}, b::BitIndex{Q}) -> BitIndex{Q}

Returns the `BitIndex` corresponding to the basis blade resulting from the geometric product of the
basis blades indexed by `a` and `b`.
"""
@inline *(a::T, b::T) where T<:BitIndex = T(signbit_of_mult(a,b), xor(UInt(a), UInt(b)))

"""
    CliffordNumbers.has_wedge(a::BitIndex{Q}, b::BitIndex{Q}) -> Bool

Returns `true` if the basis blades indexed by `a` and `b` may have a nonzero wedge product.
"""
has_wedge(a::BitIndex{Q}, b::BitIndex{Q}) where Q = iszero(UInt(a) << 1 & UInt(b) << 1)

dual(b::BitIndex{Q}) where Q = b * BitIndex{Q}(false, typemax(UInt))
undual(b::BitIndex{Q}) where Q = b * BitIndex{Q}(!iszero(dimension(Q) & 2), typemax(UInt))
