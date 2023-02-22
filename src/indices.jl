#---Struct to work with bit indices----------------------------------------------------------------#

"""
    CliffordNumbers.BitIndex{Cl<:QuadraticForm}

An index corresponding to a signed basis element of a Clifford algebra.

The basis elements of a Clifford number may be described with a binary word with length equal to the
dimension of the space. The nth bit corresponds to whether the basis element e_n is used to
construct the element with the given index.

As an example, in the algebra of physical space, the element e₁ can be represented with the binary
number `0b001`, e₂ with `0b010`, and e₃ with `0b100`. Elements of higher grades can be accessed by
using the `xor` operation: for example, the basis bivector e₁e₃ can be represented with
`xor(0b001, 0b100) == 0b101`. The scalar element is always represented with `0`.

The sign of the index is stored. A negative sign corresponds to an odd permutation of the basis
elements, and a positive sign corresponds to an even permutation.

The `QuadraticForm` is a tag that marks the Clifford algebra associated with the index. All
operations on such `BitIndex` only depend on the value modulo `elements(Cl)`. If a generic index is
desired, `QuadraticForm` (without type parameters) may be used as the tag.
"""
struct BitIndex{Cl<:QuadraticForm}
    i::Int
    # Essentially, this is a conversion to a one's complement representation,
    # since we need a signed zero to represent a negative scalar
    function BitIndex{QuadraticForm{P,Q,R}}(i::Signed) where {P,Q,R}
        sz = 2^(P+Q+R)
        return new(ifelse(signbit(i), xor(typemin(Int), mod(abs(i), sz)), i % sz))
    end
    # Used when the quadratic form dimensions are unspecified
    BitIndex{QuadraticForm}(i::Signed) = new(ifelse(signbit(i), typemin(i) | -i, i))    
end

# Useful for direct binary construction
BitIndex{Cl}(i::Unsigned) where Cl = BitIndex{Cl}(signed(i + signbit(signed(i))))

const GenericBitIndex = BitIndex{QuadraticForm}
BitIndex(i::Integer) = GenericBitIndex(i)

function Base.show(io::IO, i::BitIndex{Cl}) where Cl
    print(
        io, "BitIndex{", Cl, "}(", "~"^signbit(i.i),
        "0b", bitstring(i.i)[end-dimension(Cl)+1:end], ")"
    )
end

function Base.show(io::IO, i::GenericBitIndex)
    print(io, "BitIndex(", signbit(i.i) ? "-" : "+", "0b" * bitstring(i.i)[2:end] * ")")
end

#=
function Base.show(io::IO, ::MIME"text/plain", i::BitIndex)
    println(io, "BitIndex(0b" * bitstring(i.i) * "):")
end
=#

# Define all Boolean operators on BitIndex
for fn in (:&, :|, :~, :⊻, :⊼, :⊽, :>>>, :>>, :<<)
    @eval begin
        Base.$fn(i1::BitIndex{Cl}, i2::BitIndex{Cl}) where Cl = BitIndex{Cl}($fn(i1.i, i2.i))
    end
end

Base.signbit(i::BitIndex) = signbit(i.i)
Base.sign(i::BitIndex) = Int8(-1)^signbit(i)
Base.:-(i::BitIndex{Cl}) where Cl = BitIndex{Cl}(-i.i)
Base.:abs(i::BitIndex{Cl}) where Cl = BitIndex{Cl}(i.i & typemax(Int))

#---BitIndices represents the valid index range for a Clifford number------------------------------#

"""
    BitIndices{Cl<:QuadraticForm}

Stores the range of valid `BitIndices` for the Clifford algebra with quadratic form `Cl`.
"""
struct BitIndices{Cl<:QuadraticForm} # <: AbstractUnitRange{BitIndex{Cl}}
end

BitIndices(Cl::Type{<:QuadraticForm}) = BitIndices{Cl}()
BitIndices(::CliffordNumber{Cl}) where Cl = BitIndices{Cl}()

Base.length(::BitIndices{Cl}) where Cl = elements(Cl)
Base.first(::BitIndices{Cl}) where Cl = BitIndex{Cl}(0)
Base.last(::BitIndices{Cl}) where Cl = BitIndex{Cl}(typemax(UInt))
# TODO: Do we need to define Base.step(::BitIndices{Cl})?

function Base.iterate(::BitIndices{Cl}, state::Integer=0) where Cl
    return 0 <= state < elements(Cl) ? (BitIndex{Cl}(state), state + 1) : nothing
end

#---Extending mathematical operations on Clifford numbers to their indices-------------------------#

"""
    CliffordNumbers.sign_of_mult([Cl::Type{QuadraticForm{P,Q,R}}], a::Integer, b::Integer)

Calculates the sign associated with multiplying basis elements indexed with bit indices supplied as
integers. The sign reverses whenever the order of `a` and `b` are reverse, provided they are not the
same.
"""
function sign_of_mult(a::Integer, b::Integer)
    s = signbit(a*b)
    a = abs(a) >>> 1
    sum = 0
    while !iszero(a)
        sum += hamming_weight(a & abs(b))
        a = a >>> 1
    end
    return Int8(-1)^(!iszero(sum & 1)) * Int8(-1)^s
end

sign_of_mult(a::Integer) = sign_of_mult(a,a)

"""
    *(i1::BitIndex{Cl}, i2::BitIndex{Cl}) -> BitIndex{Cl}

Returns the basis element that results from the multiplication of basis elements `i1` and `i2`.
"""
function Base.:*(i1::BitIndex{Cl}, i2::BitIndex{Cl}) where Cl
    return BitIndex{Cl}(sign_of_mult(i1.i, i2.i) * xor(i1.i, i2.i))
end

function Base.:*(i1::BitIndex{QuadraticForm}, i2::BitIndex{QuadraticForm})
    return BitIndex(xor(i1.i, i2.i, typemin(Int) * signbit(sign_of_mult(i1.i, i2.i))))
end

#---Grades of indices-----------------------------------------------------------------------------#

grade(i::BitIndex{Cl}) where Cl = hamming_weight(i.i)

#---Indexing of Clifford numbers (zero-based)-----------------------------------------------------#

Base.getindex(m::CliffordNumber{Cl}, i::BitIndex{Cl}) where Cl = sign(i) * m.data[abs(i).i + 1]

function Base.getindex(m::CliffordNumber{Cl}, i::GenericBitIndex) where Cl
    return sign(i) * m.data[mod(abs(i).i, elements(Cl)) + 1]
end

function Base.getindex(m::CliffordNumber, i::Integer)
    return i in 0:length(m)-1 ? m.data[i+1] : throw(BoundsError(m, i))
end

Base.keys(::Type{<:CliffordNumber{Cl}}) where Cl = BitIndices{Cl}()
Base.keys(m::CliffordNumber{Cl}) where Cl = keys(typeof(m))
