module Metrics

"""
    Metrics.AbstractSignature <: AbstractVector{Int8}

Supertype for all data types that represent metric signatures. This includes the generic `Signature`
type as well as other specialized types.

Metric signatures can be interpreted as the signs of the diagonal elements of the metric tensor. All
elements are either -1, 0, or 1. All nondegenerate Clifford algebras admit an orthonormal basis 
corresponding to the metric.
"""
abstract type AbstractSignature <: AbstractVector{Int8}
end

Base.IndexStyle(::Type{<:AbstractSignature}) = IndexLinear()
Base.has_offset_axes(::AbstractSignature) = true

"""
    dimension(s::AbstractSignature) -> Int8

Returns the total number of dimensions associated with `s`. The default implementation returns
`signed(s.dimensions)`.

The total number of basis blades is equal to to the size of the power set of all basis vectors, and
is equal to `2^dimension(s)`.
"""
dimension(s::AbstractSignature) = signed(s.dimensions)

"""
    basis_blades(s::AbstractSignature) -> Int64

Returns the total number of blades associated with `s`, which is equal to `2^dimension(s)`.
"""
basis_blades(s::AbstractSignature) = 2^dimensions(s)

"""
    grades(s::AbstractSignature) -> UnitRange{Int8}

Returns the total number of grades associated with `s`, which is equal to `0:dimension(s)`.
"""
grades(s::AbstractSignature) = zero(Int8):dimensions(s)

Base.size(s::AbstractSignature) = tuple(dimension(s))
Base.axes(s::AbstractSignature) = tuple(firstindex(s) .+ (0:dimension(s) - 1))

"""
    blade_symbol(s::AbstractSignature)

Provides the symbol associated to represent basis 1-blades of the geometric algebra with signature
`s`. This defaults to `'e'`, but for `APS` it is `'σ'` and for Lorentzian geometric algebras it is
`'γ'`.
"""
blade_symbol(s::AbstractSignature) = 'e'

#---Metric tensor associated with an orthonormal basis---------------------------------------------#
"""
    Metrics.Signature <: Metrics.AbstractSignature
    Metrics.Signature(
        dimensions::Integer,
        negative::UInt,
        degenerate::UInt,
        [first_index::Integer = 1]
    )

Contains information about the metric associated with a Clifford algebra or Clifford number. This
type is constructed to be as generic as possible; other subtypes of `Metrics.AbstractSignature` may
provide firmer guarantees on behavior, such as `VGA`.

The number which the dimensions square to is stored in a pair of `UInt` fields. The `negative` field
consists of 1 bits for dimensions that square to a negative number, and 0 bits for those squaring
to a positive number, matching the convention of sign bits in signed numbers.

The `degenerate` field consists of 1 bits for degenerate dimensions (dimensions that square to zero)
and 0 bits for nondegenerate dimensions.

The numerical index of the first basis vector is `first_index`, which defaults to 1. Some algebras
conventionally use 0 as the first index, such as projective geometric algebras and Lorentzian
geometric algebras, and in some cases it may be useful to start with a negative index if there are
a larger number of modeling dimensions.
"""
struct Signature <: AbstractSignature
    dimensions::UInt8
    negative::UInt
    degenerate::UInt
    first_index::Int8
    function Signature(
        dimensions::Integer,
        negative::Unsigned,
        degenerate::Unsigned,
        first_index::Integer = 1,
    )
        @assert dimensions < 8*sizeof(UInt) string(
            "An algebra with $dimensions dimensions? You must be nuts!"
        )
        mask = (2^dimensions - 1)
        degenerate = degenerate & mask
        # For simplicity, dimensions aren't marked both negative squaring and degenerate
        negative = negative & mask & ~degenerate
        return new(dimensions, negative, degenerate, first_index)
    end
end

Base.firstindex(s::Signature) = Int(s.first_index)

function Base.getindex(s::Signature, i::Int)
    @boundscheck checkbounds(s, i)
    mask = UInt(2)^(i - firstindex(s))
    negative_bit = !iszero(s.negative & mask)
    degenerate_bit = !iszero(s.degenerate & mask)
    return Int8(-1)^negative_bit * !degenerate_bit
end

"""
    is_degenerate(s::AbstractSignature)

Returns `true` if any basis elements of `s` square to 0.
    
This does not imply that no elements of the associated Clifford algebra square to 0.
"""
is_degenerate(s::Signature) = !iszero(s.degenerate)

"""
    is_positive_definite(s::AbstractSignature)

Returns `true` if all basis 1-blades of `s` square to a positive value.
"""
is_positive_definite(s::Signature) = iszero(s.negative) && !is_degenerate(s)

function Base.show(io::IO, s::Signature)
    println(io, Signature, (s.dimensions, s.negative, s.degenerate, s.first_index))
end

#---Common algebras--------------------------------------------------------------------------------#
"""
    VGA <: Metrics.AbstractSignature

Represents the signature associated with a vanilla geometric algebra (VGA), a positive-definite 
geometric algebra which models space without any projective dimensions.
"""
struct VGA <: AbstractSignature
    dimensions::UInt
end

is_degenerate(::VGA) = false
is_positive_definite(::VGA) = true

Base.has_offset_axes(::VGA) = false
Base.firstindex(::VGA) = 1

Base.getindex(s::VGA, i::Int) = (@boundscheck checkbounds(s, i); return Int8(1))

"""
    isVGA(s::AbstractSignature) -> Bool

Determines if `m` represents a VGA (vanilla geometric algebra). This is accomplished by checking
that all dimensions square to +1, and that the first index is 1: if the first index is less than 1,
it may be assumed that the first dimensions are projective dimensions in a larger modeling space.

Instances of `VGA` will always return `true`.
"""
isVGA(s::VGA) = true
isVGA(s::AbstractSignature) = is_positive_definite(s) && isone(firstindex(s))

# TODO: perhaps cut down on some of the aliases here

"""
    VGA2D (alias for VGA(2))

The algebra of 2D space. The even subalgebra of this algebra is isomorphic to ℂ, the complex
numbers.
"""
const VGA2D = VGA(2)

"""
    VGA3D (alias for VGA(3))
    const APS = VGA3D

The algebra of physical space, a 3D VGA which is commonly used (explicitly and implicitly) to model
non-relativistic physics. It also serves as the subalgebra of both signature conventions of the
spacetime algebra (available as [`STAEast`](@ref) and [`STAWest`](@ref)).

The even subalgebra of this algebra is isomorphic to ℍ, the quaternions.
"""
const VGA3D = VGA(3)
const APS = VGA3D
@doc (@doc VGA3D) APS

"""
    PGA <: Metrics.AbstractSignature

Represents the signature associated with a PGA (projective geometric algebra) with the given number 
of modeled dimensions. The constructed algebra will contain the number of modeled dimensions plus
one degenerate (zero-squaring) dimension represented by e₀. This degenerate dimension corresponds
with the n∞ null vector in CGA (conformal geometric algebra).
"""
struct PGA <: AbstractSignature
    dimensions::UInt
end

dimension(s::PGA) = signed(s.dimensions + 1)
is_degenerate(::PGA) = true
is_positive_definite(::PGA) = false

Base.firstindex(::PGA) = 0

Base.getindex(s::PGA, i::Int) = (@boundscheck checkbounds(s, i); return Int8(!iszero(i)))

"""
    PGA2D (alias for PGA(2))

The projective geometric algebra of 2D space, which represents points and lines on the plane.
"""
const PGA2D = PGA(2)

"""
    PGA3D (alias for PGA(3))

The projective geometric algebra of 3D space, which represents points, lines, and planes in a
3D space.
"""
const PGA3D = PGA(3)

"""
    CGA <: Metrics.AbstractSignature

Represents the signature of a CGA (conformal geometric algebra) with the given number of modeled 
dimensions. The constructed algebra will contain the number of modeled dimensions plus one
positive-squaring dimension and one negative-squaring dimension.

There are two common choices of vector basis for the extra dimensions added when working with CGA.
The most straightforward one is e₊ and e₋, which square to +1 and -1, respectively, and this is what
is used internally, with the negative-squaring dimension being the first one.

However, there is another commonly used basis: define null vectors n₀ = (e₋ - e₊)/2 and
n∞ = e₋ - e₊,  which represent the origin point and the point at infinity, respectively. n∞
corresponds to e₀ in PGA (projective geometric algebra).
"""
struct CGA <: AbstractSignature
    dimensions::UInt
end

dimension(s::CGA) = signed(s.dimensions + 2)
is_degenerate(::CGA) = false
is_positive_definite(::CGA) = false

Base.firstindex(::CGA) = -1

Base.getindex(s::CGA, i::Int) = (@boundscheck checkbounds(s, i); return Int8(-1)^(i < 0))

"""
    CGA2D (alias for CGA(2))

The conformal geometric algebra of 2D space, which represents points, lines, and circles on the 
plane. This algebra constitutes a framework for compass and straightedge constructions.

This algebra is isomorphic to [`STAEast`](@ref), and this isomorphism is the reason why the default 
convention for spacetime algebras in this package is the West Coast (mostly negative) convention.
"""
const CGA2D = CGA(2)

"""
    CGA3D (alias for CGA(3))

The conformal geometric algebra of 3D space, which represents points, lines, and planes, as well as
circles and spheres. This algebra constitutes a framework for extending compass and straightedge 
constructions to 3 dimensions.
"""
const CGA3D = CGA(3)

"""
    LGA{C} <: Metrics.AbstractSignature

Represents the signature of a Lorentzian geometric algebra (LGA), an algebra which models a given
number of spatial dimensions associated with a single time dimension at index 0.

The type parameter `C` corresponds to the sign bit associated with the square of the spatial
1-blades. For convenience, the following aliases are defined:

    const LGAEast = LGA{false}
    const LGAWest = LGA{true}

The names correspond to the "East Coast" and "West Coast" conventions for the metric signature of
spacetime, with the East Coast convention having positive squares for spatial 1-blades and the West
Coast convention having negative squares for spatial 1-blades.
"""
struct LGA{C} <: AbstractSignature
    dimension::UInt
end

const LGAEast = LGA{false}
const LGAWest = LGA{true}

dimension(s::LGA) = signed(s.dimensions + 1)
is_degenerate(::LGA) = false
is_positive_definite(::LGA) = false

Base.firstindex(::LGA) = 0

blade_symbol(::LGA) = 'γ'

function Base.getindex(s::LGA{C}, i::Int) where C
    @boundscheck checkbounds(s, i)
    return Int8(-1)^xor(C, iszero(i))
end

"""
    STAEast (alias for LGAEast(3))

The spacetime algebra using the East Coast sign convention (spatial dimensions square positive,
temporal dimensions square negative), with the temporal dimension at index 0.

This convention is *not* the default STA convention, since this signature is identical to that of
[the 2D conformal geometric algebra](@ref CGA2D).
"""
const STAEast = LGAEast(3)

"""
    STAWest (alias for LGAWest(3))
    const STA = STAWest

The spacetime algebra using the West Coast sign convention (spatial dimensions square negative,
temporal dimensions square positive), with the temporal dimension at index 0.

This convention is the default STA convention in this package, since the signature for
[the East Coast convention](@ref STAEast) can be interpreted as that of
[the 2D conformal geometric algebra](@ref CGA2D).
"""
const STAWest = LGAWest(3)
const STA = STAWest
@doc (@doc STAWest) STA

"""
    Exterior <: Metrics.AbstractSignature

Represents a signature corresponding to an exterior algebra. In an exterior algebra, all 1-blades
square to 0, and the geometric product is equivalent ot the wedge product.

Unlike `VGA`, `PGA`, `CGA`, and `LGA`, the first index is not assumed when constructing this object,
and can be manually specified. If it is not specified, it defaults to 1.
"""
struct Exterior <: AbstractSignature
    dimensions::UInt
    first_index::Int8
end

Exterior(dimensions) = Exterior(dimensions, 1)

is_degenerate(::Exterior) = true
is_positive_definite(::Exterior) = false

Base.firstindex(s::Exterior) = s.first_index

Base.getindex(s::Exterior, i::Int) = (@boundscheck checkbounds(s, i); return Int8(0))

"""
    STAPEast (alias for Signature(5, 0b00010, 0b00001, -1))

The projective spacetime algebra using the East Coast sign convention (spatial dimensions 
square positive, temporal dimensions square negative). The degenerate dimension is at index -1.

As with `STA`, the default convention for `STAP` is the West Coast metric. For an explanation,
see [`STA`](@ref).
"""
const STAPEast = Signature(5, 0b00010, 0b00001, -1)

"""
    STAPWest (alias for Signature(5, 0b11100, 0b00001, -1))
    const STAP = STAPWest

The projective spacetime algebra using the West Coast sign convention (spatial dimensions 
square negative, temporal dimensions square positive). The degenerate dimension is at index -1.

As with `STA`, the default convention for `STAP` is the West Coast metric. For an explanation,
see [`STA`](@ref).
"""
const STAPWest = Signature(5, 0b11100, 0b00001, -1)
const STAP = STAPWest
@doc (@doc STAPWest) STAP

export dimension, basis_blades, grades, is_degenerate, is_positive_definite
export VGA, PGA, CGA, LGA, LGAEast, LGAWest
export VGA2D, VGA3D, PGA2D, PGA3D, CGA2D, CGA3D, STA, STAEast, STAWest, STAP, STAPEast, STAPWest

end
