module Metrics

#---Metric tensor associated with an orthonormal basis---------------------------------------------#
"""
    OrthonormalMetric <: AbstractVector{Int8}
    OrthonormalMetric(
        dimensions::Integer,
        negative::UInt,
        degenerate::UInt,
        [first_index::Integer = 1]
    )

Contains information about the metric associated with a Clifford algebra or Clifford number.

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
struct OrthonormalMetric <: AbstractVector{Int8}
    dimensions::UInt8
    negative::UInt
    degenerate::UInt
    first_index::Int8
    function OrthonormalMetric(
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

Base.has_offset_axes(::OrthonormalMetric) = true
Base.IndexStyle(::Type{<:OrthonormalMetric}) = IndexLinear()

Base.size(m::OrthonormalMetric) = tuple(m.dimensions)
Base.axes(m::OrthonormalMetric) = tuple(m.first_index .+ (0:m.dimensions-1))

function Base.getindex(m::OrthonormalMetric, i::Int)
    @boundscheck checkbounds(m, i)
    mask = UInt(2)^(i - m.first_index)
    negative_bit = !iszero(m.negative & mask)
    degenerate_bit = !iszero(m.degenerate & mask)
    return Int8(-1)^negative_bit * !degenerate_bit
end

"""
    VGA(dimensions::Integer) -> OrthonormalMetric

Constructs an `OrthonormalMetric` object representing a VGA (vanilla geometric algebra) with the 
given number of dimensions.
"""
VGA(dimensions::Integer) = OrthonormalMetric(dimensions, UInt(0), UInt(0),  1)

"""
    isVGA(m::OrthonormalMetric) -> Bool

Determines if `m` represents a VGA (vanilla geometric algebra). This is accomplished by checking
that all dimensions square to +1, and that the first index is 1: if the first index is less than 1,
it may be assumed that the first dimensions are projective dimensions in a larger modeling space.
"""
isVGA(m::OrthonormalMetric) = iszero(m.negative) && iszero(m.degenerate) && isone(firstindex(m))

"""
    PGA(modeled_dims::Integer) -> OrthonormalMetric

Constructs an `OrthonormalMetric` object representing a PGA (projective geometric algebra) with the
given number of modeled dimensions. The constructed algebra will contain the number of modeled 
dimensions plus one degenerate (zero-squaring) dimension.
"""
PGA(modeled_dims::Integer) = OrthonormalMetric(modeled_dims + 1, 0b0, 0b1, 0)

"""
    CGA(modeled_dims::Integer) -> OrthonormalMetric

Constructs an `OrthonormalMetric` object representing a CGA (conformal geometric algebra) with the
given number of modeled dimensions. The constructed algebra will contain the number of modeled 
dimensions plus one positive-squaring dimension (with index 0) and one negative-squaring dimension
(with index -1).
"""
CGA(modeled_dims::Integer) = OrthonormalMetric(modeled_dims + 2, 0b1, 0b0, -1)

"""
    LGA(spatial_dims::Integer, time_signbit::Bool) -> OrthonormalMetric

Constructs an `OrthonormalMetric` object representing an LGA (Lorentzian geometric algebra) with the
given number of positive-squaring spatial dimensions.

If `time_signbit` is `true`, the spatial dimensions will square to a positive value and the temporal
dimension will square to a negative value; this will be reversed for a `false` input.

By convention, LGAs use the symbol `γ` for their basis blades (related to the standard symbol for
the Dirac matrices).
"""
function LGA(spatial_dims::Integer, time_signbit::Bool)
    return OrthonormalMetric(spatial_dims + 1, ifelse(time_signbit, 0b1, ~0b1), 0b0, 0)
end

"""
    LGAEast(spatial_dims::Integer)

Generates a Lorentzian geometric algebra with the East Coast sign convention (spatial dimensions 
square positive, temporal dimensions square negative). This is equivalent to
[`LGA(spatial_dims, true)`](@ref LGA).

For the spacetime algebra with this convention, use [`STAEast`](@ref).

For the opposite sign convention, use [`LGAWest`](@ref).
"""
LGAEast(spatial_dims::Integer) = LGA(spatial_dims, true)

"""
    LGAWest(spatial_dims::Integer)

Generates a Lorentzian geometric algebra with the West Coast sign convention (spatial dimensions 
square negative, temporal dimensions square positive). This is equivalent to
[`LGA(spatial_dims, false)`](@ref LGA).

For the spacetime algebra with this convention, use [`STAWest`](@ref).

For the opposite sign convention, use [`LGAEast`](@ref).
"""
LGAWest(spatial_dims::Integer) = LGA(spatial_dims, false)

"""
    Exterior(dimensions::Integer) -> OrthonormalMetric

Constructs an `OrthonormalMetric` object representing an exterior algebra. Exterior algebras are 
totally degenerate Clifford algebras: all dimensions square to 0.
"""
Exterior(dimensions::Integer) = OrthonormalMetric(dimensions, 0, ~zero(UInt), 1)

#---Common algebras--------------------------------------------------------------------------------#

# TODO: perhaps cut down on some of the aliases here

"""
    VGA2D (alias for OrthonormalMetric(2, 0b00, 0b00, 1) or VGA(2))

The algebra of 2D space. The even subalgebra of this algebra is isomorphic to ℂ, the complex
numbers.
"""
const VGA2D = VGA(2)

"""
    VGA3D (alias for OrthonormalMetric(2, 0b000, 0b000, 1) or VGA(3))
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
    PGA2D (alias for OrthonormalMetric(3, 0b000, 0b001, 0) or PGA(2))

The projective geometric algebra of 2D space, which represents points and lines on the plane.
"""
const PGA2D = PGA(2)

"""
    PGA3D (alias for OrthonormalMetric(4, 0b0000, 0b0001, 0) or PGA(3))

The projective geometric algebra of 3D space, which represents points, lines, and planes in a
3D space.
"""
const PGA3D = PGA(3)

"""
    CGA2D (alias for OrthonormalMetric(4, 0b0001, 0b0000, 0) or CGA(2))

The conformal geometric algebra of 2D space, which represents points, lines, and circles on the 
plane. This algebra constitutes a framework for Euclidean geometry.

This algebra is isomorphic to [`STAEast`](@ref), and this isomorphism is the reason why the default 
convention for spacetime algebras in this package is the West Coast (mostly negative) convention.
"""
const CGA2D = CGA(2)

"""
    STAEast (alias for OrthonormalMetric(4, 0b0001, 0b0000, 0) or LGAEast(3))

The spacetime algebra using the East Coast sign convention (spatial dimensions square positive,
temporal dimensions square negative), with the temporal dimension at index 0.

This convention is *not* the default STA convention, since this object can  be interpreted as 
[the 2D conformal geometric algebra](@ref CGA2D).
"""
const STAEast = LGAEast(3)

"""
    STAWest (alias for OrthonormalMetric(4, 0b0001, 0b0000, 0) or LGAWest(3))
    const STA = STAWest

The spacetime algebra using the West Coast sign convention (spatial dimensions square negative,
temporal dimensions square positive), with the temporal dimension at index 0. This con

This convention is the default STA convention in this package, since 
[the East Coast convention](@ref STAEast) can be interpreted as 
[the 2D conformal geometric algebra](@ref CGA2D).
"""
const STAWest = LGAWest(3)
const STA = STAWest
@doc (@doc STAWest) STA

"""
    CGA3D (alias for OrthonormalMetric(5, 0b00001, 0b00000, 0) or CGA(3))

The conformal geometric algebra of 3D space, which represents points, lines, planes, circles, and
spheres in 3D space. This algebra extends Euclidean geometry to 3 dimensions.
"""
const CGA3D = CGA(3)

"""
    STAPEast (alias for OrthonormalMetric(5, 0b00010, 0b00001, -1))

The projective spacetime algebra using the East Coast sign convention (spatial dimensions 
square positive, temporal dimensions square negative). The degenerate dimension is at index -1.

As with `STA`, the default convention for `STAP` is the West Coast metric. For an explanation,
see [`STA`](@ref).
"""
const STAPEast = OrthonormalMetric(5, 0b00010, 0b00001, -1)

"""
    STAPWest (alias for OrthonormalMetric(5, 0b11100, 0b00001, -1))
    const STAP = STAPWest

The projective spacetime algebra using the West Coast sign convention (spatial dimensions 
square negative, temporal dimensions square positive). The degenerate dimension is at index -1.

As with `STA`, the default convention for `STAP` is the West Coast metric. For an explanation,
see [`STA`](@ref).
"""
const STAPWest = OrthonormalMetric(5, 0b11100, 0b00001, -1)
const STAP = STAPWest
@doc (@doc STAPWest) STAP

end
