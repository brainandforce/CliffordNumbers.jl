module CliffordNumbers

import Base: length, size, axes, getindex, to_index
import Base: zero, one, oneunit, iszero, isone
import Base: signbit, sign, copysign, flipsign
import Base: reverse, adjoint, conj, ~
import Base: ==, isequal, isapprox, +, -, *, /, //, muladd, abs, abs2, exp, ^
import Base: promote_rule, convert, similar, float, big
import Base: print, show, summary

import Base.Broadcast

"""
    CliffordNumbers.BaseNumber

Union of subtypes of `Number` provided in the Julia `Base` module: `Real` and `Complex`. This
encompasses all types that may be used to construct a `CliffordNumber`.
"""
const BaseNumber = Union{Real,Complex}

# Contains tools for working with Hamming weights of integers
include("hamming.jl")
using .Hamming
# New module for metric signatures
include("metrics.jl")
using .Metrics
import .Metrics: dimension, blade_count, grades, is_degenerate, is_positive_definite
export Metrics
export Signature, VGA, PGA, CGA, LGA, LGAEast, LGAWest, Exterior
export dimension, blade_count, grades, is_degenerate, is_positive_definite
export VGA2D, VGA3D, PGA2D, PGA3D, CGA2D, CGA3D, STA, STAEast, STAWest, STAP, STAPEast, STAPWest
# Abstract supertype for all Clifford numbers
include("abstract.jl")
export AbstractCliffordNumber
export nblades, signature, scalar_type
# Working with the grades represented by an AbstractCliffordNumber subtype
include("grades.jl")
export nonzero_grades, has_grades_of
# Indexing each graded element of an AbstractCliffordNumber
include("bitindex.jl")
export BitIndex
export grade, scalar_index, pseudoscalar_index, grade_involution, dual, undual, left_complement,
    right_complement
include("bitindices.jl")
export AbstractBitIndices, BitIndices, TransformedBitIndices, ReversedBitIndices,
    GradeInvolutedBitIndices, ConjugatedBitIndices
# Dense representation of a Clifford number
include("cliffordnumber.jl")
export CliffordNumber
# Representations of even/odd grade Clifford numbers
include("even.jl")
export EvenCliffordNumber, OddCliffordNumber
# Compact representation of k-vectors
include("kvector.jl")
export KVector
export pseudoscalar
# Conversion and promotion tools
include("convert.jl")
export scalar_convert
include("promote.jl")
export widen_grade, scalar_promote
# Mathematical operations defined for all AbstractCliffordNumber instances
# Fast multiplication kernels
include("math/multiply.jl")
# Addition, subtraction, (approximate) equality
include("math/arithmetic.jl")
# Grade automorphisms, complements, duals, etc.
include("math/duals.jl")
# Scalar products, absolute values, scalar multiplication
include("math/scalar.jl")
export isscalar, scalar, ispseudoscalar, scalar_product, normalize
# Definitions of operators for products
include("math/products.jl")
export left_contraction, right_contraction, wedge, regressive, commutator, anticommutator
export ⨼, ⨽, ∧, ∨, ×, ⨰
# Inverses, if they exist
include("math/inverse.jl")
export versor_inverse
# Exponentiation
include("math/exponential.jl")
export exppi, exptau
# Pretty printing
include("show.jl")
# Generate basis variables with a macro
include("variables.jl")
export @basis_vars

end
