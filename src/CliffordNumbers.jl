module CliffordNumbers

import Base: length, size, axes, getindex, to_index
import Base: zero, one, oneunit
import Base: signbit, sign
import Base: reverse, conj, ~
import Base: ==, isapprox, +, -, *, /, //, muladd, abs, abs2, exp, ^
import Base: promote_rule, convert, similar
import Base: summary, show

import Base.Broadcast

"""
    CliffordNumbers.BaseNumber

Union of subtypes of `Number` provided in the Julia `Base` module: `Real` and `Complex`. This
encompasses all types that may be used to construct a `CliffordNumber`.
"""
const BaseNumber = Union{Real,Complex}

# Contains tools for working with Hamming weights of integers
include("hamming.jl")
# New module for metric signatures
include("metrics.jl")
export Metrics
import .Metrics: dimension, blade_count, grades, is_degenerate, is_positive_definite
# Tools for defining quadratic forms/metric signatures (convention may not be great)
include("quadratic.jl")
export QuadraticForm, APS, STA, VGA, PGA, CGA
export dimension, blade_count, grades
# Abstract supertype for all Clifford numbers
include("abstract.jl")
export AbstractCliffordNumber
export signature, numeric_type
# Working with the grades represented by an AbstractCliffordNumber subtype
include("grades.jl")
export RepresentedGrades
export nonzero_grades, has_grades_of
# Indexing each graded element of an AbstractCliffordNumber
include("bitindex.jl")
export BitIndex
export grade, scalar_index, pseudoscalar_index, grade_involution, dual, undual
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
# Conversion and promotion tools
include("convert.jl")
export scalar_convert
include("promote.jl")
export widen_grade, scalar_promote
# Fast multiplication kernels
include("multiply.jl")
# Mathematical operations defined for all AbstractCliffordNumber instances
include("math.jl")
export isscalar, ispseudoscalar, scalar, select_grade, scalar_product, normalize, left_contraction, 
    right_contraction, dot, hestenes_product, wedge, commutator, anticommutator, versor_inverse, 
    sandwich, exppi, exptau
export ⨼, ⨽, ∧, ×, ⨰
# Pretty printing
include("show.jl")

end
