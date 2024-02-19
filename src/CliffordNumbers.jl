module CliffordNumbers

"""
    CliffordNumbers.BaseNumber

Union of subtypes of `Number` provided in the Julia `Base` module: `Real` and `Complex`. This
encompasses all types that may be used to construct a `CliffordNumber`.
"""
const BaseNumber = Union{Real,Complex}

"""
    CliffordNumbers.subscript_string(x::Number) -> String

Produces a string representation of a number in subscript format.
"""
function subscript_string(x::Number)
    str = collect(string(x))
    for (n,c) in enumerate(str)
        ('0' <= c <= '9') && (str[n] = c + 0x2050)
        (c === '-') && (str[n] = '₋')
        (c === '+') && (str[n] = '₊')
    end
    return String(str)
end

# Contains tools for working with Hamming weights of integers
include("hamming.jl")
# Tools for defining quadratic forms/metric signatures (convention may not be great)
include("quadratic.jl")
export QuadraticForm, APS, STA, VGA, PGA, CGA
export dimension, elements, grades
# Abstract supertype for all Clifford numbers
include("abstract.jl")
export AbstractCliffordNumber
export numeric_type
# Working with the grades represented by an AbstractCliffordNumber subtype
include("grades.jl")
export RepresentedGrades
export nonzero_grades
# Indexing each graded element of an AbstractCliffordNumber
include("bitindex.jl")
export BitIndex
export grade, grade_involution, dual, undual
include("bitindices.jl")
export AbstractBitIndices, BitIndices, TransformedBitIndices
# Dense representation of a Clifford number
include("cliffordnumber.jl")
export CliffordNumber
export pseudoscalar, isscalar, ispseudoscalar
include("promote.jl")
include("convert.jl")
# Mathematical operations defined for all AbstractCliffordNumber instances
include("math.jl")
export select_grade, scalar_product, normalize, left_contraction,
    right_contraction, dot, hestenes_product, wedge, ∧, versor_inverse, sandwich, exppi, exptau
# Compact representation of k-vectors
include("sparse/kvectors.jl")
export KVector

end
