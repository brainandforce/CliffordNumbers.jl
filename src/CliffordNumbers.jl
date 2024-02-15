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
# nothing to export
include("quadratic.jl")
export QuadraticForm, APS, STA, VGA, PGA, CGA
export dimension, elements, grades
include("numbers.jl")
export AbstractCliffordNumber, CliffordNumber
export pseudoscalar, isscalar, ispseudoscalar
include("bitindex.jl")
export BitIndex, BitIndices
export grade, dual, undual
include("promote.jl")
include("convert.jl")
include("math.jl")
export select_grade, grade_involution, scalar_product, normalize, left_contraction,
    right_contraction, dot, hestenes_product, wedge, ∧, versor_inverse, sandwich, exppi, exptau

include("sparse/kvectors.jl")
export KVector

end
