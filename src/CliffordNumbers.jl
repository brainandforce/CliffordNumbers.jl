module CliffordNumbers

"""
    CliffordNumbers.BaseNumber

Union of subtypes of `Number` provided in the Julia `Base` module: `Real` and `Complex`. This
encompasses all types that may be used to construct a `CliffordNumber`.
"""
const BaseNumber = Union{Real,Complex}

"""
    CliffordNumbers.isevil(i::Integer) -> Bool

Determines whether a number is evil, meaning that its Hamming weight (sum of its binary digits) is
even.
"""
isevil(i::Integer) = iseven(count_ones(i))

"""
    CliffordNumbers.isodious(i::Integer) -> Bool

Determines whether a number is odious, meaning that its Hamming weight (sum of its binary digits) is
odd.
"""
isodious(i::Integer) = !isevil(i)

"""
    CliffordNumbers.evil_number(n::Integer)

Returns the nth evil number, with the first evil number (`n == 1`) defined to be 0.

Evil numbers are numbers which have an even Hamming weight (sum of its binary digits).
"""
evil_number(n::Integer) = (x = 2*(n-1); return x + !isevil(x))

"""
    CliffordNumbers.odious_number(n::Integer)

Returns the nth odious number, with the first odious number (`n == 1`) defined to be 1.

Odious numbers are numbers which have an odd Hamming weight (sum of its binary digits).
"""
odious_number(n::Integer) = (x = 2*(n-1); return x + !isodious(x))

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
