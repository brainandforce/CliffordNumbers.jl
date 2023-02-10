module CliffordNumbers

"""
    CliffordNumbers.BaseNumber

Union of subtypes of `Number` provided in the Julia `Base` module: `Real` and `Complex`. This
encompasses all types that may be used to construct a `CliffordNumber`.
"""
const BaseNumber = Union{Real,Complex}

"""
    CliffordNumbers.hamming_weight(i::Integer) -> Int

Calculates the Hamming weight of an integer. This is used to determine the grade of a component of a
`StaticMultivector`.
"""
hamming_weight(i::Integer) = sum(!iszero(i & typeof(i)(2)^n) for n in 0:8*sizeof(i) - 1)

"""
    CliffordNumbers.isevil(i::Integer) -> Bool

Determines whether a number is evil, meaning that its Hamming weight (sum of its binary digits) is
even.
"""
isevil(i::Integer) = iseven(hamming_weight(i))

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
export QuadraticForm, VGA, PGA, APS, STA
export dimension, elements, grades
include("numbers.jl")
export CliffordNumber
export pseudoscalar, isscalar, ispseudoscalar
include("indices.jl")
export BitIndex, BitIndices
export grade
include("promote.jl")
include("convert.jl")
include("math.jl")
export dot, ⋆

end
