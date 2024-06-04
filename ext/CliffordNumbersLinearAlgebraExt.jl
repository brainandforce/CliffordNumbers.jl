module CliffordNumbersLinearAlgebraExt

using CliffordNumbers
using LinearAlgebra

LinearAlgebra.dot(x::AbstractCliffordNumber, y::AbstractCliffordNumber) = CliffordNumbers.dot(x, y)
LinearAlgebra.normalize(x::AbstractCliffordNumber) = CliffordNumbers.normalize(x)

end
