module CliffordNumbersLinearAlgebraExt

using CliffordNumbers
using LinearAlgebra

LinearAlgebra.dot(x::AbstractCliffordNumber, y::AbstractCliffordNumber) = CliffordNumbers.dot(x, y)
LinearAlgebra.normalize(x::AbstractCliffordNumber) = CliffordNumbers.normalize(x)
LinearAlgebra.det(x::AbstractCliffordNumber) = CliffordNumbers.det(x)
LinearAlgebra.tr(x::AbstractCliffordNumber) = CliffordNumbers.tr(x)

end
