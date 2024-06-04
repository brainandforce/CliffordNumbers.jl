module LinearAlgebraExt

using CliffordNumbers
using LinearAlgebra

LinearAlgebra.dot(x::AbstractCliffordNumber, y::AbstractCliffordNumber) = CliffordNumbers.dot(x, y)

end
