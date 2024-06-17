#---Geometric product, wedge product, and contractions---------------------------------------------#
for op in (:*, :∧, :⨼, :⨽, :dot)
    @eval begin
        @inline function $op(x::AbstractCliffordNumber, y::AbstractCliffordNumber)
            return @inline mul(x, y, GradeFilter{$(Expr(:quote, op))}())
        end
    end
end

# Parenthetical multiplication
(x::AbstractCliffordNumber)(y::AbstractCliffordNumber) = @inline x * y

# Scalars in wedge products
∧(x::BaseNumber, y::BaseNumber) = x * y
∧(x::BaseNumber, y::AbstractCliffordNumber) = x * y
∧(x::AbstractCliffordNumber, y::BaseNumber) = x * y

# Optimized versions for 0-blade arguments
for op in (:*, :∧)
    @eval begin
        @inline $op(k::KVector{0,Q}, x::AbstractCliffordNumber{Q}) where Q = scalar(k) * x
        @inline $op(x::AbstractCliffordNumber{Q}, k::KVector{0,Q}) where Q = x * scalar(k)
        @inline $op(k::KVector{0,Q}, l::KVector{0,Q}) where Q = KVector{0,Q}(scalar(k) * scalar(l))
    end
end

"""
    ∨(x::AbstractCliffordNumber, y::AbstractCliffordNumber)
    regressive(x::AbstractCliffordNumber, y::AbstractCliffordNumber)

Calculates the regressive product of `x` and `y`. This is accomplished by taking the wedge product
of the left complements of `x` and `y`, then taking the right complement of the result.
"""
function ∨(x::AbstractCliffordNumber, y::AbstractCliffordNumber)
    return right_complement(left_complement(x) ∧ left_complement(y))
end

@doc """
    *(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})
    (x::AbstractCliffordNumber{Q})(y::AbstractCliffordNumber{Q})

Calculates the geometric product of `x` and `y`, returning the smallest type which is able to
represent all nonzero basis blades of the result.
"""
*

@doc """
    ∧(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})
    wedge(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the wedge (outer) product of two Clifford numbers `x` and `y` with quadratic form `Q`.

Note that the wedge product, in general, is *not* equal to the commutator product (or antisymmetric
product), which may be invoked with the `commutator` function or the `×` operator.
"""
∧

@doc """
    left_contraction(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})
    ⨼(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the left contraction of `x` and `y`.

For basis blades `A` of grade `m` and `B` of grade `n`, the left contraction is zero if `n < m`,
otherwise it is `KVector{n-m,Q}(A*B)`.
"""
⨼

@doc """
    right_contraction(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})
    ⨽(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the right contraction of `x` and `y`.

For basis blades `A` of grade `m` and `B` of grade `n`, the right contraction is zero if `m < n`,
otherwise it is `KVector{m-n,Q}(A*B)`.
"""
⨽

@doc """
    CliffordNumbers.dot(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the dot product of `x` and `y`.

For basis blades `A` of grade `m` and `B` of grade `n`, the dot product is equal to the left
contraction when `m >= n` and is equal to the right contraction (up to sign) when `n >= m`.

# Why is this function not exported?

The LinearAlgebra package also defines a `dot` function, and if both packages are used together,
this will cause a name conflict if `CliffordNumbers.dot` is exported. In the future, we will try to
resolve this without requiring a LinearAlgebra dependency.

Additionally, there is reason to prefer the use of the left and right contractions over the dot
product because the contractions require fewer exceptions in their definitions and properties.
"""
dot

"""
    CliffordNumbers.hestenes_dot(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Returns the Hestenes product: this is equal to the dot product given by `dot(x, y)` but is equal to
to zero when either `x` or `y` is a scalar.

# Why is this function not exported?

In almost every case, left and right contractions are preferable - the dot product and the Hestenes
product are less regular in algebraic sense, and the conditionals present in its implementation 
slow it down relative to contractions. It is provided for the sake of exact reproducibility of
results which use it.
"""
function hestenes_dot(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    return !(isscalar(x) || isscalar(y)) * dot(x, y)
end

# Long names for operations
const wedge = ∧
const regressive = ∨
const left_contraction = ⨼
const right_contraction = ⨽

#---Commutator and anticommutator products---------------------------------------------------------#
"""
    ×(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})
    commutator(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the commutator (or antisymmetric) product, equal to `1//2 * (x*y - y*x)`.

Note that the commutator product, in general, is *not* equal to the wedge product, which may be
invoked with the `wedge` function or the `∧` operator.

# Type promotion

Because of the rational `1//2` factor in the product, inputs with scalar types subtyping `Integer`
will be promoted to `Rational` subtypes.
"""
function ×(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    return 1//2 * (x*y - y*x)
end

"""
    ⨰(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})
    anticommutator(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the anticommutator (or symmetric) product, equal to `1//2 * (x*y + y*x)`.

Note that the dot product, in general, is *not* equal to the anticommutator product, which may be
invoked with `dot`. In some cases, the preferred operators might be the left and right contractions,
which use infix operators `⨼` and `⨽` respectively.

# Type promotion

Because of the rational `1//2` factor in the product, inputs with scalar types subtyping `Integer`
will be promoted to `Rational` subtypes.
"""
function ⨰(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    return 1//2 * (x*y + y*x)
end

# Long names for operations
const commutator = ×
const anticommutator = ⨰
