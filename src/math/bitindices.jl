#---Compile-time blade data------------------------------------------------------------------------#
# The blade indices and per-slot automorphism signs of a Clifford number type are pure functions of
# its type parameters, so they resolve to literal tuples at compile time. `determinant.jl` and the
# automorphisms in `duals.jl` share these instead of re-deriving blade indices and signs per call.

"""
    CliffordNumbers.bitindices(::Type{T}) where T<:AbstractCliffordNumber
    CliffordNumbers.bitindices(x::AbstractCliffordNumber)

Returns the basis blade indices of `T` in storage order as a compile-time constant
`NTuple{nblades(T),BladeIndex}`. This is the tuple form of [`BladeIndices`](@ref), provided as a
trait that other generated functions can splice directly.
"""
bitindices(::Type{T}) where T<:AbstractCliffordNumber = Tuple(BladeIndices(T))
bitindices(x::AbstractCliffordNumber) = bitindices(typeof(x))

"""
    CliffordNumbers.blade_signs(::Type{T}, ::Val{F}) where {T<:AbstractCliffordNumber,F}

Returns, as a compile-time constant `NTuple{nblades(T),Int8}`, the per-slot sign that the sign-only
grade automorphism named by the `Symbol` `F` (`:reverse`, `:adjoint`, `:conj`, or
`:grade_involution`) applies to the stored coefficients of `T`.

Each named automorphism maps every basis blade to ±itself with no permutation, so applying it is the
branch-free `T(Tuple(x) .* blade_signs(T, Val(F)))`; see [`_sign_automorphism`](@ref).
"""
@generated function blade_signs(::Type{T}, ::Val{F}) where {T<:AbstractCliffordNumber,F}
    f = getfield(@__MODULE__, F)
    signs = map(b -> sign(f(b)), Tuple(BladeIndices(T)))
    return :($signs)
end

"""
    CliffordNumbers.grade_negate(x::T, ::Val{G}) where {T<:AbstractCliffordNumber,G}

Returns `x` with the coefficients of every grade in `G` (a `Tuple` of integer grades) negated, all
other coefficients unchanged. The per-slot signs are resolved at compile time.
"""
@generated function grade_negate(x::T, ::Val{G}) where {T<:AbstractCliffordNumber,G}
    signs = map(b -> ifelse(grade(b) in G, Int8(-1), Int8(1)), Tuple(BladeIndices(T)))
    return :(T(map(*, Tuple(x), $signs)))
end
