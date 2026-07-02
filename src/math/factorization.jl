#---Blade and versor structure---------------------------------------------------------------------#
# Structure predicates and factorizations: whether a multivector is a blade (an outer product of
# 1-vectors) or a versor (a geometric product of invertible 1-vectors), the factorizations
# themselves, and the invariant decomposition of a bivector into commuting simple planes.
#
# Blade structure is a property of the exterior algebra alone, so the blade test and the blade
# factorization reinterpret the coefficients in the positive-definite algebra of the same
# dimension, where every nonzero blade is invertible; the results transfer back unchanged because
# the wedge product never consults the metric. Versor structure does depend on the metric and is
# probed through the twisted conjugation v ↦ x v̂ x⁻¹. (Lounesto, "Clifford Algebras and Spinors",
# ch. 7; Lundholm & Svensson, "Clifford algebra, geometric algebra, and applications", §6; Dorst,
# Fontijne & Mann, "Geometric Algebra for Computer Science", ch. 21; Roelfs & De Keninck, "Graded
# Symmetry Groups: Plane and Simple".)

"""
    CliffordNumbers._structure_rtol(::Type{<:BaseNumber})

Returns the relative zero threshold for the structure tests in this file: `sqrt(eps)` for
floating-point carriers, exact zero (`false`) for exact carriers such as `Rational` or `Integer`.
"""
_structure_rtol(::Type{T}) where T<:AbstractFloat = sqrt(eps(T))
_structure_rtol(::Type{Complex{T}}) where T = _structure_rtol(T)
_structure_rtol(::Type{<:BaseNumber}) = false

"""
    CliffordNumbers._euclidean(x::KVector{K,Q,T,L}) -> KVector{K,VGA(dimension(Q)),T,L}

Reinterprets the coefficients of `x` in the positive-definite algebra of the same dimension. Blade
storage order depends only on the index bit patterns, never on the metric, so this is a plain data
copy. Exterior-algebraic (metric-free) questions about `x` can be answered on the result, where
every nonzero blade is invertible.
"""
_euclidean(x::KVector{K,Q,T,L}) where {K,Q,T,L} = KVector{K,VGA(dimension(Q)),T,L}(Tuple(x))

"""
    CliffordNumbers._parity_energy(x::AbstractCliffordNumber) -> Tuple

Returns the energies (sums of `abs2`) of the stored coefficients of even and odd grade, as a
2-tuple, with the parity mask resolved at compile time. A parity-homogeneous multivector has one of
the two energies equal to zero.
"""
@generated function _parity_energy(x::AbstractCliffordNumber)
    evens = map(b -> iseven(grade(b)), Tuple(BladeIndices(x)))
    return quote
        t = map(abs2, Tuple(x))
        return (mapreduce(*, +, $evens, t), mapreduce(*, +, $(map(!, evens)), t))
    end
end

#---Blade predicate--------------------------------------------------------------------------------#
"""
    CliffordNumbers._plucker_ok(A::KVector{K,E,T}) -> Bool

Checks the Plücker (decomposability) relations of `A` in the Euclidean reinterpretation: `A` is a
blade iff `(β ⨼ A) ∧ A` vanishes for every basis `(K-1)`-blade `β`. Contraction with a basis blade
of a positive-definite metric is coefficient extraction, so the relations are the metric-free ones.
Fully unrolled over `β`.
"""
@generated function _plucker_ok(A::KVector{K,E,T}) where {K,E,T}
    L = binomial(dimension(E), K - 1)
    checks = map(1:L) do i
        t = ntuple(isequal(i), L)
        return :(sum(abs2, Tuple((KVector{K - 1,E,T}($t) ⨼ A) ∧ A)) <= tol || return false)
    end
    return quote
        tol = (_structure_rtol(T) * sum(abs2, Tuple(A)))^2
        $(checks...)
        return true
    end
end

"""
    isblade(x::AbstractCliffordNumber) -> Bool

Returns `true` if `x` is a blade: a homogeneous multivector that factors as a wedge product of
`grade(x)` 1-vectors. Scalars, vectors, pseudovectors, pseudoscalars, and zero are always blades;
between those grades the Plücker relations `(β ⨼ x) ∧ x == 0` over the basis `(k-1)`-blades `β`
decide decomposability. Blade structure belongs to the exterior algebra, so the answer is
metric-independent and valid in degenerate and indefinite signatures. Floating-point carriers are
tested to a relative tolerance of `sqrt(eps)`.

See also: [`factor_blade`](@ref), [`isversor`](@ref), [`bivector_decomposition`](@ref).
"""
function isblade(x::KVector{K,Q}) where {K,Q}
    (K <= 1 || K >= dimension(Q) - 1) && return true
    return _plucker_ok(_euclidean(x))
end

function isblade(x::AbstractCliffordNumber{Q,T}) where {Q,T}
    tol = _structure_rtol(T)^2 * sum(abs2, Tuple(x))
    k = -1
    for g in nonzero_grades(typeof(x))
        sum(abs2, Tuple(KVector{g,Q,T}(x))) <= tol && continue
        k >= 0 && return false      # more than one nonzero grade
        k = g
    end
    k <= 1 && return true           # zero, scalar, or vector content
    return isblade(KVector{k,Q,T}(x))
end

#---Versor predicate-------------------------------------------------------------------------------#
"""
    isversor(x::AbstractCliffordNumber) -> Bool

Returns `true` if `x` is a versor: a geometric product of invertible 1-vectors, an element of the
Lipschitz group covering the orthogonal group through the twisted conjugation `v ↦ x v̂ x⁻¹`.
Checks that `x` has homogeneous grade parity, that `x x̃` is a nonzero scalar, and that conjugation
by `x` maps every basis 1-vector to a 1-vector. Nonzero scalars are versors (empty products up to
scale); null vectors are not (they are not invertible). Floating-point carriers are tested to a
relative tolerance of `sqrt(eps)`.

See also: [`factor_versor`](@ref), [`isblade`](@ref).
"""
function isversor(x::AbstractCliffordNumber{Q,T}) where {Q,T}
    scale = sum(abs2, Tuple(x))
    iszero(scale) && return false
    κ = _structure_rtol(T)
    min(_parity_energy(x)...) <= κ^2 * scale || return false
    r = x * x'
    s = scalar(r)
    sum(abs2, Tuple(r)) - s^2 <= (κ * scale)^2 || return false
    abs2(s) > (κ * scale)^2 || return false
    xgi = grade_involution(x)
    for i in 1:dimension(Q)
        v = KVector{1,Q,T}(ntuple(isequal(i), Val(Int(dimension(Q)))))
        y = xgi * v * x'                # ∝ x v̂ x⁻¹, scaled by the scalar x x̃
        sum(abs2, Tuple(y)) - sum(abs2, Tuple(KVector{1,Q}(y))) <= (κ * scale)^2 || return false
    end
    return true
end

#---Blade factorization----------------------------------------------------------------------------#
"""
    factor_blade(x::KVector{K,Q}) -> NTuple{K,KVector{1,Q}}

Factors the blade `x` into `K` 1-vectors whose wedge product reproduces `x`. The factors are
projections of the basis vectors of the largest basis blade of `x` onto `x`, orthogonalized in the
positive-definite coefficient metric — not in `Q`, where orthogonality is not generally attainable
— so the reassembly holds in every signature, including degenerate ones. The blade's weight is
carried by the first factor. Exact on exact (e.g. `Rational`) carriers; carriers that are not
closed under division are promoted. Throws a `DomainError` if `x` is zero, a scalar, or not a
blade.

See also: [`isblade`](@ref), [`factor_versor`](@ref).
"""
function factor_blade(x::KVector{K,Q,T}) where {K,Q,T}
    iszero(K) && throw(DomainError(x, "scalars have no 1-vector factors."))
    (iszero(x) || !isblade(x)) && throw(DomainError(x, "not a nonzero blade."))
    K == 1 && return (x,)
    A = _euclidean(x)
    Ainv = versor_inverse(A)
    E = signature(A)
    n = Int(dimension(E))
    p = argmax(map(abs2, Tuple(A)))
    bits = Int(BladeIndices(typeof(A))[p])
    # Project the basis vectors of the dominant basis blade onto the blade, then orthogonalize
    # (plain Gram–Schmidt, no normalization, so exact carriers stay exact).
    gs = map(filter(i -> !iszero(bits >> (i - 1) & 1), 1:n)) do i
        e = KVector{1,E,T}(ntuple(isequal(i), Val(n)))
        return (e ⨼ A) ⨼ Ainv
    end
    for i in 2:K, j in 1:(i - 1)
        d = scalar_product(gs[j], gs[j])
        iszero(d) && throw(DomainError(x, "blade factorization is degenerate."))
        gs[i] = gs[i] - (scalar_product(gs[i], gs[j]) / d) * gs[j]
    end
    # The orthogonalized wedge is proportional to the blade; rescale the first factor to match.
    w = reduce(∧, gs)
    r = Tuple(A)[p] / Tuple(w)[p]
    S = promote_type(scalar_type(eltype(gs)), typeof(r))
    return ntuple(i -> KVector{1,Q,S}(Tuple(isone(i) ? r * gs[1] : gs[i])), Val(K))
end

#---Versor factorization---------------------------------------------------------------------------#
"""
    CliffordNumbers._probe_vectors(::Type{KVector{1,Q,S}}) -> Vector{KVector{1,Q,S}}

Returns the probe vectors for the reflection peel of [`factor_versor`](@ref): the non-null basis
vectors, then the non-null pairwise sums and differences (needed when every basis vector
individually misses the moving subspace).
"""
function _probe_vectors(::Type{KVector{1,Q,S}}) where {Q,S}
    n = Int(dimension(Q))
    q = metric_tuple(Q)
    probes = KVector{1,Q,S}[]
    for i in 1:n
        iszero(q[i]) || push!(probes, KVector{1,Q,S}(ntuple(isequal(i), Val(n))))
    end
    for i in 1:n, j in (i + 1):n
        iszero(q[i] + q[j]) && continue
        push!(probes, KVector{1,Q,S}(ntuple(k -> S(k == i || k == j), Val(n))))
        push!(probes, KVector{1,Q,S}(ntuple(k -> S(k == i) - S(k == j), Val(n))))
    end
    return probes
end

"""
    factor_versor(x::AbstractCliffordNumber{Q}) -> Tuple{Vararg{KVector{1,Q}}}

Factors the versor `x` into invertible 1-vectors whose geometric product reproduces `x`
(Cartan–Dieudonné: an orthogonal transformation of an `n`-dimensional space is a product of at
most `n` reflections). Reflections are peeled off the twisted conjugation `v ↦ x v̂ x⁻¹`: a probe
vector `a` moved to `f(a)` is returned by the single reflection along `f(a) - a` when that
difference is non-null; when every displacement is null the residual is `σ(1 + N)` with `N` a
simple nilpotent bivector (a translation or null rotation) whose mirror pair is read off in closed
form. The residual scalar magnitude is folded into the last factor; a scalar `x` factors as a
parallel pair. Exact on exact carriers; carriers not closed under division are promoted. Throws a
`DomainError` if `x` is not a versor or the peel stalls.

See also: [`isversor`](@ref), [`factor_blade`](@ref).
"""
function factor_versor(x::AbstractCliffordNumber{Q,T}) where {Q,T}
    isversor(x) || throw(DomainError(x, "not a versor."))
    S = typeof(inv(one(T)))
    n = Int(dimension(Q))
    V = CliffordNumber{Q,S}(x)
    κ = _structure_rtol(S)
    mirrors = KVector{1,Q,S}[]
    probes = _probe_vectors(KVector{1,Q,S})
    for _ in 0:(n + 1)
        scale = sum(abs2, Tuple(V))
        s0 = scalar(V)
        if scale - s0^2 <= κ^2 * scale
            if isempty(mirrors)
                i = findfirst(!iszero, metric_tuple(Q))
                isnothing(i) && throw(DomainError(x, "the algebra has no invertible 1-vectors."))
                e = KVector{1,Q,S}(ntuple(isequal(i), Val(n)))
                return (e, s0 * versor_inverse(e))
            end
            mirrors[end] = s0 * mirrors[end]
            return Tuple(mirrors)
        end
        Vinv = versor_inverse(V)
        Vgi = grade_involution(V)
        # Pass 1: the most non-null single mirror f(a) - a over the probes.
        bestrel = zero(real(S))
        bestd = zero(eltype(mirrors))
        for a in probes
            d = KVector{1,Q}(Vgi * a * Vinv) - a
            nd = sum(abs2, Tuple(d))
            nd <= κ^2 * sum(abs2, Tuple(a)) && continue
            rel = abs(scalar_product(d, d)) / nd
            rel > bestrel && ((bestrel, bestd) = (rel, d))
        end
        if bestrel > κ
            push!(mirrors, bestd)
            V = CliffordNumber{Q,S}(versor_inverse(bestd) * V)
            continue
        end
        # Pass 2: every displacement is null, so the residual is σ(1 + N) with N a simple
        # nilpotent bivector — a translation or null rotation, the two-mirror case of
        # Cartan–Dieudonné. With m a non-null vector in the plane of N and n = -(m ⨼ N)/m² its
        # orthogonal partner (N = n ∧ m, n·m = 0), the pair is exact: 1 + n m = (m + m² n)(m / m²).
        abs2(s0) > κ^2 * scale || throw(DomainError(x, "the reflection peel stalled."))
        N = KVector{2,Q}(V) / s0
        isblade(N) || throw(DomainError(x, "the reflection peel stalled."))
        (g1, g2) = factor_blade(N)
        bestrel = zero(real(S))
        m = zero(eltype(mirrors))
        for c in (g1, g2, g1 + g2, g1 - g2)
            nc = sum(abs2, Tuple(c))
            iszero(nc) && continue
            rel = abs(scalar_product(c, c)) / nc
            rel > bestrel && ((bestrel, m) = (rel, c))
        end
        bestrel > κ || throw(DomainError(x, "the reflection peel stalled."))
        nv = -(m ⨼ N) / scalar_product(m, m)
        v2 = m + scalar_product(m, m) * nv
        v1 = m / scalar_product(m, m)
        push!(mirrors, KVector{1,Q,S}(Tuple(v2)))
        push!(mirrors, KVector{1,Q,S}(Tuple(v1)))
        V = CliffordNumber{Q,S}(versor_inverse(v1) * (versor_inverse(v2) * V))
    end
    throw(DomainError(x, "the reflection peel did not terminate."))
end

#---Invariant decomposition of bivectors-----------------------------------------------------------#
"""
    CliffordNumbers._isoclinic_split(b::KVector{2,Q,T}) -> NTuple{2,KVector{2,Q,T}}

Splits an isoclinic bivector (equal plane squares, where the eigen-split of
[`_bivector_split`](@ref CliffordNumbers._bivector_split) degenerates) into commuting simple
planes. An isoclinic bivector has infinitely many invariant plane pairs — every vector spans an
invariant plane with its own image `a ⨼ b` — so the basis vectors are probed and `b` is projected
onto the best-conditioned (least null) such plane; the complement is the commuting partner.
"""
function _isoclinic_split(b::KVector{2,Q,T}) where {Q,T}
    n = Int(dimension(Q))
    bestrel = zero(T)
    bestU = zero(b)
    bestm = one(T)
    for i in 1:n
        a = KVector{1,Q,T}(ntuple(isequal(i), Val(n)))
        U = a ∧ (a ⨼ b)
        w = sum(abs2, Tuple(U))
        iszero(w) && continue
        m = scalar_product(U, U')
        rel = abs(m) / w
        (rel > bestrel) && ((bestrel, bestU, bestm) = (rel, U, m))
    end
    bestrel <= sqrt(eps(T)) && return (b, zero(b))    # no usable plane: nothing to split off
    b1 = bestU * (scalar_product(b, bestU') / bestm)
    return (b1, b - b1)
end

"""
    CliffordNumbers._bivector_split(b::KVector{2,Q,T}) -> NTuple{2,KVector{2,Q,T}}

Splits a float-carrier bivector into commuting simple bivectors `(b1, b2)` with `b1 + b2 == b` and
`b2 == 0` when `b` is simple: the generated kernel behind [`bivector_decomposition`](@ref). Below
six dimensions a bivector has at most two invariant planes, whose squares `λ` are the roots of
`λ² − P λ + g/4` with `P = ⟨b²⟩₀` and `g = ⟨(⟨b²⟩₄)²⟩₀`; distinct roots give the closed-form
split, repeated roots the isoclinic probe. Shared by the rotor logarithm.
"""
@generated function _bivector_split(b::KVector{2,Q,T}) where {Q,T}
    # In fewer than four dimensions every bivector is simple (there is no grade-4 part).
    dimension(Q) < 4 && return :((b, zero(b)))
    dimension(Q) >= 6 && return :(throw(DomainError(b,
        "the two-plane invariant decomposition is only defined below six dimensions.")))
    return quote
        bb = b * b
        G4 = _grade4(bb)
        # Simplicity is decided on the coefficients of ⟨b²⟩₄, not its metric norm, which vanishes
        # spuriously in degenerate signatures.
        scale = sum(abs2, Tuple(b))
        sum(abs2, Tuple(G4)) <= eps(T) * scale^2 && return (b, zero(b))
        P = scalar(bb)
        g = scalar_product(G4, G4)
        disc = P * P - g
        sq = sqrt(max(disc, zero(disc)))
        # Numerically stable quadratic roots (λ1 λ2 = g/4): take the larger-magnitude root
        # directly and the smaller from the product to avoid cancellation when one plane nearly
        # vanishes.
        λ1 = (P + copysign(sq, P)) / 2
        λ2 = iszero(λ1) ? (P - sq) / 2 : (g / 4) / λ1
        abs(λ1 - λ2) <= sqrt(eps(T)) * (abs(P) + one(T)) && return _isoclinic_split(b)
        bg = _grade2(b * G4)
        b1 = _grade2((λ1 * b - bg / 2) * inv(λ1 - λ2))
        b2 = _grade2((λ2 * b - bg / 2) * inv(λ2 - λ1))
        return (b1, b2)
    end
end

"""
    bivector_decomposition(b::KVector{2,Q}) -> NTuple{2,KVector{2,Q}}

Splits a bivector into two commuting simple bivectors `(b1, b2)` with `b1 + b2 ≈ b` — the
invariant decomposition of Roelfs & De Keninck ("Graded Symmetry Groups: Plane and Simple"). `b2`
is zero when `b` is simple. The two plane squares are the roots of
`λ² − ⟨b²⟩₀ λ + ⟨(⟨b²⟩₄)²⟩₀ / 4`; repeated roots (isoclinic bivectors, where the split is not
unique) are resolved by probing an invariant plane through a basis vector. Only defined below six
dimensions, where at most two invariant planes exist. The carrier is converted to floating point.
The rotor logarithm and the analytic functions of versors are built on this split.

See also: [`log(::EvenCliffordNumber)`](@ref), [`isblade`](@ref).
"""
bivector_decomposition(b::KVector{2}) = _bivector_split(float(b))
