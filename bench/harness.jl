# Shared scaffolding for the benchmark suite. Every workload is measured against
# the Julia primitive a user would otherwise reach for (`Complex{T}` or
# `Quaternion{T}`) and reports a timing ratio plus an algebra-drift check, so a
# mismatch flags a regression even when the timing budget is too small to be
# useful.

using BenchmarkTools
using CliffordNumbers
using Quaternions
using Printf
using Random

# Per-workload time budget. 0.5 s is a CI smoke value; 5 s gives real numbers.
const BENCH_SECONDS = parse(Float64, get(ENV, "BENCH_SECONDS", "0.5"))

#---Algebra ↔ primitive isomorphisms---------------------------------------------------------------#

# The even subalgebra of VGA(2) is isomorphic to ℂ, coordinate-for-coordinate:
# the scalar maps to the real part and the e₁₂ coefficient to the imaginary part.
const Spinor{T} = EvenCliffordNumber{VGA(2),T,2}

tospinor(z::Complex{T}) where T = Spinor{T}((real(z), imag(z)))
ascomplex(s::Spinor) = (t = Tuple(s); complex(t[1], t[2]))

# The even subalgebra of VGA(3) is isomorphic to ℍ, but `Quaternion(::Even...)`
# is not a ring homomorphism in raw coordinates: the two trailing bivector slots
# match the quaternion j/k axes in swapped order. `as_hamilton` applies that swap
# so `as_hamilton(a*b) == as_hamilton(a) * as_hamilton(b)`, making ℍ
# multiplication a valid drift oracle for the compose workloads.
const Rotor{T} = EvenCliffordNumber{VGA(3),T,4}
const Vec3{T} = KVector{1,VGA(3),T,3}

as_hamilton(r::Rotor) = (t = Tuple(r); Quaternion(t[1], t[2], t[4], t[3]))

"""
    randrotor(rng) -> Rotor{Float64}

A unit rotor `exp(B)` for a random bivector `B`; `abs2` is exactly 1.
"""
randrotor(rng::AbstractRNG) = exp(KVector{2,VGA(3),Float64,3}(ntuple(_ -> randn(rng), 3)))

"""
    ga_rotate(r::Rotor, v::Vec3)

Sandwich a grade-1 vector with a rotor, `r v r̃`. The result is an
`OddCliffordNumber`; its grade-1 part is the rotated vector.
"""
ga_rotate(r::Rotor, v::Vec3) = r * v * reverse(r)

"""
    q_rotate(q::Quaternion, v::NTuple{3})

Rotate a 3-tuple by conjugating the corresponding pure-imaginary quaternion. Used
only as a timing baseline for `ga_rotate`.
"""
function q_rotate(q::Quaternion, v::NTuple{3})
    o = q * Quaternion(zero(v[1]), v[1], v[2], v[3]) * conj(q)
    return (o.v1, o.v2, o.v3)
end

#---Measurement and reporting----------------------------------------------------------------------#

"""Element-wise `isapprox` over two equal-length collections."""
approxall(a, b; kw...) = length(a) == length(b) && all(((x, y),) -> isapprox(x, y; kw...), zip(a, b))

# The multivector workloads have no stdlib primitive, so they time a compact
# representation against the same computation on the dense `CliffordNumber`. The
# drift check is convention-free coefficient equality between the two.
dense(x::AbstractCliffordNumber{Q}) where Q = CliffordNumber{Q}(x)
dense(x::Number) = x   # scalar-valued ops (e.g. `scalar_product`) compare directly

"""Compact-vs-dense agreement: densify each compact result and compare to the dense result."""
densematch(cs, rs; kw...) = approxall(map(dense, cs), rs; kw...)

struct WResult
    name::String
    cn_ns::Float64
    ref_ns::Float64
    ok::Bool
end

_min_ns(b::BenchmarkTools.Trial) = minimum(b).time

function _measure(cn_call, ref_call, check, divisor::Real)
    ok = check(cn_call(), ref_call())::Bool
    cn = _min_ns(@benchmark $cn_call() samples = 10_000 evals = 1 seconds = BENCH_SECONDS)
    rf = _min_ns(@benchmark $ref_call() samples = 10_000 evals = 1 seconds = BENCH_SECONDS)
    return (cn / divisor, rf / divisor, ok)
end

"""
    run_workload(name, cn_call, ref_call, check) -> WResult

`cn_call` and `ref_call` are zero-argument closures returning the Clifford and
primitive results respectively. `check(cn_result, ref_result) -> Bool` is the
algebra-drift comparison, evaluated once before timing.
"""
function run_workload(name::AbstractString, cn_call, ref_call, check)
    cn, rf, ok = _measure(cn_call, ref_call, check, 1)
    return WResult(name, cn, rf, ok)
end

"""
    run_per_op(name, n, cn_call, ref_call, check) -> WResult

Per-element microbenchmark: `cn_call`/`ref_call` each process `n` items (e.g.
`map(inv, xs)`), and the reported times are divided by `n`. Batching defeats the
constant-folding that erases a literal single call, while amortizing the closure
call overhead to nothing — so the CN/ref ratio is the per-op acceptance metric.
"""
function run_per_op(name::AbstractString, n::Integer, cn_call, ref_call, check)
    cn, rf, ok = _measure(cn_call, ref_call, check, n)
    return WResult(name, cn, rf, ok)
end

function report(title::AbstractString, results::Vector{WResult})
    println()
    println(title, "    [BENCH_SECONDS = ", BENCH_SECONDS, "]")
    @printf("  %-32s %12s %12s %9s   %s\n", "workload", "CN", "ref", "CN/ref", "check")
    println("  ", repeat('-', 76))
    for r in results
        @printf("  %-32s %12s %12s %9.2f   %s\n",
            r.name,
            BenchmarkTools.prettytime(r.cn_ns),
            BenchmarkTools.prettytime(r.ref_ns),
            r.cn_ns / r.ref_ns,
            r.ok ? "ok" : "FAIL")
    end
    nfail = count(r -> !r.ok, results)
    nfail == 0 || @warn "$nfail workload(s) failed the algebra-drift check in '$title'"
    return nfail
end
