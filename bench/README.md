# CliffordNumbers benchmark suite

A benchmark suite keyed off the Julia primitives users reach for in geometric
algebra. Each workload pairs a Clifford type with the stdlib or `Quaternions.jl`
baseline it is isomorphic to and reports a timing ratio together with an
algebra-drift check (the Clifford result compared to the primitive result), so
the suite is also a correctness regression when the timing budget is too small
for meaningful numbers.

## Running

```sh
# CI smoke (default 0.5 s per workload)
julia bench/runbenchmarks.jl

# Real numbers, threaded
BENCH_SECONDS=5 julia -t auto bench/runbenchmarks.jl
```

The runner activates and instantiates the `bench/` environment itself; the
`[sources]` entry in `bench/Project.toml` pins it to the in-repo CliffordNumbers
rather than a registered release. The process exits non-zero if any
algebra-drift check fails.

`BENCH_SECONDS` is the per-workload time budget: 0.5 s is a quick smoke value,
5 s gives stable numbers.

## Recorded results

The recorded runs are split one file per benchmark group, each with its own
environment header, current `CN/ref` table, and analysis:

- [`results_complex.md`](results_complex.md): `bench_complex` (even VGA(2) ≅ ℂ),
  including the before/after for the closed-form even-subalgebra fast paths.
- [`results_quaternion.md`](results_quaternion.md): `bench_quaternion`
  (even VGA(3) ≅ ℍ), including the remaining `compose` gap and the dropped
  closed-form Hamilton product.
- [`results_multivector.md`](results_multivector.md): the compact-vs-dense
  groups (`bench_kvector`, `bench_even`, `bench_odd`, `bench_general`,
  `bench_simd`), including the closed hot spots (fast `reverse`, generated
  complements and fast `join`, fast rotor `inv`) and the one remaining hot spot
  (PGA(3) motor `inv` validation).

All three were recorded together from a single `BENCH_SECONDS=5` run on the
consolidated branch; every row passes its algebra-drift check.

## Layout

- `harness.jl`: `BENCH_SECONDS`, the ℂ/ℍ isomorphism helpers, the
  `dense`/`densematch` compact-vs-dense helpers, and the measure-and-report
  scaffolding.
- `bench_complex.jl`: `Spinor{T} = EvenCliffordNumber{VGA(2),T,2}` (≅ ℂ):
  pointwise `*`, `z² + c`, Horner deg-8, Mandelbrot/Julia escape grids, a
  contour integral, and single-call `*` / `inv` / `versor_inverse` microbenches.
  This is the acceptance baseline for the closed-form spinor fast paths.
- `bench_quaternion.jl`: `Rotor{T} = EvenCliffordNumber{VGA(3),T,4}` (Cl⁺(3) ≅ ℍ):
  compose, single and batched sandwich (`N ∈ {100, 100k, 1M}`), slerp, and
  chain-compose. The batched-sandwich rows are the hot spots the StructArrays
  SoA extension targets.
- `bench_multivector.jl`: `KVector` / `EvenCliffordNumber` / `OddCliffordNumber`
  / `CliffordNumber` over PGA/CGA/VGA, mirroring the operations of
  ComputationalGeometricAlgebra.jl: `meet` (`∧`), `join` (`∨`), `project`,
  `move`/`rotate` and `reflect` (the sandwich `M X M̃`), `exp`-generated rotors /
  motors, contractions, complements, automorphisms and inverses. Each row times
  the compact representation against the same computation on the dense
  `CliffordNumber` (see the note below). Groups: `bench_kvector` (primitives and
  incidence), `bench_even` (rotors / motors), `bench_odd` (reflections),
  `bench_general` (full multivectors and inverse validation), and `bench_simd`,
  the same kernels broadcast over arrays at `N ∈ {1024, 65536, 1M}`.

## A note on the compact-vs-dense baseline

The general Clifford types have no single stdlib primitive to race against the
way `Complex`/`Quaternion` serve the ℂ/ℍ benches. Instead `bench_multivector.jl`
times the compact representation (`KVector`, `EvenCliffordNumber`, or
`OddCliffordNumber`) against the same operation on the dense `CliffordNumber` of
the same algebra. `CN/ref < 1` means the compact type wins, since it never
computes the zero coefficients dense storage carries; `CN/ref > 1` flags a
compact path slower than the dense one, an optimization target. The drift check
(`densematch`) embeds the compact result into the dense type and compares
coefficients, so it is convention-free.

## A note on the ℍ isomorphism

The vector-space conversion `Quaternion(::EvenCliffordNumber)` shipped by
CliffordNumbers is not a ring homomorphism in raw coordinates: the two trailing
bivector slots line up with the quaternion j/k axes in swapped order. The
compose/chain workloads therefore go through `as_hamilton` (in `harness.jl`),
which applies that swap so `as_hamilton(a*b) == as_hamilton(a) * as_hamilton(b)`
and ℍ multiplication is a valid drift oracle. The sandwich workloads instead
check norm preservation, which is convention-free.
