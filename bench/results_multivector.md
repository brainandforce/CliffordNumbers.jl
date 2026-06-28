# Benchmark results — `bench_multivector` (compact vs dense)

Recorded run of the multivector groups of the [`bench/`](README.md) suite on the
consolidated `clifford-improvements` branch: `bench_kvector`, `bench_even`,
`bench_odd`, `bench_general`, and `bench_simd`, covering `KVector`,
`EvenCliffordNumber`, `OddCliffordNumber`, and dense `CliffordNumber` over
PGA / CGA / VGA. The operations and grade conventions mirror how the types are
used in applications (ComputationalGeometricAlgebra.jl): `meet` (`∧`), `join`
(`∨`), `project`, `move` / `rotate` / `reflect` (the sandwich `M X M̃`),
`exp`-generated rotors and motors, contractions, complements, automorphisms and
inverses.

There is no stdlib oracle for the general types, so `ref` is the **identical
computation on the dense `CliffordNumber`** of the same algebra: `CN/ref < 1`
means the compact type wins (it never computes the zero coefficients dense
storage carries); **`CN/ref > 1` flags a compact path that is slower than
brute-forcing the dense multivector** — a concrete optimization target. Every row
also passed its drift check (compact result ≈ dense result, embedded
coordinate-for-coordinate via `densematch`), so the numbers double as a
correctness snapshot.

## Environment

| | |
|---|---|
| Commit | `2d99643` (branch `clifford-improvements`) |
| Date | 2026-06-28 |
| Julia | 1.12.6 |
| CPU | Intel Core i7-10700KF @ 3.80 GHz, 16 threads |
| OS | Windows 11 (`x86_64-w64-mingw32`) |
| `BENCH_SECONDS` | 5.0 |
| Invocation | `BENCH_SECONDS=5 julia -t auto bench/runbenchmarks.jl` |

Per-op rows (`bench_kvector` / `bench_even` / `bench_odd` / `bench_general`) are
per-element times from a batched microbench of 1000 items, so the op cannot be
constant-folded. The `bench_simd` rows are wall-clock for the whole broadcast.

## bench_kvector — geometric primitives & incidence (compact vs dense)

`KVector` of the application grade (PGA(3): point=3, line=2, plane=1; PGA(2):
point=2, line=1) vs the dense `CliffordNumber` of the same algebra.

| workload | CN | ref | CN/ref | check |
|---|---:|---:|---:|:---:|
| meet: plane ∧ plane → line (PGA3) | 10.7 ns | 34.9 ns | **0.31** | ok |
| meet: line ∧ plane → point (PGA3) | 15.0 ns | 34.5 ns | **0.43** | ok |
| meet: 3 planes → point (PGA3) | 10.5 ns | 80.0 ns | **0.13** | ok |
| join: point ∨ point → line (PGA3) | 10.2 ns | 52.5 ns | **0.19** | ok |
| join: line ∨ point → plane (PGA3) | 8.4 ns | 52.3 ns | **0.16** | ok |
| project: (π⌋P)∧π → point (PGA3) | 8.9 ns | 47.3 ns | **0.19** | ok |
| meet: line ∧ line → point (PGA2) | 2.0 ns | 13.2 ns | **0.15** | ok |
| join: point ∨ point → line (PGA2) | 2.1 ns | 18.9 ns | **0.11** | ok |
| wedge: v ∧ w → bivector (VGA3) | 2.3 ns | 15.4 ns | **0.15** | ok |
| scalar_product ⟨v,w⟩ (VGA3) | 0.9 ns | 2.3 ns | **0.39** | ok |
| left_contraction v⌋B (VGA3) | 1.8 ns | 16.7 ns | **0.11** | ok |
| right_complement bivector (VGA3) | 1.5 ns | 5.3 ns | **0.28** | ok |
| commutator B₁×B₂ (VGA3) | 9.4 ns | 29.2 ns | **0.32** | ok |
| reverse line bivector (PGA3) | 2.8 ns | 7.8 ns | **0.36** | ok |

## bench_even — rotors & motors (EvenCliffordNumber vs dense)

| workload | CN | ref | CN/ref | check |
|---|---:|---:|---:|:---:|
| rotor compose r₁r₂ (VGA3) | 11.0 ns | 11.5 ns | **0.96** | ok |
| rotor sandwich r·v·r̃ (VGA3) | 9.1 ns | 21.1 ns | **0.43** | ok |
| rotor from bivector exp(B) (VGA3) | 44.7 ns | 187.8 ns | **0.24** | ok |
| relative rotor r₁/r₂ (VGA3) | 19.0 ns | 70.0 ns | **0.27** | ok |
| motor compose m₁m₂ (PGA3) | 14.0 ns | 39.5 ns | **0.35** | ok |
| move point M·P·M̃ (PGA3) | 19.3 ns | 84.7 ns | **0.23** | ok |
| move line M·L·M̃ (PGA3) | 20.5 ns | 84.4 ns | **0.24** | ok |
| move plane M·π·M̃ (PGA3) | 19.5 ns | 84.3 ns | **0.23** | ok |
| motor from bivector exp(B) (PGA3) | 327.0 ns | 564.8 ns | **0.58** | ok |
| CGA move point M·P·M̃ (CGA3) | 51.0 ns | 62.257 µs | **≈0.001** | ok |
| rotor compose (VGA4) | 12.0 ns | 31.9 ns | **0.38** | ok |

The CGA row reports `0.00` in the raw table because the dense
`CliffordNumber{CGA(3)}` carries 32 coefficients: its geometric product is
~1200× slower than the compact even motor (16 coefficients) — the clearest
illustration of why applications must not fall back to the dense type in
higher-dimensional algebras.

## bench_odd — reflections (OddCliffordNumber vs dense)

| workload | CN | ref | CN/ref | check |
|---|---:|---:|---:|:---:|
| reflect point π·P·π (PGA3) | 10.5 ns | 80.9 ns | **0.13** | ok |
| reflect plane π·π′·π (PGA3) | 17.4 ns | 81.1 ns | **0.21** | ok |
| reflect vector n·v·n (VGA3) | 6.7 ns | 17.5 ns | **0.38** | ok |
| odd versor sandwich u·v·u⁻¹ (VGA3) | 9.2 ns | 24.6 ns | **0.37** | ok |
| odd∘odd → rotor u₂u₁ (VGA3) | 10.4 ns | 11.3 ns | **0.92** | ok |
| rotor∘reflection r·u → odd (VGA3) | 10.3 ns | 11.3 ns | **0.91** | ok |
| grade_involution odd versor (VGA3) | 1.8 ns | 3.9 ns | **0.46** | ok |

## bench_general — full multivectors & inverses (CliffordNumber)

Here the dense `CliffordNumber` *is* the natural type, so `ref` is a second,
algebraically-equivalent path (block decomposition / the validation-free
inverse), and `CN/ref` measures structural or validation overhead rather than
representation size.

| workload | CN | ref | CN/ref | check |
|---|---:|---:|---:|:---:|
| geometric product full vs even/odd block (VGA4) | 30.2 ns | 52.7 ns | **0.57** | ok |
| inv vs versor_inverse: rotor (VGA3) | 2.3 ns | 2.3 ns | **1.00** | ok |
| inv vs versor_inverse: motor (PGA3) | 82.6 ns | 8.0 ns | **10.32** ⚠ | ok |

`geometric product full vs block` at **0.57** says the dense full product beats
summing the four even/odd sub-products — splitting a general multivector into
graded blocks is *not* a win when you genuinely need all grades.

The `motor (PGA3)` row is the one remaining hot spot — see below.

## bench_simd — batched broadcast at increasing N (compact vs dense)

One fixed transform broadcast over an array of primitives at
`N ∈ {1024, 65536, 1M}` — the vertex-buffer / point-cloud workload. A flat ratio
across N means the kernel scales cleanly.

| workload | CN | ref | CN/ref | check |
|---|---:|---:|---:|:---:|
| move point M·P·M̃ (PGA3) N=1024 | 18.9 µs | 88.1 µs | **0.21** | ok |
| move point M·P·M̃ (PGA3) N=65536 | 1.523 ms | 6.630 ms | **0.23** | ok |
| move point M·P·M̃ (PGA3) N=1M | 25.150 ms | 107.580 ms | **0.23** | ok |
| reflect point π·P·π (PGA3) N=1024 | 10.9 µs | 87.6 µs | **0.12** | ok |
| reflect point π·P·π (PGA3) N=65536 | 1.015 ms | 6.584 ms | **0.15** | ok |
| reflect point π·P·π (PGA3) N=1M | 16.911 ms | 107.730 ms | **0.16** | ok |
| rotor sandwich r·v·r̃ (VGA3) N=1024 | 8.9 µs | 21.1 µs | **0.42** | ok |
| rotor sandwich r·v·r̃ (VGA3) N=65536 | 731.6 µs | 1.702 ms | **0.43** | ok |
| rotor sandwich r·v·r̃ (VGA3) N=1M | 12.227 ms | 28.335 ms | **0.43** | ok |
| meet plane ∧ plane (PGA3) N=1024 | 3.7 µs | 24.4 µs | **0.15** | ok |
| meet plane ∧ plane (PGA3) N=65536 | 488.6 µs | 2.368 ms | **0.21** | ok |
| meet plane ∧ plane (PGA3) N=1M | 9.235 ms | 39.296 ms | **0.23** | ok |
| join point ∨ point (PGA3) N=1024 | 11.4 µs | 53.2 µs | **0.21** | ok |
| join point ∨ point (PGA3) N=65536 | 973.5 µs | 4.280 ms | **0.23** | ok |
| join point ∨ point (PGA3) N=1M | 16.497 ms | 70.246 ms | **0.23** | ok |
| rotor compose r₁r₂ (VGA3) N=1024 | 2.2 µs | 7.8 µs | **0.28** | ok |
| rotor compose r₁r₂ (VGA3) N=65536 | 338.0 µs | 869.8 µs | **0.39** | ok |
| rotor compose r₁r₂ (VGA3) N=1M | 7.021 ms | 15.659 ms | **0.45** | ok |

SIMD reading: `move point`, `reflect point`, and `rotor sandwich` hold a flat
~0.12–0.43 ratio from 1k to 1M — the sandwich kernels scale. `join` is now flat
at **0.21–0.23** across all N (it was the slow incidence path at ~0.64 before the
generated-complement fix; see below). `rotor compose` drifts 0.28 → 0.45 as the
dense product's extra arithmetic amortizes away once both sides are memory-bound.

## What the consolidated PRs fixed

Three hot spots flagged by the original enlarged run (on `main`) are now closed.
Before/after `CN/ref` (lower is better):

| symptom | before | after | fix |
|---|---:|---:|---|
| `reverse line bivector (PGA3)` | 5.10 ⚠ | **0.36** | closed-form `reverse(::KVector)` (PR-GPUEnablement) |
| `join: point ∨ point → line (PGA3)` | 0.66 | **0.19** | `@generated` complements (PR-EnlargedBenchFixes) |
| `join: line ∨ point → plane (PGA3)` | 0.67 | **0.16** | same |
| `join point ∨ point (PGA3)` N=1M | 0.64 | **0.23** | same, under broadcast |
| `inv vs versor_inverse: rotor (VGA3)` | 5.40 ⚠ | **1.00** | fast `inv` on the ℍ even subalgebra (PR-SpinorFastPaths) |

- **`reverse(::KVector)`** previously fell back to the generic
  `T(x[reverse.(BladeIndices(T))])` indexing path (a runtime `to_index` blade
  search); it now uses the closed form `ifelse(iszero(K & 2), x, -x)`, identical
  to `adjoint` for real scalars. This was also required for GPU-safety.
- **The regressive product `∨` / `join`** is
  `right_complement(left_complement(x) ∧ left_complement(y))`; the two
  complements ran the same runtime `to_index` search. Making
  `left_complement` / `right_complement` `@generated` (resolving every output
  coefficient's position and sign at compile time) lifted the whole incidence
  layer below the dense path.
- **`inv` on the rotor (VGA(3) ≅ ℍ)** now delegates to the validation-free
  closed-form `versor_inverse`, so the row collapses to **1.00×** (they are the
  same call). The `exp`-of-mixed-signature-bivector crash the original run had to
  work around is also fixed (PR-RotorPrimitives' Taylor-scaling clamp); the suite
  keeps `randeven2` for CGA only to stay on the closed conformal-versor path, not
  to avoid a throw.

## Remaining hot spot

### `inv` validation on the PGA(3) motor — `inv vs versor_inverse: motor` 10.32×

The motor row is now the only `CN/ref > 1` in the suite. It got *relatively*
larger than the original 5.32× — not a regression but the opposite: the
`@generated` `adjoint` fix dropped the **denominator** (`versor_inverse`) from
~52 ns to ~8 ns, while motor `inv` still computes `versor_inverse` and then
*validates* `x·inv ≈ 1 && inv·x ≈ 1` (~83 ns). PGA(3)'s even subalgebra is
**degenerate** (the metric has a null direction), so it is not isomorphic to
ℝ/ℂ/ℍ and the fast-`inv` specialization was deliberately scoped to exclude it —
for a degenerate metric the existence of the inverse is not guaranteed by
`abs2 ≠ 0`, so the validation is doing real work. Specializing a validation-free
motor inverse (with a degenerate-aware existence check) would close this; it is
out of scope for the shipped fast-`inv` work and is the next candidate if the
incidence-layer cost matters in practice.

### Where the compact types already win big

No fix needed — the payoff rows that justify reaching for the tight type:

- **sandwich kernels** (`move`/`reflect`/`rotate`): **0.12–0.24×**, flat from
  N=1k to N=1M — the rigid-motion and reflection workloads scale cleanly.
- **`exp`-generated rotors/motors**: **0.24–0.58×**.
- **CGA move point**: **≈0.001×** (≈1200×) — dense storage is catastrophic in
  CGA(3)'s 32-dimensional algebra.
- **contractions / wedge / project / join**: **0.11–0.28×**.

## Reproducing

```sh
BENCH_SECONDS=5 julia -t auto bench/runbenchmarks.jl
```

Absolute timings are machine-specific; the `CN/ref` ratios are the portable
signal. For the multivector groups, a ratio `> 1` is the actionable output: it
names a compact-type path that loses to brute-forcing the dense multivector.
