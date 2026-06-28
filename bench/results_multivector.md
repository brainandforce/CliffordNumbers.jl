# Benchmark results — `bench_multivector` (compact vs dense)

Recorded run of the multivector groups of the [`bench/`](README.md) suite:
`bench_kvector`, `bench_even`, `bench_odd`, `bench_general`, and `bench_simd`,
over `KVector`, `EvenCliffordNumber`, `OddCliffordNumber`, and dense
`CliffordNumber` in PGA / CGA / VGA. The operations and grade conventions mirror
application use (ComputationalGeometricAlgebra.jl): `meet` (`∧`), `join` (`∨`),
`project`, `move` / `rotate` / `reflect` (the sandwich `M X M̃`), `exp`-generated
rotors and motors, contractions, complements, automorphisms, and inverses.

There is no stdlib oracle for the general types, so `ref` is the identical
computation on the dense `CliffordNumber` of the same algebra. `CN/ref < 1` means
the compact type wins (it never computes the zero coefficients dense storage
carries); `CN/ref > 1` flags a compact path slower than the dense one, which is an
optimization target. Every row passed its drift check (compact result ≈ dense
result via `densematch`).

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

The `bench_kvector` / `bench_even` / `bench_odd` / `bench_general` rows are
per-element times from a batched microbench of 1000 items; the `bench_simd` rows
are wall-clock for the whole broadcast.

## bench_kvector — geometric primitives & incidence (compact vs dense)

`KVector` of the application grade (PGA(3): point=3, line=2, plane=1; PGA(2):
point=2, line=1) vs the dense `CliffordNumber`.

| workload | CN | ref | CN/ref | check |
|---|---:|---:|---:|:---:|
| meet: plane ∧ plane → line (PGA3) | 10.7 ns | 34.9 ns | 0.31 | ok |
| meet: line ∧ plane → point (PGA3) | 15.0 ns | 34.5 ns | 0.43 | ok |
| meet: 3 planes → point (PGA3) | 10.5 ns | 80.0 ns | 0.13 | ok |
| join: point ∨ point → line (PGA3) | 10.2 ns | 52.5 ns | 0.19 | ok |
| join: line ∨ point → plane (PGA3) | 8.4 ns | 52.3 ns | 0.16 | ok |
| project: (π⌋P)∧π → point (PGA3) | 8.9 ns | 47.3 ns | 0.19 | ok |
| meet: line ∧ line → point (PGA2) | 2.0 ns | 13.2 ns | 0.15 | ok |
| join: point ∨ point → line (PGA2) | 2.1 ns | 18.9 ns | 0.11 | ok |
| wedge: v ∧ w → bivector (VGA3) | 2.3 ns | 15.4 ns | 0.15 | ok |
| scalar_product ⟨v,w⟩ (VGA3) | 0.9 ns | 2.3 ns | 0.39 | ok |
| left_contraction v⌋B (VGA3) | 1.8 ns | 16.7 ns | 0.11 | ok |
| right_complement bivector (VGA3) | 1.5 ns | 5.3 ns | 0.28 | ok |
| commutator B₁×B₂ (VGA3) | 9.4 ns | 29.2 ns | 0.32 | ok |
| reverse line bivector (PGA3) | 2.8 ns | 7.8 ns | 0.36 | ok |

## bench_even — rotors & motors (EvenCliffordNumber vs dense)

| workload | CN | ref | CN/ref | check |
|---|---:|---:|---:|:---:|
| rotor compose r₁r₂ (VGA3) | 11.0 ns | 11.5 ns | 0.96 | ok |
| rotor sandwich r·v·r̃ (VGA3) | 9.1 ns | 21.1 ns | 0.43 | ok |
| rotor from bivector exp(B) (VGA3) | 44.7 ns | 187.8 ns | 0.24 | ok |
| relative rotor r₁/r₂ (VGA3) | 19.0 ns | 70.0 ns | 0.27 | ok |
| motor compose m₁m₂ (PGA3) | 14.0 ns | 39.5 ns | 0.35 | ok |
| move point M·P·M̃ (PGA3) | 19.3 ns | 84.7 ns | 0.23 | ok |
| move line M·L·M̃ (PGA3) | 20.5 ns | 84.4 ns | 0.24 | ok |
| move plane M·π·M̃ (PGA3) | 19.5 ns | 84.3 ns | 0.23 | ok |
| motor from bivector exp(B) (PGA3) | 327.0 ns | 564.8 ns | 0.58 | ok |
| CGA move point M·P·M̃ (CGA3) | 51.0 ns | 62.257 µs | ≈0.001 | ok |
| rotor compose (VGA4) | 12.0 ns | 31.9 ns | 0.38 | ok |

The CGA row reports `0.00` in the raw table because the dense
`CliffordNumber{CGA(3)}` carries 32 coefficients, so its geometric product is
~1200× slower than the compact 16-coefficient motor.

## bench_odd — reflections (OddCliffordNumber vs dense)

| workload | CN | ref | CN/ref | check |
|---|---:|---:|---:|:---:|
| reflect point π·P·π (PGA3) | 10.5 ns | 80.9 ns | 0.13 | ok |
| reflect plane π·π′·π (PGA3) | 17.4 ns | 81.1 ns | 0.21 | ok |
| reflect vector n·v·n (VGA3) | 6.7 ns | 17.5 ns | 0.38 | ok |
| odd versor sandwich u·v·u⁻¹ (VGA3) | 9.2 ns | 24.6 ns | 0.37 | ok |
| odd∘odd → rotor u₂u₁ (VGA3) | 10.4 ns | 11.3 ns | 0.92 | ok |
| rotor∘reflection r·u → odd (VGA3) | 10.3 ns | 11.3 ns | 0.91 | ok |
| grade_involution odd versor (VGA3) | 1.8 ns | 3.9 ns | 0.46 | ok |

## bench_general — full multivectors & inverses (CliffordNumber)

The dense `CliffordNumber` is the natural type here, so `ref` is a second
algebraically-equivalent path (block decomposition, or the validation-free
inverse), and `CN/ref` measures structural or validation overhead.

| workload | CN | ref | CN/ref | check |
|---|---:|---:|---:|:---:|
| geometric product full vs even/odd block (VGA4) | 30.2 ns | 52.7 ns | 0.57 | ok |
| inv vs versor_inverse: rotor (VGA3) | 2.3 ns | 2.3 ns | 1.00 | ok |
| inv vs versor_inverse: motor (PGA3) | 82.6 ns | 8.0 ns | 10.32 | ok |

The full product beating the summed even/odd blocks (0.57) means splitting a
general multivector into graded blocks is not a win when all grades are needed.
The motor row is the one remaining hot spot, discussed below.

## bench_simd — batched broadcast at increasing N (compact vs dense)

One transform broadcast over an array at `N ∈ {1024, 65536, 1M}` — the
vertex-buffer / point-cloud workload. A flat ratio across N means the kernel
scales cleanly.

| workload | CN | ref | CN/ref | check |
|---|---:|---:|---:|:---:|
| move point M·P·M̃ (PGA3) N=1024 | 18.9 µs | 88.1 µs | 0.21 | ok |
| move point M·P·M̃ (PGA3) N=65536 | 1.523 ms | 6.630 ms | 0.23 | ok |
| move point M·P·M̃ (PGA3) N=1M | 25.150 ms | 107.580 ms | 0.23 | ok |
| reflect point π·P·π (PGA3) N=1024 | 10.9 µs | 87.6 µs | 0.12 | ok |
| reflect point π·P·π (PGA3) N=65536 | 1.015 ms | 6.584 ms | 0.15 | ok |
| reflect point π·P·π (PGA3) N=1M | 16.911 ms | 107.730 ms | 0.16 | ok |
| rotor sandwich r·v·r̃ (VGA3) N=1024 | 8.9 µs | 21.1 µs | 0.42 | ok |
| rotor sandwich r·v·r̃ (VGA3) N=65536 | 731.6 µs | 1.702 ms | 0.43 | ok |
| rotor sandwich r·v·r̃ (VGA3) N=1M | 12.227 ms | 28.335 ms | 0.43 | ok |
| meet plane ∧ plane (PGA3) N=1024 | 3.7 µs | 24.4 µs | 0.15 | ok |
| meet plane ∧ plane (PGA3) N=65536 | 488.6 µs | 2.368 ms | 0.21 | ok |
| meet plane ∧ plane (PGA3) N=1M | 9.235 ms | 39.296 ms | 0.23 | ok |
| join point ∨ point (PGA3) N=1024 | 11.4 µs | 53.2 µs | 0.21 | ok |
| join point ∨ point (PGA3) N=65536 | 973.5 µs | 4.280 ms | 0.23 | ok |
| join point ∨ point (PGA3) N=1M | 16.497 ms | 70.246 ms | 0.23 | ok |
| rotor compose r₁r₂ (VGA3) N=1024 | 2.2 µs | 7.8 µs | 0.28 | ok |
| rotor compose r₁r₂ (VGA3) N=65536 | 338.0 µs | 869.8 µs | 0.39 | ok |
| rotor compose r₁r₂ (VGA3) N=1M | 7.021 ms | 15.659 ms | 0.45 | ok |

The sandwich kernels (`move`, `reflect`, `rotor sandwich`) hold a flat 0.12–0.43
from 1k to 1M. `join` is now flat at 0.21–0.23 (it was ~0.64 before the
generated-complement fix). `rotor compose` drifts 0.28 → 0.45 as the dense
product's extra arithmetic amortizes away once both sides are memory-bound.

## What the consolidated PRs fixed

Three hot spots from the original run on `main` are now closed:

| symptom | before | after | fix |
|---|---:|---:|---|
| reverse line bivector (PGA3) | 5.10 | 0.36 | closed-form `reverse(::KVector)` (PR-GPUEnablement) |
| join: point ∨ point (PGA3) | 0.66 | 0.19 | `@generated` complements (PR-EnlargedBenchFixes) |
| join: line ∨ point (PGA3) | 0.67 | 0.16 | same |
| inv vs versor_inverse: rotor (VGA3) | 5.40 | 1.00 | fast `inv` on the ℍ even subalgebra (PR-SpinorFastPaths) |

- `reverse(::KVector)` previously fell back to the generic
  `T(x[reverse.(BladeIndices(T))])` indexing path (a runtime `to_index` blade
  search); it now uses the closed form `ifelse(iszero(K & 2), x, -x)`.
- The regressive product `∨` is
  `right_complement(left_complement(x) ∧ left_complement(y))`; the two
  complements ran the same runtime search. Making `left_complement` /
  `right_complement` `@generated` resolves every position and sign at compile
  time.
- `inv` on the VGA(3) rotor now delegates to the validation-free closed-form
  `versor_inverse`, so the row is 1.00 (the same call). The
  `exp`-of-mixed-signature-bivector crash the original run worked around is also
  fixed (PR-RotorPrimitives' Taylor-scaling clamp); the suite keeps `randeven2`
  for CGA only to stay on the closed conformal-versor path.

## Remaining hot spot — motor `inv` validation (10.32×)

The PGA(3) motor row is the only `CN/ref > 1` left. It grew from the original
5.32× because the `@generated` `adjoint` fix dropped the denominator
(`versor_inverse`) from ~52 ns to ~8 ns, while motor `inv` still computes
`versor_inverse` and then validates `x·inv ≈ 1 && inv·x ≈ 1` (~83 ns). PGA(3)'s
even subalgebra is degenerate (the metric has a null direction), so it is not
isomorphic to ℝ/ℂ/ℍ and the fast-`inv` specialization deliberately excludes it:
for a degenerate metric the inverse is not guaranteed by `abs2 ≠ 0`, so the
validation does real work. A validation-free motor inverse with a
degenerate-aware existence check would close this; it is out of scope for the
shipped work.

## Where the compact types win

- sandwich kernels (`move` / `reflect` / `rotate`): 0.12–0.24, flat from 1k to 1M.
- `exp`-generated rotors / motors: 0.24–0.58.
- CGA move point: ≈0.001 (≈1200×); dense storage is impractical in CGA(3).
- contractions / wedge / project / join: 0.11–0.28.

## Reproducing

```sh
BENCH_SECONDS=5 julia -t auto bench/runbenchmarks.jl
```

Absolute timings are machine-specific; the `CN/ref` ratios are the portable
signal. A ratio `> 1` is the actionable output: a compact path that loses to the
dense one.
