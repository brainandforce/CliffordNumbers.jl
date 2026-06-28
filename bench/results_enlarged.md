# Benchmark results — enlarged suite

Recorded run of the full [`bench/`](README.md) suite, including the new
`bench_multivector.jl` workloads that cover `KVector`, `EvenCliffordNumber`,
`OddCliffordNumber`, and dense `CliffordNumber` over PGA / CGA / VGA. The
operations and grade conventions mirror how the types are used in applications
(ComputationalGeometricAlgebra.jl): `meet` (`∧`), `join` (`∨`), `project`,
`move` / `rotate` / `reflect` (the sandwich `M X M̃`), `exp`-generated rotors and
motors, contractions, complements, automorphisms and inverses.

For the ℂ / ℍ groups, `CN/ref` is the Clifford time over the stdlib /
`Quaternions.jl` primitive (lower is better; `<1` means CliffordNumbers is
faster). For the multivector groups there is no stdlib oracle, so `ref` is the
**identical computation on the dense `CliffordNumber`** of the same algebra:
`CN/ref < 1` means the compact type wins (it never computes the zero
coefficients dense storage carries); **`CN/ref > 1` flags a compact path that is
slower than brute-forcing the dense multivector** — a concrete optimization
target. Every row also passed its drift check (compact result ≈ dense result,
embedded coordinate-for-coordinate via `densematch`), so the numbers double as a
correctness snapshot.

## Environment

| | |
|---|---|
| Commit | `1da0d6e` (branch `pr-benchsuite`) |
| Date | 2026-06-04 |
| Julia | 1.12.5 |
| CPU | x86_64 skylake, 16 threads |
| OS | Windows 11 (`x86_64-w64-mingw32`) |
| `BENCH_SECONDS` | 5.0 |
| Invocation | `BENCH_SECONDS=5 julia -t auto bench/runbenchmarks.jl` |

Per-op rows (`(per op)`, and every row in the `bench_kvector` / `bench_even` /
`bench_odd` / `bench_general` groups) are per-element times from a batched
microbench of 1000 items, so the op cannot be constant-folded; the ratio is the
per-call acceptance metric. The `bench_simd` rows and the named array workloads
are wall-clock for the whole broadcast.

## bench_complex — even VGA(2) ≅ ℂ

`Spinor{T} = EvenCliffordNumber{VGA(2),T,2}` vs `Complex{T}`.

| workload | CN | ref | CN/ref | check |
|---|---:|---:|---:|:---:|
| pointwise `*` (N=65536) | 161.5 µs | 160.7 µs | **1.00** | ok |
| `z² + c` (N=65536) | 162.7 µs | 161.8 µs | **1.01** | ok |
| Horner deg-8 (N=65536) | 591.7 µs | 653.3 µs | **0.91** | ok |
| Mandelbrot 512² (1 thread) | 18.946 ms | 9.126 ms | **2.08** | ok |
| Mandelbrot 512² (threaded) | 2.711 ms | 1.256 ms | **2.16** | ok |
| Julia set 512² | 13.922 ms | 6.644 ms | **2.10** | ok |
| contour integral (N=4096) | 688.0 µs | 38.3 µs | **17.96** | ok |
| microbench: `*` (per op) | 1.0 ns | 0.9 ns | **1.11** | ok |
| microbench: `inv` (per op) | 56.2 ns | 2.7 ns | **20.81** | ok |
| microbench: `versor_inverse` (per op) | 9.8 ns | 2.7 ns | **3.63** | ok |

## bench_quaternion — even VGA(3) ≅ ℍ

`Rotor{T} = EvenCliffordNumber{VGA(3),T,4}` vs `Quaternion{T}`.

| workload | CN | ref | CN/ref | check |
|---|---:|---:|---:|:---:|
| compose (per op) | 11.0 ns | 3.8 ns | **2.89** | ok |
| sandwich (per op) | 13.5 ns | 18.6 ns | **0.73** | ok |
| batched sandwich (N=100) | 1.5 µs | 1.9 µs | **0.79** | ok |
| batched sandwich (N=100k) | 1.737 ms | 2.085 ms | **0.83** | ok |
| batched sandwich (N=1M) | 17.802 ms | 21.184 ms | **0.84** | ok |
| slerp (per op) | 57.3 ns | 57.9 ns | **0.99** | ok |
| chain-compose (N=10k) | 53.0 µs | 70.3 µs | **0.75** | ok |

## bench_kvector — geometric primitives & incidence (compact vs dense)

`KVector` of the application grade (PGA(3): point=3, line=2, plane=1; PGA(2):
point=2, line=1) vs the dense `CliffordNumber` of the same algebra.

| workload | CN | ref | CN/ref | check |
|---|---:|---:|---:|:---:|
| meet: plane ∧ plane → line (PGA3) | 10.7 ns | 34.9 ns | **0.31** | ok |
| meet: line ∧ plane → point (PGA3) | 14.7 ns | 34.9 ns | **0.42** | ok |
| meet: 3 planes → point (PGA3) | 10.2 ns | 79.5 ns | **0.13** | ok |
| join: point ∨ point → line (PGA3) | 156.0 ns | 237.2 ns | **0.66** | ok |
| join: line ∨ point → plane (PGA3) | 160.3 ns | 237.9 ns | **0.67** | ok |
| project: (π⌋P)∧π → point (PGA3) | 8.8 ns | 47.2 ns | **0.19** | ok |
| meet: line ∧ line → point (PGA2) | 1.9 ns | 13.3 ns | **0.14** | ok |
| join: point ∨ point → line (PGA2) | 44.4 ns | 116.0 ns | **0.38** | ok |
| wedge: v ∧ w → bivector (VGA3) | 2.0 ns | 13.9 ns | **0.14** | ok |
| scalar_product ⟨v,w⟩ (VGA3) | 1.0 ns | 2.2 ns | **0.45** | ok |
| left_contraction v⌋B (VGA3) | 1.8 ns | 15.2 ns | **0.12** | ok |
| right_complement bivector (VGA3) | 12.7 ns | 31.6 ns | **0.40** | ok |
| commutator B₁×B₂ (VGA3) | 9.4 ns | 29.2 ns | **0.32** | ok |
| reverse line bivector (PGA3) | 162.6 ns | 31.9 ns | **5.10** ⚠ | ok |

## bench_even — rotors & motors (EvenCliffordNumber vs dense)

| workload | CN | ref | CN/ref | check |
|---|---:|---:|---:|:---:|
| rotor compose r₁r₂ (VGA3) | 10.7 ns | 11.4 ns | **0.94** | ok |
| rotor sandwich r·v·r̃ (VGA3) | 13.8 ns | 32.6 ns | **0.42** | ok |
| rotor from bivector exp(B) (VGA3) | 42.5 ns | 210.9 ns | **0.20** | ok |
| relative rotor r₁/r₂ (VGA3) | 116.7 ns | 260.3 ns | **0.45** | ok |
| motor compose m₁m₂ (PGA3) | 14.3 ns | 38.8 ns | **0.37** | ok |
| move point M·P·M̃ (PGA3) | 37.6 ns | 147.2 ns | **0.26** | ok |
| move line M·L·M̃ (PGA3) | 36.4 ns | 148.4 ns | **0.25** | ok |
| move plane M·π·M̃ (PGA3) | 37.4 ns | 148.4 ns | **0.25** | ok |
| motor from bivector exp(B) (PGA3) | 328.1 ns | 592.5 ns | **0.55** | ok |
| CGA move point M·P·M̃ (CGA3) | 85.6 ns | 63.346 µs | **≈0.001** | ok |
| rotor compose (VGA4) | 12.1 ns | 30.6 ns | **0.40** | ok |

The CGA row reports `0.00` in the raw table because the dense `CliffordNumber{CGA(3)}`
carries 32 coefficients: its geometric product is ~740× slower than the compact
even motor (16 coefficients) — the clearest illustration of why applications must
not fall back to the dense type in higher-dimensional algebras.

## bench_odd — reflections (OddCliffordNumber vs dense)

| workload | CN | ref | CN/ref | check |
|---|---:|---:|---:|:---:|
| reflect point π·P·π (PGA3) | 11.0 ns | 80.9 ns | **0.14** | ok |
| reflect plane π·π′·π (PGA3) | 17.2 ns | 80.0 ns | **0.21** | ok |
| reflect vector n·v·n (VGA3) | 6.9 ns | 17.2 ns | **0.40** | ok |
| odd versor sandwich u·v·u⁻¹ (VGA3) | 28.8 ns | 59.4 ns | **0.48** | ok |
| odd∘odd → rotor u₂u₁ (VGA3) | 10.6 ns | 11.5 ns | **0.92** | ok |
| rotor∘reflection r·u → odd (VGA3) | 10.6 ns | 12.2 ns | **0.87** | ok |
| grade_involution odd versor (VGA3) | 9.6 ns | 18.2 ns | **0.53** | ok |

## bench_general — full multivectors & inverses (CliffordNumber)

Here the dense `CliffordNumber` *is* the natural type, so `ref` is a second,
algebraically-equivalent path (block decomposition / the validation-free
inverse), and `CN/ref` measures structural or validation overhead rather than
representation size.

| workload | CN | ref | CN/ref | check |
|---|---:|---:|---:|:---:|
| geometric product full vs even/odd block (VGA4) | 30.4 ns | 52.8 ns | **0.58** | ok |
| inv vs versor_inverse: rotor (VGA3) | 109.1 ns | 20.2 ns | **5.40** ⚠ | ok |
| inv vs versor_inverse: motor (PGA3) | 276.3 ns | 51.9 ns | **5.32** ⚠ | ok |

`geometric product full vs block` at **0.58** says the dense full product beats
summing the four even/odd sub-products — splitting a general multivector into
graded blocks is *not* a win when you genuinely need all grades.

## bench_simd — batched broadcast at increasing N (compact vs dense)

One fixed transform / operation broadcast over an array of primitives at
`N ∈ {1024, 65536, 1M}` — the vertex-buffer / point-cloud workload. A flat ratio
across N means the kernel scales cleanly; a ratio that drifts toward 1 as N grows
means the compact advantage erodes once the workload turns memory-bound.

| workload | CN | ref | CN/ref | check |
|---|---:|---:|---:|:---:|
| move point M·P·M̃ (PGA3) N=1024 | 37.8 µs | 153.3 µs | **0.25** | ok |
| move point M·P·M̃ (PGA3) N=65536 | 2.733 ms | 10.965 ms | **0.25** | ok |
| move point M·P·M̃ (PGA3) N=1M | 44.730 ms | 182.667 ms | **0.24** | ok |
| reflect point π·P·π (PGA3) N=1024 | 33.8 µs | 153.3 µs | **0.22** | ok |
| reflect point π·P·π (PGA3) N=65536 | 2.472 ms | 11.175 ms | **0.22** | ok |
| reflect point π·P·π (PGA3) N=1M | 42.259 ms | 187.505 ms | **0.23** | ok |
| rotor sandwich r·v·r̃ (VGA3) N=1024 | 14.4 µs | 33.3 µs | **0.43** | ok |
| rotor sandwich r·v·r̃ (VGA3) N=65536 | 1.076 ms | 2.461 ms | **0.44** | ok |
| rotor sandwich r·v·r̃ (VGA3) N=1M | 17.651 ms | 40.554 ms | **0.44** | ok |
| meet plane ∧ plane (PGA3) N=1024 | 3.6 µs | 24.1 µs | **0.15** | ok |
| meet plane ∧ plane (PGA3) N=65536 | 478.9 µs | 2.336 ms | **0.20** | ok |
| meet plane ∧ plane (PGA3) N=1M | 9.031 ms | 38.859 ms | **0.23** | ok |
| join point ∨ point (PGA3) N=1024 | 160.7 µs | 243.6 µs | **0.66** | ok |
| join point ∨ point (PGA3) N=65536 | 10.591 ms | 16.479 ms | **0.64** | ok |
| join point ∨ point (PGA3) N=1M | 169.846 ms | 267.139 ms | **0.64** | ok |
| rotor compose r₁r₂ (VGA3) N=1024 | 2.1 µs | 7.7 µs | **0.27** | ok |
| rotor compose r₁r₂ (VGA3) N=65536 | 323.1 µs | 824.7 µs | **0.39** | ok |
| rotor compose r₁r₂ (VGA3) N=1M | 6.595 ms | 15.680 ms | **0.42** | ok |

SIMD reading: `move point`, `reflect point`, and `rotor sandwich` hold a flat
~0.22–0.44 ratio from 1k to 1M — the sandwich kernels scale. `meet plane ∧ plane`
drifts 0.15 → 0.23 (the wedge is so cheap per element that the 1M case becomes
memory-bound and the compact lead narrows). `rotor compose` drifts 0.27 → 0.42:
the small-N case is dominated by the dense product's extra arithmetic, which
amortizes away once both sides are memory-bound. `join` is flat at 0.64 — the
regressive product is the limiter (see hot spots below), not the memory traffic.

## Hot spots

Ranked by how far CliffordNumbers trails its own faster alternative. The
multivector groups surface **two clusters where the compact path is slower than
it should be** — both with a known, validated one-line fix — plus a re-confirmation
of the two clusters already on the roadmap from `results.md`.

### 1. `reverse(::KVector)` has no fast sign-flip — drives PR-KVectorReverse (new)

| symptom | ratio | root cause |
|---|---:|---|
| `reverse line bivector (PGA3)` | **5.10×** | `reverse(::KVector)` is the only grade automorphism with no `KVector` fast path |

`src/math/duals.jl` defines closed-form `KVector` methods for `adjoint`, `conj`,
and `grade_involution` (`ifelse(iszero(K & 2), x, -x)` and friends), but **not**
for `reverse`. So `reverse(::KVector{2,PGA(3)})` falls back to the generic
`T(x[reverse.(BladeIndices(T))])` indexing path, which for a `KVector` runs the
`@generated` per-blade `is_same_blade` linear search in `to_index` — `158 ns`,
versus `2.1 ns` for the equivalent `-x`. (The dense `CliffordNumber` is *faster*
here, `31.9 ns`, only because its `to_index` is a single `Int % blade_count + 1`.)

Because for real scalars `reverse == adjoint`, the fix is the identical method:

```julia
# src/math/duals.jl, alongside the existing adjoint/conj/grade_involution methods
reverse(x::KVector{K}) where K = ifelse(iszero(K & 2), x, -x)
```

Verified: this matches the generic `reverse` for every grade `K` across VGA(3),
PGA(3), and VGA(4). Expected to collapse the 5.10× to ~1× and is the single
highest-leverage one-liner the enlarged suite found.

### 2. `inv` validation — re-confirms PR-SpinorFastPaths (fast `inv`)

| symptom | ratio | root cause |
|---|---:|---|
| `inv vs versor_inverse: rotor (VGA3)` | **5.40×** | `inv` computes `versor_inverse` then *validates* `x·inv≈1 && inv·x≈1` |
| `inv vs versor_inverse: motor (PGA3)` | **5.32×** | same validation, on the wider PGA(3) even type |

This is the same cluster the ℂ `microbench inv` (20.81×) and `contour integral`
(17.96×) rows flag, now confirmed for general even versors: rotors and motors are
versors by construction, so the `x·inv ≈ 1` validation is wasted work. Specialising
`inv = versor_inverse` for `EvenCliffordNumber` (and `KVector{1}`, already done)
collapses the ratio.

### 3. regressive product (`∨` / `join`) is the slow incidence path — drives PR-Perf-B

| symptom | ratio | root cause |
|---|---:|---|
| `join: point ∨ point → line (PGA3)` | **0.66** | `∨` = `right_complement(left_complement(x) ∧ left_complement(y))` — two complements + a wedge |
| `join: line ∨ point → plane (PGA3)` | **0.67** | same triple-hop pipeline |
| `join point ∨ point (PGA3)` (N=1M) | **0.64** | flat under broadcast — it is arithmetic-bound, not memory-bound |

Every other primitive op is **0.12–0.42×** against dense; `join` lags at ~0.66
because the regressive product is implemented as three passes (left-complement
both args, wedge, right-complement). ComputationalGeometricAlgebra.jl already
hand-rolls closed-form `Point∨Point` and `Line∨Point` for exactly this reason
(its `wiki/concepts/performance.md` PR-Perf-B notes 208 B / 8 allocs per generic
`∨`). A closed-form or allocation-free regressive in CliffordNumbers would lift
the whole incidence layer.

### Where the compact types already win big

No fix needed — these are the payoff rows that justify ever reaching for the
tight type instead of `CliffordNumber`:

- **sandwich kernels** (`move`/`reflect`/`rotate`): **0.14–0.26×**, flat from N=1k
  to N=1M — the rigid-motion and reflection workloads scale cleanly.
- **`exp`-generated rotors/motors**: **0.20–0.55×** (the closed-form `cos/sin`
  path on the compact even type vs the dense one).
- **CGA move point**: **≈0.001×** (≈740×) — dense storage is catastrophic in
  CGA(3)'s 32-dimensional algebra.
- **contractions / wedge / project**: **0.12–0.19×**.

## Fixable issues surfaced by this run (plan)

Two items beyond the performance hot spots above — the first is a correctness
crash the suite had to work around:

### PR-ExpTaylorClamp (new) — `exp` of a CGA bivector throws

`exp` of an `EvenCliffordNumber`/`KVector` bivector in a **mixed-signature**
algebra (CGA, STA, …) whose `abs2` lands in `(-1, 0)` throws
`DomainError: Cannot raise an integer x to a negative power`:

```
exp_taylor → r = div(exponent(2*abs2(x) + 1), 2)   # r < 0 when abs2 ∈ (-1, 0)
            y = x / (2^r)                            # 2^r with r::Int < 0 → DomainError
```

The scaling step (`src/math/exponential.jl:103-105`) only ever wants to scale a
*large* argument **down** before the Taylor expansion; a negative `r` would scale
it *up*, defeating the purpose. The one-line fix is to clamp:

```julia
r = max(0, div(exponent(2*abs2(x) + 1), 2))
```

Verified: random CGA(3) bivectors hit `abs2 ≈ -0.51 → r = -3` (throws today);
with `max(0, …)` the value takes the safe scaled-Taylor path. Until this lands,
`bench_multivector.jl` builds CGA motors as the geometric product of two
1-vectors (`randeven2`) rather than `exp(bivector)` to keep the suite runnable —
see the comment on `randeven2`.

### Harness note — `dense(::Number)`

The compact-vs-dense check needed `dense` to be a no-op on plain scalars so that
scalar-valued ops (`scalar_product`) compare directly. Added to `harness.jl`;
not a library change.

## Reproducing

```sh
BENCH_SECONDS=5 julia -t auto bench/runbenchmarks.jl
```

Absolute timings are machine-specific; the `CN/ref` ratios are the portable
signal. For the multivector groups, a ratio `> 1` is the actionable output: it
names a compact-type path that loses to brute-forcing the dense multivector.
