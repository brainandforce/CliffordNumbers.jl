# Benchmark results

Recorded run of the [`bench/`](README.md) suite. Each workload is timed against
the Julia primitive its even subalgebra is isomorphic to; `CN/ref` is the
Clifford time divided by the primitive time (lower is better; `<1` means
CliffordNumbers is faster). Every row also passed its algebra-drift check
(Clifford result ≈ primitive result), so these numbers double as a correctness
snapshot.

## Environment

| | |
|---|---|
| Commit | `d3eb750` (branch `pr-benchsuite`) |
| Date | 2026-06-04 |
| Julia | 1.12.5 |
| CPU | x86_64 skylake, 16 threads |
| OS | Windows 11 (`x86_64-w64-mingw32`) |
| `BENCH_SECONDS` | 5.0 |
| Invocation | `BENCH_SECONDS=5 julia -t auto bench/runbenchmarks.jl` |

`* (per op)` / `inv (per op)` / etc. are per-element times from a batched
microbench (so the op can't be constant-folded); the ratio is the per-call
acceptance metric. All other rows are wall-clock for the whole workload.

## bench_complex — even VGA(2) ≅ ℂ

`Spinor{T} = EvenCliffordNumber{VGA(2),T,2}` vs `Complex{T}`.

| workload | CN | ref | CN/ref | check |
|---|---:|---:|---:|:---:|
| pointwise `*` (N=65536) | 160.5 µs | 160.8 µs | **1.00** | ok |
| `z² + c` (N=65536) | 161.5 µs | 161.3 µs | **1.00** | ok |
| Horner deg-8 (N=65536) | 590.4 µs | 653.0 µs | **0.90** | ok |
| Mandelbrot 512² (1 thread) | 18.359 ms | 9.040 ms | **2.03** | ok |
| Mandelbrot 512² (threaded) | 2.717 ms | 1.268 ms | **2.14** | ok |
| Julia set 512² | 13.899 ms | 6.611 ms | **2.10** | ok |
| contour integral (N=4096) | 675.3 µs | 38.3 µs | **17.63** | ok |
| microbench: `*` (per op) | 0.9 ns | 0.9 ns | **1.00** | ok |
| microbench: `inv` (per op) | 56.2 ns | 2.7 ns | **20.81** | ok |
| microbench: `versor_inverse` (per op) | 9.8 ns | 2.7 ns | **3.63** | ok |

## bench_quaternion — even VGA(3) ≅ ℍ

`Rotor{T} = EvenCliffordNumber{VGA(3),T,4}` vs `Quaternion{T}`.

| workload | CN | ref | CN/ref | check |
|---|---:|---:|---:|:---:|
| compose (per op) | 11.0 ns | 3.8 ns | **2.89** | ok |
| sandwich (per op) | 13.6 ns | 18.4 ns | **0.74** | ok |
| batched sandwich (N=100) | 1.4 µs | 1.9 µs | **0.74** | ok |
| batched sandwich (N=100k) | 1.566 ms | 2.091 ms | **0.75** | ok |
| batched sandwich (N=1M) | 16.799 ms | 21.137 ms | **0.79** | ok |
| slerp (per op) | 59.1 ns | 56.9 ns | **1.04** | ok |
| chain-compose (N=10k) | 52.8 µs | 70.3 µs | **0.75** | ok |

## Hot spots

Ranked by how far CliffordNumbers trails the primitive. The two clusters map
cleanly onto the two algebraic-win PRs.

### 1. `inv` validation — drives PR-SpinorFastPaths (fast `inv`)

| symptom | ratio | root cause |
|---|---:|---|
| microbench `inv` | **20.81×** | generic `inv` computes `versor_inverse` then *validates* `x·inv≈1 && inv·x≈1` |
| contour integral (N=4096) | **17.63×** | the `1/(z−pₖ)` in the integrand is exactly that `inv` path, per sample |
| microbench `versor_inverse` | **3.63×** | the validation-free path is already 5.7× faster than `inv`, but still trails `Complex.inv` |

`versor_inverse` (3.63×) being 5.7× faster than `inv` (20.81×) is direct
evidence that the validation is the dominant cost. The even VGA(2) subalgebra is
isomorphic to ℂ, where the inverse always exists when `abs2 ≠ 0`, so the
validation is wasted work — specialising `inv = versor_inverse` for these types
should collapse the 20.81× toward ~1× and pull the contour integral down with
it.

### 2. Scalar `*` / `²` in tight loops — drives PR-SpinorFastPaths (closed-form `*`/`^2`)

| symptom | ratio | root cause |
|---|---:|---|
| Mandelbrot 512² (threaded) | **2.14×** | scalar `z = z² + c` per pixel can't vectorize; the generic `map(muladd,…)` carries surplus `vxorpd; vaddsd` from the leading `muladd(0,0,x)` |
| Julia set 512² | **2.10×** | same scalar squaring path |
| Mandelbrot 512² (1 thread) | **2.03×** | same |
| rotor compose (per op) | **2.89×** | the VGA(3) even × even geometric product vs a hand-rolled Hamilton product |

Note the contrast: **vectorized** `pointwise *` and `z²+c` over an array are at
parity (1.00×) — LLVM hides the surplus instructions across the SIMD batch — but
the **scalar** iteration loops (Mandelbrot/Julia) pay ~2× per element. Emitting
the closed form `(a,b)² = (a²−b², 2ab)` directly should bring those to ~1×.

### Where CliffordNumbers already wins

The rotor/array side needs no algebraic fix — it is already faster than
`Quaternions.jl`:

- **sandwich** (single and batched): **0.74–0.79×** across N ∈ {1, 100, 100k, 1M}
- **chain-compose** (N=10k): **0.75×**
- **Horner deg-8**: **0.90×**

These are the throughput rows PR-ArrayVectorization builds on: the batched
sandwich already beats the quaternion baseline on the `isbits` scalar path, and
the ratio is flat from N=100 to N=1M, so SoA + SIMD has clean headroom to widen
the lead rather than recover a deficit.

## Reproducing

```sh
BENCH_SECONDS=5 julia -t auto bench/runbenchmarks.jl
```

Absolute timings are machine-specific; the `CN/ref` ratios are the portable
signal and what later PRs should move.
