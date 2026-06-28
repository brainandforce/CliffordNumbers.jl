# Benchmark results — `bench_complex` (even VGA(2) ≅ ℂ)

Recorded run of the `bench_complex` group of the [`bench/`](README.md) suite on
the consolidated `clifford-improvements` branch (all topic-PR work applied).
`Spinor{T} = EvenCliffordNumber{VGA(2),T,2}` is timed against `Complex{T}`;
`CN/ref` is the Clifford time over the primitive time (lower is better; `<1`
means CliffordNumbers is faster). Every row passed its algebra-drift check
(Clifford result ≈ `Complex` result), so the numbers double as a correctness
snapshot.

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

`* (per op)` / `inv (per op)` / `versor_inverse (per op)` are per-element times
from a batched microbench (so the op can't be constant-folded); the ratio is the
per-call acceptance metric. All other rows are wall-clock for the whole workload.

## Current run

| workload | CN | ref | CN/ref | check |
|---|---:|---:|---:|:---:|
| pointwise `*` (N=65536) | 163.6 µs | 164.8 µs | **0.99** | ok |
| `z² + c` (N=65536) | 163.8 µs | 164.9 µs | **0.99** | ok |
| Horner deg-8 (N=65536) | 635.2 µs | 657.8 µs | **0.97** | ok |
| Mandelbrot 512² (1 thread) | 8.811 ms | 9.160 ms | **0.96** | ok |
| Mandelbrot 512² (threaded) | 1.269 ms | 1.325 ms | **0.96** | ok |
| Julia set 512² | 6.331 ms | 6.662 ms | **0.95** | ok |
| contour integral (N=4096) | 18.3 µs | 38.4 µs | **0.48** | ok |
| microbench: `*` (per op) | 0.9 ns | 0.9 ns | **1.00** | ok |
| microbench: `inv` (per op) | 1.1 ns | 2.8 ns | **0.39** | ok |
| microbench: `versor_inverse` (per op) | 1.1 ns | 2.8 ns | **0.39** | ok |

The spinor is now at parity-or-better than `Complex{T}` on every ℂ workload, and
`inv` / `versor_inverse` / the contour integral are well below it.

## Before/after — what PR-SpinorFastPaths moved

**Before** is the original baseline on `main` (`c5ca350`, recorded before any of
the topic-PR work); **after** is this consolidated branch. Both runs used
`BENCH_SECONDS=5` on a comparable machine. Every workload still passes its
algebra-drift check — the fast paths are pure speedups, not behaviour changes.

| workload | before CN/ref | after CN/ref | note |
|---|---:|---:|---|
| microbench: `inv` (per op) | 20.81 | **0.39** | closed-form reciprocal, no validation |
| contour integral (N=4096) | 17.63 | **0.48** | integrand is `1/(z−pₖ)`, i.e. `inv` per sample |
| microbench: `versor_inverse` (per op) | 3.63 | **0.39** | shares the closed form `inv` delegates to |
| Mandelbrot 512² (1 thread) | 2.03 | **0.96** | closed-form `abs2` + `*`/`²` |
| Mandelbrot 512² (threaded) | 2.14 | **0.96** | same |
| Julia set 512² | 2.10 | **0.95** | same |
| Horner deg-8 (N=65536) | 0.90 | **0.97** | already vectorized; noise |
| pointwise `*` (N=65536) | 1.00 | **0.99** | already vectorized; noise |
| `z² + c` (N=65536) | 1.00 | **0.99** | already vectorized; noise |
| microbench: `*` (per op) | 1.00 | **1.00** | unchanged |

Three changes hit three bottlenecks:

- **Closed-form reciprocal.** Beyond dropping the validating `x·inv ≈ 1` check,
  the inverse is emitted as one fused expression — reverse folded into the signs
  (no intermediate `x'`), `abs2` inlined, modulus reciprocated once and multiplied
  through (1 divide vs 2). This takes `inv` from **20.81× → 0.39×** and the contour
  integral from **17.63× → 0.48×**. The closed form lives in `versor_inverse`
  (so the public versor inverse is fast too, **3.63× → 0.39×**); `inv` delegates
  to it.
- **Closed-form `abs2`.** `abs2(x) = scalar_product(x, x')` built the reverse and
  ran a three-tuple `mapreduce`; for the positive-definite even subalgebras every
  blade has `sign_of_square · reverse_sign = +1`, so it is a plain sum of squares.
  The `abs2(z) > 4` escape test runs every iteration, so Mandelbrot/Julia went
  **~2.1× → ~0.95×** (now faster than `Complex{T}`).
- **Closed-form `*`/`²`.** Removes the un-foldable `muladd(0,0,x)` the generic
  kernel seeds with. It matters only where the multiply can't vectorize: the
  scalar escape loops improve, while the already-vectorized array rows
  (`pointwise *`, `z² + c`) sit at parity throughout.

## Test coverage

The specialized methods (`inv`, `*`, `abs2` on the even subalgebras) are covered
by `@testset "Even subalgebra fast paths"` in `test/operations.jl`: exact integer
cases pin the product formula across all three σ branches (ℂ / split / dual), and
randomized cases cross-check the closed forms against the generic kernel,
`versor_inverse`, `scalar_product`, and `Complex`.

## Reproducing

```sh
BENCH_SECONDS=5 julia -t auto bench/runbenchmarks.jl
```

Absolute timings are machine-specific; the `CN/ref` ratios are the portable
signal.
