# Benchmark results — `bench_complex` (even VGA(2) ≅ ℂ)

Recorded run of the `bench_complex` group of the [`bench/`](README.md) suite.
`Spinor{T} = EvenCliffordNumber{VGA(2),T,2}` is timed against `Complex{T}`;
`CN/ref` is the Clifford time over the primitive time, so a value below 1 means
CliffordNumbers is faster. Every row passed its algebra-drift check, so the
numbers double as a correctness snapshot.

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

The `(per op)` rows are per-element times from a batched microbench, so the op
cannot be constant-folded; the others are wall-clock for the whole workload.

## Current run

| workload | CN | ref | CN/ref | check |
|---|---:|---:|---:|:---:|
| pointwise `*` (N=65536) | 163.6 µs | 164.8 µs | 0.99 | ok |
| `z² + c` (N=65536) | 163.8 µs | 164.9 µs | 0.99 | ok |
| Horner deg-8 (N=65536) | 635.2 µs | 657.8 µs | 0.97 | ok |
| Mandelbrot 512² (1 thread) | 8.811 ms | 9.160 ms | 0.96 | ok |
| Mandelbrot 512² (threaded) | 1.269 ms | 1.325 ms | 0.96 | ok |
| Julia set 512² | 6.331 ms | 6.662 ms | 0.95 | ok |
| contour integral (N=4096) | 18.3 µs | 38.4 µs | 0.48 | ok |
| microbench: `*` (per op) | 0.9 ns | 0.9 ns | 1.00 | ok |
| microbench: `inv` (per op) | 1.1 ns | 2.8 ns | 0.39 | ok |
| microbench: `versor_inverse` (per op) | 1.1 ns | 2.8 ns | 0.39 | ok |

The spinor is at parity or better than `Complex{T}` on every workload.

## Before/after — PR-SpinorFastPaths

The baseline is the original `main` (`c5ca350`); the "after" column is this
branch. Both used `BENCH_SECONDS=5`, and every workload still passes its drift
check, so the fast paths are speedups, not behaviour changes.

| workload | before CN/ref | after CN/ref |
|---|---:|---:|
| microbench: `inv` (per op) | 20.81 | 0.39 |
| contour integral (N=4096) | 17.63 | 0.48 |
| microbench: `versor_inverse` (per op) | 3.63 | 0.39 |
| Mandelbrot 512² (1 thread) | 2.03 | 0.96 |
| Mandelbrot 512² (threaded) | 2.14 | 0.96 |
| Julia set 512² | 2.10 | 0.95 |

Three changes account for this:

- Closed-form reciprocal. `inv`/`versor_inverse` skip the validating `x·inv ≈ 1`
  check and emit the inverse as one fused expression (reverse folded into the
  signs, `abs2` inlined, one reciprocal multiplied through). This drives the
  `inv` and contour-integral rows. `inv` delegates to `versor_inverse`.
- Closed-form `abs2`. For the positive-definite even subalgebras every blade has
  `sign_of_square · reverse_sign = +1`, so `abs2` is a plain sum of squares. The
  escape loops test `abs2(z) > 4` every iteration, so this drives Mandelbrot and
  Julia.
- Closed-form `*`/`²`. Drops the un-foldable `muladd(0,0,x)` the generic kernel
  seeds with. It only helps where the multiply cannot vectorize; the already
  vectorized array rows stay at parity.

## Test coverage

`@testset "Even subalgebra fast paths"` in `test/operations.jl` covers the
specialized `inv`, `*`, and `abs2`: integer cases pin the product formula across
the ℂ / split / dual branches, and random cases cross-check against the generic
kernel, `versor_inverse`, `scalar_product`, and `Complex`.

## Reproducing

```sh
BENCH_SECONDS=5 julia -t auto bench/runbenchmarks.jl
```

Absolute timings are machine-specific; the `CN/ref` ratios are the portable
signal.
