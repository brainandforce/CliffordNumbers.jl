# Benchmark results: `bench_quaternion` (even VGA(3) ≅ ℍ)

Recorded run of the `bench_quaternion` group of the [`bench/`](README.md) suite.
`Rotor{T} = EvenCliffordNumber{VGA(3),T,4}` is timed against `Quaternion{T}` from
`Quaternions.jl`; `CN/ref` is the Clifford time over the primitive time, so a
value below 1 means CliffordNumbers is faster. Every row passed its drift check.

The compose and chain workloads compare against ℍ multiplication through
`as_hamilton` (see `harness.jl`), which swaps the two trailing bivector slots so
that `as_hamilton(a*b) == as_hamilton(a) * as_hamilton(b)`. The sandwich
workloads check norm preservation instead. See the README's note on the ℍ
isomorphism.

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

The `(per op)` rows are per-element times from a batched microbench; the others
are wall-clock.

## Current run

| workload | CN | ref | CN/ref | check |
|---|---:|---:|---:|:---:|
| compose (per op) | 11.0 ns | 4.1 ns | 2.68 | ok |
| sandwich (per op) | 8.7 ns | 18.4 ns | 0.47 | ok |
| batched sandwich (N=100) | 900.0 ns | 1.9 µs | 0.47 | ok |
| batched sandwich (N=100k) | 1.152 ms | 2.137 ms | 0.54 | ok |
| batched sandwich (N=1M) | 11.900 ms | 21.641 ms | 0.55 | ok |
| slerp (per op) | 58.5 ns | 58.8 ns | 0.99 | ok |
| chain-compose (N=10k) | 53.0 µs | 70.2 µs | 0.75 | ok |

The rotor matches or beats `Quaternions.jl` everywhere except the single-product
compose: the sandwich kernels are 0.47–0.55 across N, chain-compose is 0.75, and
slerp is at parity.

## The compose gap (2.68×)

The single rotor × rotor geometric product trails the hand-tuned `Quaternion`
`*`. A closed-form Hamilton product was tried (a scalar `muladd` chain and a
SIMD-friendly broadcast form) and dropped. Unlike the 2-component spinor, the
generic `@generated` kernel already SIMD-vectorizes the 4-component product over
4-tuples, so removing its `zero_tuple` seed gained nothing, and the closed forms
measured slower in the compute-bound reduce case (chain-compose ≈ 5.3 ns/mul
generic vs ≈ 7.4–7.9 ns/mul closed). The gap is the cost of the
promotion/`muladd` path, not the product itself, which is why the batched compose
and chain-compose still come out ahead. See the NOTE in `src/math/multiply.jl`.

## `inv` on the rotor subalgebra

The closed-form reciprocal from PR-SpinorFastPaths also covers the VGA(3) rotor.
A focused per-op microbench (`map(inv, xs)` over 1000 elements):

| `inv` per op | previous (generic versor) | closed form | primitive ref |
|---|---:|---:|---:|
| rotor (VGA(3) ≅ ℍ) | 11.1 ns | 2.1 ns | `Quaternion.inv` 29.2 ns |

This collapses 4 component divides to 1 reciprocal and 4 muls. Both Base and
`Quaternions.jl` `inv` do overflow-scaling the versor formula skips. The
validation-overhead view of the same effect is the `inv vs versor_inverse: rotor`
row in [`results_multivector.md`](results_multivector.md), now 1.00.

## Reproducing

```sh
BENCH_SECONDS=5 julia -t auto bench/runbenchmarks.jl
```

Absolute timings are machine-specific; the `CN/ref` ratios are the portable
signal.
