# Benchmark results — `bench_quaternion` (even VGA(3) ≅ ℍ)

Recorded run of the `bench_quaternion` group of the [`bench/`](README.md) suite
on the consolidated `clifford-improvements` branch. `Rotor{T} =
EvenCliffordNumber{VGA(3),T,4}` (Cl⁺(3) ≅ ℍ) is timed against `Quaternion{T}`
from `Quaternions.jl`; `CN/ref` is the Clifford time over the primitive time
(lower is better; `<1` means CliffordNumbers is faster). Every row passed its
algebra-drift check, so the numbers double as a correctness snapshot.

The compose/chain workloads compare against ℍ multiplication through
`as_hamilton` (defined in `harness.jl`), which swaps the two trailing bivector
slots so that `as_hamilton(a*b) == as_hamilton(a) * as_hamilton(b)`; the sandwich
workloads check norm preservation instead, which is convention-free. See the
README's "note on the ℍ isomorphism".

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

`(per op)` rows are per-element times from a batched microbench; all other rows
are wall-clock for the whole workload.

## Current run

| workload | CN | ref | CN/ref | check |
|---|---:|---:|---:|:---:|
| compose (per op) | 11.0 ns | 4.1 ns | **2.68** | ok |
| sandwich (per op) | 8.7 ns | 18.4 ns | **0.47** | ok |
| batched sandwich (N=100) | 900.0 ns | 1.9 µs | **0.47** | ok |
| batched sandwich (N=100k) | 1.152 ms | 2.137 ms | **0.54** | ok |
| batched sandwich (N=1M) | 11.900 ms | 21.641 ms | **0.55** | ok |
| slerp (per op) | 58.5 ns | 58.8 ns | **0.99** | ok |
| chain-compose (N=10k) | 53.0 µs | 70.2 µs | **0.75** | ok |

## Reading

The rotor side needs no algebraic fix — it already beats or matches
`Quaternions.jl` on every workload except the single-product compose:

- **sandwich** (single and batched): **0.47–0.55×** across N ∈ {1, 100, 100k, 1M}
  — the rotation kernel scales cleanly (the `bench_simd` group in
  [`results_multivector.md`](results_multivector.md) confirms the broadcast path).
- **chain-compose** (N=10k): **0.75×**.
- **slerp**: **0.99×** (parity).

### The one remaining gap — `compose` at 2.68×

The single rotor×rotor geometric product trails the hand-tuned `Quaternion` `*`.
A closed-form Hamilton product was tried (both a scalar `muladd` chain and a
SIMD-friendly broadcast form) and **dropped**: unlike the 2-component spinor, the
generic `@generated` kernel already SIMD-vectorizes the 4-component product over
4-tuples, so removing its `zero_tuple` seed gave nothing, and the closed forms
measured *slower* in the compute-bound reduce case (`chain-compose` ≈ 5.3 ns/mul
generic vs ≈ 7.4–7.9 ns/mul closed). The gap is not the geometric product itself
— `Quaternions.jl` `*` is a hand-tuned 4-term form, while the GA value passes
through promotion/`muladd` machinery that is already well vectorized, which is why
the *batched* compose (in `bench_simd`) and `chain-compose` both come out ahead.
See the NOTE in `src/math/multiply.jl` and the negative result recorded here.

## `inv` on the rotor subalgebra

The ℍ workloads above exercise the geometric product, not `inv`. The closed-form
reciprocal that PR-SpinorFastPaths added also covers the VGA(3) rotor: a focused
per-op microbench (`map(inv, xs)` over 1000 elements) gives

| `inv` per op | previous (generic versor) | closed form | primitive ref |
|---|---:|---:|---:|
| rotor (VGA(3) ≅ ℍ) | 11.1 ns | **2.1 ns** | `Quaternion.inv` 29.2 ns |

Collapsing **4 component divides → 1 reciprocal + 4 muls** is ~5× over the old
path and ~14× over `Quaternion.inv`; both Base and `Quaternions.jl` `inv` do
overflow-scaling the versor formula skips. (The dense-vs-validation view of the
same effect is the `inv vs versor_inverse: rotor (VGA3)` row in
[`results_multivector.md`](results_multivector.md), now **1.00×**.)

## Reproducing

```sh
BENCH_SECONDS=5 julia -t auto bench/runbenchmarks.jl
```

Absolute timings are machine-specific; the `CN/ref` ratios are the portable
signal.
