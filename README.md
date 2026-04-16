# Structured Poisson: LDU vs DIA vs RB-DIA Gauss-Seidel Benchmark

A standalone C/C++ benchmark comparing three storage/ordering schemes for the
Gauss-Seidel smoother on a 3D Poisson problem discretised on a structured
hexahedral mesh:

- **LDU** — indirect addressing via (lower, diag, upper) face lists, mirroring
  OpenFOAM's native `lduMatrix` layout
- **DIA** — direct diagonal storage, using the fixed offset stencil {1, N, N²}
  of a structured hex mesh to eliminate indirect addressing
- **RB-DIA** — Red-Black ordered DIA variant with OpenMP parallelism across
  cells of the same colour, enabling shared-memory multi-core execution

The problem is a 3D manufactured-solution Poisson equation
`-∇²u = 3π² sin(πx) sin(πy) sin(πz)` with homogeneous Dirichlet BCs, discretised
on an N×N×N unit cube. All three solvers produce bit-for-bit identical
converged solutions, verified by matching L2 error against the manufactured
exact solution with second-order mesh convergence.

## Why this benchmark exists

This repo is the standalone testbed backing a larger OpenFOAM plugin project
([DIAGaussSeidel-Smoother-OpenFOAM](https://github.com/amartyadav/DIAGaussSeidel-Smoother-OpenFOAM)),
which replaces OpenFOAM's stock Gauss-Seidel smoother with a DIA-based variant
for structured hex meshes. Inside OpenFOAM the plugin yields ~50% reduction in
cache misses at all levels and a ~28% speedup on the GAMG smoother kernel.

The benchmark here serves two purposes:

1. **Clean, reproducible measurement of the core kernel** outside the OpenFOAM
   framework, useful for profiling under LIKWID and for studying algorithmic
   properties without the surrounding simulation machinery.

2. **A portable platform** for running the same kernel on non-x86 hardware —
   specifically Calligo Technologies' **Uttunga** RISC-V processor with native
   posit arithmetic. On Uttunga, the goal shifts from cache efficiency (which
   motivated the OpenFOAM plugin) to studying whether the higher precision of
   posit arithmetic reduces iteration counts to convergence on ill-conditioned
   problems — analogous to results previously observed in HYPRE BoomerAMG on
   multimaterial thermal analysis (aluminium + CFRP).

## Modes

The binary has three modes:

| Mode | Command | Purpose |
|---|---|---|
| Correctness | `./poisson_benchmark` | Runs all three solvers at N ∈ {4, 6, 8, 12, 16, 20}, verifies second-order convergence of the manufactured-solution L2 error |

| Kernel profiling | `./poisson_benchmark profile [N] [sweeps]` | Times a fixed number of sweeps (default N=100, 1000 sweeps) with optional LIKWID markers |
| Convergence study | `./poisson_benchmark convergence [N] [tol] [max_sweeps]` | Iterates until residual drops below tolerance, reports sweep count and final L2 error |

## Build

Requires a C++14 compiler and OpenMP. For LIKWID profiling, pass
`-DLIKWID_PERFMON` and link against `-llikwid`.

## Results so far (x86, AMD Zen4, 8 cores)

| Solver | 1000-sweep time at N=100 (s) | OMP scaling (8 threads) |
|---|---|---|
| DIA (serial) | ~7.5 | n/a |
| RB-DIA (1 thread) | ~8.2 | baseline |
| RB-DIA (8 threads) | ~2.1 | 3.9× speedup |

DIA is inherently serial due to the forward-push data dependency. RB-DIA breaks
that dependency via Red-Black ordering and scales well under OpenMP.

## Uttunga (in progress)

The same kernel compiles and runs on Uttunga v1.0 (RISC-V with posit arithmetic
co-processor, 4 cores at 400 MHz). The relevant experiment is not runtime
(which is expected to be slower at 400 MHz vs Zen4) but **iteration count to
convergence**. Preliminary results suggest that posit arithmetic lowers the
starting residual and reaches a given tolerance in fewer sweeps than IEEE 754
doubles on x86 — consistent with the BoomerAMG multimaterial study. Full
results pending further runs.

If this holds, it motivates a paper comparing iteration-count-to-convergence
across the two number systems for a set of structured-mesh iterative solvers.

## Author

Amartya Yadav
HPC Software Engineer, Calligo Technologies

## License

MIT
