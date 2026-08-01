## Overview

This repository compares Column-and-Constraint Generation (CCG) approaches to solve robust optimization problems with mixed-integer recourse. Two problems are implemented:

- Unit Commitment (UC)
- Rostering

The job entry points are defined in [run.jl](run.jl): `unit_commitment_job()` and `rostering_job()` for single runs, and `run_unit_commitment_parallel()` / `run_rostering_parallel()` for batch/parallel runs.

## Parameters

### Unit Commitment (UC)
- budget: 1, 2, 3
- size: 1 = small, 2 = medium, 3 = large
- method: `CCGL` (ours) or `CCGM` (literature)
- time_limit: in seconds

### Rostering
- seed: fixes the random seed (try 1–100 for broader coverage)
- scale: 1 → horizon 21; 2 → horizon 42
- budget: 3, 6, 9 (as in the literature)
- method: one of the five formulations from the companion article:
  - `PdM`: componentwise dualization after fixing the binary recourse decision
  - `PdPrimeM`: componentwise dualization with a continuous auxiliary uncertainty copy
  - `HatPdPrimeM`: componentwise dualization with a binary auxiliary copy
  - `PdDoublePrimeM`: scalar dualization with a continuous auxiliary copy
  - `HatPdDoublePrimeM`: scalar dualization with a binary auxiliary copy
- time_limit: in seconds (2 hours used in the literature, but more is allowed)

The instance generator follows the staff-rostering experiment of Subramanyam
(2022): ten instances of size `(12, 3, 21)`, ten scaled instances of size
`(24, 6, 42)`, and budgets `3`, `6`, and `9`. Costs are left unchanged when
scaling, while staffing limits and demand parameters are doubled.

## Parallel Runs

The number of workers is set in [run.jl](run.jl) via `Nbworkers = 10`. Adjust that value to change parallelism.

Example:

### UC: 18 runs over budgets × sizes × methods with a 24h time limit
run_unit_commitment_parallel([1,2,3], [1,2,3], [CCGL, CCGM], 24*3600.0)

### Rostering: up to 1200 runs across seeds × budgets × scales × methods
run_rostering_parallel(1:100, [3,6,9], [1,2], [CCGL, CCGM], 2*3600.0)

## Five-formulation staff-rostering experiment

The complete experiment contains 300 runs (10 seeds × 2 sizes × 3 budgets ×
5 formulations). By default, the script follows the paper's two-hour limit and
uses eight Gurobi threads for one run at a time:

```bash
julia --project=. experiments/rostering_five_approaches.jl
```

The settings can be changed through environment variables. This command is a
small smoke benchmark over one seed, one size, and one budget:

```bash
ROSTERING_SEEDS=1 ROSTERING_SCALES=1 ROSTERING_BUDGETS=3 \
ROSTERING_TIME_LIMIT=60 ROSTERING_WORKERS=1 ROSTERING_THREADS=1 \
julia --project=. experiments/rostering_five_approaches.jl
```

For parallel runs, avoid CPU oversubscription. For example, use ten workers and
one Gurobi thread per run:

```bash
ROSTERING_WORKERS=10 ROSTERING_THREADS=1 \
julia --project=. experiments/rostering_five_approaches.jl
```

Each run writes its detailed metrics to `results/Rostering/`. The experiment
script additionally writes a timestamped `five_approaches_summary_*.csv` with
the solution time, gap, outer iterations, and total inner iterations.


## Results

Each run writes a CSV containing computational statistics (times, gaps, etc.) under the `results/` folder.
