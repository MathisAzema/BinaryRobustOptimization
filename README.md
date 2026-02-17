## Overview

This repository compares Column-and-Constraint Generation (CCG) approaches to solve robust optimization problems with mixed-integer recourse. Two problems are implemented:

- Unit Commitment (UC)
- Rostering

The job entry points are defined in [run.jl](run.jl): `unit_commitment_job()` and `rostering_job()` for single runs, and `run_unit_commitment_parallel()` / `run_rostering_parallel()` for batch/parallel runs.

## Parameters

### Unit Commitment (UC)
- budget: typical values 1, 3, 5
- size: 1 = small, 2 = medium, 3 = large
- method: `CCGL` (ours) or `CCGM` (literature)
- time_limit: in seconds

### Rostering
- seed: fixes the random seed (try 1–100 for broader coverage)
- scale: 1 → horizon 21; 2 → horizon 42
- budget: 3, 6, 9 (as in the literature)
- method: `CCGL` (ours) or `CCGM` (literature)
- time_limit: in seconds (2 hours used in the literature, but more is allowed)

## Parallel Runs

Example:

# UC: 18 runs over budgets × sizes × methods with a 24h time limit
run_unit_commitment_parallel([1,3,5], [1,2,3], [CCGL, CCGM], 24*3600.0)

# Rostering: up to 1200 runs across seeds × budgets × scales × methods
run_rostering_parallel(1:100, [3,6,9], [1,2], [CCGL, CCGM], 2*3600.0)
```

The number of workers is set in [run.jl](run.jl) via `Nbworkers = 10`. Adjust that value to change parallelism.

## Results

Each run writes a CSV containing computational statistics (times, gaps, etc.) under the `results/` folder.