# run_unit_commitment_parallel([1,2], [1, 2, 3], [CCGL, CCGM], 10.0)
# run_rostering_parallel(1:10, [1], [1, 2], [CCGL, CCGM], 10.0)

import Pkg
Pkg.activate(".")

using Distributed

Nbworkers = 10
if nworkers() >= Nbworkers+1
    rmprocs(workers())
    addprocs(Nbworkers)
else
    addprocs(Nbworkers - nworkers())
end

# Ensure package availability on workers
@everywhere using BinaryRobustOptimization
@everywhere set_solver_Gurobi()
@everywhere set_num_threads(1)

@everywhere function unit_commitment_job(budget, size, method, time_limit)
    if size == 1
        problem = UnitCommitment("6bus_JEAS", budget, 1)
    elseif size == 2
        problem = UnitCommitment("medium", budget, 5)
    elseif size == 3
        problem = UnitCommitment("118_syst_JEAS", budget, 35)
    else
        error("Unknown size value: $(size)")
    end
    run_ccg(problem, method, time_limit)
end

function run_unit_commitment_parallel(budget_list::Vector{Int}, size_list, method_list, time_limit)
    combos = [(b, s, m) for b in budget_list for s in size_list for m in method_list]

    results = pmap(combos) do (b, s, m)
        unit_commitment_job(b, s, m, time_limit)
    end
    return 
end

@everywhere function rostering_job(seed::Int, scale, budget, method, time_limit)
    problem = Rostering(budget, scale, scale, seed)
    run_ccg(problem, method, time_limit)
end

function run_rostering_parallel(seed_list, budget_list, scale_list, method_list, time_limit)
    combos = [(seed, b, s, m) for seed in seed_list for b in budget_list for s in scale_list for m in method_list]

    results = pmap(combos) do (seed, b, s, m)
        rostering_job(seed, s, b, m, time_limit)
    end
    return 
end