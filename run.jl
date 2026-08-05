# run_unit_commitment_parallel([1,2], [1, 2, 3], [CCGL, CCGM], 10.0)
# run_rostering_parallel(1:10, [1], [1, 2], [CCGL, CCGM], 10.0)

import Pkg
Pkg.activate(".")

using Distributed

Nbworkers = 10
println(nworkers())
if nworkers() >= Nbworkers+1
    rmprocs(workers())
    addprocs(Nbworkers)
elseif nworkers() ==1
    addprocs(Nbworkers - nworkers()+1)
else
    addprocs(Nbworkers - nworkers())
end
println(nworkers())

# Ensure package availability on workers
@everywhere using BinaryRobustOptimization
@everywhere set_solver_Gurobi()
@everywhere set_num_threads(1)

@everywhere function rostering_job(
    seed::Int,
    scale,
    budget,
    method,
    time_limit;
    inner_presolve::Bool=true,
)
    problem = Rostering(budget, scale, scale, seed)
    run_ccg(problem, method, time_limit; inner_presolve = inner_presolve)
end

function run_rostering_parallel(
    seed_list,
    budget_list,
    scale_list,
    method_list,
    time_limit;
    inner_presolve::Bool=true,
)
    combos = [(seed, b, s, m) for seed in seed_list for b in budget_list for s in scale_list for m in method_list]

    results = pmap(combos) do (seed, b, s, m)
        rostering_job(
            seed,
            s,
            b,
            m,
            time_limit;
            inner_presolve = inner_presolve,
        )
    end
    return 
end
