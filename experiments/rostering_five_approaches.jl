import Pkg
using Distributed
using Dates

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT)
cd(PROJECT_ROOT)

function parse_integer_spec(specification::AbstractString)
    if occursin(':', specification)
        endpoints = parse.(Int, split(specification, ':'))
        length(endpoints) == 2 || error("Expected a range such as 1:10")
        return collect(endpoints[1]:endpoints[2])
    end
    return parse.(Int, split(specification, ','))
end

function parse_boolean_spec(specification::AbstractString)
    normalized = lowercase(strip(specification))
    normalized in ("true", "1", "yes", "on") && return true
    normalized in ("false", "0", "no", "off") && return false
    error("Expected true/false for ROSTERING_INNER_PRESOLVE, got: $(specification)")
end

time_limit = parse(Float64, get(ENV, "ROSTERING_TIME_LIMIT", "7200"))
requested_workers = parse(Int, get(ENV, "ROSTERING_WORKERS", "1"))
threads_per_run = parse(Int, get(ENV, "ROSTERING_THREADS", "8"))
inner_presolve = parse_boolean_spec(get(ENV, "ROSTERING_INNER_PRESOLVE", "true"))
seeds = parse_integer_spec(get(ENV, "ROSTERING_SEEDS", "1:10"))
budgets = parse_integer_spec(get(ENV, "ROSTERING_BUDGETS", "3,6,9"))
scales = parse_integer_spec(get(ENV, "ROSTERING_SCALES", "1,2"))

requested_workers >= 1 || error("ROSTERING_WORKERS must be positive")
threads_per_run >= 1 || error("ROSTERING_THREADS must be positive")

if nworkers() < requested_workers
    addprocs(
        requested_workers - nworkers();
        exeflags = "--project=$(PROJECT_ROOT)",
    )
end

@everywhere begin
    using BinaryRobustOptimization
    using DataFrames
    using CSV
    set_solver_Gurobi()
end

methods = collect(five_lagrangian_formulations())
configurations = [
    (seed = seed, scale = scale, budget = budget, method = method)
    for seed in seeds
    for scale in scales
    for budget in budgets
    for method in methods
]

@everywhere function run_rostering_configuration(
    config,
    time_limit,
    threads_per_run,
    inner_presolve,
)
    set_num_threads(threads_per_run)
    try
        problem = Rostering(config.budget, config.scale, config.scale, config.seed)
        raw_result = run_ccg(
            problem,
            config.method,
            time_limit;
            inner_presolve = inner_presolve,
        )
        inner_times = raw_result[7]
        relative_gap = raw_result[6]
        solved = isfinite(relative_gap) && relative_gap <= 1e-3
        status = if solved
            "solved"
        elseif raw_result[4] >= 0.99 * time_limit
            "time_limit"
        else
            "stopped"
        end
        return (
            seed = config.seed,
            scale = config.scale,
            I = problem.I,
            J = problem.J,
            T = problem.T,
            budget = config.budget,
            method = string(config.method),
            inner_presolve = inner_presolve,
            time_limit = time_limit,
            elapsed_time = raw_result[4],
            lower_bound = raw_result[5],
            gap_percent = 100 * relative_gap,
            outer_iterations = length(inner_times),
            inner_iterations = sum(length, inner_times),
            solved = solved,
            status = status,
            error = "",
            result_file = raw_result[1] * ".csv",
        )
    catch exception
        return (
            seed = config.seed,
            scale = config.scale,
            I = 12 * config.scale,
            J = 3 * config.scale,
            T = 21 * config.scale,
            budget = config.budget,
            method = string(config.method),
            inner_presolve = inner_presolve,
            time_limit = time_limit,
            elapsed_time = NaN,
            lower_bound = NaN,
            gap_percent = NaN,
            outer_iterations = 0,
            inner_iterations = 0,
            solved = false,
            status = "error",
            error = sprint(showerror, exception),
            result_file = "",
        )
    end
end

println(
    "Running $(length(configurations)) configurations with " *
    "$(requested_workers) worker(s), $(threads_per_run) Gurobi thread(s) per run, " *
    "MP_inner presolve=$(inner_presolve), and a $(time_limit)-second time limit.",
)

results = pmap(configurations) do config
    run_rostering_configuration(
        config,
        time_limit,
        threads_per_run,
        inner_presolve,
    )
end

summary = DataFrame(results)
results_directory = joinpath(PROJECT_ROOT, "results", "Rostering")
mkpath(results_directory)
timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
summary_path = joinpath(
    results_directory,
    "five_approaches_summary_$(timestamp).csv",
)
CSV.write(summary_path, summary)

println("Summary written to $(summary_path)")
println("Finished without errors: $(count(!=("error"), summary.status))/$(nrow(summary))")
println("Solved to 0.1%: $(count(summary.solved))/$(nrow(summary))")
println("Time limit reached: $(count(==("time_limit"), summary.status))/$(nrow(summary))")
