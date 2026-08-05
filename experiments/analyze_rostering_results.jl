import Pkg

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT)

using CSV
using Printf
using Statistics

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

const DEFAULT_RESULTS_DIRECTORY = joinpath(PROJECT_ROOT, "results", "Rostering")
const DEFAULT_OUTPUT_FILE = joinpath(
    PROJECT_ROOT, "results",
    "rostering_numerical_table.txt",
)
const DEFAULT_TIME_LIMIT = 7200.0

# A run is considered solved when its final relative gap, stored as a percentage
# in the CSV file, is at most 0.1%.
const SOLVED_GAP_TOLERANCE_PERCENT = 0.1

# Worst-case statistics use every available instance, whether or not the
# complete rostering instance reached the target optimality gap. The final
# timing observation of each run is excluded because it may be incomplete;
# every retained time t[i,j,k] is normalized by its inner-iteration index k.

const METHODS = [
    "PdM",
    "PdPrimeM",
    "PdDoublePrimeM",
    "PdDoublePrimeUL",
    "HatPdPrimeM",
    "HatPdDoublePrimeM",
]

const METHOD_LATEX_LABELS = Dict(
    "PdM" => "\$P_d^M\$",
    "PdPrimeM" => "\$P_d^{\\prime M}\$",
    "HatPdPrimeM" => "\$\\widehat{P}_d^{\\prime M}\$",
    "PdDoublePrimeM" => "\$P_d^{\\prime\\prime M}\$",
    "HatPdDoublePrimeM" => "\$\\widehat{P}_d^{\\prime\\prime M}\$",
    "PdDoublePrimeUL" => "\$P_d^{\\prime\\prime UL}\$",
)

struct RunRecord
    seed::Int
    T::Int
    budget::Int
    method::String
    presolve::Bool
    time_limit::Float64
    elapsed_time::Float64
    gap_percent::Float64
    normalized_worst_case_times::Vector{Float64}
    path::String
end

struct AggregateStatistics
    solved::Int
    available::Int
    mean_solved_elapsed_time::Union{Missing, Float64}
    mean_worst_cases::Union{Missing, Float64}
    mean_normalized_worst_case_time::Union{Missing, Float64}
end

function usage()
    println(
        "Usage: julia --project=. experiments/analyze_rostering_results.jl " *
        "[results_directory] [output_file] [time_limit]",
    )
    println("Defaults:")
    println("  results_directory = $(DEFAULT_RESULTS_DIRECTORY)")
    println("  output_file       = $(DEFAULT_OUTPUT_FILE)")
    println("  time_limit        = $(Int(DEFAULT_TIME_LIMIT))")
end

function parse_arguments(arguments)
    if any(argument -> argument in ("-h", "--help"), arguments)
        usage()
        exit(0)
    end
    length(arguments) <= 3 || error("Expected at most three arguments. Use --help.")

    results_directory = length(arguments) >= 1 ? abspath(arguments[1]) : DEFAULT_RESULTS_DIRECTORY
    output_file = length(arguments) >= 2 ? abspath(arguments[2]) : DEFAULT_OUTPUT_FILE
    time_limit = if length(arguments) >= 3
        parse(Float64, arguments[3])
    else
        DEFAULT_TIME_LIMIT
    end
    return results_directory, output_file, time_limit
end

function require_scalar(scalars, metric, path)
    haskey(scalars, metric) || error("Missing metric '$metric' in $path")
    return scalars[metric]
end

function parse_number(::Type{T}, value, metric, path) where {T<:Real}
    parsed = tryparse(T, strip(value))
    parsed === nothing && error("Invalid $metric='$value' in $path")
    return parsed
end

function parse_boolean(value, metric, path)
    normalized = lowercase(strip(value))
    normalized in ("true", "1", "yes", "on") && return true
    normalized in ("false", "0", "no", "off") && return false
    error("Invalid $metric='$value' in $path")
end

function parse_result_file(path)
    scalars = Dict{String, String}()
    normalized_worst_case_times = Float64[]

    for row in CSV.File(path; types = Dict(:metric => String, :value => String))
        metric = strip(row.metric)
        value = strip(row.value)
        if metric == "Time_per_iteration"
            ismissing(row.j) && error(
                "Missing inner-iteration index for Time_per_iteration in $path",
            )
            inner_iteration = parse_number(
                Int,
                string(row.j),
                "inner iteration",
                path,
            )
            inner_iteration >= 1 || error(
                "Inner-iteration indices must be positive in $path",
            )
            worst_case_time = parse_number(
                Float64,
                value,
                "Time_per_iteration",
                path,
            )
            push!(
                normalized_worst_case_times,
                worst_case_time / inner_iteration,
            )
        elseif !isempty(metric)
            scalars[metric] = value
        end
    end

    # The final inner-master solve of a run may have been interrupted by the
    # global time limit. Exclude that last timing observation from every
    # statistic, while retaining the unfiltered count for the consistency
    # check against `inner_iterations` below.
    recorded_worst_case_count = length(normalized_worst_case_times)
    isempty(normalized_worst_case_times) || pop!(normalized_worst_case_times)

    filename_match = match(r"^(\d+)_", basename(path))
    filename_match === nothing && error("Cannot read the seed from $(basename(path))")
    seed = parse(Int, filename_match.captures[1])

    record = RunRecord(
        seed,
        parse_number(Int, require_scalar(scalars, "T", path), "T", path),
        parse_number(Int, require_scalar(scalars, "budget", path), "budget", path),
        require_scalar(scalars, "method", path),
        parse_boolean(
            require_scalar(scalars, "inner_presolve", path),
            "inner_presolve",
            path,
        ),
        parse_number(
            Float64,
            require_scalar(scalars, "time_limit", path),
            "time_limit",
            path,
        ),
        parse_number(Float64, require_scalar(scalars, "Time", path), "Time", path),
        parse_number(Float64, require_scalar(scalars, "gap", path), "gap", path),
        normalized_worst_case_times,
        path,
    )

    if haskey(scalars, "inner_iterations")
        reported_iterations = parse_number(
            Int,
            scalars["inner_iterations"],
            "inner_iterations",
            path,
        )
        reported_iterations == recorded_worst_case_count || @warn(
            "The number of timing rows does not match inner_iterations",
            file = path,
            reported = reported_iterations,
            timing_rows = recorded_worst_case_count,
        )
    end

    return record
end

function load_records(results_directory, selected_time_limit)
    isdir(results_directory) || error("Results directory does not exist: $results_directory")

    # Summary CSV files do not follow this result-file naming convention and are
    # intentionally ignored.
    result_filename = r"^\d+_\d+_[A-Za-z]+_\d+_presolve_(?:on|off)_\d+\.csv$"
    paths = sort(filter(
        path -> occursin(result_filename, basename(path)),
        readdir(results_directory; join = true),
    ))
    isempty(paths) && error("No rostering result CSV found in $results_directory")

    records = RunRecord[]
    ignored_time_limit = 0
    for path in paths
        record = parse_result_file(path)
        if isapprox(record.time_limit, selected_time_limit; atol = 1e-6, rtol = 0.0)
            record.method in METHODS && push!(records, record)
        else
            ignored_time_limit += 1
        end
    end
    isempty(records) && error(
        "No result uses the selected time limit of $(selected_time_limit) seconds",
    )

    keys_seen = Set{Tuple{Int, Int, Int, String, Bool}}()
    for record in records
        key = (record.seed, record.T, record.budget, record.method, record.presolve)
        key in keys_seen && error("Duplicate experimental configuration: $key")
        push!(keys_seen, key)
    end

    return records, ignored_time_limit
end

function is_solved(record)
    return isfinite(record.gap_percent) &&
           record.gap_percent <= SOLVED_GAP_TOLERANCE_PERCENT + 1e-10
end

function aggregate(records)
    solved_records = filter(is_solved, records)
    if isempty(records)
        return AggregateStatistics(0, 0, missing, missing, missing)
    end

    worst_case_counts = Int[]
    normalized_worst_case_times = Float64[]
    for record in records
        push!(worst_case_counts, length(record.normalized_worst_case_times))
        append!(
            normalized_worst_case_times,
            record.normalized_worst_case_times,
        )
    end

    mean_worst_cases = mean(worst_case_counts)
    mean_solved_elapsed_time = if isempty(solved_records)
        missing
    else
        mean(record.elapsed_time for record in solved_records)
    end
    mean_time = if isempty(normalized_worst_case_times)
        missing
    else
        mean(normalized_worst_case_times)
    end
    return AggregateStatistics(
        length(solved_records),
        length(records),
        mean_solved_elapsed_time,
        mean_worst_cases,
        mean_time,
    )
end

function configuration_statistics(records, T, budget, method, presolve)
    selected = filter(
        record -> record.T == T &&
                  record.budget == budget &&
                  record.method == method &&
                  record.presolve == presolve,
        records,
    )
    return aggregate(selected)
end

solved_cell(statistics) = statistics.available == 0 ? "--" :
                          "$(statistics.solved)/$(statistics.available)"
mean_solved_time_cell(statistics) =
    ismissing(statistics.mean_solved_elapsed_time) ? "--" :
    (@sprintf("%.1f", statistics.mean_solved_elapsed_time))
worst_case_cell(statistics) = statistics.available == 0 ? "--" :
                              (@sprintf("%.1f", statistics.mean_worst_cases))
mean_time_cell(statistics) =
    ismissing(statistics.mean_normalized_worst_case_time) ? "--" :
    (@sprintf("%.4f", statistics.mean_normalized_worst_case_time))

function write_presolve_table(io, records, selected_time_limit, presolve)
    selected_records = filter(record -> record.presolve == presolve, records)
    problem_classes = sort(unique([
        (record.T, record.budget) for record in selected_records
    ]))
    presolve_text = presolve ? "enabled" : "disabled"
    presolve_label = presolve ? "on" : "off"

    println(io, raw"\begin{table*}[htbp]")
    println(io, raw"\centering")
    println(
        io,
        "\\caption{Computational comparison of the six worst-case formulations " *
        "for the staff rostering problem when presolve is $presolve_text for the inner CCG master problems.}",
    )
    println(io, "\\label{tab:rostering-worst-case-presolve-$presolve_label}")
    println(io, raw"\setlength{\tabcolsep}{3pt}")
    println(io, raw"\renewcommand{\arraystretch}{0.9}")
    println(io, raw"\begin{tabular}{cc*{6}{c}}")
    println(io, raw"\toprule")

    header_cells = String[raw"$(T,k)$", "Statistic"]
    append!(header_cells, [METHOD_LATEX_LABELS[method] for method in METHODS])
    println(io, join(header_cells, " & "), raw" \\\\")
    println(io, raw"\midrule")

    for (class_index, (T, budget)) in enumerate(problem_classes)
        class_index > 1 && println(io, raw"\midrule")
        statistics_by_method = [
            configuration_statistics(
                selected_records,
                T,
                budget,
                method,
                presolve,
            )
            for method in METHODS
        ]
        class_cell = "\$(" * string(T) * ", " * string(budget) * ")\$"

        solved_cells = String[class_cell, "Solved"]
        append!(solved_cells, solved_cell.(statistics_by_method))
        println(io, join(solved_cells, " & "), raw" \\\\")

        solved_time_cells = String["", raw"Time (s)"]
        append!(solved_time_cells, mean_solved_time_cell.(statistics_by_method))
        println(io, join(solved_time_cells, " & "), raw" \\\\")

        # worst_case_cells = String["", raw"\# WC"]
        # append!(worst_case_cells, worst_case_cell.(statistics_by_method))
        # println(io, join(worst_case_cells, " & "), raw" \\\\")

        time_cells = String["", raw"$\overline{\tau}$ (s)"]
        append!(time_cells, mean_time_cell.(statistics_by_method))
        println(io, join(time_cells, " & "), raw" \\\\")
    end

    println(io, raw"\botrule")
    println(io, raw"\end{tabular}")
    println(io, raw"\end{table*}")
end

function write_latex_table(output_file, records, selected_time_limit)
    mkpath(dirname(output_file))
    open(output_file, "w") do io
        write_presolve_table(io, records, selected_time_limit, true)
        println(io)
        write_presolve_table(io, records, selected_time_limit, false)
    end
end

function report_coverage(records)
    println("Available files by (T, k):")
    problem_classes = sort(unique([
        (record.T, record.budget) for record in records
    ]))
    for (T, budget) in problem_classes
        class_records = filter(
            record -> record.T == T && record.budget == budget,
            records,
        )
        counts = [
            count(
                record -> record.method == method && record.presolve == presolve,
                class_records,
            )
            for method in METHODS
            for presolve in (false, true)
        ]
        println(
            "  (T=$T, k=$budget): $(minimum(counts))-$(maximum(counts)) " *
            "file(s) per method/presolve configuration",
        )
    end
end

function main(arguments = ARGS)
    results_directory, output_file, selected_time_limit = parse_arguments(arguments)
    records, ignored_time_limit = load_records(results_directory, selected_time_limit)
    write_latex_table(output_file, records, selected_time_limit)

    println("Analyzed $(length(records)) result files.")
    ignored_time_limit > 0 && println(
        "Ignored $ignored_time_limit file(s) with another time limit.",
    )
    report_coverage(records)
    println("LaTeX table written to $output_file")
end

main()
