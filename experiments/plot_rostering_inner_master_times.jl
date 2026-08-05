import Pkg

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT)

# Render plots without requiring a graphical display.
ENV["GKSwstype"] = "100"

using CSV
using LaTeXStrings
using Plots

const DEFAULT_RESULTS_DIRECTORY = joinpath(PROJECT_ROOT, "results", "Rostering")
const DEFAULT_OUTPUT_DIRECTORY = joinpath(
    PROJECT_ROOT,
    "results",
    "rostering_inner_master_plots",
)
const DEFAULT_TIME_LIMIT = 7200.0

const METHODS = [
    "PdM",
    "PdPrimeM",
    "PdDoublePrimeM",
    "PdDoublePrimeUL",
]

const METHOD_LABELS = Dict(
    "PdM" => L"P_d^M",
    "PdPrimeM" => L"P_d^{\prime M}",
    "PdDoublePrimeM" => L"P_d^{\prime\prime M}",
    "PdDoublePrimeUL" => L"P_d^{\prime\prime UL}",
)

# Color-blind-friendly colors distinguish the four continuous curves.
const METHOD_COLORS = ["#0072B2", "#E69F00", "#009E73", "#CC79A7"]
const MIN_INSTANCES_PER_POINT = 10

function usage()
    println(
        "Usage: julia --project=. experiments/plot_rostering_inner_master_times.jl " *
        "[results_directory] [output_directory] [time_limit]",
    )
    println("Defaults:")
    println("  results_directory = $(DEFAULT_RESULTS_DIRECTORY)")
    println("  output_directory  = $(DEFAULT_OUTPUT_DIRECTORY)")
    println("  time_limit        = $(Int(DEFAULT_TIME_LIMIT))")
end

function parse_arguments(arguments)
    if any(argument -> argument in ("-h", "--help"), arguments)
        usage()
        exit(0)
    end
    length(arguments) <= 3 || error("Expected at most three arguments. Use --help.")

    results_directory = if length(arguments) >= 1
        abspath(arguments[1])
    else
        DEFAULT_RESULTS_DIRECTORY
    end
    output_directory = if length(arguments) >= 2
        abspath(arguments[2])
    else
        DEFAULT_OUTPUT_DIRECTORY
    end
    time_limit = if length(arguments) >= 3
        parse(Float64, arguments[3])
    else
        DEFAULT_TIME_LIMIT
    end
    return results_directory, output_directory, time_limit
end

function parse_boolean(value, metric, path)
    normalized = lowercase(strip(value))
    normalized in ("true", "1", "yes", "on") && return true
    normalized in ("false", "0", "no", "off") && return false
    error("Invalid $metric='$value' in $path")
end

function parse_number(::Type{T}, value, metric, path) where {T<:Real}
    parsed = tryparse(T, strip(value))
    parsed === nothing && error("Invalid $metric='$value' in $path")
    return parsed
end

function load_inner_time_averages(results_directory, selected_time_limit)
    isdir(results_directory) || error("Results directory does not exist: $results_directory")

    result_filename = r"^\d+_\d+_[A-Za-z]+_\d+_presolve_(?:on|off)_\d+\.csv$"
    paths = sort(filter(
        path -> occursin(result_filename, basename(path)),
        readdir(results_directory; join = true),
    ))
    isempty(paths) && error("No rostering result CSV found in $results_directory")

    # Key: (T, budget, method, presolve, inner iteration).
    time_sums = Dict{Tuple{Int, Int, String, Bool, Int}, Float64}()
    observation_counts = Dict{Tuple{Int, Int, String, Bool, Int}, Int}()
    instance_counts = Dict{Tuple{Int, Int, String, Bool, Int}, Int}()
    files_used = 0

    for path in paths
        horizon = nothing
        budget = nothing
        method = nothing
        presolve = nothing
        time_limit = nothing
        used_this_file = false
        timing_rows = Tuple{Int, Float64}[]
        iterations_seen_in_file = Set{Tuple{Int, Int, String, Bool, Int}}()

        for row in CSV.File(path; types = Dict(:metric => String, :value => String))
            metric = strip(row.metric)
            value = strip(row.value)

            if metric == "T"
                horizon = parse_number(Int, value, "T", path)
            elseif metric == "budget"
                budget = parse_number(Int, value, "budget", path)
            elseif metric == "method"
                method = value
            elseif metric == "inner_presolve"
                presolve = parse_boolean(value, "inner_presolve", path)
            elseif metric == "time_limit"
                time_limit = parse_number(Float64, value, "time_limit", path)
            elseif metric == "Time_per_iteration"
                method in METHODS || continue
                time_limit === nothing && error(
                    "Missing time_limit before timing rows in $path",
                )
                isapprox(time_limit, selected_time_limit; atol = 1e-6, rtol = 0.0) ||
                    continue
                horizon === nothing && error("Missing T before timing rows in $path")
                budget === nothing && error("Missing budget before timing rows in $path")
                presolve === nothing && error("Missing presolve flag before timing rows in $path")
                ismissing(row.j) && error("Missing inner-iteration index in $path")

                inner_iteration = parse_number(
                    Int,
                    string(row.j),
                    "inner iteration",
                    path,
                )
                inner_iteration >= 1 || error(
                    "Inner-iteration indices must be positive in $path",
                )
                elapsed_time = parse_number(
                    Float64,
                    value,
                    "Time_per_iteration",
                    path,
                )

                push!(timing_rows, (inner_iteration, elapsed_time))
            end
        end

        # The last inner-master solve in a run may have been stopped by the
        # global time limit, so omit that final timing row before aggregating.
        isempty(timing_rows) || pop!(timing_rows)
        for (inner_iteration, elapsed_time) in timing_rows
            key = (horizon, budget, method, presolve, inner_iteration)
            time_sums[key] = get(time_sums, key, 0.0) + elapsed_time
            observation_counts[key] = get(observation_counts, key, 0) + 1
            push!(iterations_seen_in_file, key)
            used_this_file = true
        end
        for key in iterations_seen_in_file
            instance_counts[key] = get(instance_counts, key, 0) + 1
        end
        files_used += used_this_file
    end

    isempty(time_sums) && error(
        "No timing observation found for the requested methods and time limit",
    )
    return time_sums, observation_counts, instance_counts, files_used
end

function problem_classes(time_sums, presolve)
    classes = unique([
        (key[1], key[2]) for key in keys(time_sums) if key[4] == presolve
    ])
    # Arrange each budget on one row: T=21 on the left and T=42 on the right.
    return sort(classes; by = problem_class -> (problem_class[2], problem_class[1]))
end

function integer_ticks(max_iteration; target_intervals = 6)
    max_iteration <= 10 && return collect(1:max_iteration)

    raw_step = max_iteration / target_intervals
    magnitude = 10.0 ^ floor(log10(raw_step))
    normalized_step = raw_step / magnitude
    nice_multiplier = if normalized_step <= 1
        1
    elseif normalized_step <= 2
        2
    elseif normalized_step <= 5
        5
    else
        10
    end
    step = max(1, round(Int, nice_multiplier * magnitude))
    ticks = collect(step:step:max_iteration)
    first(ticks) == 1 || pushfirst!(ticks, 1)
    return ticks
end

function method_series(
    time_sums,
    observation_counts,
    instance_counts,
    horizon,
    budget,
    method,
    presolve,
)
    candidate_iterations = sort(unique([
        key[5] for key in keys(time_sums)
        if key[1] == horizon &&
           key[2] == budget &&
           key[3] == method &&
           key[4] == presolve
    ]))
    iterations = Int[]
    for iteration in candidate_iterations
        key = (horizon, budget, method, presolve, iteration)
        get(instance_counts, key, 0) < MIN_INSTANCES_PER_POINT && break
        push!(iterations, iteration)
    end
    means = [
        time_sums[(horizon, budget, method, presolve, iteration)] /
        observation_counts[(horizon, budget, method, presolve, iteration)]
        for iteration in iterations
    ]
    return iterations, means
end

function create_panel(
    time_sums,
    observation_counts,
    instance_counts,
    presolve,
)
    classes = problem_classes(time_sums, presolve)
    isempty(classes) && error("No data for presolve=$(presolve)")

    panels = Plots.Plot[]
    for (horizon, budget) in classes
        panel_max_iteration = 1
        panel = plot(
            ;
            title = "(T, k) = ($horizon, $budget)",
            xlabel = "Inner iteration",
            ylabel = "Mean time (s)",
            legend = false,
            grid = true,
            minorgrid = true,
            minorticks = 5,
            framestyle = :box,
            titlefontsize = 22,
            guidefontsize = 22,
            tickfontsize = 22,
            left_margin = 7Plots.mm,
            right_margin = 5Plots.mm,
            top_margin = 4Plots.mm,
            bottom_margin = 8Plots.mm,
        )

        for (method_index, method) in enumerate(METHODS)
            iterations, means = method_series(
                time_sums,
                observation_counts,
                instance_counts,
                horizon,
                budget,
                method,
                presolve,
            )
            isempty(iterations) && continue
            panel_max_iteration = max(panel_max_iteration, maximum(iterations))
            plot!(
                panel,
                iterations,
                means;
                label = METHOD_LABELS[method],
                color = METHOD_COLORS[method_index],
                linestyle = :solid,
                linewidth = 2.3,
            )
        end
        plot!(panel; xticks = integer_ticks(panel_max_iteration))
        push!(panels, panel)
    end

    columns = min(2, length(panels))
    rows = cld(length(panels), columns)

    # Fill a potentially incomplete final row so the legend always occupies
    # its own full-width row underneath the complete panel grid.
    while length(panels) < rows * columns
        push!(panels, plot(; framestyle = :none, axis = false, legend = false))
    end
    legend_panel = plot(
        ;
        framestyle = :none,
        axis = false,
        grid = false,
        legend = :bottom,
        legend_columns = -1,
        legendfontsize = 22,
        foreground_color_legend = :transparent,
        background_color_legend = :transparent,
    )
    for (method_index, method) in enumerate(METHODS)
        plot!(
            legend_panel,
            [NaN],
            [NaN];
            label = METHOD_LABELS[method],
            color = METHOD_COLORS[method_index],
            linestyle = :solid,
            linewidth = 2.3,
        )
    end

    panel_layout = @layout [grid(rows, columns); legend_area{0.08h}]
    return plot(
        panels...,
        legend_panel;
        layout = panel_layout,
        size = (1500, 450 * rows + 140),
        dpi = 300,
    )
end

function save_panel(panel, output_directory, presolve)
    mkpath(output_directory)
    presolve_text = presolve ? "on" : "off"
    base_path = joinpath(
        output_directory,
        "rostering_inner_master_time_presolve_$presolve_text",
    )
    savefig(panel, base_path * ".pdf")
    println("Wrote $(base_path).pdf")
end


function main(arguments = ARGS)
    results_directory, output_directory, selected_time_limit =
        parse_arguments(arguments)
    time_sums, observation_counts, instance_counts, files_used =
        load_inner_time_averages(
            results_directory,
            selected_time_limit,
        )

    gr()
    save_panel(
        create_panel(time_sums, observation_counts, instance_counts, true),
        output_directory,
        true,
    )
    save_panel(
        create_panel(time_sums, observation_counts, instance_counts, false),
        output_directory,
        false,
    )
    println("Used $files_used result files with a $(selected_time_limit)-second limit.")
end

main()
