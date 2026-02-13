"""
    run_ccg(problem::AbstractProblem, subproblemtype::SubproblemType, time_limit::Float64, opt_tol::Float64 = 1e-4, feas_tol::Float64 = 1e-5)

Solve `problem` using the column-and-constraint generation algorithm using the
`subproblemtype` subproblem with specified `time_limit` (in seconds) and
optimality and feasibility tolerances of `opt_tol` and `feas_tol`, respectively.

Returns `(num_iter, lb, ub, total_time, iter_inner)`
"""
function run_ccg(
    problem::AbstractProblem,
    subproblemtype::SubproblemType,
    time_limit::Float64,
    opt_tol::Float64=1e-3,
    feas_tol::Float64=1e-5,
)
    if mixed_integer_recourse(problem)
        if !complete_recourse(problem)
            @error("to do - CCG algorithm for mixed-integer instances without complete recourse")
            return nothing
        end
        if subproblemtype == LagrangianDual && indicator_uncertainty(problem)
            @error("to do - CCG algorithm for mixed-integer instances with indicator uncertainties")
            return nothing
        end
        return run_ccg_mixed_integer_recourse(problem, subproblemtype, time_limit, opt_tol, feas_tol)
    end
    return run_iterative_continuous_recourse(problem, CCG, subproblemtype, time_limit, opt_tol, feas_tol)
end

function run_ccg_mixed_integer_recourse(
    problem::AbstractProblem,
    subproblemtype::SubproblemType,
    time_limit::Float64,
    opt_tol::Float64,
    feas_tol::Float64,
)
    @assert complete_recourse(problem)
    @assert subproblemtype != LagrangianDual || !indicator_uncertainty(problem)

    masterproblemtype = CCG
    algname = "$(subproblemtype)-$(masterproblemtype)"
    LB = -Inf
    UB = +Inf
    iter = 0
    start_t = time()
    iter_inner = Vector{Int}()

    MP_outer = init_master(problem)
    scenario_list = Dict()

    λ = nothing
    if subproblemtype == LagrangianDual
        λ = 1.0
    end

    previous_scenario = zeros(21)

    Time_MP_inner = Vector{Vector{Float64}}()

    @timeit algname begin
        while time() - start_t <= time_limit && iter <= 10
            iter += 1
            push!(iter_inner, 0)

            lb = solve_MP(problem, MP_outer, time_limit - (time() - start_t))
            !isnan(lb) && (LB = lb) # normal termination (optimal or infeasible)
            !isfinite(lb) && break  # infeasible or non-normal termination

            @timeit "$algname-InnerLevel" begin
                LB_inner = -Inf
                UB_inner = +Inf
                MP_inner = init_master_inner_level(problem, subproblemtype)
                MP_inner_opt = nothing
                discrete_decision_list = Dict()

                if subproblemtype == LagrangianDual
                    @timeit "compute_lagrangian_parameter" begin
                        λ = max(λ, compute_lagrangian_coefficient(problem, MP_outer))
                    end
                end

                normal_termination = true

                # scenario_list_inner = Dict()

                Time_MP_inner_iter = Float64[]
                if subproblemtype == Enumeration
                    start_inner = time()
                    ub = solve_MP_inner_enumeration(problem, MP_outer, MP_inner)
                    UB = min(UB, ub)
                    MP_inner_opt = init_master_inner_level(problem, MP_inner)
                    push!(Time_MP_inner, [time() - start_inner])
                else
                    while time() - start_t <= time_limit #&& iter_inner[end]<=60
                        iter_inner[end] += 1
                        # if subproblemtype == LagrangianDualbis
                        #     solve_MP_FW(problem, MP_inner, previous_scenario, scenario_list_inner)
                        # end 
                        start = time()
                        ub = solve_MP(problem, MP_inner, time_limit - (time() - start_t))
                        push!(Time_MP_inner_iter, time()-start)

                        if isnan(ub) # non-normal termination
                            normal_termination = false
                            break
                        end
                        UB_inner = ub
                        UB = min(UB, UB_inner)

                        SP = build_second_stage_problem(problem, MP_outer, MP_inner)
                        lb = solve_SP(problem, SP, time_limit - (time() - start_t))
                        if !isfinite(lb) # non-normal termination
                            normal_termination = false
                            break
                        end
                        if lb > LB_inner
                            # previous_scenario = JuMP.value.(MP_inner[:g])
                            LB_inner = lb
                            MP_inner_opt = init_master_inner_level(problem, MP_inner)
                            # if gap(LB_inner, LB) > 10.0*opt_tol
                            #     break
                            # end
                        end
                        if subproblemtype == LagrangianDual
                            record_discrete_second_stage_decision(problem, SP, discrete_decision_list)
                            while true
                                @timeit "solve_second_stage_problem_lagrangian" begin
                                    step = solve_second_stage_problem_lagrangian(problem, MP_outer, MP_inner, λ)
                                end
                                if step <= feas_tol
                                    break
                                end
                                λ *= 2.0
                            end
                        end

                        # Print progress
                        print_progress(iter_inner[end], LB_inner, UB_inner, time() - start_t, λ, true)

                        # Terminate or update master
                        if gap(UB_inner, LB_inner) > opt_tol
                            update_master_inner_level(problem, MP_outer, MP_inner, SP, subproblemtype, λ)
                        else
                            break
                        end
                    end
                end
                push!(Time_MP_inner, Time_MP_inner_iter)
                normal_termination || break
            end

            # Print progress
            print_progress(iter, LB, UB, time() - start_t, λ)

            # Terminate or update master
            if gap(UB, LB) > opt_tol
                if MP_inner_opt == nothing
                    computational_time = round(time() - start_t)
                    return return_solution(problem, computational_time, LB, UB, Time_MP_inner, subproblemtype)
                end
                # if debug_repeated_scenarios(scenario_list, masterproblemtype, problem, MP_inner_opt, LB, UB, false)
                #     return string(subproblemtype), problem.T, problem.budget, time() - start_t, 100*(UB-LB)/UB, sum([Time_MP_inner[i][t] for i in 1:length(Time_MP_inner) for t in 1:length(Time_MP_inner[i])]), compute_time_worstcase(Time_MP_inner), LB
                # end
                update_master_mixed_integer(problem, MP_outer, MP_inner_opt, masterproblemtype)
            else
                break
            end
        end
    end
    # return Time_MP_inner
    # println(UB, LB)
    computational_time = round(time() - start_t)
    return return_solution(problem, computational_time, LB, UB, Time_MP_inner, subproblemtype)

    nb_worst_case = sum([1 for i in 1:length(Time_MP_inner) for t in 1:length(Time_MP_inner[i])]; init = 0)
    # return string(subproblemtype), problem.T, problem.budget, time() - start_t, 100*(UB-LB)/UB, sum([Time_MP_inner[i][t] for i in 1:length(Time_MP_inner) for t in 1:length(Time_MP_inner[i])]; init = 0), compute_time_worstcase(Time_MP_inner), LB, nb_worst_case
    # return iter, LB, UB, time() - start_t, iter_inner
end