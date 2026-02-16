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
        if subproblemtype == CCGM && indicator_uncertainty(problem)
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
    @assert subproblemtype != CCGM || !indicator_uncertainty(problem)

    masterproblemtype = CCG
    LB = -Inf
    UB = +Inf
    iter = 0
    start_t = time()
    iter_inner = Vector{Int}()
    scenario_list = Dict()

    MP_outer = init_master(problem)

    λ = nothing
    if subproblemtype == CCGM
        λ = 1.0
    end

    Time_MP_inner = Vector{Vector{Float64}}()

    while time() - start_t <= time_limit
        iter += 1
        push!(iter_inner, 0)

        # unset_silent(MP_outer)
        lb = solve_MP(problem, MP_outer, time_limit - (time() - start_t))
        if !isfinite(lb) 
            break  # infeasible or non-normal termination
        end
        LB = max(LB, lb)

        LB_inner = -Inf
        UB_inner = +Inf
        ξ = [0 for t in 1:problem.T]
        MP_inner = init_master_inner_level(problem, subproblemtype)
        discrete_decision_list = Dict()

        if subproblemtype == CCGM
            λ = max(λ, compute_lagrangian_coefficient(problem, MP_outer))
        end

        normal_termination = true

        Time_MP_inner_iter = Float64[]
        if subproblemtype == Enumeration
            start_inner = time()
            ub = solve_MP_inner_enumeration(problem, MP_outer, MP_inner, time_limit - (time() - start_t))
            push!(Time_MP_inner, [time() - start_inner])
            if isnan(ub) # non-normal termination
                normal_termination = false
            else 
                UB = min(UB, ub)
            end
        else
            while time() - start_t <= time_limit
                ub = -Inf
                DCheuristic = false
                if subproblemtype == CCGLDC
                    ub_fw = solve_MP_FW(problem, MP_inner, ξ)
                    println((ub_fw, LB_inner, UB_inner))
                    if ub_fw <= 1.001* LB_inner
                        println("HHH")
                        JuMP.unfix.(MP_inner[:ξ])
                        start = time()
                        ub = solve_MP(problem, MP_inner, time_limit - (time() - start_t))
                        push!(Time_MP_inner_iter, time()-start)

                        if isnan(ub) # non-normal termination
                            normal_termination = false
                            break
                        end
                        UB_inner = min(UB_inner, ub)
                    else
                        DCheuristic = true
                        ub = ub_fw
                    end
                else
                    start = time()
                    ub = solve_MP(problem, MP_inner, time_limit - (time() - start_t))
                    push!(Time_MP_inner_iter, time()-start)

                    if isnan(ub) # non-normal termination
                        normal_termination = false
                        break
                    end
                    UB_inner = min(UB_inner, ub)
                end
                iter_inner[end] += 1

                SP = build_second_stage_problem(problem, MP_outer, MP_inner, subproblemtype, λ)
                lb = solve_SP(problem, SP, time_limit - (time() - start_t))

                if !isfinite(lb) # non-normal termination
                    normal_termination = false
                    break
                end
                # LB_inner = max(lb, LB_inner)
                if lb > LB_inner
                    LB_inner = lb
                    for t in 1:problem.T
                        ξ[t] = Float64.(Int.(round(value(MP_inner[:ξ][t]))))
                    end
                end

                if subproblemtype == CCGM
                    while true
                        step = compute_grad_lagrangian(problem, SP, MP_inner)
                        if step <= feas_tol
                            break
                        end
                        λ *= 2.0
                        SP = build_second_stage_problem(problem, MP_outer, MP_inner, subproblemtype, λ)
                        lb = solve_SP(problem, SP, time_limit - (time() - start_t))
                    end
                end

                # Print progress
                print_progress(iter_inner[end], LB_inner, UB_inner, time() - start_t, λ, true)

                if record_discrete_second_stage_decision(problem, SP, discrete_decision_list)
                    if DCheuristic
                        JuMP.unfix.(MP_inner[:ξ])
                        start = time()
                        ub = solve_MP(problem, MP_inner, time_limit - (time() - start_t))
                        push!(Time_MP_inner_iter, time()-start)

                        if isnan(ub) # non-normal termination
                            normal_termination = false
                            break
                        end
                        UB_inner = min(UB_inner, ub)

                        SP = build_second_stage_problem(problem, MP_outer, MP_inner, subproblemtype, λ)
                        lb = solve_SP(problem, SP, time_limit - (time() - start_t))

                        if !isfinite(lb) # non-normal termination
                            normal_termination = false
                            break
                        end
                        if lb > LB_inner
                            LB_inner = lb
                            for t in 1:problem.T
                                ξ[t] = Float64.(Int.(round(value(MP_inner[:ξ][t]))))
                            end
                        end
                    else
                        break
                    end
                end

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
        UB = min(UB, UB_inner)

        # Print progress
        print_progress(iter, LB, UB, time() - start_t, λ)

        if record_scenario(problem, ξ, scenario_list)
            break
        end

        # Terminate or update master
        if gap(UB, LB) > opt_tol
            update_master_mixed_integer(problem, MP_outer, ξ, masterproblemtype)
        else
            break
        end
    end

    computational_time = round(time() - start_t)
    return return_solution(problem, computational_time, LB, UB, Time_MP_inner, subproblemtype)
end