
function gap(ub, lb)
    return isinf(ub) ? Inf : ((ub - lb) / ub)
end

function print_progress(iter, lb, ub, elapsed_time, λ=nothing, inner=false)
    elapsed_time = round(elapsed_time, digits=2)
    str = inner ? "\tinner-iter" : "iter"
    LB=round(lb, digits=2)
    UB=round(ub, digits=2)
    if λ === nothing
        println("$(str) $(iter): LB = $(LB) UB = $(UB) gap = $(round(gap(ub, lb) * 100.0, digits=2))% time=$(elapsed_time)sec")
    else
        println("$(str) $(iter): LB = $(LB) UB = $(UB) gap = $(round(gap(ub, lb) * 100.0, digits=2))% time=$(elapsed_time)sec λ = $(λ)")
    end
end

function set_optimizer_time_limit(m::JuMP.Model, time_limit::Float64)
    if SOLVER == "Mosek"
        JuMP.set_optimizer_attribute(m, "MSK_DPAR_OPTIMIZER_MAX_TIME", max(time_limit, 0.01))
    end
    if SOLVER == "Gurobi"
        JuMP.set_optimizer_attribute(m, "TimeLimit", max(time_limit, 0.01))
    end
    if SOLVER == "CPLEX"
        JuMP.set_optimizer_attribute(m, "CPXPARAM_TimeLimit", max(time_limit, 0.01))
    end
    if SOLVER == "SCIP"
        JuMP.set_optimizer_attribute(m, "limits/time", max(time_limit, 0.01))
    end
end

function initializeJuMPModel()
    if SOLVER == "Mosek"
        return Model(optimizer_with_attributes(
            Mosek.Optimizer,
            "MSK_IPAR_LOG" => 0,
            "MSK_DPAR_MIO_TOL_REL_GAP" => 0,
            "MSK_IPAR_NUM_THREADS" => THREADLIM
        ))
    end
    if SOLVER == "Gurobi"
        m = Model(optimizer_with_attributes(
            () -> Gurobi.Optimizer(GUROBI_ENV),
            "OutputFlag" => 0,
            "MIPGap" => 0.001,
            "Threads" => THREADLIM,
            "Presolve" => 1, #Presolve deactivated : 0
        ))
        JuMP.set_silent(m)
        # JuMP.unset_silent(m)
        return m
    end
    if SOLVER == "CPLEX"
        return Model(optimizer_with_attributes(
            CPLEX.Optimizer,
            "CPXPARAM_ScreenOutput" => 0,
            "CPXPARAM_MIP_Tolerances_MIPGap" => 0.001,
            "CPXPARAM_Threads" => THREADLIM
        ))
    end
    if SOLVER == "SCIP"
        return Model(optimizer_with_attributes(
            SCIP.Optimizer,
            "display/verblevel" => 0,
            "limits/gap" => 0,
            "parallel/maxnthreads" => THREADLIM
        ))
    end
end

function solve_MP(problem::AbstractProblem, MP::JuMP.Model, time_limit::Float64)
    set_optimizer_time_limit(MP, time_limit)
    @timeit "Optimize Master" optimize!(MP)
    status = termination_status(MP)
    if complete_recourse(problem)
        if status != MOI.OPTIMAL
            # unset_silent(MP)
            optimize!(MP)
            @warn("could not solve master problem to optimality. status = $(status)")
            return NaN
        end
    else
        if status ∉ [MOI.OPTIMAL, MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
            unset_silent(MP)
            optimize!(MP)
            @warn("could not solve master problem to optimality or infeasibility. status = $(status)")
            return NaN
        end
    end
    if status ∈ [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
        @warn("potentially infeasible problem instance")
        return +Inf
    end
    return objective_scale(problem) * objective_value(MP)
end

function solve_SP(problem::AbstractProblem, SP::JuMP.Model, time_limit::Float64)
    set_optimizer_time_limit(SP, time_limit)
    @timeit "Optimize Subproblem" optimize!(SP)
    status = termination_status(SP)
    if complete_recourse(problem)
        if status != MOI.OPTIMAL
            @warn("could not solve subproblem to optimality. status = $(status)")
            return NaN
        end
    else
        if status ∉ [MOI.OPTIMAL, MOI.SOLUTION_LIMIT, MOI.OBJECTIVE_LIMIT]
            @warn("could not solve subproblem to optimality, solution_limit, or objective_limit. status = $(status)")
            return NaN
        end
    end
    if status == MOI.OBJECTIVE_LIMIT
        return -Inf
    end
    return objective_scale(problem) * objective_value(SP)
end

function debug_repeated_scenarios(scenario_list::Dict, masterproblemtype::MasterType, problem::AbstractProblem, SP::JuMP.Model, LB::Float64, UB::Float64, found_infeasible_scenario::Bool)
    in_list = record_scenario(problem, SP, scenario_list)
    if in_list
        return true
    else
        return false
    end
    # if in_list && masterproblemtype == CCG
    #     found_infeasible_scenario && @warn("repeated scenario - claimed infeasible", objective_value(SP))
    #     found_infeasible_scenario || @warn("repeated scenario - claimed better", gap(UB, LB))
    #     @assert false
    # end
end

function compute_time_worstcase(res)
    N = maximum([length(res[i]) for i in 1:length(res)])
    Time_worstcase = [0.0 for j in 1:N]
    Nb = [0 for j in 1:N]
    for i in 1:length(res)
        for j in 1:length(res[i])
            Time_worstcase[j] += res[i][j]
            Nb[j] += 1
        end
    end
    for j in 1:N
        Time_worstcase[j] /= Nb[j]
    end
    return Time_worstcase
end

function generate_all_Γ_tuple(n::Int64, Γ::Int64)
    D_set=[]
    function generate_Γ_tuple(n::Int64, Γ::Int64, current_tuple::Vector{Any}=[])
        if length(current_tuple) == Γ
            push!(D_set, current_tuple)
        else
            start_val = isempty(current_tuple) ? 1 : current_tuple[end] + 1
            for i in start_val:n
                generate_Γ_tuple(n, Γ, vcat(current_tuple, i))
            end
        end
    end
    generate_Γ_tuple(n, Γ)
    uncertainty_set = [zeros(Int64, n) for _ in 1:length(D_set)]
    for (idx, tuple) in enumerate(D_set)
        for t in tuple
            uncertainty_set[idx][t] = 1
        end
    end
    return uncertainty_set
end
