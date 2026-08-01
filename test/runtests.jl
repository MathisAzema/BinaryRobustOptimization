using Test
using JuMP
using BinaryRobustOptimization

function tiny_rostering_instance()
    R = Rostering(2, 1, 1, 1)
    R.name = "tiny"
    R.I = 2
    R.J = 1
    R.T = 4
    R.N = 2
    R.FixedCostRegular = [1.0 1.0 1.0 1.0; 1.5 1.5 1.5 1.5]
    R.FixedCostPartTime = reshape([2.0, 3.0, 2.0, 3.0], 1, :)
    R.HourlyCostPartTime = reshape([1.0, 1.5, 1.0, 1.5], 1, :)
    R.MinShiftsRegular = [0, 0]
    R.MaxShiftsRegular = [4, 4]
    R.MinShiftsPartTime = [0]
    R.MaxShiftsPartTime = [2]
    R.PenaltyCost = [9.0, 10.0, 11.0, 12.0]
    R.DemandAvg = [2.0, 3.0, 2.0, 3.0]
    R.DemandDev = [1.0, 1.0, 2.0, 2.0]
    R.budget = 2
    empty!(R.β)
    empty!(R.γ)
    return R
end

function binary_vectors(n; budget = n)
    return [
        Float64.([(mask >> (t - 1)) & 1 for t in 1:n])
        for mask in 0:(2^n - 1)
        if count_ones(mask) <= budget
    ]
end

function feasible_part_time_schedules(R)
    schedules = Matrix{Float64}[]
    for values in binary_vectors(R.T; budget = R.MaxShiftsPartTime[1])
        all(values[t] + values[t + 1] <= 1 for t in 1:(R.T - 1)) || continue
        push!(schedules, reshape(values, 1, :))
    end
    return schedules
end

function direct_recourse_value(R, x, ξ)
    model = BinaryRobustOptimization.initializeJuMPModel()
    @variable(model, y[1:R.J, 1:R.T], Bin)
    @variable(model, z[1:R.J, 1:R.T] >= 0)
    @variable(model, w[1:R.T] >= 0)
    @constraint(model, [j in 1:R.J, t in 1:(R.T - 1)], y[j,t] + y[j,t + 1] <= 1)
    @constraint(model, [j in 1:R.J],
        R.MinShiftsPartTime[j] <= sum(y[j,t] for t in 1:R.T) <= R.MaxShiftsPartTime[j]
    )
    @constraint(model, [j in 1:R.J, t in 1:R.T], z[j,t] <= R.N * y[j,t])
    @constraint(model, [t in 1:R.T],
        sum(R.N * x[i,t] for i in 1:R.I) + sum(z[j,t] for j in 1:R.J) + w[t] >=
        R.DemandAvg[t] + R.DemandDev[t] * ξ[t]
    )
    @objective(model, Min,
        sum(R.FixedCostRegular[i,t] * x[i,t] for i in 1:R.I, t in 1:R.T) +
        sum(R.FixedCostPartTime[j,t] * y[j,t] for j in 1:R.J, t in 1:R.T) +
        sum(R.HourlyCostPartTime[j,t] * z[j,t] for j in 1:R.J, t in 1:R.T) +
        sum(R.PenaltyCost[t] * w[t] for t in 1:R.T)
    )
    optimize!(model)
    @test termination_status(model) == JuMP.MOI.OPTIMAL
    return objective_value(model)
end

function exhaustive_formulation_value(R, x, method)
    inner = BinaryRobustOptimization.init_master_inner_level(R, method)
    schedules = feasible_part_time_schedules(R)
    if method in (HatPdPrimeM, HatPdDoublePrimeM)
        for y in schedules, auxiliary_copy in binary_vectors(R.T)
            BinaryRobustOptimization.update_master_inner_level(
                R,
                inner,
                x,
                y,
                method;
                auxiliary_copy = auxiliary_copy,
            )
        end
    else
        for y in schedules
            BinaryRobustOptimization.update_master_inner_level(R, inner, x, y, method)
        end
    end
    optimize!(inner)
    @test termination_status(inner) == JuMP.MOI.OPTIMAL
    return objective_value(inner)
end

@testset "Five rostering reformulations" begin
    set_num_threads(1)
    set_solver_Gurobi()
    R = tiny_rostering_instance()
    x = zeros(R.I, R.T)

    expected = maximum(
        direct_recourse_value(R, x, ξ)
        for ξ in binary_vectors(R.T; budget = R.budget)
    )

    formulation_values = Dict(
        method => exhaustive_formulation_value(R, x, method)
        for method in five_lagrangian_formulations()
    )

    @test all(
        isapprox(value, expected; atol = 1e-6)
        for value in Base.values(formulation_values)
    )
    @test formulation_values[PdM] ≈ formulation_values[PdPrimeM] atol = 1e-6
end

@testset "MP_inner presolve switch" begin
    set_solver_Gurobi()
    R = tiny_rostering_instance()
    outer = BinaryRobustOptimization.init_master(R)
    inner = BinaryRobustOptimization.init_master_inner_level(R, PdM)

    BinaryRobustOptimization.set_optimizer_presolve(inner, false)
    @test JuMP.get_optimizer_attribute(inner, "Presolve") == 0
    @test JuMP.get_optimizer_attribute(outer, "Presolve") == 1

    BinaryRobustOptimization.set_optimizer_presolve(inner, true)
    @test JuMP.get_optimizer_attribute(inner, "Presolve") == 1
end
