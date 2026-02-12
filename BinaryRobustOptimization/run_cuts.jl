# using Revise; using Distributed
# t = @async addprocs(2)
# Distributed.nworkers()
# include("run_cuts.jl")

@everywhere import Pkg
@everywhere Pkg.activate(".")
@everywhere include("package.jl")
# run_rostering(10, [1], [1], [LagrangianDual], 30.0)
# run_rostering(10, [1,2], [3,6,9], [LagrangianDual, LagrangianDualbis], 3600.0)


function run_rostering(N::Int, scale_list::Vector{Int}, budget_list::Vector{Int}, method_list, end_time)
    futures_results = Vector{Future}()
    for i in 1:N
        for scale in scale_list
            for budget in budget_list
                for method in method_list
                    problem = Rostering(budget, scale, scale, i)
                    println("Instance ", i, " Scale ", scale, " Budget ", budget, " Method ", method)
                    push!(futures_results, @spawn RobustOptLagrangianDual.run_ccg(problem, method, end_time))
                end
            end
        end
    end
    results = map(fetch, futures_results)
    println("Writing excel files")
    XLSX.openxlsx("Rostering.xlsx", mode="w") do xf
        sheet = xf[1]
        XLSX.rename!(sheet, "results")
        sheet["A1"]="method"
        sheet["B1"]="objective"
        sheet["C1"]="scale"
        sheet["D1"]="budget"
        sheet["E1"]="Time"
        sheet["F1"]="gap"
        sheet["G1"]="Time tot worstcase"
        sheet["H1"]="Time per iteration"
        for i in 1:length(results)
            if results[i]!=nothing
                sheet["A"*string(i+1)]= results[i][1]
                sheet["B"*string(i+1)]= results[i][8]
                sheet["C"*string(i+1)]= results[i][2]
                sheet["D"*string(i+1)]= results[i][3]
                sheet["E"*string(i+1)]= results[i][4]
                sheet["F"*string(i+1)]= results[i][5]
                sheet["G"*string(i+1)]= results[i][6]
                sheet["H"*string(i+1)]= string(round.(results[i][7], digits=2))
            end
        end
    end
end

function run_unit_commitment(budget_list::Vector{Int}, method_list, end_time)
    futures_results = Vector{Future}()
    for budget in budget_list
        for method in method_list
            # problem1 = UnitCommitment("6bus_JEAS", budget, 1)
            # push!(futures_results, @spawn BinaryRobustOptimization.run_ccg(problem1, method, end_time))
            problem2 = UnitCommitment("medium", budget, 5)
            push!(futures_results, @spawn BinaryRobustOptimization.run_ccg(problem2, method, end_time))
            # problem3 = UnitCommitment("118_syst_JEAS", budget, 35)
            # push!(futures_results, @spawn BinaryRobustOptimization.run_ccg(problem3, method, end_time))
        end
    end
    results = map(fetch, futures_results)
    println("Writing excel files")
    N = budget_list[end]
    XLSX.openxlsx("UnitCommitment_$N.xlsx", mode="w") do xf
        sheet = xf[1]
        XLSX.rename!(sheet, "results")
        sheet["A1"]="method"
        sheet["B1"]="objective"
        sheet["C1"]="scale"
        sheet["D1"]="budget"
        sheet["E1"]="Time"
        sheet["F1"]="gap"
        sheet["G1"]="Time tot worstcase"
        sheet["H1"]="Time per iteration"
        for i in 1:length(results)
            if results[i]!=nothing
                sheet["A"*string(i+1)]= results[i][1]
                sheet["B"*string(i+1)]= results[i][8]
                sheet["C"*string(i+1)]= results[i][2]
                sheet["D"*string(i+1)]= results[i][3]
                sheet["E"*string(i+1)]= results[i][4]
                sheet["F"*string(i+1)]= results[i][5]
                sheet["G"*string(i+1)]= results[i][6]
                sheet["H"*string(i+1)]= string(round.(results[i][7], digits=2))
            end
        end
    end
end