"""
    UnitCommitment <: AbstractProblem

Default implementation of `UnitCommitment` problem.
"""

struct ThermalUnit
    name::Int64
    Bus::Int64
    MinPower::Float64  
    MaxPower::Float64  
    DeltaRampUp::Float64  
    DeltaRampDown::Float64 
    QuadTerm::Float64  
    StartUpCost::Float64
    StartDownCost::Float64  
    LinearTerm::Float64
    ConstTerm::Float64 
    InitialPower::Float64
    InitUpDownTime::Int64
    MinUpTime::Int64  
    MinDownTime::Int64
end

struct Line
    id::Int64
    b1::Int64
    b2::Int64 
    Fmax::Float64  
    B12::Float64 
end

struct UnitCommitment <: AbstractProblem
    name::String
    T::Int64
    N::Int64
    Thermalunits::Vector{ThermalUnit}
    Buses::Int64
    Lines::Vector{Line}
    Demandbus::Vector{Vector{Float64}}
    BusWind::Vector{Int64}
    DemandDev::Vector{Vector{Float64}}
    PenaltyCost::Float64
    budget::Int64

    # function UnitCommitment(filename::String, budget::Int)
    #     data=Dataset(filename)
    #     Blocks=keys(data.group)
    #     data_block=data.group[Blocks[1]]
    #     TimeHorizon= 24
    #     demand=round.(data_block["ActivePowerDemand"].var[:])
    #     N=size(keys(data_block.group))[1]-1
    #     Thermal_units_v1=Vector{ThermalUnit}(undef, N)
    #     Thermal_units_v2=Vector{ThermalUnit}(undef, N)
    #     k=0
    #     for unit in data_block.group
    #         if k < N
    #             if occursin("Block", first(unit))
    #                 k += 1
    #                 unit_name=parse(Int64, split(first(unit), "_")[end])+1
    #                 type=last(unit).attrib["type"]
    #                 if type == "ThermalUnitBlock" 
    #                     Bus = 1
    #                     MinPower=round.(last(unit)["MinPower"].var[1])
    #                     MaxPower =round.(last(unit)["MaxPower"].var[1]) 
    #                     DeltaRampUp  =round.(last(unit)["DeltaRampUp"].var[1])
    #                     DeltaRampDown  =round.(last(unit)["DeltaRampDown"].var[1])
    #                     QuadTerm =round(last(unit)["QuadTerm"].var[1]) 
    #                     StartUpCost=round(last(unit)["StartUpCost"].var[1])
    #                     StartDownCost=0.0 
    #                     LinearTerm=round(last(unit)["LinearTerm"].var[1])
    #                     ConstTerm=round(last(unit)["ConstTerm"].var[1])
    #                     InitialPower=round.(last(unit)["InitialPower"].var[1])  
    #                     InitUpDownTime =round.(last(unit)["InitUpDownTime"].var[1]) 
    #                     MinUpTime=round.(last(unit)["MinUpTime"].var[1])  
    #                     MinDownTime=round.(last(unit)["MinDownTime"].var[1])

    #                     unit=ThermalUnit(unit_name, Bus, MinPower, MaxPower, DeltaRampUp, DeltaRampDown, QuadTerm, StartUpCost, StartDownCost, LinearTerm, ConstTerm, InitialPower, InitUpDownTime, Int64(MinUpTime), Int64(MinDownTime))
    #                     Thermal_units_v1[unit_name]=unit
    #                 end
    #             end
    #         end
    #     end
    #     H = [Thermal_units_v1[i].LinearTerm for i in 1:N]
    #     idx_sort = sortperm(H; rev=false)
    #     for (k,i) in enumerate(idx_sort)
    #         unit = Thermal_units_v1[i]
    #         unit2=ThermalUnit(k, unit.Bus, unit.MinPower, unit.MaxPower, unit.DeltaRampUp, unit.DeltaRampDown, unit.QuadTerm, unit.StartUpCost, unit.StartDownCost, unit.LinearTerm, unit.ConstTerm, unit.InitialPower, unit.InitUpDownTime, Int64(unit.MinUpTime), Int64(unit.MinDownTime))
    #         Thermal_units_v2[k]=unit2
    #     end
    #     Lines = Line[]
    #     Demandbus=[demand]
    #     DemandDev = [[demand[t]*1.96*0.025 for t in 1:TimeHorizon]]
    #     BusWind=[1]
    #     PenaltyCost = 300.0
    #     Buses = 1

    #     new(filename, 
    #         TimeHorizon, 
    #         N, 
    #         Thermal_units_v2, 
    #         Buses,
    #         Lines, 
    #         Demandbus, 
    #         BusWind, 
    #         DemandDev, 
    #         PenaltyCost,
    #         budget)
    # end

    function UnitCommitment(folder::String, budget::Int)
        """
        Parse the IEEE 118-bus instnace
        """
        TimeHorizon= 24
        syst = "data/UC/"*folder
        generators = CSV.read(joinpath(pwd(), syst, "generators.csv"), DataFrame; header=false)
        NumberUnits= parse(Int64,generators[end,1]) 
        N=NumberUnits
        name_instance="IEEE"*string(NumberUnits)
        Thermal_units=Vector{ThermalUnit}(undef, N)
        for i in 2:NumberUnits+1
            unit_name=i-1
            Bus = parse(Int64,generators[i,2])
            ConstTerm = parse(Float64,generators[i,3])
            LinearTerm = parse(Float64,generators[i,4])
            MaxPower = parse(Float64,generators[i,6])
            MinPower = parse(Float64,generators[i,7])
            DeltaRampUp = parse(Float64,generators[i,14])
            DeltaRampDown = parse(Float64,generators[i,14])
            StartUpCost = parse(Float64,generators[i,15])
            StartDownCost = 0.0*parse(Float64,generators[i,15])
            MinUpTime=parse(Int64,generators[i,13]) 
            MinDownTime=parse(Int64,generators[i,12])
            QuadTerm =0.0        
            InitialPower=parse(Float64,generators[i,11])
            InitUpDownTime =parse(Int64,generators[i,10])
            unit=ThermalUnit(unit_name, Bus, MinPower, MaxPower, DeltaRampUp, DeltaRampDown, QuadTerm, StartUpCost, StartDownCost, LinearTerm, ConstTerm, InitialPower, InitUpDownTime, MinUpTime, MinDownTime)
            Thermal_units[unit_name]=unit
        end

        maximum_load = CSV.read(joinpath(pwd(), syst, "maximum_load.csv"), DataFrame; header=false)
        Numbus=maximum_load[end,1]
        Buses=1:Numbus
        load_distribution_profile = CSV.read(joinpath(pwd(), syst, "load_distribution_profile.csv"), DataFrame; header=false)[:,2]/100
        Demandbus=[maximum_load[b,2]*load_distribution_profile for b in Buses]
        df_lines = CSV.read(joinpath(pwd(), syst, "lines.csv"), DataFrame; header=false)
        Numlines=parse(Int64,df_lines[end,1])
        Lines=Vector{Line}(undef, Numlines)
        for i in 2:Numlines+1
            id = parse(Int64, df_lines[i,1])
            b1 = parse(Int64, df_lines[i,2])
            b2 = parse(Int64, df_lines[i,3])
            fmax = parse(Float64, df_lines[i,6])
            X = 1/parse(Float64, df_lines[i,5])
            Lines[i-1]=Line(id, b1, b2, fmax, X)
        end

        BusWind=[b for b in Buses if sum(Demandbus[b])>=1]

        DemandDev = [[Demandbus[b][t]*1.96*0.025 for t in 1:TimeHorizon] for b in Buses]

        PenaltyCost = 300.0

        new(name_instance, 
            TimeHorizon, 
            N, 
            Thermal_units, 
            Buses[end],
            Lines, 
            Demandbus, 
            BusWind, 
            DemandDev, 
            PenaltyCost,
            budget)
    end
end

mixed_integer_recourse(UC::UnitCommitment) = false

complete_recourse(UC::UnitCommitment) = true

objective_scale(UC::UnitCommitment) = 1.0

indicator_uncertainty(UC::UnitCommitment) = false

function solve_deterministic_problem(UC::UnitCommitment)
    m = initializeJuMPModel()

    newgap=0.1/100
    set_optimizer_attribute(m, "MIPGap", newgap)

    @variable(m, is_on[i in 1:UC.N, t in 0:UC.T], Bin)
    @variable(m, start_up[i in 1:UC.N, t in 1:UC.T], Bin)
    @variable(m, start_down[i in 1:UC.N, t in 1:UC.T], Bin)

    @variable(m, thermal_cost>=0)
    @variable(m, thermal_fixed_cost>=0)

    @constraint(m, thermal_fixed_cost>=sum(unit.ConstTerm*is_on[unit.name, t]+unit.StartUpCost*start_up[unit.name, t]+unit.StartDownCost*start_down[unit.name, t] for unit in UC.Thermalunits for t in 1:UC.T))

    @objective(m, Min, thermal_cost + thermal_fixed_cost)

    @constraint(m,  [unit in UC.Thermalunits, t in 1:UC.T], is_on[unit.name, t]-is_on[unit.name, t-1]==start_up[unit.name, t]-start_down[unit.name, t])
    @constraint(m,  [unit in UC.Thermalunits, t in 1:UC.T], start_up[unit.name, t]<=is_on[unit.name, t])
    @constraint(m,  [unit in UC.Thermalunits, t in 1:UC.T], start_down[unit.name, t]<=1-is_on[unit.name, t])
    @constraint(m,  [unit in UC.Thermalunits, t in 1:UC.T], sum(start_up[unit.name, τ] for τ in max(1, t-unit.MinUpTime+1):t)+1*(t<unit.MinUpTime-unit.InitUpDownTime+1)*(unit.InitUpDownTime>0) <= is_on[unit.name, t])
    @constraint(m,  [unit in UC.Thermalunits, t in 1:UC.T], sum(start_down[unit.name, τ] for τ in max(1, t-unit.MinDownTime+1):t)+1*(t<unit.MinDownTime+unit.InitUpDownTime+1)*(unit.InitUpDownTime<0) <= 1-is_on[unit.name, t])
    @constraint(m,  [unit in UC.Thermalunits], is_on[unit.name, 0]==(unit.InitUpDownTime>=0))

    for unit in UC.Thermalunits
        if (unit.InitialPower-unit.MinPower)/unit.DeltaRampDown <0
            limit=-1
        else
            limit=Int64(ceil((unit.InitialPower-unit.MinPower)/unit.DeltaRampDown))
        end
        for t in 0:limit
            @constraint(m, is_on[unit.name, t]==1)
        end
    end

    @variable(m, power[i in 1:UC.N, t in 0:UC.T] >= 0)
    @variable(m, power_shedding[b in 1:UC.Buses, t in 1:UC.T] >= 0)
    @variable(m, power_curtailement[b in 1:UC.Buses, t in 1:UC.T] >= 0)

    @constraint(m,  [unit in UC.Thermalunits], power[unit.name, 0]==unit.InitialPower)
    @constraint(m,  [unit in UC.Thermalunits, t in 1:UC.T], power[unit.name, t]<=unit.MaxPower*is_on[unit.name, t])
    @constraint(m,  [unit in UC.Thermalunits, t in 1:UC.T], power[unit.name, t]>=unit.MinPower*is_on[unit.name, t])
    @constraint(m,  [unit in UC.Thermalunits, t in 1:UC.T], power[unit.name, t]-power[unit.name, t-1]<=-unit.DeltaRampUp*start_up[unit.name, t]-unit.MinPower*is_on[unit.name, t-1] + (unit.MinPower+unit.DeltaRampUp)*is_on[unit.name, t])
    @constraint(m,  [unit in UC.Thermalunits, t in 1:UC.T], power[unit.name, t-1]-power[unit.name, t]<=-unit.DeltaRampDown*start_down[unit.name, t] - unit.MinPower*is_on[unit.name, t] + (unit.MinPower+unit.DeltaRampDown)*is_on[unit.name, t-1])

    @constraint(m, thermal_cost >= sum(unit.LinearTerm*power[unit.name, t] for unit in UC.Thermalunits for t in 1:UC.T) + UC.PenaltyCost*sum(power_shedding[b,t] + power_curtailement[b,t] for b in 1:UC.Buses for t in 1:UC.T))

    @variable(m, flow[l in 1:length(UC.Lines), t in 1:UC.T])
    @variable(m, angle[b in 1:UC.Buses, t in 1:UC.T])
    @constraint(m, [line in UC.Lines, t in 1:UC.T], flow[line.id,t]<=line.Fmax)
    @constraint(m, [line in UC.Lines, t in 1:UC.T], flow[line.id,t]>=-line.Fmax)
    @constraint(m, [line in UC.Lines, t in 1:UC.T], flow[line.id,t] == line.B12*(angle[line.b1,t]-angle[line.b2,t]))

    @constraint(m, [t in 1:UC.T, b in 1:UC.Buses], sum(power[unit.name, t] for unit in UC.Thermalunits if unit.Bus==b) + power_shedding[b,t]- power_curtailement[b,t] + sum(flow[line.id, t] for line in UC.Lines if line.b2==b) - sum(flow[line.id, t] for line in UC.Lines if line.b1==b) == UC.Demandbus[b][t])

    optimize!(m)

    println((JuMP.value(thermal_fixed_cost), JuMP.value(thermal_cost)))

    return objective_value(m), m
end

function init_master(UC::UnitCommitment)
    m = initializeJuMPModel()

    @variable(m, is_on[i in 1:UC.N, t in 0:UC.T], Bin)
    @variable(m, start_up[i in 1:UC.N, t in 1:UC.T], Bin)
    @variable(m, start_down[i in 1:UC.N, t in 1:UC.T], Bin)

    @constraint(m,  [unit in UC.Thermalunits, t in 1:UC.T], is_on[unit.name, t]-is_on[unit.name, t-1]==start_up[unit.name, t]-start_down[unit.name, t])
    @constraint(m,  [unit in UC.Thermalunits, t in 1:UC.T], start_up[unit.name, t]<=is_on[unit.name, t])
    @constraint(m,  [unit in UC.Thermalunits, t in 1:UC.T], start_down[unit.name, t]<=1-is_on[unit.name, t])
    @constraint(m,  [unit in UC.Thermalunits, t in 1:UC.T], sum(start_up[unit.name, τ] for τ in max(1, t-unit.MinUpTime+1):t)+1*(t<unit.MinUpTime-unit.InitUpDownTime+1)*(unit.InitUpDownTime>0) <= is_on[unit.name, t])
    @constraint(m,  [unit in UC.Thermalunits, t in 1:UC.T], sum(start_down[unit.name, τ] for τ in max(1, t-unit.MinDownTime+1):t)+1*(t<unit.MinDownTime+unit.InitUpDownTime+1)*(unit.InitUpDownTime<0) <= 1-is_on[unit.name, t])
    @constraint(m,  [unit in UC.Thermalunits], is_on[unit.name, 0]==(unit.InitUpDownTime>=0))

    for unit in UC.Thermalunits
        if (unit.InitialPower-unit.MinPower)/unit.DeltaRampDown <0
            limit=-1
        else
            limit=Int64(ceil((unit.InitialPower-unit.MinPower)/unit.DeltaRampDown))
        end
        for t in 0:limit
            @constraint(m, is_on[unit.name, t]==1)
        end
    end

    @variable(m, thermal_fixed_cost>=0)

    @constraint(m, thermal_fixed_cost>=sum(unit.ConstTerm*is_on[unit.name, t]+unit.StartUpCost*start_up[unit.name, t]+unit.StartDownCost*start_down[unit.name, t] for unit in UC.Thermalunits for t in 1:UC.T))

    @variable(m, s >= 0)

    @objective(m, Min, s + thermal_fixed_cost)

    return m
end

function update_master_continuous(UC::UnitCommitment, MP::JuMP.Model, SP::JuMP.Model, master::MasterType, subproblem::SubproblemType)
    s = MP[:s]
    is_on = MP[:is_on]
    start_up = MP[:start_up]
    start_down = MP[:start_down]

    if master == CCG
        ξ = JuMP.value.(SP[:ξ])
        ξ = Float64.(Int.(round.(ξ))) # should be 0-1

        power = @variable(MP, [i in 1:UC.N, t in 0:UC.T])
        power_shedding = @variable(MP, [b in 1:UC.Buses, t in 1:UC.T], lower_bound = 0)
        power_curtailement = @variable(MP, [b in 1:UC.Buses, t in 1:UC.T], lower_bound = 0)

        @constraint(MP,  [unit in UC.Thermalunits], power[unit.name, 0]==unit.InitialPower)
        @constraint(MP,  [unit in UC.Thermalunits, t in 1:UC.T], power[unit.name, t]<=unit.MaxPower*is_on[unit.name, t])
        @constraint(MP,  [unit in UC.Thermalunits, t in 1:UC.T], power[unit.name, t]>=unit.MinPower*is_on[unit.name, t])
        @constraint(MP,  [unit in UC.Thermalunits, t in 1:UC.T], power[unit.name, t]-power[unit.name, t-1]<=-unit.DeltaRampUp*start_up[unit.name, t]-unit.MinPower*is_on[unit.name, t-1] + (unit.MinPower+unit.DeltaRampUp)*is_on[unit.name, t])
        @constraint(MP,  [unit in UC.Thermalunits, t in 1:UC.T], power[unit.name, t-1]-power[unit.name, t]<=-unit.DeltaRampDown*start_down[unit.name, t] - unit.MinPower*is_on[unit.name, t] + (unit.MinPower+unit.DeltaRampDown)*is_on[unit.name, t-1])

        @constraint(MP, s >= sum(unit.LinearTerm*power[unit.name, t] for unit in UC.Thermalunits for t in 1:UC.T) + UC.PenaltyCost*sum(power_shedding[b,t] + power_curtailement[b,t] for b in 1:UC.Buses for t in 1:UC.T))

        flow = @variable(MP, [l in 1:length(UC.Lines), t in 1:UC.T])
        angle = @variable(MP, [b in 1:UC.Buses, t in 1:UC.T])
        @constraint(MP, [line in UC.Lines, t in 1:UC.T], flow[line.id,t]<=line.Fmax)
        @constraint(MP, [line in UC.Lines, t in 1:UC.T], flow[line.id,t]>=-line.Fmax)
        @constraint(MP, [line in UC.Lines, t in 1:UC.T], flow[line.id,t] == line.B12*(angle[line.b1,t]-angle[line.b2,t]))
        @constraint(MP, [t in 1:UC.T, b in 1:UC.Buses], sum(power[unit.name, t] for unit in UC.Thermalunits if unit.Bus==b) + power_shedding[b,t]- power_curtailement[b,t] + sum(flow[line.id, t] for line in UC.Lines if line.b2==b) - sum(flow[line.id, t] for line in UC.Lines if line.b1==b) == UC.Demandbus[b][t] + UC.DemandDev[b][t]*ξ[t])
    end
end

function build_sp(UC::UnitCommitment, MP::JuMP.Model, subproblem::SubproblemType, λ::Float64 = 1.0)
    is_on = JuMP.value.(MP[:is_on])
    is_on = Float64.(Int.(round.(is_on))) # should be 0-1
    start_up = JuMP.value.(MP[:start_up])
    start_up = Float64.(Int.(round.(start_up)))
    start_down = JuMP.value.(MP[:start_down])
    start_down = Float64.(Int.(round.(start_down)))
    SP = initializeJuMPModel()

    # dual variables (common to all)
    @variable(SP, μinit[1:UC.N])
    @variable(SP, μmax[1:UC.N, 1:UC.T] >= 0)
    @variable(SP, μmin[1:UC.N, 1:UC.T] >= 0)
    @variable(SP, μup[1:UC.N, 1:UC.T] >= 0)
    @variable(SP, μdown[1:UC.N, 1:UC.T] >= 0)

    @variable(SP, ν[t in 1:UC.T, b in 1:UC.Buses])

    @variable(SP, γmax[l in 1:length(UC.Lines), t in 1:UC.T]>=0)
    @variable(SP, γmin[l in 1:length(UC.Lines), t in 1:UC.T]>=0)
    @variable(SP, γangle[l in 1:length(UC.Lines), t in 1:UC.T])

    @constraint(SP, [t in 1:1:UC.T, b in 1:UC.Buses], ν[t, b] >= - UC.PenaltyCost)
    @constraint(SP, [t in 1:1:UC.T, b in 1:UC.Buses], ν[t, b] <= UC.PenaltyCost)

    @constraint(SP, [t in 1:UC.T-1, unit in UC.Thermalunits],μmin[unit.name, t] - μmax[unit.name, t] + μup[unit.name, t+1] - μup[unit.name, t] + μdown[unit.name, t] - μdown[unit.name, t+1] + ν[t, unit.Bus] == unit.LinearTerm)
    @constraint(SP, [unit in UC.Thermalunits],μmin[unit.name, UC.T] - μmax[unit.name, UC.T] - μup[unit.name, UC.T] + μdown[unit.name, UC.T] + ν[UC.T, unit.Bus] == unit.LinearTerm)
    @constraint(SP, [unit in UC.Thermalunits], μup[unit.name, 1] - μdown[unit.name, 1] + μinit[unit.name] == 0)

    @constraint(SP, [line in UC.Lines, t in 1:UC.T], γmin[line.id, t] - γmax[line.id, t] + ν[t, line.b2] - ν[t, line.b1] + γangle[line.id, t] == 0)

    @constraint(SP, [b in 1:UC.Buses, t in 1:UC.T], sum(γangle[line.id, t]*line.B12 for line in UC.Lines if line.b2==b) - sum(γangle[line.id, t]*line.B12 for line in UC.Lines if line.b1==b) == 0)
    # uncertainty set
    @variable(SP, ξ[1:UC.T], Bin)
    @constraint(SP, sum(ξ[t] for t in 1:UC.T) <= UC.budget)

    thermal_fixed_cost = sum(unit.ConstTerm*is_on[unit.name, t]+unit.StartUpCost*start_up[unit.name, t]+unit.StartDownCost*start_down[unit.name, t] for unit in UC.Thermalunits for t in 1:UC.T)

    hexpr = @expression(SP, thermal_fixed_cost
                            + sum(μinit[unit.name] * unit.InitialPower for unit in UC.Thermalunits)
                            + sum(μmin[unit.name,t] * unit.MinPower*is_on[unit.name, t] for unit in UC.Thermalunits, t in 1:UC.T)
                            - sum(μmax[unit.name,t] * unit.MaxPower*is_on[unit.name, t] for unit in UC.Thermalunits, t in 1:UC.T)
                            - sum(μup[unit.name, t] * (-unit.DeltaRampUp*start_up[unit.name, t]-unit.MinPower*is_on[unit.name, t-1] + (unit.MinPower+unit.DeltaRampUp)*is_on[unit.name, t]) for unit in UC.Thermalunits, t in 1:UC.T)
                            - sum(μdown[unit.name, t] * (-unit.DeltaRampDown*start_down[unit.name, t] - unit.MinPower*is_on[unit.name, t] + (unit.MinPower+unit.DeltaRampDown)*is_on[unit.name, t-1]) for unit in UC.Thermalunits, t in 1:UC.T)
                            - sum(line.Fmax*(γmin[line.id,t]+γmax[line.id, t]) for line in UC.Lines, t in 1:UC.T)
                            + sum(ν[t, b]*UC.Demandbus[b][t] for b in 1:UC.Buses, t in 1:UC.T))

    if subproblem ∈ [LinearizedKKT]
        @variable(SP, α[1:UC.T])
        @constraint(SP, [t in 1:UC.T], α[t]- sum(UC.DemandDev[b][t]*ν[t, b] for b in 1:UC.Buses)== 0)

        @variable(SP, δ[1:UC.T])
        @constraint(SP, [t in 1:UC.T], δ[t] <= sum(UC.DemandDev[b][t] for b in 1:UC.Buses) * UC.PenaltyCost * ξ[t])
        # @constraint(SP, [t in 1:UC.T], δ[t] >= -sum(UC.DemandDev[b,t] for b in Buses) * ξ[t])
        @constraint(SP, [t in 1:UC.T], δ[t] <= α[t]+sum(UC.DemandDev[b][t] for b in 1:UC.Buses) * UC.PenaltyCost * (1 - ξ[t]))

        @objective(SP, Max, hexpr + sum(δ[t] for t in 1:UC.T))

    end

    if subproblem ∈ [LagrangianDual]

        # @variable(SP, σ[1:UC.T]>=0)

        # @constraint(SP, [t in 1:UC.T], - sum(UC.DemandDev[b][t]*ν[t, b] for b in 1:UC.Buses) - σ[t] <= λ*(1-2*ξ[t]))
        # @objective(SP, Max, hexpr
        #     + λ * sum(ξ[t] for t in 1:UC.T)
        #     - sum(σ[t] for t in 1:UC.T)
        #     )

        #sclae sigma
        @variable(SP, σ[1:UC.T])

        @constraint(SP, [t in 1:UC.T], σ[t] <= sum(UC.DemandDev[b][t]*ν[t, b] for b in 1:UC.Buses)/λ + 1-ξ[t])
        @constraint(SP, [t in 1:UC.T], σ[t] <= ξ[t])

        @objective(SP, Max, hexpr
                    + λ * sum(σ[t] for t in 1:UC.T)
                    )
    end

    if subproblem ∈ [LagrangianDualbis]
        @variable(SP, α[1:UC.T])
        @variable(SP, σ[1:UC.T]>=0)
        @constraint(SP, [t in 1:UC.T], α[t] <= sum(UC.DemandDev[b][t]*ν[t, b] for b in 1:UC.Buses) + σ[t])

        @variable(SP, δ[1:UC.T])
        @constraint(SP, [t in 1:UC.T], δ[t] <= sum(UC.DemandDev[b][t] for b in 1:UC.Buses) * UC.PenaltyCost * ξ[t])
        # @constraint(SP, [t in 1:UC.T], δ[t] >= -sum(UC.DemandDev[b,t] for b in Buses) * ξ[t])
        @constraint(SP, [t in 1:UC.T], δ[t] <= α[t])

        @objective(SP, Max, hexpr + sum(δ[t] for t in 1:UC.T) - sum(σ[t] for t in 1:UC.T)
                            )

    end

    return SP
end

# function build_checking_sp(FL::FacilityLocation, MP::JuMP.Model)
#     y = JuMP.value.(MP[:y])
#     y = Float64.(Int.(round.(y))) # should be 0-1
#     SP = initializeJuMPModel()

#     # dual variables (common to all)
#     @variable(SP, α[FL.Customers] >= 0)
#     @variable(SP, β[FL.Facilities] >= 0)
#     @variable(SP, z[FL.Facilities], Bin)
#     @variable(SP, ρ[FL.Facilities] >= 0)

#     # uncertainty set
#     @constraint(SP, sum(z[j] for j in FL.Facilities) <= FL.budget)

#     # dual feasibility
#     @constraint(SP, [i in FL.Customers, j in FL.Facilities],
#         α[i] - β[j] <= FL.Distance[(i,j)] + (ρ[j])
#     )
#     @constraint(SP, [i in FL.Customers],
#         α[i] <= FL.PenaltyCost[i]
#     )
#     # indicators
#     @constraint(SP, [j in FL.Facilities],
#         !z[j] => {ρ[j] == 0}
#     )

#     # objective
#     @objective(SP, Max,
#         +sum(FL.FixedCost[j]*y[j] for j in FL.Facilities)
#         +sum(FL.Demand[i]*α[i] for i in FL.Customers)
#         -sum(FL.Capacity[j]*y[j]*β[j] for j in FL.Facilities)
#     )

#     return SP
# end

# function init_lagrangian_coefficient(FL::FacilityLocation, SP::JuMP.Model)
#     return maximum(value.(SP[:ρ]))
# end

function solve_second_stage_problem_lagrangian(UC::UnitCommitment, MP::JuMP.Model, SP::JuMP.Model, λ::Float64)

    is_on = JuMP.value.(MP[:is_on])
    is_on = Float64.(Int.(round.(is_on))) # should be 0-1
    start_up = JuMP.value.(MP[:start_up])
    start_up = Float64.(Int.(round.(start_up)))
    start_down = JuMP.value.(MP[:start_down])
    start_down = Float64.(Int.(round.(start_down)))
    ξ = JuMP.value.(SP[:ξ])
    ξ = Float64.(Int.(round.(ξ)))
    m = initializeJuMPModel()

    @variable(m, power[i in 1:UC.N, t in 0:UC.T] >= 0)
    @variable(m, power_shedding[b in 1:UC.Buses, t in 1:UC.T] >= 0)
    @variable(m, power_curtailement[b in 1:UC.Buses, t in 1:UC.T] >= 0)

    @constraint(m,  [unit in UC.Thermalunits], power[unit.name, 0]==unit.InitialPower)
    @constraint(m,  [unit in UC.Thermalunits, t in 1:UC.T], power[unit.name, t]<=unit.MaxPower*is_on[unit.name, t])
    @constraint(m,  [unit in UC.Thermalunits, t in 1:UC.T], power[unit.name, t]>=unit.MinPower*is_on[unit.name, t])
    @constraint(m,  [unit in UC.Thermalunits, t in 1:UC.T], power[unit.name, t]-power[unit.name, t-1]<=-unit.DeltaRampUp*start_up[unit.name, t]-unit.MinPower*is_on[unit.name, t-1] + (unit.MinPower+unit.DeltaRampUp)*is_on[unit.name, t])
    @constraint(m,  [unit in UC.Thermalunits, t in 1:UC.T], power[unit.name, t-1]-power[unit.name, t]<=-unit.DeltaRampDown*start_down[unit.name, t] - unit.MinPower*is_on[unit.name, t] + (unit.MinPower+unit.DeltaRampDown)*is_on[unit.name, t-1])

    @variable(m, flow[l in 1:length(UC.Lines), t in 1:UC.T])
    @variable(m, angle[b in 1:UC.Buses, t in 1:UC.T])
    @constraint(m, [line in UC.Lines, t in 1:UC.T], flow[line.id,t]<=line.Fmax)
    @constraint(m, [line in UC.Lines, t in 1:UC.T], flow[line.id,t]>=-line.Fmax)
    @constraint(m, [line in UC.Lines, t in 1:UC.T], flow[line.id,t] == line.B12*(angle[line.b1,t]-angle[line.b2,t]))

    @variable(m, u[1:UC.T] >= 0)
    @constraint(m, [t in 1:UC.T], u[t] <= 1)

    @constraint(m, [t in 1:UC.T, b in 1:UC.Buses], sum(power[unit.name, t] for unit in UC.Thermalunits if unit.Bus==b) + power_shedding[b,t]- power_curtailement[b,t] + sum(flow[line.id, t] for line in UC.Lines if line.b2==b) - sum(flow[line.id, t] for line in UC.Lines if line.b1==b) == UC.Demandbus[b][t]+UC.DemandDev[b][t]*u[t])

    @objective(m, Min,
        sum(unit.LinearTerm*power[unit.name, t] for unit in UC.Thermalunits for t in 1:UC.T) 
        + UC.PenaltyCost*sum(power_shedding[b,t] + power_curtailement[b,t] for b in 1:UC.Buses for t in 1:UC.T)
        + λ*sum(u[t]+ξ[t] - 2*ξ[t]*u[t] for t in 1:UC.T)
    )

    optimize!(m)

    step = sum(value(u[t])+ξ[t] - 2*ξ[t]*value(u[t]) for t in 1:UC.T)

    return step
end

function record_scenario(UC::UnitCommitment, SP::JuMP.Model, scenario_list::Dict)
    ξ = JuMP.value.(SP[:ξ])
    ξ = Int.(round.(ξ)) # should be 0-1
    deviation = Vector{Int}()
    for t in 1:UC.T
        if ξ[t] == 1
            push!(deviation, t)
        end
    end
    println(deviation)
    in_list = true
    if !haskey(scenario_list, deviation)
        in_list = false
        scenario_list[deviation] = true
    end

    return in_list
end

function solve_MP_FW(UC::UnitCommitment, SP::JuMP.Model, previous_scenario::Vector{Float64})
    println(previous_scenario)
    ξ_k = ones(UC.T)
    ξ_k1 = previous_scenario
    k = 0
    ub = 0.0
    while 0 <= k <= 50 && sum(abs.(ξ_k1 .- ξ_k)) >= 1
        ξ_k = copy(ξ_k1)
        k += 1
        JuMP.fix.(SP[:ξ], ξ_k; force = true)
        optimize!(SP)
        ub = JuMP.objective_value(SP)
        αval = JuMP.value.(SP[:α])
        ξ_k1 = zeros(UC.T)
        budget_indices = partialsort(1:UC.T, 1:min(UC.budget, UC.T), by=i->-αval[i])
        ξ_k1[budget_indices] .= 1
    end
    deviation = Vector{Int}()
    for t in 1:UC.T
        if ξ_k[t] == 1
            push!(deviation, t)
        end
    end
    return ub, ξ_k
end

function unfix_uncertainty_variables(UC::UnitCommitment, SP::JuMP.Model)
    JuMP.unfix.(SP[:ξ])
end

function compute_lagrangian_coefficient(UC::UnitCommitment, SP::JuMP.Model)
    return sum(sum(UC.DemandDev))*UC.PenaltyCost/15
end

