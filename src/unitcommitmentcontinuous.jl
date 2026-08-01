struct ThermalUnitContinuous
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

struct LineContinuous
    id::Int64
    b1::Int64
    b2::Int64 
    Fmax::Float64  
    B12::Float64 
end

mutable struct UnitCommitmentContinuous <: AbstractProblem
    name::String
    T::Int64
    N::Int64
    units::Vector{ThermalUnitContinuous}
    Buses::Int64
    Lines::Vector{LineContinuous}
    Demandbus::Vector{Vector{Float64}}
    BusWind::Vector{Int64}
    DemandDev::Vector{Vector{Float64}}
    PenaltyCost::Float64
    budget::Int64

    function UnitCommitmentContinuous(folder::String, budget::Int)
        """
        Parse the IEEE 118-bus instnace
        """
        TimeHorizon= 24
        syst = "data/UC/"*folder
        generators = CSV.read(joinpath(pwd(), syst, "generators.csv"), DataFrame; header=false)
        NumberUnits= parse(Int64,generators[end,1]) 
        N=NumberUnits
        name_instance="IEEE"*string(NumberUnits)
        units=Vector{ThermalUnitContinuous}(undef, N)
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
            unit=ThermalUnitContinuous(unit_name, Bus, MinPower, MaxPower, DeltaRampUp, DeltaRampDown, QuadTerm, StartUpCost, StartDownCost, LinearTerm, ConstTerm, InitialPower, InitUpDownTime, MinUpTime, MinDownTime)
            units[unit_name]=unit
        end

        maximum_load = CSV.read(joinpath(pwd(), syst, "maximum_load.csv"), DataFrame; header=false)
        Numbus=maximum_load[end,1]
        Buses=1:Numbus
        load_distribution_profile = CSV.read(joinpath(pwd(), syst, "load_distribution_profile.csv"), DataFrame; header=false)[:,2]/100
        Demandbus=[maximum_load[b,2]*load_distribution_profile for b in Buses]
        df_lines = CSV.read(joinpath(pwd(), syst, "lines.csv"), DataFrame; header=false)
        Numlines=parse(Int64,df_lines[end,1])
        Lines=Vector{LineContinuous}(undef, Numlines)
        for i in 2:Numlines+1
            id = parse(Int64, df_lines[i,1])
            b1 = parse(Int64, df_lines[i,2])
            b2 = parse(Int64, df_lines[i,3])
            fmax = parse(Float64, df_lines[i,6])
            X = 1/parse(Float64, df_lines[i,5])
            Lines[i-1]=LineContinuous(id, b1, b2, fmax, X)
        end

        BusWind=[b for b in Buses if sum(Demandbus[b])>=1]

        DemandDev = [[Demandbus[b][t]*1.96*0.025 for t in 1:TimeHorizon] for b in Buses]

        PenaltyCost = 300.0

        new(name_instance, 
            TimeHorizon, 
            N,
            units,
            Buses[end],
            Lines, 
            Demandbus, 
            BusWind, 
            DemandDev, 
            PenaltyCost,
            budget)
    end
end

mixed_integer_recourse(UC::UnitCommitmentContinuous) = false

complete_recourse(UC::UnitCommitmentContinuous) = true

objective_scale(UC::UnitCommitmentContinuous) = 1.0

indicator_uncertainty(UC::UnitCommitmentContinuous) = false

function solve_deterministic_problem(UC::UnitCommitmentContinuous)
    m = initializeJuMPModel()

    newgap=0.1/100
    set_optimizer_attribute(m, "MIPGap", newgap)

    @variable(m, is_on[i in 1:UC.N, t in 0:UC.T], Bin)
    @variable(m, start_up[i in 1:UC.N, t in 1:UC.T], Bin)
    @variable(m, start_down[i in 1:UC.N, t in 1:UC.T], Bin)

    @variable(m, thermal_cost>=0)
    @variable(m, thermal_fixed_cost>=0)
    @constraint(m, thermal_fixed_cost>=sum(unit.ConstTerm*is_on[unit.name, t]+unit.StartUpCost*start_up[unit.name, t]+unit.StartDownCost*start_down[unit.name, t] for unit in UC.units for t in 1:UC.T))
    
    @objective(m, Min, thermal_cost + thermal_fixed_cost)

    @constraint(m,  [unit in UC.units, t in 1:UC.T], is_on[unit.name, t]-is_on[unit.name, t-1]==start_up[unit.name, t]-start_down[unit.name, t])
    @constraint(m,  [unit in UC.units, t in 1:UC.T], start_up[unit.name, t]<=is_on[unit.name, t])
    @constraint(m,  [unit in UC.units, t in 1:UC.T], start_down[unit.name, t]<=1-is_on[unit.name, t])
    @constraint(m,  [unit in UC.units, t in 1:UC.T], sum(start_up[unit.name, τ] for τ in max(1, t-unit.MinUpTime+1):t)+1*(t<unit.MinUpTime-unit.InitUpDownTime+1)*(unit.InitUpDownTime>0) <= is_on[unit.name, t])
    @constraint(m,  [unit in UC.units, t in 1:UC.T], sum(start_down[unit.name, τ] for τ in max(1, t-unit.MinDownTime+1):t)+1*(t<unit.MinDownTime+unit.InitUpDownTime+1)*(unit.InitUpDownTime<0) <= 1-is_on[unit.name, t])
    @constraint(m,  [unit in UC.units], is_on[unit.name, 0]==(unit.InitUpDownTime>=0))

    for unit in UC.units
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

    @constraint(m,  [unit in UC.units], power[unit.name, 0]==unit.InitialPower)
    @constraint(m,  [unit in UC.units, t in 1:UC.T], power[unit.name, t]<=unit.MaxPower*is_on[unit.name, t])
    @constraint(m,  [unit in UC.units, t in 1:UC.T], power[unit.name, t]>=unit.MinPower*is_on[unit.name, t])
    @constraint(m,  [unit in UC.units, t in 1:UC.T], power[unit.name, t]-power[unit.name, t-1]<=-unit.DeltaRampUp*start_up[unit.name, t]-unit.MinPower*is_on[unit.name, t-1] + (unit.MinPower+unit.DeltaRampUp)*is_on[unit.name, t])
    @constraint(m,  [unit in UC.units, t in 1:UC.T], power[unit.name, t-1]-power[unit.name, t]<=-unit.DeltaRampDown*start_down[unit.name, t] - unit.MinPower*is_on[unit.name, t] + (unit.MinPower+unit.DeltaRampDown)*is_on[unit.name, t-1])


    @constraint(m, thermal_cost >= sum(unit.LinearTerm*power[unit.name, t] for unit in UC.units for t in 1:UC.T) + UC.PenaltyCost*sum(power_shedding[b,t] + power_curtailement[b,t] for b in 1:UC.Buses for t in 1:UC.T))

    @variable(m, flow[l in 1:length(UC.Lines), t in 1:UC.T])
    @variable(m, angle[b in 1:UC.Buses, t in 1:UC.T])
    @constraint(m, [line in UC.Lines, t in 1:UC.T], flow[line.id,t]<=line.Fmax)
    @constraint(m, [line in UC.Lines, t in 1:UC.T], flow[line.id,t]>=-line.Fmax)
    @constraint(m, [line in UC.Lines, t in 1:UC.T], flow[line.id,t] == line.B12*(angle[line.b1,t]-angle[line.b2,t]))

    @constraint(m, [t in 1:UC.T, b in 1:UC.Buses], sum(power[unit.name, t] for unit in UC.units if unit.Bus==b) + power_shedding[b,t]- power_curtailement[b,t] + sum(flow[line.id, t] for line in UC.Lines if line.b2==b) - sum(flow[line.id, t] for line in UC.Lines if line.b1==b) == UC.Demandbus[b][t])

    optimize!(m)

    return objective_value(m)
end

function init_master(UC::UnitCommitmentContinuous)
    m = initializeJuMPModel()

    @variable(m, is_on[i in 1:UC.N, t in 0:UC.T], Bin)
    @variable(m, start_up[i in 1:UC.N, t in 1:UC.T], Bin)
    @variable(m, start_down[i in 1:UC.N, t in 1:UC.T], Bin)

    @constraint(m,  [unit in UC.units, t in 1:UC.T], is_on[unit.name, t]-is_on[unit.name, t-1]==start_up[unit.name, t]-start_down[unit.name, t])
    @constraint(m,  [unit in UC.units, t in 1:UC.T], start_up[unit.name, t]<=is_on[unit.name, t])
    @constraint(m,  [unit in UC.units, t in 1:UC.T], start_down[unit.name, t]<=1-is_on[unit.name, t])
    @constraint(m,  [unit in UC.units, t in 1:UC.T], sum(start_up[unit.name, τ] for τ in max(1, t-unit.MinUpTime+1):t)+1*(t<unit.MinUpTime-unit.InitUpDownTime+1)*(unit.InitUpDownTime>0) <= is_on[unit.name, t])
    @constraint(m,  [unit in UC.units, t in 1:UC.T], sum(start_down[unit.name, τ] for τ in max(1, t-unit.MinDownTime+1):t)+1*(t<unit.MinDownTime+unit.InitUpDownTime+1)*(unit.InitUpDownTime<0) <= 1-is_on[unit.name, t])
    @constraint(m,  [unit in UC.units], is_on[unit.name, 0]==(unit.InitUpDownTime>=0))

    for unit in UC.units
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
    @constraint(m, thermal_fixed_cost>=sum(unit.ConstTerm*is_on[unit.name, t]+unit.StartUpCost*start_up[unit.name, t]+unit.StartDownCost*start_down[unit.name, t] for unit in UC.units for t in 1:UC.T))

    @variable(m, s >= 0)

    @objective(m, Min, s + thermal_fixed_cost)

    return m
end

function update_master_continous(UC::UnitCommitmentContinuous, MP::JuMP.Model, ξ::Vector{Int64}, master::MasterType)
    s = MP[:s]
    is_on = MP[:is_on]
    start_up = MP[:start_up]
    start_down = MP[:start_down]

    if master == CCG

        power_shedding = @variable(MP, [b in 1:UC.Buses, t in 1:UC.T], lower_bound = 0)
        power_curtailement = @variable(MP, [b in 1:UC.Buses, t in 1:UC.T], lower_bound = 0)

        power = @variable(MP, [i in 1:UC.N, t in 0:UC.T])
        @constraint(MP,  [unit in UC.units], power[unit.name, 0]==unit.InitialPower)
        @constraint(MP,  [unit in UC.units, t in 1:UC.T], power[unit.name, t]<=unit.MaxPower*is_on[unit.name, t])
        @constraint(MP,  [unit in UC.units, t in 1:UC.T], power[unit.name, t]>=unit.MinPower*is_on[unit.name, t])
        @constraint(MP,  [unit in UC.units, t in 1:UC.T], power[unit.name, t]-power[unit.name, t-1]<=-unit.DeltaRampUp*start_up[unit.name, t]-unit.MinPower*is_on[unit.name, t-1] + (unit.MinPower+unit.DeltaRampUp)*is_on[unit.name, t])
        @constraint(MP,  [unit in UC.units, t in 1:UC.T], power[unit.name, t-1]-power[unit.name, t]<=-unit.DeltaRampDown*start_down[unit.name, t] - unit.MinPower*is_on[unit.name, t] + (unit.MinPower+unit.DeltaRampDown)*is_on[unit.name, t-1])

        @constraint(MP, s >= sum(unit.LinearTerm*power[unit.name, t] for unit in UC.units for t in 1:UC.T) + UC.PenaltyCost*sum(power_shedding[b,t] + power_curtailement[b,t] for b in 1:UC.Buses for t in 1:UC.T))

        flow = @variable(MP, [l in 1:length(UC.Lines), t in 1:UC.T])
        angle = @variable(MP, [b in 1:UC.Buses, t in 1:UC.T])
        @constraint(MP, [line in UC.Lines, t in 1:UC.T], flow[line.id,t]<=line.Fmax)
        @constraint(MP, [line in UC.Lines, t in 1:UC.T], flow[line.id,t]>=-line.Fmax)
        @constraint(MP, [line in UC.Lines, t in 1:UC.T], flow[line.id,t] == line.B12*(angle[line.b1,t]-angle[line.b2,t]))
        @constraint(MP, [t in 1:UC.T, b in 1:UC.Buses], sum(power[unit.name, t] for unit in UC.units if unit.Bus==b) + power_shedding[b,t]- power_curtailement[b,t] + sum(flow[line.id, t] for line in UC.Lines if line.b2==b) - sum(flow[line.id, t] for line in UC.Lines if line.b1==b) == UC.Demandbus[b][t] + UC.DemandDev[b][t]*ξ[t])
    end
end

function build_sp(UC::UnitCommitmentContinuous, MP::JuMP.Model, subproblem::SubproblemType, λ=nothing)
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

    @constraint(SP, [t in 1:UC.T-1, unit in UC.units],μmin[unit.name, t] - μmax[unit.name, t] + μup[unit.name, t+1] - μup[unit.name, t] + μdown[unit.name, t] - μdown[unit.name, t+1] + ν[t, unit.Bus] == unit.LinearTerm)
    @constraint(SP, [unit in UC.units],μmin[unit.name, UC.T] - μmax[unit.name, UC.T] - μup[unit.name, UC.T] + μdown[unit.name, UC.T] + ν[UC.T, unit.Bus] == unit.LinearTerm)
    @constraint(SP, [unit in UC.units], μup[unit.name, 1] - μdown[unit.name, 1] + μinit[unit.name] == 0)

    @constraint(SP, [line in UC.Lines, t in 1:UC.T], γmin[line.id, t] - γmax[line.id, t] + ν[t, line.b2] - ν[t, line.b1] + γangle[line.id, t] == 0)

    @constraint(SP, [b in 1:UC.Buses, t in 1:UC.T], sum(γangle[line.id, t]*line.B12 for line in UC.Lines if line.b2==b) - sum(γangle[line.id, t]*line.B12 for line in UC.Lines if line.b1==b) == 0)
    # uncertainty set
    @variable(SP, ξ[1:UC.T], Bin)
    @constraint(SP, sum(ξ[t] for t in 1:UC.T) <= UC.budget)

    thermal_fixed_cost = sum(unit.ConstTerm*is_on[unit.name, t]+unit.StartUpCost*start_up[unit.name, t]+unit.StartDownCost*start_down[unit.name, t] for unit in UC.units for t in 1:UC.T)

    @variable(SP, αr[1:UC.T])
    @constraint(SP, [t in 1:UC.T], αr[t] == sum(UC.DemandDev[b][t]*ν[t, b] for b in 1:UC.Buses))

    hexpr = @expression(SP, thermal_fixed_cost
                            + sum(μinit[unit.name] * unit.InitialPower for unit in UC.units)
                            + sum(μmin[unit.name,t] * unit.MinPower*is_on[unit.name, t] for unit in UC.units, t in 1:UC.T)
                            - sum(μmax[unit.name,t] * unit.MaxPower*is_on[unit.name, t] for unit in UC.units, t in 1:UC.T)
                            - sum(μup[unit.name, t] * (-unit.DeltaRampUp*start_up[unit.name, t]-unit.MinPower*is_on[unit.name, t-1] + (unit.MinPower+unit.DeltaRampUp)*is_on[unit.name, t]) for unit in UC.units, t in 1:UC.T)
                            - sum(μdown[unit.name, t] * (-unit.DeltaRampDown*start_down[unit.name, t] - unit.MinPower*is_on[unit.name, t] + (unit.MinPower+unit.DeltaRampDown)*is_on[unit.name, t-1]) for unit in UC.units, t in 1:UC.T)
                            - sum(line.Fmax*(γmin[line.id,t]+γmax[line.id, t]) for line in UC.Lines, t in 1:UC.T)
                            + sum(ν[t, b]*UC.Demandbus[b][t] for b in 1:UC.Buses, t in 1:UC.T))

    if subproblem ∈ [CCGL]
        @variable(SP, α[1:UC.T])
        @constraint(SP, [t in 1:UC.T], α[t]- sum(UC.DemandDev[b][t]*ν[t, b] for b in 1:UC.Buses)== 0)

        @variable(SP, δ[1:UC.T])
        @constraint(SP, [t in 1:UC.T], δ[t] <= sum(UC.DemandDev[b][t] for b in 1:UC.Buses) * UC.PenaltyCost * ξ[t])
        # @constraint(SP, [t in 1:UC.T], δ[t] >= -sum(UC.DemandDev[b,t] for b in Buses) * ξ[t])
        @constraint(SP, [t in 1:UC.T], δ[t] <= α[t]+sum(UC.DemandDev[b][t] for b in 1:UC.Buses) * UC.PenaltyCost * (1 - ξ[t]))

        @objective(SP, Max, hexpr + sum(δ[t] for t in 1:UC.T))

    end

    if subproblem ∈ [CCGM]

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

    if subproblem ∈ [CCGM2]
        # println(round.([UC.PenaltyCost * sum(UC.DemandDev[b][t] for b in 1:UC.Buses) for t in 1:UC.T], digits=2))

        @variable(SP, σ[1:UC.T] >= 0)

        @constraint(SP, [t in 1:UC.T], - UC.PenaltyCost * sum(UC.DemandDev[b][t] for b in 1:UC.Buses)  * (1-2*ξ[t])<=sum(UC.DemandDev[b][t]*ν[t, b] for b in 1:UC.Buses) + σ[t] )
        @objective(SP, Max, hexpr + sum(UC.DemandDev[b][t]*UC.PenaltyCost*ξ[t] for b in 1:UC.Buses for t in 1:UC.T) - sum(σ[t] for t in 1:UC.T))

    end

    return SP
end

function compute_lagrangian_coefficient(UC::UnitCommitmentContinuous, MP::JuMP.Model)
    # return sum(sum(UC.DemandDev))*UC.PenaltyCost/15
    return 4020.0
end

function solve_second_stage_problem_lagrangian(UC::UnitCommitmentContinuous, MP::JuMP.Model, SP::JuMP.Model, λ::Float64)

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

    @constraint(m,  [unit in UC.units], power[unit.name, 0]==unit.InitialPower)
    @constraint(m,  [unit in UC.units, t in 1:UC.T], power[unit.name, t]<=unit.MaxPower*is_on[unit.name, t])
    @constraint(m,  [unit in UC.units, t in 1:UC.T], power[unit.name, t]>=unit.MinPower*is_on[unit.name, t])
    @constraint(m,  [unit in UC.units, t in 1:UC.T], power[unit.name, t]-power[unit.name, t-1]<=-unit.DeltaRampUp*start_up[unit.name, t]-unit.MinPower*is_on[unit.name, t-1] + (unit.MinPower+unit.DeltaRampUp)*is_on[unit.name, t])
    @constraint(m,  [unit in UC.units, t in 1:UC.T], power[unit.name, t-1]-power[unit.name, t]<=-unit.DeltaRampDown*start_down[unit.name, t] - unit.MinPower*is_on[unit.name, t] + (unit.MinPower+unit.DeltaRampDown)*is_on[unit.name, t-1])

    @variable(m, flow[l in 1:length(UC.Lines), t in 1:UC.T])
    @variable(m, angle[b in 1:UC.Buses, t in 1:UC.T])
    @constraint(m, [line in UC.Lines, t in 1:UC.T], flow[line.id,t]<=line.Fmax)
    @constraint(m, [line in UC.Lines, t in 1:UC.T], flow[line.id,t]>=-line.Fmax)
    @constraint(m, [line in UC.Lines, t in 1:UC.T], flow[line.id,t] == line.B12*(angle[line.b1,t]-angle[line.b2,t]))

    @variable(m, u[1:UC.T] >= 0)
    @constraint(m, [t in 1:UC.T], u[t] <= 1)

    @constraint(m, [t in 1:UC.T, b in 1:UC.Buses], sum(power[unit.name, t] for unit in UC.units if unit.Bus==b) + power_shedding[b,t]- power_curtailement[b,t] + sum(flow[line.id, t] for line in UC.Lines if line.b2==b) - sum(flow[line.id, t] for line in UC.Lines if line.b1==b) == UC.Demandbus[b][t]+UC.DemandDev[b][t]*u[t])

    @objective(m, Min,
        sum(unit.LinearTerm*power[unit.name, t] for unit in UC.units for t in 1:UC.T) 
        + UC.PenaltyCost*sum(power_shedding[b,t] + power_curtailement[b,t] for b in 1:UC.Buses for t in 1:UC.T)
        + λ*sum(u[t]+ξ[t] - 2*ξ[t]*u[t] for t in 1:UC.T)
    )

    optimize!(m)

    step = sum(value(u[t])+ξ[t] - 2*ξ[t]*value(u[t]) for t in 1:UC.T)

    return step
end

function record_scenario(UC::UnitCommitmentContinuous, ξ::Vector{Int64}, scenario_list::Dict)
    demandDeviations = Vector{Int}()
    for t in 1:UC.T
        if ξ[t] == 1
            push!(demandDeviations, t)
        end
    end
    in_list = true
    if !haskey(scenario_list, demandDeviations)
        in_list = false
        scenario_list[demandDeviations] = true
    end

    return in_list
end

function return_solution(UC::UnitCommitmentContinuous, computational_time::Float64, timelimit::Float64, LB::Float64, UB::Float64, Time_W::Vector{Any}, subproblemtype::SubproblemType)
    name_csv = "$(UC.name)_$(subproblemtype)_"*string(UC.budget)*"_"*string(Int(timelimit))
    # write results to CSV (long format: one datum per line)
    try
        df = DataFrame(
            metric = String[],
            value = String[],
            i = Vector{Union{Missing, Int}}(),
        )
        # scalars
        push!(df, ("name", "$(UC.name)_$(subproblemtype)", missing))
        push!(df, ("T", string(UC.T), missing))
        push!(df, ("budget", string(UC.budget), missing))
        push!(df, ("Time", string(computational_time), missing))
        push!(df, ("LB", string(LB), missing))
        push!(df, ("UB", string(UB), missing))
        push!(df, ("gap", string(gap(UB, LB)*100), missing))
        # arrays
        for (iter_idx, time_w) in enumerate(Time_W)
            push!(df, ("Time_worstcase_pb", string(round(time_w, digits=2)), iter_idx))
        end
        filepath = joinpath(pwd(), "results/UC", name_csv*".csv")
        CSV.write(filepath, df)
    catch e
        @warn("Failed to write results CSV", error = e)
    end
    return name_csv, UC.T, UC.budget, computational_time, round(LB, digits=2), round(gap(UB, LB), digits=2), Time_W
end