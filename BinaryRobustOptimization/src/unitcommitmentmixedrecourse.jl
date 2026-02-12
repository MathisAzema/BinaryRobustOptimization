struct ThermalUnitMixedRecourse
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

struct LineMixedRecourse
    id::Int64
    b1::Int64
    b2::Int64 
    Fmax::Float64  
    B12::Float64 
end

struct UnitCommitmentMixedRecourse <: AbstractProblem
    name::String
    T::Int64
    Nslow::Int64
    Nfast::Int64
    slowunits::Vector{ThermalUnitMixedRecourse}
    fastunits::Vector{ThermalUnitMixedRecourse}
    Buses::Int64
    Lines::Vector{LineMixedRecourse}
    Demandbus::Vector{Vector{Float64}}
    BusWind::Vector{Int64}
    DemandDev::Vector{Vector{Float64}}
    PenaltyCost::Float64
    budget::Int64

    # function UnitCommitmentMixedRecourse(filename::String, budget::Int, N1::Int)
    #     data=Dataset(filename)
    #     Blocks=keys(data.group)
    #     data_block=data.group[Blocks[1]]
    #     TimeHorizon= 24
    #     demand=round.(data_block["ActivePowerDemand"].var[:])
    #     N=size(keys(data_block.group))[1]-1
    #     slowunits=Vector{ThermalUnitMixedRecourse}(undef, N1)
    #     fastunits=Vector{ThermalUnitMixedRecourse}(undef, N-N1)
    #     k=0
    #     for unit in data_block.group
    #         if k < N
    #             if occursin("Block", first(unit))
    #                 unit_name=parse(Int64, split(first(unit), "_")[end])+1
    #                 k += 1
    #                 if k <= N1
    #                     unit_name=k
    #                 else
    #                     unit_name = k - N1
    #                 end
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

    #                     unit=ThermalUnitMixedRecourse(unit_name, Bus, MinPower, MaxPower, DeltaRampUp, DeltaRampDown, QuadTerm, StartUpCost, StartDownCost, LinearTerm, ConstTerm, InitialPower, InitUpDownTime, Int64(MinUpTime), Int64(MinDownTime))
    #                     if k <= N1
    #                         slowunits[unit_name]=unit
    #                     else
    #                         fastunits[unit_name]=unit
    #                     end
    #                 end
    #             end
    #         end
    #     end
    #     Lines = LineMixedRecourse[]
    #     Demandbus=[demand]
    #     DemandDev = [[demand[t]*1.96*0.025 for t in 1:TimeHorizon]]
    #     BusWind=[1]
    #     PenaltyCost = 300.0
    #     Buses = 1

    #     new(filename, 
    #         TimeHorizon, 
    #         N1,
    #         N-N1, 
    #         slowunits,
    #         fastunits,
    #         Buses,
    #         Lines, 
    #         Demandbus, 
    #         BusWind, 
    #         DemandDev, 
    #         PenaltyCost,
    #         budget)
    # end

    function UnitCommitmentMixedRecourse(folder::String, budget::Int, N1::Int)
        """
        Parse the IEEE 118-bus instnace
        """
        TimeHorizon= 24
        syst = "data/UC/"*folder
        generators = CSV.read(joinpath(pwd(), syst, "generators.csv"), DataFrame; header=false)
        NumberUnits= parse(Int64,generators[end,1]) 
        N=NumberUnits
        name_instance="IEEE"*string(NumberUnits)
        slowunits=Vector{ThermalUnitMixedRecourse}(undef, N1)
        fastunits=Vector{ThermalUnitMixedRecourse}(undef, N-N1)
        for i in 2:NumberUnits+1
            if i-1 <= N1
                unit_name=i-1
            else
                unit_name = i-1 - N1
            end
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
            unit=ThermalUnitMixedRecourse(unit_name, Bus, MinPower, MaxPower, DeltaRampUp, DeltaRampDown, QuadTerm, StartUpCost, StartDownCost, LinearTerm, ConstTerm, InitialPower, InitUpDownTime, MinUpTime, MinDownTime)
            if i-1 <= N1
                slowunits[unit_name]=unit
            else
                fastunits[unit_name]=unit
            end
        end

        maximum_load = CSV.read(joinpath(pwd(), syst, "maximum_load.csv"), DataFrame; header=false)
        Numbus=maximum_load[end,1]
        Buses=1:Numbus
        load_distribution_profile = CSV.read(joinpath(pwd(), syst, "load_distribution_profile.csv"), DataFrame; header=false)[:,2]/100
        Demandbus=[maximum_load[b,2]*load_distribution_profile for b in Buses]
        df_lines = CSV.read(joinpath(pwd(), syst, "lines.csv"), DataFrame; header=false)
        Numlines=parse(Int64,df_lines[end,1])
        Lines=Vector{LineMixedRecourse}(undef, Numlines)
        for i in 2:Numlines+1
            id = parse(Int64, df_lines[i,1])
            b1 = parse(Int64, df_lines[i,2])
            b2 = parse(Int64, df_lines[i,3])
            fmax = parse(Float64, df_lines[i,6])
            X = 1/parse(Float64, df_lines[i,5])
            Lines[i-1]=LineMixedRecourse(id, b1, b2, fmax, X)
        end

        BusWind=[b for b in Buses if sum(Demandbus[b])>=1]

        DemandDev = [[Demandbus[b][t]*1.96*0.025 for t in 1:TimeHorizon] for b in Buses]

        PenaltyCost = 300.0

        new(name_instance, 
            TimeHorizon, 
            N1,
            N-N1, 
            slowunits,
            fastunits,
            Buses[end],
            Lines, 
            Demandbus, 
            BusWind, 
            DemandDev, 
            PenaltyCost,
            budget)
    end
end

mixed_integer_recourse(UC::UnitCommitmentMixedRecourse) = true

complete_recourse(UC::UnitCommitmentMixedRecourse) = true

objective_scale(UC::UnitCommitmentMixedRecourse) = 1.0

indicator_uncertainty(UC::UnitCommitmentMixedRecourse) = false

function solve_deterministic_problem(UC::UnitCommitmentMixedRecourse)
    m = initializeJuMPModel()

    newgap=0.1/100
    set_optimizer_attribute(m, "MIPGap", newgap)

    @variable(m, is_on_slow[i in 1:UC.Nslow, t in 0:UC.T], Bin)
    @variable(m, start_up_slow[i in 1:UC.Nslow, t in 1:UC.T], Bin)
    @variable(m, start_down_slow[i in 1:UC.Nslow, t in 1:UC.T], Bin)

    @variable(m, is_on_fast[i in 1:UC.Nfast, t in 0:UC.T], Bin)
    @variable(m, start_up_fast[i in 1:UC.Nfast, t in 1:UC.T], Bin)
    @variable(m, start_down_fast[i in 1:UC.Nfast, t in 1:UC.T], Bin)

    @variable(m, thermal_cost>=0)
    @variable(m, thermal_fixed_cost_slow>=0)
    @variable(m, thermal_fixed_cost_fast>=0)
    @constraint(m, thermal_fixed_cost_slow>=sum(unit.ConstTerm*is_on_slow[unit.name, t]+unit.StartUpCost*start_up_slow[unit.name, t]+unit.StartDownCost*start_down_slow[unit.name, t] for unit in UC.slowunits for t in 1:UC.T))
    @constraint(m, thermal_fixed_cost_fast>=sum(unit.ConstTerm*is_on_fast[unit.name, t]+unit.StartUpCost*start_up_fast[unit.name, t]+unit.StartDownCost*start_down_fast[unit.name, t] for unit in UC.fastunits for t in 1:UC.T))

    @objective(m, Min, thermal_cost + thermal_fixed_cost_slow + thermal_fixed_cost_fast)

    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], is_on_slow[unit.name, t]-is_on_slow[unit.name, t-1]==start_up_slow[unit.name, t]-start_down_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], start_up_slow[unit.name, t]<=is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], start_down_slow[unit.name, t]<=1-is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], sum(start_up_slow[unit.name, τ] for τ in max(1, t-unit.MinUpTime+1):t)+1*(t<unit.MinUpTime-unit.InitUpDownTime+1)*(unit.InitUpDownTime>0) <= is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], sum(start_down_slow[unit.name, τ] for τ in max(1, t-unit.MinDownTime+1):t)+1*(t<unit.MinDownTime+unit.InitUpDownTime+1)*(unit.InitUpDownTime<0) <= 1-is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits], is_on_slow[unit.name, 0]==(unit.InitUpDownTime>=0))

    for unit in UC.slowunits
        if (unit.InitialPower-unit.MinPower)/unit.DeltaRampDown <0
            limit=-1
        else
            limit=Int64(ceil((unit.InitialPower-unit.MinPower)/unit.DeltaRampDown))
        end
        for t in 0:limit
            @constraint(m, is_on_slow[unit.name, t]==1)
        end
    end

    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], is_on_fast[unit.name, t]-is_on_fast[unit.name, t-1]==start_up_fast[unit.name, t]-start_down_fast[unit.name, t])
    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], start_up_fast[unit.name, t]<=is_on_fast[unit.name, t])
    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], start_down_fast[unit.name, t]<=1-is_on_fast[unit.name, t])
    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], sum(start_up_fast[unit.name, τ] for τ in max(1, t-unit.MinUpTime+1):t)+1*(t<unit.MinUpTime-unit.InitUpDownTime+1)*(unit.InitUpDownTime>0) <= is_on_fast[unit.name, t])
    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], sum(start_down_fast[unit.name, τ] for τ in max(1, t-unit.MinDownTime+1):t)+1*(t<unit.MinDownTime+unit.InitUpDownTime+1)*(unit.InitUpDownTime<0) <= 1-is_on_fast[unit.name, t])
    @constraint(m,  [unit in UC.fastunits], is_on_fast[unit.name, 0]==(unit.InitUpDownTime>=0))
    for unit in UC.fastunits
        if (unit.InitialPower-unit.MinPower)/unit.DeltaRampDown <0
            limit=-1
        else
            limit=Int64(ceil((unit.InitialPower-unit.MinPower)/unit.DeltaRampDown))
        end
        for t in 0:limit
            @constraint(m, is_on_fast[unit.name, t]==1)
        end
    end

    @variable(m, power_slow[i in 1:UC.Nslow, t in 0:UC.T] >= 0)
    @variable(m, power_fast[i in 1:UC.Nfast, t in 0:UC.T] >= 0)
    @variable(m, power_shedding[b in 1:UC.Buses, t in 1:UC.T] >= 0)
    @variable(m, power_curtailement[b in 1:UC.Buses, t in 1:UC.T] >= 0)

    @constraint(m,  [unit in UC.slowunits], power_slow[unit.name, 0]==unit.InitialPower)
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], power_slow[unit.name, t]<=unit.MaxPower*is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], power_slow[unit.name, t]>=unit.MinPower*is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], power_slow[unit.name, t]-power_slow[unit.name, t-1]<=-unit.DeltaRampUp*start_up_slow[unit.name, t]-unit.MinPower*is_on_slow[unit.name, t-1] + (unit.MinPower+unit.DeltaRampUp)*is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], power_slow[unit.name, t-1]-power_slow[unit.name, t]<=-unit.DeltaRampDown*start_down_slow[unit.name, t] - unit.MinPower*is_on_slow[unit.name, t] + (unit.MinPower+unit.DeltaRampDown)*is_on_slow[unit.name, t-1])

    @constraint(m,  [unit in UC.fastunits], power_fast[unit.name, 0]==unit.InitialPower)
    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], power_fast[unit.name, t]<=unit.MaxPower*is_on_fast[unit.name, t])
    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], power_fast[unit.name, t]>=unit.MinPower*is_on_fast[unit.name, t])
    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], power_fast[unit.name, t]-power_fast[unit.name, t-1]<=-unit.DeltaRampUp*start_up_fast[unit.name, t]-unit.MinPower*is_on_fast[unit.name, t-1] + (unit.MinPower+unit.DeltaRampUp)*is_on_fast[unit.name, t])
    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], power_fast[unit.name, t-1]-power_fast[unit.name, t]<=-unit.DeltaRampDown*start_down_fast[unit.name, t] - unit.MinPower*is_on_fast[unit.name, t] + (unit.MinPower+unit.DeltaRampDown)*is_on_fast[unit.name, t-1])


    @constraint(m, thermal_cost >= sum(unit.LinearTerm*power_slow[unit.name, t] for unit in UC.slowunits for t in 1:UC.T) + sum(unit.LinearTerm*power_fast[unit.name, t] for unit in UC.fastunits for t in 1:UC.T) + UC.PenaltyCost*sum(power_shedding[b,t] + power_curtailement[b,t] for b in 1:UC.Buses for t in 1:UC.T))

    @variable(m, flow[l in 1:length(UC.Lines), t in 1:UC.T])
    @variable(m, angle[b in 1:UC.Buses, t in 1:UC.T])
    @constraint(m, [line in UC.Lines, t in 1:UC.T], flow[line.id,t]<=line.Fmax)
    @constraint(m, [line in UC.Lines, t in 1:UC.T], flow[line.id,t]>=-line.Fmax)
    @constraint(m, [line in UC.Lines, t in 1:UC.T], flow[line.id,t] == line.B12*(angle[line.b1,t]-angle[line.b2,t]))

    @constraint(m, [t in 1:UC.T, b in 1:UC.Buses], sum(power_slow[unit.name, t] for unit in UC.slowunits if unit.Bus==b) + sum(power_fast[unit.name, t] for unit in UC.fastunits if unit.Bus==b) + power_shedding[b,t]- power_curtailement[b,t] + sum(flow[line.id, t] for line in UC.Lines if line.b2==b) - sum(flow[line.id, t] for line in UC.Lines if line.b1==b) == UC.Demandbus[b][t])

    optimize!(m)

    println((JuMP.value(thermal_fixed_cost_slow), JuMP.value(thermal_fixed_cost_fast), JuMP.value(thermal_cost)))

    return objective_value(m), m
end

function init_master(UC::UnitCommitmentMixedRecourse)
    m = initializeJuMPModel()

    @variable(m, is_on_slow[i in 1:UC.Nslow, t in 0:UC.T], Bin)
    @variable(m, start_up_slow[i in 1:UC.Nslow, t in 1:UC.T], Bin)
    @variable(m, start_down_slow[i in 1:UC.Nslow, t in 1:UC.T], Bin)

    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], is_on_slow[unit.name, t]-is_on_slow[unit.name, t-1]==start_up_slow[unit.name, t]-start_down_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], start_up_slow[unit.name, t]<=is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], start_down_slow[unit.name, t]<=1-is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], sum(start_up_slow[unit.name, τ] for τ in max(1, t-unit.MinUpTime+1):t)+1*(t<unit.MinUpTime-unit.InitUpDownTime+1)*(unit.InitUpDownTime>0) <= is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], sum(start_down_slow[unit.name, τ] for τ in max(1, t-unit.MinDownTime+1):t)+1*(t<unit.MinDownTime+unit.InitUpDownTime+1)*(unit.InitUpDownTime<0) <= 1-is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits], is_on_slow[unit.name, 0]==(unit.InitUpDownTime>=0))

    for unit in UC.slowunits
        if (unit.InitialPower-unit.MinPower)/unit.DeltaRampDown <0
            limit=-1
        else
            limit=Int64(ceil((unit.InitialPower-unit.MinPower)/unit.DeltaRampDown))
        end
        for t in 0:limit
            @constraint(m, is_on_slow[unit.name, t]==1)
        end
    end

    @variable(m, thermal_fixed_cost_slow>=0)

    @constraint(m, thermal_fixed_cost_slow>=sum(unit.ConstTerm*is_on_slow[unit.name, t]+unit.StartUpCost*start_up_slow[unit.name, t]+unit.StartDownCost*start_down_slow[unit.name, t] for unit in UC.slowunits for t in 1:UC.T))

    @variable(m, s >= 0)

    @objective(m, Min, s + thermal_fixed_cost_slow)

    return m
end

function update_master_mixed_integer(UC::UnitCommitmentMixedRecourse, MP_outer::JuMP.Model, MP_inner::JuMP.Model, master::MasterType)
    s = MP_outer[:s]
    is_on_slow = MP_outer[:is_on_slow]
    start_up_slow = MP_outer[:start_up_slow]
    start_down_slow = MP_outer[:start_down_slow]

    if master == CCG
        ξ = JuMP.value.(MP_inner[:ξ])
        ξ = Float64.(Int.(round.(ξ))) # should be 0-1

        is_on_fast = @variable(MP_outer, [i in 1:UC.Nfast, t in 0:UC.T], Bin)
        start_up_fast = @variable(MP_outer, [i in 1:UC.Nfast, t in 1:UC.T], Bin)
        start_down_fast = @variable(MP_outer, [i in 1:UC.Nfast, t in 1:UC.T], Bin)

        @constraint(MP_outer,  [unit in UC.fastunits, t in 1:UC.T], is_on_fast[unit.name, t]-is_on_fast[unit.name, t-1]==start_up_fast[unit.name, t]-start_down_fast[unit.name, t])
        @constraint(MP_outer,  [unit in UC.fastunits, t in 1:UC.T], start_up_fast[unit.name, t]<=is_on_fast[unit.name, t])
        @constraint(MP_outer,  [unit in UC.fastunits, t in 1:UC.T], start_down_fast[unit.name, t]<=1-is_on_fast[unit.name, t])
        @constraint(MP_outer,  [unit in UC.fastunits, t in 1:UC.T], sum(start_up_fast[unit.name, τ] for τ in max(1, t-unit.MinUpTime+1):t)+1*(t<unit.MinUpTime-unit.InitUpDownTime+1)*(unit.InitUpDownTime>0) <= is_on_fast[unit.name, t])
        @constraint(MP_outer,  [unit in UC.fastunits, t in 1:UC.T], sum(start_down_fast[unit.name, τ] for τ in max(1, t-unit.MinDownTime+1):t)+1*(t<unit.MinDownTime+unit.InitUpDownTime+1)*(unit.InitUpDownTime<0) <= 1-is_on_fast[unit.name, t])
        @constraint(MP_outer,  [unit in UC.fastunits], is_on_fast[unit.name, 0]==(unit.InitUpDownTime>=0))
        for unit in UC.fastunits
            if (unit.InitialPower-unit.MinPower)/unit.DeltaRampDown <0
                limit=-1
            else
                limit=Int64(ceil((unit.InitialPower-unit.MinPower)/unit.DeltaRampDown))
            end
            for t in 0:limit
                @constraint(MP_outer, is_on_fast[unit.name, t]==1)
            end
        end

        thermal_fixed_cost_fast = @variable(MP_outer, lower_bound=0)

        @constraint(MP_outer, thermal_fixed_cost_fast>=sum(unit.ConstTerm*is_on_fast[unit.name, t]+unit.StartUpCost*start_up_fast[unit.name, t]+unit.StartDownCost*start_down_fast[unit.name, t] for unit in UC.fastunits for t in 1:UC.T))

        power_slow = @variable(MP_outer, [i in 1:UC.Nslow, t in 0:UC.T])
        power_shedding = @variable(MP_outer, [b in 1:UC.Buses, t in 1:UC.T], lower_bound = 0)
        power_curtailement = @variable(MP_outer, [b in 1:UC.Buses, t in 1:UC.T], lower_bound = 0)

        @constraint(MP_outer,  [unit in UC.slowunits], power_slow[unit.name, 0]==unit.InitialPower)
        @constraint(MP_outer,  [unit in UC.slowunits, t in 1:UC.T], power_slow[unit.name, t]<=unit.MaxPower*is_on_slow[unit.name, t])
        @constraint(MP_outer,  [unit in UC.slowunits, t in 1:UC.T], power_slow[unit.name, t]>=unit.MinPower*is_on_slow[unit.name, t])
        @constraint(MP_outer,  [unit in UC.slowunits, t in 1:UC.T], power_slow[unit.name, t]-power_slow[unit.name, t-1]<=-unit.DeltaRampUp*start_up_slow[unit.name, t]-unit.MinPower*is_on_slow[unit.name, t-1] + (unit.MinPower+unit.DeltaRampUp)*is_on_slow[unit.name, t])
        @constraint(MP_outer,  [unit in UC.slowunits, t in 1:UC.T], power_slow[unit.name, t-1]-power_slow[unit.name, t]<=-unit.DeltaRampDown*start_down_slow[unit.name, t] - unit.MinPower*is_on_slow[unit.name, t] + (unit.MinPower+unit.DeltaRampDown)*is_on_slow[unit.name, t-1])

        power_fast = @variable(MP_outer, [i in 1:UC.Nfast, t in 0:UC.T])
        @constraint(MP_outer,  [unit in UC.fastunits], power_fast[unit.name, 0]==unit.InitialPower)
        @constraint(MP_outer,  [unit in UC.fastunits, t in 1:UC.T], power_fast[unit.name, t]<=unit.MaxPower*is_on_fast[unit.name, t])
        @constraint(MP_outer,  [unit in UC.fastunits, t in 1:UC.T], power_fast[unit.name, t]>=unit.MinPower*is_on_fast[unit.name, t])
        @constraint(MP_outer,  [unit in UC.fastunits, t in 1:UC.T], power_fast[unit.name, t]-power_fast[unit.name, t-1]<=-unit.DeltaRampUp*start_up_fast[unit.name, t]-unit.MinPower*is_on_fast[unit.name, t-1] + (unit.MinPower+unit.DeltaRampUp)*is_on_fast[unit.name, t])
        @constraint(MP_outer,  [unit in UC.fastunits, t in 1:UC.T], power_fast[unit.name, t-1]-power_fast[unit.name, t]<=-unit.DeltaRampDown*start_down_fast[unit.name, t] - unit.MinPower*is_on_fast[unit.name, t] + (unit.MinPower+unit.DeltaRampDown)*is_on_fast[unit.name, t-1])

        @constraint(MP_outer, s >= thermal_fixed_cost_fast + sum(unit.LinearTerm*power_slow[unit.name, t] for unit in UC.slowunits for t in 1:UC.T) + sum(unit.LinearTerm*power_fast[unit.name, t] for unit in UC.fastunits for t in 1:UC.T) + UC.PenaltyCost*sum(power_shedding[b,t] + power_curtailement[b,t] for b in 1:UC.Buses for t in 1:UC.T))

        flow = @variable(MP_outer, [l in 1:length(UC.Lines), t in 1:UC.T])
        angle = @variable(MP_outer, [b in 1:UC.Buses, t in 1:UC.T])
        @constraint(MP_outer, [line in UC.Lines, t in 1:UC.T], flow[line.id,t]<=line.Fmax)
        @constraint(MP_outer, [line in UC.Lines, t in 1:UC.T], flow[line.id,t]>=-line.Fmax)
        @constraint(MP_outer, [line in UC.Lines, t in 1:UC.T], flow[line.id,t] == line.B12*(angle[line.b1,t]-angle[line.b2,t]))
        @constraint(MP_outer, [t in 1:UC.T, b in 1:UC.Buses], sum(power_slow[unit.name, t] for unit in UC.slowunits if unit.Bus==b) + sum(power_fast[unit.name, t] for unit in UC.fastunits if unit.Bus==b) + power_shedding[b,t]- power_curtailement[b,t] + sum(flow[line.id, t] for line in UC.Lines if line.b2==b) - sum(flow[line.id, t] for line in UC.Lines if line.b1==b) == UC.Demandbus[b][t] + UC.DemandDev[b][t]*ξ[t])
    end
end

function init_master_inner_level(UC::UnitCommitmentMixedRecourse, master_inner::SubproblemType)
    m = initializeJuMPModel()
    # uncertainty set
    @variable(m, ξ[1:UC.T], Bin)
    @constraint(m, sum(ξ[t] for t in 1:UC.T) <= UC.budget)

    @variable(m, s,
        upper_bound = UC.PenaltyCost * sum(UC.Demandbus[b][t]+UC.DemandDev[b][t] for b in 1:UC.Buses, t in 1:UC.T)
    )

    @objective(m, Max, s)

    if master_inner ∈ [LagrangianDualbis]
        @variable(m, α[1:UC.T])
        @variable(m, δ[1:UC.T])

        @constraint(m, [t in 1:UC.T], δ[t] <= sum(UC.DemandDev[b][t] for b in 1:UC.Buses) * UC.PenaltyCost * ξ[t])
        @constraint(m, [t in 1:UC.T], δ[t] <= α[t]+sum(UC.DemandDev[b][t] for b in 1:UC.Buses) * UC.PenaltyCost * (1 - ξ[t]))
        # @constraint(m, [t in 1:UC.T], δ[t] >= α[t] - sum(UC.DemandDev[b][t] for b in 1:UC.Buses) * UC.PenaltyCost * (1 - ξ[t]))
    end

    return m
end

function init_master_inner_level(UC::UnitCommitmentMixedRecourse, MP_inner::JuMP.Model)
    ξval = JuMP.value.(MP_inner[:ξ])
    ξval = Float64.(Int.(round.(ξval))) # should be 0-1

    m = initializeJuMPModel()
    @variable(m, ξ[1:UC.T])
    @objective(m, Max, 0)
    fix.(ξ, ξval)
    optimize!(m)

    return m
end

function update_master_inner_level(UC::UnitCommitmentMixedRecourse, MP_outer::JuMP.Model, MP_inner::JuMP.Model, SP_inner::JuMP.Model, master_inner::SubproblemType, λ = nothing)
    is_on_slow = JuMP.value.(MP_outer[:is_on_slow])
    is_on_slow = Float64.(Int.(round.(is_on_slow)))
    start_up_slow = JuMP.value.(MP_outer[:start_up_slow])
    start_up_slow = Float64.(Int.(round.(start_up_slow)))
    start_down_slow = JuMP.value.(MP_outer[:start_down_slow])
    start_down_slow = Float64.(Int.(round.(start_down_slow)))

    is_on_fast = JuMP.value.(SP_inner[:is_on_fast])
    is_on_fast = Float64.(Int.(round.(is_on_fast)))
    start_up_fast = JuMP.value.(SP_inner[:start_up_fast])
    start_up_fast = Float64.(Int.(round.(start_up_fast)))
    start_down_fast = JuMP.value.(SP_inner[:start_down_fast])
    start_down_fast = Float64.(Int.(round.(start_down_fast)))

    update_master_inner_level(UC, MP_inner, is_on_slow, start_up_slow, start_down_slow, is_on_fast, start_up_fast, start_down_fast, master_inner, λ)
end

function update_master_inner_level(UC::UnitCommitmentMixedRecourse, MP_inner::JuMP.Model, is_on_slow, start_up_slow, start_down_slow, is_on_fast, start_up_fast, start_down_fast, master_inner::SubproblemType, λ = nothing)

    s = MP_inner[:s]
    ξ = MP_inner[:ξ]

    # dual variables (common to all)

    ν = @variable(MP_inner, [t in 1:UC.T, b in 1:UC.Buses])

    γmax = @variable(MP_inner, [l in 1:length(UC.Lines), t in 1:UC.T], lower_bound = 0)
    γmin = @variable(MP_inner, [l in 1:length(UC.Lines), t in 1:UC.T], lower_bound = 0)
    γangle = @variable(MP_inner, [l in 1:length(UC.Lines), t in 1:UC.T])

    @constraint(MP_inner, [t in 1:1:UC.T, b in 1:UC.Buses], ν[t, b] >= - UC.PenaltyCost)
    @constraint(MP_inner, [t in 1:1:UC.T, b in 1:UC.Buses], ν[t, b] <= UC.PenaltyCost)

    @constraint(MP_inner, [line in UC.Lines, t in 1:UC.T], γmin[line.id, t] - γmax[line.id, t] + ν[t, line.b2] - ν[t, line.b1] + γangle[line.id, t] == 0)

    @constraint(MP_inner, [b in 1:UC.Buses, t in 1:UC.T], sum(γangle[line.id, t]*line.B12 for line in UC.Lines if line.b2==b) - sum(γangle[line.id, t]*line.B12 for line in UC.Lines if line.b1==b) == 0)

    μinit_slow = @variable(MP_inner, [1:UC.Nslow])
    μmax_slow = @variable(MP_inner, [1:UC.Nslow, 1:UC.T], lower_bound = 0)
    μmin_slow = @variable(MP_inner, [1:UC.Nslow, 1:UC.T], lower_bound = 0)
    μup_slow = @variable(MP_inner, [1:UC.Nslow, 1:UC.T], lower_bound = 0)
    μdown_slow = @variable(MP_inner, [1:UC.Nslow, 1:UC.T], lower_bound = 0)
    @constraint(MP_inner, [t in 1:UC.T-1, unit in UC.slowunits],μmin_slow[unit.name, t] - μmax_slow[unit.name, t] + μup_slow[unit.name, t+1] - μup_slow[unit.name, t] + μdown_slow[unit.name, t] - μdown_slow[unit.name, t+1] + ν[t, unit.Bus] == unit.LinearTerm)
    @constraint(MP_inner, [unit in UC.slowunits],μmin_slow[unit.name, UC.T] - μmax_slow[unit.name, UC.T] - μup_slow[unit.name, UC.T] + μdown_slow[unit.name, UC.T] + ν[UC.T, unit.Bus] == unit.LinearTerm)
    @constraint(MP_inner, [unit in UC.slowunits], μup_slow[unit.name, 1] - μdown_slow[unit.name, 1] + μinit_slow[unit.name] == 0)

    μinit_fast = @variable(MP_inner, [1:UC.Nfast])
    μmax_fast = @variable(MP_inner, [1:UC.Nfast, 1:UC.T], lower_bound = 0)
    μmin_fast = @variable(MP_inner, [1:UC.Nfast, 1:UC.T], lower_bound = 0)
    μup_fast = @variable(MP_inner, [1:UC.Nfast, 1:UC.T], lower_bound = 0)
    μdown_fast = @variable(MP_inner, [1:UC.Nfast, 1:UC.T], lower_bound = 0)
    @constraint(MP_inner, [t in 1:UC.T-1, unit in UC.fastunits],μmin_fast[unit.name, t] - μmax_fast[unit.name, t] + μup_fast[unit.name, t+1] - μup_fast[unit.name, t] + μdown_fast[unit.name, t] - μdown_fast[unit.name, t+1] + ν[t, unit.Bus] == unit.LinearTerm)
    @constraint(MP_inner, [unit in UC.fastunits],μmin_fast[unit.name, UC.T] - μmax_fast[unit.name, UC.T] - μup_fast[unit.name, UC.T] + μdown_fast[unit.name, UC.T] + ν[UC.T, unit.Bus] == unit.LinearTerm)
    @constraint(MP_inner, [unit in UC.fastunits], μup_fast[unit.name, 1] - μdown_fast[unit.name, 1] + μinit_fast[unit.name] == 0)

    thermal_fixed_cost = sum(unit.ConstTerm*is_on_slow[unit.name, t]+unit.StartUpCost*start_up_slow[unit.name, t]+unit.StartDownCost*start_down_slow[unit.name, t] for unit in UC.slowunits for t in 1:UC.T; init = 0) + sum(unit.ConstTerm*is_on_fast[unit.name, t]+unit.StartUpCost*start_up_fast[unit.name, t]+unit.StartDownCost*start_down_fast[unit.name, t] for unit in UC.fastunits for t in 1:UC.T; init = 0)

    hexpr = @expression(MP_inner, thermal_fixed_cost
                    + sum(μinit_slow[unit.name] * unit.InitialPower for unit in UC.slowunits)
                    + sum(μmin_slow[unit.name,t] * unit.MinPower*is_on_slow[unit.name, t] for unit in UC.slowunits, t in 1:UC.T)
                    - sum(μmax_slow[unit.name,t] * unit.MaxPower*is_on_slow[unit.name, t] for unit in UC.slowunits, t in 1:UC.T)
                    - sum(μup_slow[unit.name, t] * (-unit.DeltaRampUp*start_up_slow[unit.name, t]-unit.MinPower*is_on_slow[unit.name, t-1] + (unit.MinPower+unit.DeltaRampUp)*is_on_slow[unit.name, t]) for unit in UC.slowunits, t in 1:UC.T)
                    - sum(μdown_slow[unit.name, t] * (-unit.DeltaRampDown*start_down_slow[unit.name, t] - unit.MinPower*is_on_slow[unit.name, t] + (unit.MinPower+unit.DeltaRampDown)*is_on_slow[unit.name, t-1]) for unit in UC.slowunits, t in 1:UC.T)
                    + sum(μinit_fast[unit.name] * unit.InitialPower for unit in UC.fastunits)
                    + sum(μmin_fast[unit.name,t] * unit.MinPower*is_on_fast[unit.name, t] for unit in UC.fastunits, t in 1:UC.T)
                    - sum(μmax_fast[unit.name,t] * unit.MaxPower*is_on_fast[unit.name, t] for unit in UC.fastunits, t in 1:UC.T)
                    - sum(μup_fast[unit.name, t] * (-unit.DeltaRampUp*start_up_fast[unit.name, t]-unit.MinPower*is_on_fast[unit.name, t-1] + (unit.MinPower+unit.DeltaRampUp)*is_on_fast[unit.name, t]) for unit in UC.fastunits, t in 1:UC.T)
                    - sum(μdown_fast[unit.name, t] * (-unit.DeltaRampDown*start_down_fast[unit.name, t] - unit.MinPower*is_on_fast[unit.name, t] + (unit.MinPower+unit.DeltaRampDown)*is_on_fast[unit.name, t-1]) for unit in UC.fastunits, t in 1:UC.T)
                    - sum(line.Fmax*(γmin[line.id,t]+γmax[line.id, t]) for line in UC.Lines, t in 1:UC.T)
                    + sum(ν[t, b]*UC.Demandbus[b][t] for b in 1:UC.Buses, t in 1:UC.T)
                    )

    if master_inner ∈ [LagrangianDualbis]
        α = MP_inner[:α]
        δ = MP_inner[:δ]
        
        σ = @variable(MP_inner, [1:UC.T], lower_bound = 0)

        @constraint(MP_inner, [t in 1:UC.T], α[t] <= sum(UC.DemandDev[b][t]*ν[t, b] for b in 1:UC.Buses) + σ[t])

        @constraint(MP_inner, s <= hexpr + sum(δ[t] for t in 1:UC.T) - sum(σ[t] for t in 1:UC.T))
    end

    if master_inner ∈ [LagrangianDual]

        # σ = @variable(MP_inner, [1:UC.T], lower_bound = 0)

        # @constraint(MP_inner, [t in 1:UC.T], - sum(UC.DemandDev[b][t]*ν[t, b] for b in 1:UC.Buses) - σ[t] <= λ*(1-2*ξ[t]))
        # @constraint(MP_inner, s <=hexpr + λ * sum(ξ[t] for t in 1:UC.T) - sum(σ[t] for t in 1:UC.T))

        #scale sigma
        σ = @variable(MP_inner, [1:UC.T])

        @constraint(MP_inner, [t in 1:UC.T], σ[t] <= sum(UC.DemandDev[b][t]*ν[t, b] for b in 1:UC.Buses)/λ + 1-ξ[t])
        @constraint(MP_inner, [t in 1:UC.T], σ[t] <= ξ[t])
        @constraint(MP_inner, s <= hexpr + λ * sum(σ[t] for t in 1:UC.T))
    end
end

function compute_lagrangian_coefficient(UC::UnitCommitmentMixedRecourse, MP_outer::JuMP.Model)
    return sum(sum(UC.DemandDev))*UC.PenaltyCost/15
end

function build_second_stage_problem(UC::UnitCommitmentMixedRecourse, MP_outer::JuMP.Model, MP_inner::JuMP.Model)
    is_on_slow = JuMP.value.(MP_outer[:is_on_slow])
    is_on_slow = Float64.(Int.(round.(is_on_slow)))
    start_up_slow = JuMP.value.(MP_outer[:start_up_slow])
    start_up_slow = Float64.(Int.(round.(start_up_slow)))
    start_down_slow = JuMP.value.(MP_outer[:start_down_slow])
    start_down_slow = Float64.(Int.(round.(start_down_slow)))

    ξ = JuMP.value.(MP_inner[:ξ])
    ξ = Float64.(Int.(round.(ξ))) # should be 0-1

    m = initializeJuMPModel()

    @variable(m, is_on_fast[i in 1:UC.Nfast, t in 0:UC.T], Bin)
    @variable(m, start_up_fast[i in 1:UC.Nfast, t in 1:UC.T], Bin)
    @variable(m, start_down_fast[i in 1:UC.Nfast, t in 1:UC.T], Bin)

    @variable(m, thermal_cost>=0)
    @variable(m, thermal_fixed_cost_slow>=0)
    @variable(m, thermal_fixed_cost_fast>=0)
    @constraint(m, thermal_fixed_cost_slow>=sum(unit.ConstTerm*is_on_slow[unit.name, t]+unit.StartUpCost*start_up_slow[unit.name, t]+unit.StartDownCost*start_down_slow[unit.name, t] for unit in UC.slowunits for t in 1:UC.T))
    @constraint(m, thermal_fixed_cost_fast>=sum(unit.ConstTerm*is_on_fast[unit.name, t]+unit.StartUpCost*start_up_fast[unit.name, t]+unit.StartDownCost*start_down_fast[unit.name, t] for unit in UC.fastunits for t in 1:UC.T))

    @objective(m, Min, thermal_cost + thermal_fixed_cost_slow + thermal_fixed_cost_fast)

    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], is_on_slow[unit.name, t]-is_on_slow[unit.name, t-1]==start_up_slow[unit.name, t]-start_down_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], start_up_slow[unit.name, t]<=is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], start_down_slow[unit.name, t]<=1-is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], sum(start_up_slow[unit.name, τ] for τ in max(1, t-unit.MinUpTime+1):t)+1*(t<unit.MinUpTime-unit.InitUpDownTime+1)*(unit.InitUpDownTime>0) <= is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], sum(start_down_slow[unit.name, τ] for τ in max(1, t-unit.MinDownTime+1):t)+1*(t<unit.MinDownTime+unit.InitUpDownTime+1)*(unit.InitUpDownTime<0) <= 1-is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits], is_on_slow[unit.name, 0]==(unit.InitUpDownTime>=0))

    for unit in UC.slowunits
        if (unit.InitialPower-unit.MinPower)/unit.DeltaRampDown <0
            limit=-1
        else
            limit=Int64(ceil((unit.InitialPower-unit.MinPower)/unit.DeltaRampDown))
        end
        for t in 0:limit
            @constraint(m, is_on_slow[unit.name, t]==1)
        end
    end

    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], is_on_fast[unit.name, t]-is_on_fast[unit.name, t-1]==start_up_fast[unit.name, t]-start_down_fast[unit.name, t])
    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], start_up_fast[unit.name, t]<=is_on_fast[unit.name, t])
    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], start_down_fast[unit.name, t]<=1-is_on_fast[unit.name, t])
    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], sum(start_up_fast[unit.name, τ] for τ in max(1, t-unit.MinUpTime+1):t)+1*(t<unit.MinUpTime-unit.InitUpDownTime+1)*(unit.InitUpDownTime>0) <= is_on_fast[unit.name, t])
    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], sum(start_down_fast[unit.name, τ] for τ in max(1, t-unit.MinDownTime+1):t)+1*(t<unit.MinDownTime+unit.InitUpDownTime+1)*(unit.InitUpDownTime<0) <= 1-is_on_fast[unit.name, t])
    @constraint(m,  [unit in UC.fastunits], is_on_fast[unit.name, 0]==(unit.InitUpDownTime>=0))
    for unit in UC.fastunits
        if (unit.InitialPower-unit.MinPower)/unit.DeltaRampDown <0
            limit=-1
        else
            limit=Int64(ceil((unit.InitialPower-unit.MinPower)/unit.DeltaRampDown))
        end
        for t in 0:limit
            @constraint(m, is_on_fast[unit.name, t]==1)
        end
    end

    @variable(m, power_slow[i in 1:UC.Nslow, t in 0:UC.T] >= 0)
    @variable(m, power_fast[i in 1:UC.Nfast, t in 0:UC.T] >= 0)
    @variable(m, power_shedding[b in 1:UC.Buses, t in 1:UC.T] >= 0)
    @variable(m, power_curtailement[b in 1:UC.Buses, t in 1:UC.T] >= 0)

    @constraint(m,  [unit in UC.slowunits], power_slow[unit.name, 0]==unit.InitialPower)
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], power_slow[unit.name, t]<=unit.MaxPower*is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], power_slow[unit.name, t]>=unit.MinPower*is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], power_slow[unit.name, t]-power_slow[unit.name, t-1]<=-unit.DeltaRampUp*start_up_slow[unit.name, t]-unit.MinPower*is_on_slow[unit.name, t-1] + (unit.MinPower+unit.DeltaRampUp)*is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], power_slow[unit.name, t-1]-power_slow[unit.name, t]<=-unit.DeltaRampDown*start_down_slow[unit.name, t] - unit.MinPower*is_on_slow[unit.name, t] + (unit.MinPower+unit.DeltaRampDown)*is_on_slow[unit.name, t-1])

    @constraint(m,  [unit in UC.fastunits], power_fast[unit.name, 0]==unit.InitialPower)
    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], power_fast[unit.name, t]<=unit.MaxPower*is_on_fast[unit.name, t])
    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], power_fast[unit.name, t]>=unit.MinPower*is_on_fast[unit.name, t])
    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], power_fast[unit.name, t]-power_fast[unit.name, t-1]<=-unit.DeltaRampUp*start_up_fast[unit.name, t]-unit.MinPower*is_on_fast[unit.name, t-1] + (unit.MinPower+unit.DeltaRampUp)*is_on_fast[unit.name, t])
    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], power_fast[unit.name, t-1]-power_fast[unit.name, t]<=-unit.DeltaRampDown*start_down_fast[unit.name, t] - unit.MinPower*is_on_fast[unit.name, t] + (unit.MinPower+unit.DeltaRampDown)*is_on_fast[unit.name, t-1])


    @constraint(m, thermal_cost >= sum(unit.LinearTerm*power_slow[unit.name, t] for unit in UC.slowunits for t in 1:UC.T) + sum(unit.LinearTerm*power_fast[unit.name, t] for unit in UC.fastunits for t in 1:UC.T) + UC.PenaltyCost*sum(power_shedding[b,t] + power_curtailement[b,t] for b in 1:UC.Buses for t in 1:UC.T))

    @variable(m, flow[l in 1:length(UC.Lines), t in 1:UC.T])
    @variable(m, angle[b in 1:UC.Buses, t in 1:UC.T])
    @constraint(m, [line in UC.Lines, t in 1:UC.T], flow[line.id,t]<=line.Fmax)
    @constraint(m, [line in UC.Lines, t in 1:UC.T], flow[line.id,t]>=-line.Fmax)
    @constraint(m, [line in UC.Lines, t in 1:UC.T], flow[line.id,t] == line.B12*(angle[line.b1,t]-angle[line.b2,t]))

    @constraint(m, [t in 1:UC.T, b in 1:UC.Buses], sum(power_slow[unit.name, t] for unit in UC.slowunits if unit.Bus==b) + sum(power_fast[unit.name, t] for unit in UC.fastunits if unit.Bus==b) + power_shedding[b,t]- power_curtailement[b,t] + sum(flow[line.id, t] for line in UC.Lines if line.b2==b) - sum(flow[line.id, t] for line in UC.Lines if line.b1==b) == UC.Demandbus[b][t]+ UC.DemandDev[b][t]*ξ[t])

    return m
end

function solve_second_stage_problem_lagrangian(UC::UnitCommitmentMixedRecourse, MP_outer::JuMP.Model, MP_inner::JuMP.Model, λ::Float64)
    is_on_slow = JuMP.value.(MP_outer[:is_on_slow])
    is_on_slow = Float64.(Int.(round.(is_on_slow)))
    start_up_slow = JuMP.value.(MP_outer[:start_up_slow])
    start_up_slow = Float64.(Int.(round.(start_up_slow)))
    start_down_slow = JuMP.value.(MP_outer[:start_down_slow])
    start_down_slow = Float64.(Int.(round.(start_down_slow)))

    ξ = JuMP.value.(MP_inner[:ξ])
    ξ = Float64.(Int.(round.(ξ))) # should be 0-1

    m = initializeJuMPModel()

    @variable(m, is_on_fast[i in 1:UC.Nfast, t in 0:UC.T], Bin)
    @variable(m, start_up_fast[i in 1:UC.Nfast, t in 1:UC.T], Bin)
    @variable(m, start_down_fast[i in 1:UC.Nfast, t in 1:UC.T], Bin)

    @variable(m, thermal_cost>=0)
    @variable(m, thermal_fixed_cost_slow>=0)
    @variable(m, thermal_fixed_cost_fast>=0)
    @constraint(m, thermal_fixed_cost_slow>=sum(unit.ConstTerm*is_on_slow[unit.name, t]+unit.StartUpCost*start_up_slow[unit.name, t]+unit.StartDownCost*start_down_slow[unit.name, t] for unit in UC.slowunits for t in 1:UC.T))
    @constraint(m, thermal_fixed_cost_fast>=sum(unit.ConstTerm*is_on_fast[unit.name, t]+unit.StartUpCost*start_up_fast[unit.name, t]+unit.StartDownCost*start_down_fast[unit.name, t] for unit in UC.fastunits for t in 1:UC.T))

    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], is_on_slow[unit.name, t]-is_on_slow[unit.name, t-1]==start_up_slow[unit.name, t]-start_down_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], start_up_slow[unit.name, t]<=is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], start_down_slow[unit.name, t]<=1-is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], sum(start_up_slow[unit.name, τ] for τ in max(1, t-unit.MinUpTime+1):t)+1*(t<unit.MinUpTime-unit.InitUpDownTime+1)*(unit.InitUpDownTime>0) <= is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], sum(start_down_slow[unit.name, τ] for τ in max(1, t-unit.MinDownTime+1):t)+1*(t<unit.MinDownTime+unit.InitUpDownTime+1)*(unit.InitUpDownTime<0) <= 1-is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits], is_on_slow[unit.name, 0]==(unit.InitUpDownTime>=0))

    for unit in UC.slowunits
        if (unit.InitialPower-unit.MinPower)/unit.DeltaRampDown <0
            limit=-1
        else
            limit=Int64(ceil((unit.InitialPower-unit.MinPower)/unit.DeltaRampDown))
        end
        for t in 0:limit
            @constraint(m, is_on_slow[unit.name, t]==1)
        end
    end

    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], is_on_fast[unit.name, t]-is_on_fast[unit.name, t-1]==start_up_fast[unit.name, t]-start_down_fast[unit.name, t])
    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], start_up_fast[unit.name, t]<=is_on_fast[unit.name, t])
    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], start_down_fast[unit.name, t]<=1-is_on_fast[unit.name, t])
    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], sum(start_up_fast[unit.name, τ] for τ in max(1, t-unit.MinUpTime+1):t)+1*(t<unit.MinUpTime-unit.InitUpDownTime+1)*(unit.InitUpDownTime>0) <= is_on_fast[unit.name, t])
    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], sum(start_down_fast[unit.name, τ] for τ in max(1, t-unit.MinDownTime+1):t)+1*(t<unit.MinDownTime+unit.InitUpDownTime+1)*(unit.InitUpDownTime<0) <= 1-is_on_fast[unit.name, t])
    @constraint(m,  [unit in UC.fastunits], is_on_fast[unit.name, 0]==(unit.InitUpDownTime>=0))
    for unit in UC.fastunits
        if (unit.InitialPower-unit.MinPower)/unit.DeltaRampDown <0
            limit=-1
        else
            limit=Int64(ceil((unit.InitialPower-unit.MinPower)/unit.DeltaRampDown))
        end
        for t in 0:limit
            @constraint(m, is_on_fast[unit.name, t]==1)
        end
    end

    @variable(m, power_slow[i in 1:UC.Nslow, t in 0:UC.T] >= 0)
    @variable(m, power_fast[i in 1:UC.Nfast, t in 0:UC.T] >= 0)
    @variable(m, power_shedding[b in 1:UC.Buses, t in 1:UC.T] >= 0)
    @variable(m, power_curtailement[b in 1:UC.Buses, t in 1:UC.T] >= 0)

    @constraint(m,  [unit in UC.slowunits], power_slow[unit.name, 0]==unit.InitialPower)
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], power_slow[unit.name, t]<=unit.MaxPower*is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], power_slow[unit.name, t]>=unit.MinPower*is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], power_slow[unit.name, t]-power_slow[unit.name, t-1]<=-unit.DeltaRampUp*start_up_slow[unit.name, t]-unit.MinPower*is_on_slow[unit.name, t-1] + (unit.MinPower+unit.DeltaRampUp)*is_on_slow[unit.name, t])
    @constraint(m,  [unit in UC.slowunits, t in 1:UC.T], power_slow[unit.name, t-1]-power_slow[unit.name, t]<=-unit.DeltaRampDown*start_down_slow[unit.name, t] - unit.MinPower*is_on_slow[unit.name, t] + (unit.MinPower+unit.DeltaRampDown)*is_on_slow[unit.name, t-1])

    @constraint(m,  [unit in UC.fastunits], power_fast[unit.name, 0]==unit.InitialPower)
    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], power_fast[unit.name, t]<=unit.MaxPower*is_on_fast[unit.name, t])
    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], power_fast[unit.name, t]>=unit.MinPower*is_on_fast[unit.name, t])
    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], power_fast[unit.name, t]-power_fast[unit.name, t-1]<=-unit.DeltaRampUp*start_up_fast[unit.name, t]-unit.MinPower*is_on_fast[unit.name, t-1] + (unit.MinPower+unit.DeltaRampUp)*is_on_fast[unit.name, t])
    @constraint(m,  [unit in UC.fastunits, t in 1:UC.T], power_fast[unit.name, t-1]-power_fast[unit.name, t]<=-unit.DeltaRampDown*start_down_fast[unit.name, t] - unit.MinPower*is_on_fast[unit.name, t] + (unit.MinPower+unit.DeltaRampDown)*is_on_fast[unit.name, t-1])


    @constraint(m, thermal_cost >= sum(unit.LinearTerm*power_slow[unit.name, t] for unit in UC.slowunits for t in 1:UC.T) + sum(unit.LinearTerm*power_fast[unit.name, t] for unit in UC.fastunits for t in 1:UC.T) + UC.PenaltyCost*sum(power_shedding[b,t] + power_curtailement[b,t] for b in 1:UC.Buses for t in 1:UC.T))

    @variable(m, flow[l in 1:length(UC.Lines), t in 1:UC.T])
    @variable(m, angle[b in 1:UC.Buses, t in 1:UC.T])
    @constraint(m, [line in UC.Lines, t in 1:UC.T], flow[line.id,t]<=line.Fmax)
    @constraint(m, [line in UC.Lines, t in 1:UC.T], flow[line.id,t]>=-line.Fmax)
    @constraint(m, [line in UC.Lines, t in 1:UC.T], flow[line.id,t] == line.B12*(angle[line.b1,t]-angle[line.b2,t]))

    @variable(m, u[1:UC.T] >= 0)
    @constraint(m, [t in 1:UC.T], u[t] <= 1)

    @constraint(m, [t in 1:UC.T, b in 1:UC.Buses], sum(power_slow[unit.name, t] for unit in UC.slowunits if unit.Bus==b) + sum(power_fast[unit.name, t] for unit in UC.fastunits if unit.Bus==b) + power_shedding[b,t]- power_curtailement[b,t] + sum(flow[line.id, t] for line in UC.Lines if line.b2==b) - sum(flow[line.id, t] for line in UC.Lines if line.b1==b) == UC.Demandbus[b][t]+ UC.DemandDev[b][t]*u[t])

    @objective(m, Min, thermal_cost + thermal_fixed_cost_slow + thermal_fixed_cost_fast +sum(λ*(((1 - 2ξ[t])*u[t]) + ξ[t]) for t in 1:UC.T))

    optimize!(m)

    step = sum(((1 - 2ξ[t])*value(u[t])) + ξ[t] for t in 1:UC.T)

    # println((λ, step))

    return step
end

function record_discrete_second_stage_decision(UC::UnitCommitmentMixedRecourse, SP_inner::JuMP.Model, discrete_decision_list::Dict)
    is_on_fast = JuMP.value.(SP_inner[:is_on_fast])
    is_on_fast = Float64.(Int.(round.(is_on_fast)))

    schedule = Vector{Tuple{Int, Int}}()
    for j in 1:UC.Nfast, t in 1:UC.T
        if is_on_fast[j,t] == 1
            push!(schedule, (j,t))
        end
    end
    in_list = true
    if !haskey(discrete_decision_list, schedule)
        in_list = false
        discrete_decision_list[schedule] = true
    end

    return in_list
end

function record_scenario(UC::UnitCommitmentMixedRecourse, MP_inner::JuMP.Model, scenario_list::Dict)
    ξ = JuMP.value.(MP_inner[:ξ])
    ξ = Int.(round.(ξ)) # should be 0-1
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