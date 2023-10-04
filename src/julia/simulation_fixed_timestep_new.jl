begin
    using Statistics          # to use mean()
    using Printf              # to use @printf()
    using Plots               # to use histogram() and histogram!()
    using DelimitedFiles      # to use writedlm() (and readdlm())
    using Distributions       # to use Exponential() distribution for sampling
    using Random              # to use rand!() inplace operation and seed!()
end

begin
    # constant parameters
    ## simulation
    Random.seed!(10)
    const ibegin = 3
    const iend = 20

    ## system 
    const NUM_SYSTEM = 1000
    # const NUM_SYSTEM = 2000
    # NUM_NODE = 10                   # this should also be optimized later
    const TIME_STEP = 1               # 1 hour
    # const LIFE_LIMIT = 60_0000      # suitable life limit (roughly 5% of bug value)
    const LIFE_LIMIT = 20_0000
    const STATE_NUM_NODE = 6
    const k = 3
    const w = 3_0000
    ## switch A
    const λA = 1 / 5.90e4             # hour
    const PA0 = exp(-λA * TIME_STEP)
    const PEA1 = 0.20 * (1 - PA0)
    const PEA2 = 0.15 * (1 - PA0)
    const PEA3 = 0.65 * (1 - PA0)
    ## switch B
    const λB = 1 / 2.20e5             # hour
    const PB0 = exp(-λB * TIME_STEP)
    const PEB1 = 0.45 * (1 - PB0)
    const PEB2 = 0.55 * (1 - PB0)
    nothing
end

mutable struct Result
    averagelife_max::Float64
    averagelife_idx::Int8
    reliability_max::Float64
    reliability_idx::Int8
end

function update_result(NUM_NODE, result, averagelife, reliability)

    if max(result.averagelife_max, averagelife) == averagelife
        result.averagelife_max = averagelife
        result.averagelife_idx = NUM_NODE
    end

    if max(result.reliability_max, reliability) == reliability
        result.reliability_max = reliability
        result.reliability_idx = NUM_NODE
    end
    return result
end

function julia_main()
    result = Result(0, 0, 0, 0)
    @time for i::Int8 = ibegin:iend
        result = simulate(i, result)
    end
    @printf("MTTF: %3d%16.4f\n", result.averagelife_idx, result.averagelife_max)
    @printf("R(w): %3d%15.2f%%\n", result.reliability_idx, result.reliability_max * 100)
end

function simulate(NUM_NODE, result)
    # variables definition for each simulation
    gA = zeros(Int8, NUM_NODE)
    gB = zeros(Int8, NUM_NODE)
    gN = zeros(Int8, NUM_NODE)

    system_life = zeros(Int32, NUM_SYSTEM)
    reliability_counter = 0
    @inbounds for i = 1:NUM_SYSTEM
        system_life[i] = simulate_fixed_timestep!(NUM_NODE, gA, gB, gN)
        system_life[i] >= w && (reliability_counter += 1)
    end

    averagelife = mean(system_life)
    reliability = reliability_counter / NUM_SYSTEM
    @printf("NUM_NODE:%3d\tMTTF: %12.6f\t\tReliability: %7.3f%%\n", NUM_NODE, averagelife, reliability * 100)

    # update return result
    return update_result(NUM_NODE, result, averagelife, reliability)
end

function simulate_fixed_timestep!(NUM_NODE, gA, gB, gN)
    master_node::Int8 = initialize!(NUM_NODE, gA, gB, gN)

    state_system = false
    life_counter::Int32 = 0
    while !state_system && life_counter < LIFE_LIMIT
        compute_switchstate!(NUM_NODE, gA, gB)
        compute_nodestate!(NUM_NODE, gA, gB, gN)
        master_node = reselect_master!(NUM_NODE, gN, master_node)

        Gsys::Int8 = compute_systemstate!(NUM_NODE, gN, master_node)
        if Gsys == 2 || Gsys == 3
            life_counter += 1
        else
            life_counter += 1
            state_system = true
        end
    end
    life_counter
end

function initialize!(NUM_NODE, gA, gB, gN)
    fill!(gA, 0)  # 状态全清0，代表正常
    fill!(gB, 0)
    # 节点状态为完好
    fill!(gN, 0)
    # 随机选取1个节点为主节点
    master_node::Int8 = 1 # rand(1:NUM_NODE)
    return master_node
end

function compute_switchstate!(NUM_NODE, gA, gB)
    @inbounds for i = 1:NUM_NODE
        sample = rand()
        if gA[i] != 0 || sample < PA0 # 当前切换器已经坏了，不需要计算||仍然正常工作
        # do nothing
        else
            tol = sample - PA0
            gA[i] = tol < PEA1 ? 1 : tol < PEA1 + PEA2 ? 2 : 3
        end
        sample = rand()
        if gB[i] != 0 || sample < PB0
            # do nothing
        else
            gB[i] = sample < PB0 + PEB1 ? 1 : 2
        end
    end
    nothing
end

function compute_nodestate!(NUM_NODE, gA, gB, gN)
    @inbounds for i = 1:NUM_NODE
        if gA[i] == 0
            gB[i] == 0 && (gN[i] = 0; continue)
            gB[i] == 1 && (gN[i] = 3; continue)
            gB[i] == 2 && (gN[i] = 1; continue)
        elseif gA[i] == 1
            gB[i] == 0 && (gN[i] = 1; continue)
            gB[i] == 1 && (gN[i] = 5; continue)
            gB[i] == 2 && (gN[i] = 1; continue)
        elseif gA[i] == 2
            gB[i] == 0 && (gN[i] = 2; continue)
            gB[i] == 1 && (gN[i] = 3; continue)
            gB[i] == 2 && (gN[i] = 4; continue)
        elseif gA[i] == 3
            gB[i] == 0 && (gN[i] = 4; continue)
            gB[i] == 1 && (gN[i] = 4; continue)
            gB[i] == 2 && (gN[i] = 4; continue)
        end
    end
    nothing
end

function ok_for_master(gNi)
    gNi == 0 && return true
    gNi == 1 && return false
    gNi == 2 && return true
    gNi == 3 && return true # MO
    gNi == 4 && return false
    gNi == 5 && return false
end

function reselect_master!(NUM_NODE, gN, master_node)
    if !ok_for_master(gN[master_node])

        alert_mintime = 1.0 # rand() [0.0,1)
        mintime_index = 0
        alert_counter = rand(NUM_NODE)

        @inbounds for i = 1:NUM_NODE
            if ok_for_master(gN[i])
                tmp = alert_counter[i]

                alert_mintime = min(alert_mintime, tmp)
                mintime_index = tmp == alert_mintime ? i : mintime_index
            end
        end
        mintime_index == 0 && return master_node

        master_node = mintime_index
    end
    cnt = 0
    idx = 0
    flag = false # flag for FB
    @inbounds for i = 1:NUM_NODE
        gN[i] == 3 && (cnt = cnt + 1; idx = i; continue)
    end
    if cnt == 1
        master_node = idx
        return master_node
        # elseif cnt >= 2
        #     master_node = 0 # 这里随意设置，因为系统必然失效
        # elseif cnt == 0 && flag
        #     master_node = 0 # 这里随意设置，因为系统必然失效
        # else
        #     # 也就是说现在信号是好的，不需要重选
    end
    master_node
end

function compute_systemstate!(NUM_NODE, gN, master_node)
    QPF::Int8 = QSO::Int8 = QDM::Int8 = QMO::Int8 = QDN::Int8 = QFB::Int8 = 0
    # Gsys::Int8 = 2

    @inbounds for elem in gN
        elem == 0 && (QPF += 1; continue)
        elem == 1 && (QSO += 1; continue)
        elem == 2 && (QDM += 1; continue)
        elem == 3 && (QMO += 1; continue)
        elem == 4 && (QDN += 1; continue)
        elem == 5 && (QFB += 1; continue)
    end

    C1::Bool = (QFB >= 1)
    C2::Bool = (QMO >= 2)
    C3::Bool = (QPF + QMO + QDM == 0)
    C4::Bool = ((QPF + QSO + 1((QMO + QDM) > 0)) < k)
    C5::Bool = (QFB == 0)
    C6::Bool = (QMO == 1 && QPF + QSO >= k - 1)
    C7::Bool = ((QMO == 0 && QPF >= 1 && QPF + QSO >= k) || (QMO == 0 && QPF == 0 && QDM >= 1 && QSO >= k - 1))
    C8::Bool = (QFB + QMO == 0)
    C9::Bool = (QPF >= 1 && (QPF + QSO == k - 1) && QDM >= 1)

    if C1 || C2 || C3 || C4
        Gsys = 1
    elseif C5 && (C6 || C7)
        Gsys = 2
    elseif C8 && C9
        # cond = QDM / (QDM + QPF)
        # if rand() < cond
        if gN[master_node] == 2
            Gsys = 3
        elseif gN[master_node] == 0
            Gsys = 4
        end
    end
    return Gsys
end

@time julia_main()