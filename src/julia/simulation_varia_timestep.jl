
begin
    using Statistics        # to use mean()
    using Printf            # to use @printf()
    using Plots             # to use histogram() and histogram!()
    using DelimitedFiles    # to use writedlm() (and readdlm())
    using Distributions     # to use Exponential() distribution for sampling
    using PlutoUI           # to use @with_terminal
    using Random            # to use rand!() inplace operation and seed!()
end

begin
    # constant parameters
    ## simulation
    # Random.seed!(10)
    const ibegin = 3
    const iend = 30
    ## system
    const NUM_SYSTEM = 100_0000
    # const NUM_SYSTEM = 10_0000
    # const NUM_NODE = 10
    const TIME_STEP = 1               # 1 hour
    const LIFE_LIMIT = 20_0000
    # const LIFE_LIMIT = 20_0000
    const STATE_NUM_NODE = 6
    const k = 3
    const w = 3_0000
    ## switch A
    const λA = 1 / 5.90e4             # hour
    const PA0 = exp(-λA * TIME_STEP)
    const PA1 = 0.20 * (1 - PA0)
    const PA2 = 0.15 * (1 - PA0)
    const PA3 = 0.65 * (1 - PA0)
    ## switch B
    const λB = 1 / 2.20e5             # hour
    const PB0 = exp(-λB * TIME_STEP)
    const PB1 = 0.45 * (1 - PB0)
    const PB2 = 0.55 * (1 - PB0)

    # counter
    QPF_counter = 0
    QSO_counter = 0
    QDM_counter = 0
    QMO_counter = 0
    QDN_counter = 0
    QFB_counter = 0
    # QPF::Int8 = QSO::Int8 = QDM::Int8 = QMO::Int8 = QDN::Int8 = QFB::Int8 = 0

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

function simulate(NUM_NODE::Int8, result::Result)
    global QPF_counter
    global QSO_counter
    global QDM_counter
    global QMO_counter
    global QDN_counter
    global QFB_counter

    gA = zeros(Int8, NUM_NODE)
    gB = zeros(Int8, NUM_NODE)
    gN = zeros(Int8, NUM_NODE)
    lifeA = zeros(NUM_NODE)
    lifeB = zeros(NUM_NODE)

    system_life = zeros(Float64, NUM_SYSTEM)
    reliability_counter = 0
    QPF_counter = 0
    QSO_counter = 0
    QDM_counter = 0
    QMO_counter = 0
    QDN_counter = 0
    QFB_counter = 0
    @inbounds for i = 1:NUM_SYSTEM
        system_life[i] = simulate_variable_timestep!(NUM_NODE, gA, gB, gN, lifeA, lifeB)
        system_life[i] >= w && (reliability_counter += 1)
    end

    averagelife = mean(system_life)
    reliability = reliability_counter / NUM_SYSTEM
    @printf("NUM_NODE:%3d\tMTTF: %12.6f\tReliability: %7.3f%%\t", NUM_NODE, averagelife, reliability * 100)
    @printf("QPF: %9.6f  QSO: %9.6f  QDM: %9.6f  QMO: %9.6f  QDN: %9.6f  QFB: %9.6f\n",
        QPF_counter / NUM_SYSTEM / NUM_NODE,
        QSO_counter / NUM_SYSTEM / NUM_NODE,
        QDM_counter / NUM_SYSTEM / NUM_NODE,
        QMO_counter / NUM_SYSTEM / NUM_NODE,
        QDN_counter / NUM_SYSTEM / NUM_NODE,
        QFB_counter / NUM_SYSTEM / NUM_NODE)
    # update return result
    return update_result(NUM_NODE, result, averagelife, reliability)
end

function simulate_variable_timestep!(NUM_NODE, gA, gB, gN, lifeA, lifeB)
    global QPF_counter
    global QSO_counter
    global QDM_counter
    global QMO_counter
    global QDN_counter
    global QFB_counter
    master_node::Int8 = initialize!(NUM_NODE, gA, gB, gN)

    life_counter::Float64 = 0
    compute_switchstate!(NUM_NODE, gA, gB, lifeA, lifeB)

    @inbounds for i = 1:2*NUM_NODE
        minA, idxA = findmin(lifeA)
        minB, idxB = findmin(lifeB)
        min_life = min(minA, minB)

        if min_life >= LIFE_LIMIT
            life_counter = LIFE_LIMIT
            break
        end

        if min_life == minA
            idx = idxA
            lifeA[idxA] = +Inf
        else
            idx = idxB
            lifeB[idxB] = +Inf
        end
        compute_node_perfstate!(NUM_NODE, gA, gB, gN, lifeA, lifeB, idx)

        master_node = compute_node_rolestate!(NUM_NODE, gN, master_node)

        Gsys::Int8 = compute_systemstate!(NUM_NODE, gN, master_node)
        if Gsys == 2 || Gsys == 3
            life_counter = min_life
        else
            life_counter = min_life
            break
        end
    end
    QPF::Int8 = QSO::Int8 = QDM::Int8 = QMO::Int8 = QDN::Int8 = QFB::Int8 = 0

    @inbounds for elem in gN
        elem == 0 && (QPF += 1; continue)
        elem == 1 && (QSO += 1; continue)
        elem == 2 && (QDM += 1; continue)
        elem == 3 && (QMO += 1; continue)
        elem == 4 && (QDN += 1; continue)
        elem == 5 && (QFB += 1; continue)
    end

    QPF > 0 && (QPF_counter += QPF)
    QSO > 0 && (QSO_counter += QSO)
    QDM > 0 && (QDM_counter += QDM)
    QMO > 0 && (QMO_counter += QMO)
    QDN > 0 && (QDN_counter += QDN)
    QFB > 0 && (QFB_counter += QFB)

    life_counter = min(life_counter, LIFE_LIMIT)
end

function initialize!(NUM_NODE, gA, gB, gN)
    fill!(gA, 0)
    fill!(gB, 0)

    fill!(gN, 0)

    master_node::Int8 = 1 # rand(1:NUM_NODE)
    return master_node
end

function compute_switchstate!(NUM_NODE, gA, gB, lifeA, lifeB)
    rand!(Exponential(1 / λA), lifeA)
    rand!(Exponential(1 / λB), lifeB)
    @inbounds for i = 1:NUM_NODE
        tolA = rand() * (1 - PA0)
        gA[i] = tolA < PA1 ? 1 : tolA < PA1 + PA2 ? 2 : 3
        tolB = rand() * (1 - PB0)
        gB[i] = tolB < PB1 ? 1 : 2
    end
    nothing
end

function compute_node_perfstate!(NUM_NODE, gA, gB, gN, lifeA, lifeB, idx)
    if lifeA[idx] != +Inf
        lifeB[idx] != +Inf && (gN[idx] = 0; return nothing)
        gB[idx] == 1 && (gN[idx] = 3; return nothing)
        gB[idx] == 2 && (gN[idx] = 1; return nothing)
    elseif gA[idx] == 1
        lifeB[idx] != +Inf && (gN[idx] = 1; return nothing)
        gB[idx] == 1 && (gN[idx] = 5; return nothing)
        gB[idx] == 2 && (gN[idx] = 1; return nothing)
    elseif gA[idx] == 2
        lifeB[idx] != +Inf && (gN[idx] = 2; return nothing)
        gB[idx] == 1 && (gN[idx] = 3; return nothing)
        gB[idx] == 2 && (gN[idx] = 4; return nothing)
    elseif gA[idx] == 3
        lifeB[idx] != +Inf && (gN[idx] = 4; return nothing)
        gB[idx] == 1 && (gN[idx] = 4; return nothing)
        gB[idx] == 2 && (gN[idx] = 4; return nothing)
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

function compute_node_rolestate!(NUM_NODE, gN, master_node)
    # role state transition process is based on the formula given above
    # node (performance) vector: gN
    if !ok_for_master(gN[master_node])

        alert_mintime = 1.0 # rand() [0.0, 1.0)
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

julia_main()
