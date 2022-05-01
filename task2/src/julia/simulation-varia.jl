### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 5ed627ad-1c4f-4050-a0d8-297222e52e54
begin
    using Statistics        # to use mean()
    using Printf            # to use @printf()
    using Plots             # to use histogram() and histogram!()
    using DelimitedFiles    # to use writedlm() (and readdlm())
    using Distributions     # to use Exponential() distribution for sampling
    using PlutoUI           # to use @with_terminal
    using Random            # to use rand!() inplace operation
end

# ╔═╡ a995f5d3-4a32-41ea-913f-f551af205197
md"## Code Implementation
In this section, I will give you the exact code implementation for this task.
"

# ╔═╡ 93f9ab9b-35f9-48ea-ac96-40e2f525cb35
md"**Notice**: In `julia` language, the `begin-end` blocks do not introduce a new scope block."

# ╔═╡ 359d6bd1-78c3-43e6-8555-2d811903eed9
md"### Packages"

# ╔═╡ cea03e40-0119-476e-92c5-7742e675b5ea
md"### Constant Parameters"

# ╔═╡ c78092ff-73c2-40c2-b8e2-3f1178db4152
begin
    # constant parameters
    ## system 
    const NUM_SYSTEM = 1000_0000
    # const NUM_SYSTEM = 2000
    const NUM_NODE = 10               # this should also be optimized later
    const TIME_STEP = 1               # 1 hour
    const LIFE_LIMIT = 200_0000        # suitable life limit (roughly 5% of bug value)
    # const LIFE_LIMIT = 20_0000
    const STATE_NUM_NODE = 6
    const k = 3
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

# ╔═╡ d3313a3a-5252-478e-9220-ef41b9921b7f
md"### Functions Defination"

# ╔═╡ 59f49b80-7861-477d-b9d5-270d12ae0fbb
md"#### Fixed Time step"

# ╔═╡ 6e802270-f77d-4d5b-a941-59b14f4acd57
function estimate_switch_state!(gA, gB)
    # 各switch有一定概率正常工作，一定概率出现异常
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


# ╔═╡ f5d6b344-cb0d-41f8-a1bf-158d7e7f7b75
function estimate_node_state!(gA, gB, gN)
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


# ╔═╡ a554951b-5dc3-40a2-90bc-d08164c63fc7
function ok_for_master(gNi)
    gNi == 0 && return true
    gNi == 1 && return false
    gNi == 2 && return true
    gNi == 3 && return true # MO
    gNi == 4 && return false
    gNi == 5 && return false
end

# ╔═╡ e22abba0-b15d-4ee4-9f0a-4bea779d2177

# 节点重选 
function reselect_master!(gN, master_node)  # 1
    if !ok_for_master(gN[master_node])
        # 说明主节点坏掉了
        alert_mintime = 1.0 # rand() [0.0,1)
        mintime_index = 0
        alert_count = rand(NUM_NODE)
        @inbounds for i = 1:NUM_NODE
            if ok_for_master(gN[i])
                tmp = alert_mintime
                alert_mintime = min(alert_mintime, alert_count[i])
                mintime_index = tmp != alert_mintime ? i : mintime_index
            end
        end
        if mintime_index == 0
            return
        end  #! 如果没有找到合适的节点，直接返回，反正后面讨论系统的状态的时候会计算到
        master_node = mintime_index
    end
    # 如果原来的主节点可用，通常情况下是不需要进行节点重选的，这个时候需要把需要进行节点重选的情况筛选出来
    # 信号出现异常，这个时候要么是MO，要么是FB
    # 我们只需要讨论一种情况，如果系统有一个MO出现
    cnt = 0
    idx = 0
    flag = false # flag for FB
    @inbounds for i = 1:NUM_NODE
        gN[i] == 3 && (cnt = cnt + 1; idx = i; continue)
        gN[i] == 5 && (flag = true; continue)
    end
    if cnt == 1
        master_node = idx
        # elseif cnt >= 2
        #     master_node = 0 # 这里随意设置，因为系统必然失效
        # elseif cnt == 0 && flag
        #     master_node = 0 # 这里随意设置，因为系统必然失效
        # else
        #     # 也就是说现在信号是好的，不需要重选
    end
    nothing
end

# ╔═╡ f806096c-8086-4e52-9e30-30d5d1087c70
# @with_terminal @time julia_main()

# ╔═╡ 86e5acd8-7dac-4cae-8e58-3fe80d21ec09
md"#### Variable Time step
Here I use variable time step methology to estimate component life, which can significantly lower the time cost during simulation

The differences between the fixed time step approach are listed as follows:
- Two vectors to store components' life are needed, but each simulation only run such sample once
- The function `update_node_state()` should be explicitly separated into function calls, due to the skip of `switch_state_estimation()` and the changes of `node_state_estimation()`, and the only need turns to be `reselect_master()`
- Actually a better choice is to rewrite code in another `Pluto` notebook, and integrate both into single one after both run correctly"

# ╔═╡ 5ec48cc2-cb7d-4e42-b1cf-56afacccc768
function estimate_switch_state!(gA, gB, lifeA, lifeB)
    distrA = Exponential(1 / λA)
    distrB = Exponential(1 / λB)
    rand!(distrA, lifeA)
    rand!(distrB, lifeB)
    @inbounds for i = 1:NUM_NODE
        tolA = rand() * (1 - PA0)
        gA[i] = tolA < PEA1 ? 1 : tolA < PEA1 + PEA2 ? 2 : 3
        tolB = rand() * (1 - PB0)
        gB[i] = tolB < PEB1 ? 1 : 2
    end
    nothing
end

# ╔═╡ b2313a28-3e93-4233-96aa-226c13c9bfdb
function estimate_node_state!(gA, gB, gN, lifeA, lifeB, switch_tag, idx)
    if switch_tag
        # 刚刚坏掉的是switchA[idx]，需要重新计算当前的节点状态
        if gA[idx] == 0  # 这条语句不会被执行
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
    else
        if gB[idx] == 0
            # 什么也不做，因为这条语句不可能被执行
        elseif gB[idx] == 1
            lifeA[idx] != +Inf && (gN[idx] == 3; return nothing)
            gA[idx] == 1 && (gN[idx] = 5; return nothing)
            gA[idx] == 2 && (gN[idx] = 3; return nothing)
            gA[idx] == 3 && (gN[idx] = 4; return nothing)
        elseif gB[idx] == 2
            lifeA[idx] != +Inf && (gN[idx] == 1; return nothing)
            gA[idx] == 1 && (gN[idx] = 1; return nothing)
            gA[idx] == 2 && (gN[idx] = 4; return nothing)
            gA[idx] == 3 && (gN[idx] = 4; return nothing)
        end
    end
    nothing
end


# ╔═╡ 137b8727-37aa-4909-9308-1785d84361fd
function initialize!(gA, gB, gN)
    fill!(gA, 0)  # 状态全清0，代表正常
    fill!(gB, 0)
    # 计算节点状态
    estimate_node_state!(gA, gB, gN)
    # 随机选取1个节点为主节点
    master_node::Int8 = rand(1:10)
    return master_node
end


# ╔═╡ ddd2e917-2bd4-4a88-a670-c37c82265dd4
function update_node_state!(gA, gB, gN, master_node) # 2 当准备开始写了，发现有需要的前置条件，
    estimate_switch_state!(gA, gB)
    estimate_node_state!(gA, gB, gN)
    # 在这两步骤之后，应该是已经可以知道系统是否能正常运作了
    # 问题在于，用什么来衡量系统时钟出错了
    reselect_master!(gN, master_node)
    nothing
end


# ╔═╡ 379c5bd4-d2ee-48e4-bbdb-b2b7c5f1619e
# 我发现文本的解释好像有点问题，或者是我的理解有误

# > 有且仅有一个节点处于 gMO（注：按前文描述的系统工作机制，该节点必然担当主节点，虽然因为随机因素，过程可能曲折）

# 按理说，MO 节点应当是始终处于主模式的，无法退出到从模式，所以必然直接占据了主模式的信号，从而不需要重选，重选的情况应该只会发生在，主模式因故失效的情况下
function simulate!(gA, gB, gN)
    master_node::Int8 = initialize!(gA, gB, gN)
    Gsys::Int8 = 0  # 初始化后的系统可以正常工作
    state_system = false
    life_counter::Int32 = 0
    QPF::Int8 = QSO::Int8 = QDM::Int8 = QMO::Int8 = QDN::Int8 = QFB::Int8 = 0
    while !state_system && life_counter < LIFE_LIMIT
        update_node_state!(gA, gB, gN, master_node)
        # 计算系统的状态，我是先写的这一部分，没有按着顺序从头写到尾
        QPF = QSO = QDM = QMO = QDN = QFB = 0
        @inbounds for elem in gN
            elem == 0 && (QPF += 1; continue)
            elem == 1 && (QSO += 1; continue)
            elem == 2 && (QDM += 1; continue)
            elem == 3 && (QMO += 1; continue)
            elem == 4 && (QDN += 1; continue)
            elem == 5 && (QFB += 1; continue)
        end
        # 这之前 Gsys=0，初始化0不是一个好的选择，而且以下判断可能需要放到前面去
        if QFB >= 1 || QMO >= 2 || QPF + QMO + QDM == 0 || (QPF + QSO + 1((QMO + QDM) > 0)) < k
            Gsys = 1
        elseif QFB == 0 && ((QMO == 1 && QPF + QSO >= k - 1) || ((QMO == 0 && QPF >= 1 && QPF + QSO >= k) || (QMO == 0 && QPF == 0 && QDM >= 1 && QSO >= k - 1)))
            Gsys = 2
            # 这里的条件文本没有给清楚，但大致能猜到C5 && (C6 || C7)
        elseif QFB + QMO == 0 && (QPF >= 1 && QPF + QSO == k - 1 && QDM >= 1)
            # Todo: 明确节点的选择是在哪里发生的
            # Todo: 明确指导书中写在这里的条件概率有什么用，我是在做仿真而不是理论分析
            if gN[master_node] == 2
                Gsys = 3  # if one of gDM is selected as master
            elseif gN[master_node] == 0
                Gsys = 4  # if one of gPF is selected as master, fail to satisfy valid node limit k
            end
        end
        if Gsys == 2 || Gsys == 3
            life_counter += 1
        else
            life_counter += 1
            state_system = true
        end
    end
    @printf("\ngA:")
    for i = 1:NUM_NODE
        @printf("%3d", gA[i])
    end
    @printf("\ngB:")
    for i = 1:NUM_NODE
        @printf("%3d", gB[i])
    end
    @printf("\ngN:")
    for i = 1:NUM_NODE
        @printf("%3d", gN[i])
    end
    @printf("\nGsys:%3d\tQPF%3d\tQSO%3d\tQDM%3d\tQMO%3d\tQDN%3d\tQFB%3d\n", Gsys, QPF, QSO, QDM, QMO, QDN, QFB)
    life_counter
end


# ╔═╡ 6ce1bbb7-51a5-41b1-92e9-347533b0c64b
function julia_main()
    # variables definition for each simulation
    gA = zeros(Int8, NUM_NODE)
    gB = zeros(Int8, NUM_NODE)
    gN = zeros(Int8, NUM_NODE)
    system_life = zeros(Int32, NUM_SYSTEM)
    # run simulation
    @time @inbounds for i = 1:NUM_SYSTEM
        @printf("\nsystem number:%8d\n", i)
        system_life[i] = simulate!(gA, gB, gN)
        @printf("system life:%8d\n", system_life[i])
    end
    bugval = count(i -> (i == LIFE_LIMIT), system_life)
    writedlm("out/csv/sys$NUM_SYSTEM-lim$LIFE_LIMIT-fixed-has-bugval.csv", system_life, ',')
    @printf("Invalid simulations: %8d (%5.2f%%)\n", bugval, bugval * 100.0 / NUM_SYSTEM)
    deleteat!(system_life, system_life .== LIFE_LIMIT)
    writedlm("out/csv/sys$NUM_SYSTEM-lim$LIFE_LIMIT-fixed-no-bugval.csv", system_life, ',')
    # readdlm("out/csv/FileName.csv",',',Int32)
    @printf("Avg Life: %12.6f\n", mean(system_life))
    p = histogram(system_life, bins=min(NUM_SYSTEM, 256), label="$(NUM_SYSTEM-bugval) valid samples")
    p = histogram!(legend=:topright, bar_edges=true)
    p = histogram!(title="Fixed Time Step Simulation")
    savefig(p, "out/figure/sys$NUM_SYSTEM-lim$LIFE_LIMIT-fixed.png")
end

# ╔═╡ 0fa93817-eb9b-4699-aa3b-a0d2e75d6573
md"#### Variable time step implementation
Two approaches are suitable here for implement such behavior:
- sample all life time once, and set following function code should be rewrite to fit current need.
- sample 1 life time each time, the following function code will not change significantly, but initialization process would set all to `Inf` instead of using `zeros()`"

# ╔═╡ 7e5a32c8-da2d-4e4f-8a72-cf9012b922cb
function simulate_variabletimestep!(gA, gB, gN, lifeA, lifeB)
    master_node::Int8 = initialize!(gA, gB, gN)
    Gsys::Int8 = 0  # 初始化后的系统可以正常工作
    state_system = false
    life_counter::Float64 = 0
    QPF::Int8 = QSO::Int8 = QDM::Int8 = QMO::Int8 = QDN::Int8 = QFB::Int8 = 0
    estimate_switch_state!(gA, gB, lifeA, lifeB)
    @inbounds for i = 1:2*NUM_NODE
        minA, idxA = findmin(lifeA)
        minB, idxB = findmin(lifeB)
        min_life = min(minA, minB) # 按理说，这里我应该已经知道是哪个switch坏了

        min_life > LIFE_LIMIT && break

        switch_tag = (min_life == minA) # switch_tag=true => minA, switch_tag=false =>minB
        idx = switch_tag ? idxA : idxB
        estimate_node_state!(gA, gB, gN, lifeA, lifeB, switch_tag, idx)
        if switch_tag
            lifeA[idxA] = +Inf
        else
            lifeB[idxB] = +Inf
        end
        # start from here, 2 methods for variable time step
        reselect_master!(gN, master_node)
        # 计算系统的状态，我是先写的这一部分，没有按着顺序从头写到尾
        QPF = QSO = QDM = QMO = QDN = QFB = 0
        @inbounds for elem in gN
            elem == 0 && (QPF += 1; continue)
            elem == 1 && (QSO += 1; continue)
            elem == 2 && (QDM += 1; continue)
            elem == 3 && (QMO += 1; continue)
            elem == 4 && (QDN += 1; continue)
            elem == 5 && (QFB += 1; continue)
        end
        # 这之前 Gsys=0，初始化0不是一个好的选择，而且以下判断可能需要放到前面去
        if QFB >= 1 || QMO >= 2 || QPF + QMO + QDM == 0 || (QPF + QSO + 1((QMO + QDM) > 0)) < k
            Gsys = 1
        elseif QFB == 0 && ((QMO == 1 && QPF + QSO >= k - 1) || ((QMO == 0 && QPF >= 1 && QPF + QSO >= k) || (QMO == 0 && QPF == 0 && QDM >= 1 && QSO >= k - 1)))
            Gsys = 2
            # 这里的条件文本没有给清楚，但大致能猜到C5 && (C6 || C7)
        elseif QFB + QMO == 0 && (QPF >= 1 && QPF + QSO == k - 1 && QDM >= 1)
            if gN[master_node] == 2
                Gsys = 3  # if one of gDM is selected as master
            elseif gN[master_node] == 0
                Gsys = 4  # if one of gPF is selected as master, fail to satisfy valid node limit k
            end
        end
        if Gsys == 2 || Gsys == 3
            life_counter = max(min_life, life_counter)
        else
            life_counter = max(min_life, life_counter)
            state_system = true
        end
    end
    # @printf("\ngA:")
    # for i = 1:NUM_NODE
    #     @printf("%3d", gA[i])
    # end
    # @printf("\ngB:")
    # for i = 1:NUM_NODE
    #     @printf("%3d", gB[i])
    # end
    # @printf("\ngN:")
    # for i = 1:NUM_NODE
    #     @printf("%3d", gN[i])
    # end
    # @printf("\nGsys:%3d\tQPF%3d\tQSO%3d\tQDM%3d\tQMO%3d\tQDN%3d\tQFB%3d\n", Gsys, QPF, QSO, QDM, QMO, QDN, QFB)
    # @printf("%d,%d,%d,%d,%d,%d\n", QPF, QSO, QDM, QMO, QDN, QFB)
    life_counter
end

# ╔═╡ 680b3e22-6d04-41b1-83ec-2e6754d6016f
function julia_main_new()
    # variables definition for each simulation
    gA = zeros(Int8, NUM_NODE)
    gB = zeros(Int8, NUM_NODE)
    gN = zeros(Int8, NUM_NODE)
    system_life = zeros(Float64, NUM_SYSTEM)
    lifeA = zeros(NUM_NODE)
    lifeB = zeros(NUM_NODE)
    # run simulation
    @time @inbounds for i = 1:NUM_SYSTEM
        # @printf("\nsystem number:%8d\n", i)
        system_life[i] = simulate_variabletimestep!(gA, gB, gN, lifeA, lifeB)
        # @printf("system life:%16.4f\n", system_life[i])
    end
    # bugval = count(i -> (i >= LIFE_LIMIT), system_life)
    # writedlm("out/csv/sys$NUM_SYSTEM-lim$LIFE_LIMIT-variable-has-bugval.csv", system_life, ',')
    # @printf("Invalid simulations: %8d (%5.2f%%)\n", bugval, bugval * 100.0 / NUM_SYSTEM)
    # deleteat!(system_life, system_life .>= LIFE_LIMIT)
    # writedlm("out/csv/sys$NUM_SYSTEM-lim$LIFE_LIMIT-variable-no-bugval.csv", system_life, ',')
    # readdlm("out/csv/FileName.csv",',',Int32)
    @printf("Avg Life: %12.6f\n", mean(system_life))
    # p = histogram(system_life, bins=min(NUM_SYSTEM, 256), label="$(NUM_SYSTEM-bugval) valid samples")
    p = histogram(system_life, bins=min(NUM_SYSTEM, 256), label="$(NUM_SYSTEM) samples")
    # p = histogram(system_life, bins=NUM_SYSTEM, label="$(NUM_SYSTEM) samples")
    p = histogram!(legend=:topright, bar_edges=true)
    p = histogram!(title="Variable Time Step Simulation")
    savefig(p, "out/figure/sys$NUM_SYSTEM-lim$LIFE_LIMIT-variable.png")
end

@time julia_main_new()