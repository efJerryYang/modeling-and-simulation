### A Pluto.jl notebook ###
# v0.17.5

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

# ╔═╡ 86e5acd8-7dac-4cae-8e58-3fe80d21ec09
md"## 时间按变化步长推进
变化步长最主要的思想就是，从仿真一开始就通过从分布取样(sample from the distribution)，确定各个元件的使用寿命，从而将仿真过程的内层循环次数减少。

具体而言，就是在时间轴上首先取样(sample)出切换器A、B的故障发生时间，充分利用故障只能发生一次的特点，每次内循环只需要按照发生故障的先后顺序更新节点的状态即可。

因固定步长的代码略去不讲，这里不再详细对比二者的区别，仅在必要的时候引用固定步长代码以帮助理解。
"

# ╔═╡ a995f5d3-4a32-41ea-913f-f551af205197
md"## 案例2代码实现
本部分将详细解释各部分代码的具体含义
"

# ╔═╡ 93f9ab9b-35f9-48ea-ac96-40e2f525cb35
md"**注意**：在`julia`语言中, `begin-end`语句并不会引入一个新的命名空间。(命名空间为namespace，但原本的用词是scope block，这里为方便理解和查询相关词汇含义而改写用词，有不恰当之处)"

# ╔═╡ 359d6bd1-78c3-43e6-8555-2d811903eed9
md"### 导入标准库
各个标准库主要使用到的函数已经列在行末注释，无需特别关注本部分，函数名均较为直观，容易理解。"

# ╔═╡ cea03e40-0119-476e-92c5-7742e675b5ea
md"### 常量参数
本部分为所需的常量参数，命名与指导书中一致。

**注意**：使用了部分语言不支持的Unicode变量名λA,λB。"

# ╔═╡ c78092ff-73c2-40c2-b8e2-3f1178db4152
begin
    # constant parameters
    ## system 
    const NUM_SYSTEM = 10_0000
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
md"### 函数定义
函数末尾带有`!`的含义是该函数将参数的引用传入内部，函数内部修改参数值会直接作用到原变量，以及不拷贝参数副本。这样做是为了提高运行效率。

**注意**：非常不建议按顺序逐个阅读函数，推荐的方式是打开2份本文档，跟随`julia_main_varia()`和`simulate_variable_timestep!()`的代码执行过程，查看对应函数的具体实现。"

# ╔═╡ 908b0932-85a8-4463-95e8-d9850b4bf3d2
md"#### 初始化运行
本函数清空上一轮的仿真状态，用于新一轮仿真的启动"

# ╔═╡ 137b8727-37aa-4909-9308-1785d84361fd
function initialize!(NUM_NODE, gA, gB, gN, gR)
    fill!(gA, 0)
    fill!(gB, 0)

    fill!(gN, 0)
    fill!(gR, 0)

    master_node::Int8 = 1#rand(1:NUM_NODE)
    gR[master_node] = 1
    return master_node
end


# ╔═╡ 5ec48cc2-cb7d-4e42-b1cf-56afacccc768
function estimate_switch_state!(NUM_NODE, gA, gB, lifeA, lifeB)
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

# ╔═╡ 190eed91-a3c6-440f-b042-e80c70b1bc3f
md"#### Estimate Node Performance State
The code here is not obvious, you may refer to the `fixed-timestep` version of this part for better understanding."

# ╔═╡ b2313a28-3e93-4233-96aa-226c13c9bfdb
function estimate_node_performance_state!(NUM_NODE, gA, gB, gN, lifeA, lifeB, switch_tag, idx)
    if switch_tag
        if gA[idx] == 1
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
        if gB[idx] == 1
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

# ╔═╡ a554951b-5dc3-40a2-90bc-d08164c63fc7
function ok_for_master(gNi)
    gNi == 0 && return true
    gNi == 1 && return false
    gNi == 2 && return true
    gNi == 3 && return true # MO
    gNi == 4 && return false
    gNi == 5 && return false
end

# ╔═╡ 1c0f6a74-d178-4fe9-81cf-e936018d827d
md"#### Available Role States
Considering the role states to be chosen:
- `0` slave node
- `1` master node
- `2` disable node
- `3` block node
Both slave nodes and master node are connected to bus, and only disable nodes are disconnected. The block node may block bus or be disable, it depends on

The classes here (slave, master, disable) are so heuristic that I cannot guarantee the correctness. You may want to discuss with your friends or teacher for a better classification strategy.

**Warning**: not sure of state coverage
"

# ╔═╡ a2592eff-4763-4a8f-ac38-46b1b34bb0fa
md"#### Role State Transition Equation
```math
\{G_{role}^{(m)}\}=F_{G}(\{G_{role}^{(m-1)}\},\{G_{N}^{(m)}\},[random])
```
"

# ╔═╡ e44a2cf9-9fd3-44e6-828c-ba969ff632e1
md"#### Reselect Master
The master reselection should based on both node preformance state and role state.

Here I discard the `reselect_master!()` function, and integrate its code into `estimate_role_state!()` for better compatibility considering the possible use of `gR` vector."

# ╔═╡ 2af9a147-3785-4ff2-a940-4079be53d6cd
function estimate_node_role_state!(NUM_NODE, gN, gR, master_node)
    # role state transition process is based on the formula given above
    # role vector: gR
    # node (performance) vector: gN
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
        mintime_index == 0 && return master_node
        master_node = mintime_index
    end
    cnt = 0
    idx = 0
    # flag = false # flag for FB
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
    master_node

end

# ╔═╡ 029505db-edcb-46e8-accd-a0a20439efbe
function estimate_system_state!(NUM_NODE, gN, master_node)
    QPF::Int8 = QSO::Int8 = QDM::Int8 = QMO::Int8 = QDN::Int8 = QFB::Int8 = 0
    Gsys::Int8 = 2
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
        cond = QDM / (QDM + QPF)
        if rand() < cond
            Gsys = 3
        else
            Gsys = 4
        end
    end
    return Gsys
end

# ╔═╡ 926fdff5-3cec-4bc3-9937-f81f71dd25c7
md"#### 仿真执行
单次仿真的完整执行过程
- `initialize!()`
- `estimate_switch_state!()`
- loop: iter switches, select current `min_life`
  - find `min_life`, find its `index` and switch type `X`
  - if `min_life` > `LIFE_LIMIT`, break loop, stop simulation, return previous `life_counter`
  - else, `estimate_node_performance_state!()`
  - set corresponding `lifeX[index]` to `+Inf`, avoid repeat selection in the following simulation iterations
  - `estimate_node_role_state!()`, include the master reselection process
  - `estimate_system_state!()`
  - update `life_counter`
- return `life_counter`
"

# ╔═╡ 7e5a32c8-da2d-4e4f-8a72-cf9012b922cb
function simulate_variable_timestep!(NUM_NODE, gA, gB, gN, gR, lifeA, lifeB)
    master_node::Int8 = initialize!(NUM_NODE, gA, gB, gN, gR)
    state_system = false
    life_counter::Float64 = 0
    estimate_switch_state!(NUM_NODE, gA, gB, lifeA, lifeB)
    @inbounds for i = 1:2*NUM_NODE
        minA, idxA = findmin(lifeA)
        minB, idxB = findmin(lifeB)
        min_life = min(minA, minB)

        min_life > LIFE_LIMIT && (life_counter = LIFE_LIMIT; break)

        switch_tag, idx = (min_life == minA) ? (true, idxA) : (false, idxB)
        estimate_node_performance_state!(NUM_NODE, gA, gB, gN, lifeA, lifeB, switch_tag, idx)
        if switch_tag
            lifeA[idxA] = +Inf
        else
            lifeB[idxB] = +Inf
        end
        Gsys::Int8 = estimate_system_state!(NUM_NODE, gN, master_node)
        if Gsys == 2 || Gsys == 3
            life_counter = min_life
        else
            break
        end
    end
    life_counter
end

# ╔═╡ 6a2c4add-cf6e-41f5-abb3-a5cc2470c7d4
md"#### Simulation Executor

`julia_main_varia()` is not exactly the **real main** function. 

Typically, a function called `real_main()` would possibly be chosen as the **real main** function, but here for convenience, I just run the code in global scope out of a `real_main()` function."

# ╔═╡ 680b3e22-6d04-41b1-83ec-2e6754d6016f
function julia_main_varia(NUM_NODE::Int8, avg_life_max, avg_life_idx, reliability_max, reliability_idx)
    gA = zeros(Int8, NUM_NODE)
    gB = zeros(Int8, NUM_NODE)
    gN = zeros(Int8, NUM_NODE)
    gR = zeros(Int8, NUM_NODE)
    system_life = zeros(Float64, NUM_SYSTEM)
    lifeA = zeros(NUM_NODE)
    lifeB = zeros(NUM_NODE)
    reliability_counter = 0

    @inbounds for i = 1:NUM_SYSTEM
        # @printf("\nsystem number:%8d\n", i)
        system_life[i] = simulate_variable_timestep!(NUM_NODE, gA, gB, gN, gR, lifeA, lifeB)
        system_life[i] > w && (reliability_counter += 1)
        # @printf("system life:%16.4f\n", system_life[i])
    end

    avg_life = mean(system_life)
    reliability = reliability_counter / NUM_SYSTEM
    @printf("NUM_NODE:%3d\tAvg: %12.6f\tReliability: %7.3f%%\n", NUM_NODE, avg_life, reliability * 100)
    avg_life_max = max(avg_life_max, avg_life)
    avg_life_max == avg_life && (avg_life_idx = NUM_NODE)
    reliability_max = max(reliability_max, reliability)
    reliability_max == reliability && (reliability_idx = NUM_NODE)
    return avg_life_max, avg_life_idx, reliability_max, reliability_idx
end

# ╔═╡ 54ab68ae-a358-43c2-91b2-618178ea227b
md"### Execute with Fixed NUM_NODE
The fixed `NUM_NODE` is 10, simply for testing the correctness."

# ╔═╡ fb90255d-b546-4bfd-99f3-a5e3d0f529c1
# avg_life_max = 0.0
# avg_life_idx = 0
# reliability_max = 0.0
# reliability_idx = 0
# avg_life_max, avg_life_idx, reliability_max, reliability_idx = julia_main_varia(convert(Int8, 10), avg_life_max, avg_life_idx, reliability_max, reliability_idx)

# ╔═╡ d3d04a6f-bed6-413d-9f44-7547f7487765
md"### Execute with A Range of NUM_NODE
This part is designated for collecting problem results, including the most average MTTF, and the greatest reliability.

**Warning**: The results may possibly not correct, some modifications are needed within the function `estimate_node_role_state!()`. 

**Notice**: Currently, the function `estimate_node_role_state!()` is just a rename of the discarded `reselect_master!()` function, and you may need to add code for solving the node role transition problem according to the *Role State Transition Equation* given above."

# ╔═╡ 0662362e-2115-4158-b7b9-81802b00d8cf
function real_main()
    avg_life_max = 0.0
    avg_life_idx = 0
    reliability_max = 0.0
    reliability_idx = 0
    @time for i::Int8 = 3:20
        avg_life_max, avg_life_idx, reliability_max, reliability_idx = julia_main_varia(i, avg_life_max, avg_life_idx, reliability_max, reliability_idx)
    end

    # ╔═╡ 53eea4fb-418f-4bcb-bb51-45a84e55eb76
    println("MTTF: $avg_life_idx\t$avg_life_max")

    # ╔═╡ 37f916dd-7059-45fc-80d0-28b9caf1ecfe
    println("R(w): $reliability_idx\t$(reliability_max*100)%")
end
real_main()