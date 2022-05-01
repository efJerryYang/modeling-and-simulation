using Statistics
using Printf
using LinearAlgebra
using Plots
using DelimitedFiles
# using StaticArrays # 不能使用fill!()方法, is designed to be immutable
# using SharedArrays
# using Distributed

# constant parameters
## system 
# NUM_SYSTEM = 10_0000
const NUM_SYSTEM = 25_0001
const NUM_NODE = 10               # this should also be optimized later
const TIME_STEP = 1               # 1 hour
const LIFE_LIMIT = 60_0000
# const LIFE_LIMIT = 2_000_0000
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


function estimate_switch_state!(gA, gB) # 考虑到节点的状态完全由切换器的状态决定，又回想起来前面写过了什么，然后把最后一块拼图拼上了
    # 各switch有一定概率正常工作，一定概率出现异常
    @inbounds for i = 1:NUM_NODE
        sample = rand()
        if gA[i] != 0 || sample < PA0 #当前切换器已经坏了，不需要计算||仍然正常工作
            continue
        else
            tol = sample - PA0
            gA[i] = tol < PEA1 ? 1 : tol < PEA1 + PEA2 ? 2 : 3
        end
        # sample = rand() # 为什么这里出现了问题，原本几乎无法产生合理的可以使B改变的rand()
        if gB[i] != 0 || sample < PB0
            continue
        else
            gB[i] = sample < PB0 + PEB1 ? 1 : 2
        end
        # Todo: 以上代码应该可以改成位运算
    end
    nothing
end
function estimate_node_state!(gA, gB, gN)
    # 这之前已经修改过switch state了
    @inbounds for i = 1:NUM_NODE
        if gA[i] == 0
            # 为了好的观感稍稍的降低了效率，相当于if并列而没有elseif
            gB[i] == 0 && (gN[i] = 0) # dont change masterflag
            gB[i] == 1 && (gN[i] = 3)
            gB[i] == 2 && (gN[i] = 1) #, (master_flag[i] = false))
        elseif gA[i] == 1
            gB[i] == 0 && (gN[i] = 1) #, master_flag[i]=false)
            gB[i] == 1 && (gN[i] = 5)
            gB[i] == 2 && (gN[i] = 1) #, master_flag[i]=false)
        elseif gA[i] == 2
            gB[i] == 0 && (gN[i] = 2)
            gB[i] == 1 && (gN[i] = 3)
            gB[i] == 2 && (gN[i] = 4)
        elseif gA[i] == 3
            gB[i] == 0 && (gN[i] = 4)
            gB[i] == 1 && (gN[i] = 4)
            gB[i] == 2 && (gN[i] = 4)
        end
    end
    nothing
end
function initialize!(gA, gB, gN, master_flag, master_node)
    # 最开始，不妨设第1个节点为主节点，其它节点为从节点
    # Todo: 这是可以修改的，也可以用随机算法，但意义不大
    fill!(gA, 0)  # 状态全清0，代表正常
    fill!(gB, 0)

    # 设置第1个节点为主节点
    fill!(master_flag, false) # 这里应该是可以优化的，我不需要散列，因为只有一个master
    master_node = 1
    master_flag[master_node] = 1
    # 计算节点状态
    estimate_node_state!(gA, gB, gN)
    nothing
end
function count_master_num(master_flag)
    num = sum(master_flag)
end
function update_node_state!(gA, gB, gN, master_flag, master_node) # 2 当准备开始写了，发现有需要的前置条件，
    estimate_switch_state!(gA, gB)
    estimate_node_state!(gA, gB, gN)
    # 在这两步骤之后，应该是已经可以知道系统是否能正常运作了
    # 问题在于，用什么来衡量系统时钟出错了
    reselect_master!(gA, gB, gN, master_flag, master_node)
    nothing
end

#! 这部分应该是整个系统里最麻烦的地方
# 注意，从描述中可知，如果系统处于正常工作状态，哪怕有异常节点，也不会触发主
# 节点重选，重选主节点一定发生在系统处于异常时，但问题在于，这个时候不应该认
# 为系统失效，所以不能完全应用Gsys的状态（否则就要统计两次Gsys的状态）
# 除此以外，这里的“戒备状态”事实上并不是真的时间，而是随机产生一组随机数，选
# 择最小数下标的作为重选的主节点很容易实现： 
# _, new_master_index = findmin(rand(NUM_NODE))
# 这里再根据具体的要求改改就可以了
# 问题在于，如何判断总线时钟失效，以及如果当前有多个主节点怎么办？
# 先停这里，写update_node_state
# if clock_invalid
#     master_flag
# end
function ok_for_master(gNi)
    gNi == 0 && return true
    gNi == 1 && return false
    gNi == 2 && return true
    gNi == 3 && return true # MO
    gNi == 4 && return false
    gNi == 5 && return false
end
# 节点重选 
# 这里节点重选的地方还没能想清楚
# 从最初始的地方开始，是主节点在 1 号，然后剩余节点都是从节点
# 
# 在 update_node_state 之后，由于对于 switch 状态的重新估计，导致了 node 的状态可能发生改变。
# 
# 这时候，需要对接下来可能发生变化的系统状态进行讨论，
# 
# 我们要进行节点重选，是因为更新节点状态之后可能会导致系统出现异常，那么，这里能不能套用后面计算系统状态的结论？
# 如果不行的话，这个分析过程将会变得相当复杂，等价于分析后面系统状态变化的工作量。
# 
# 重选节点的算法本身很简单，但是如何选择可以用于重选的节点？这个限制条件要怎么加？
# 
# 现在我已经有了什么：所有节点的状态，以及一个最初设定的主节点 1
# 
# 现在遇到的困难是，如果主节点不是 1 了，我要如何找出所选的主节点，以及当存在了多个主节点的时候，要如何选择去除多余的主节点。
# 
# 那么，现在为了简化问题，我先强行要求只能有一个主节点，否则系统异常，查找可以去除的主节点，并完成节点重选
# 
# 我认为，主节点可以从 MO 里面添加，但是不能从节点更新就开始就添加会比较好处理
# 即，在节点更新函数里面，可以使主节点个数减少到 0，但不能主动设置主节点

function reselect_master!(gA, gB, gN, master_flag, master_node)  # 1
    # index = findfirst(master_flag)
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
        master_flag[master_node] = 0
        master_node = mintime_index  #! 对于没有找到合适的节点如何处理，当时没有考虑
        master_flag[master_node] = 1  #! bug: master_flag[1, 1, 0, 0, 0, 0, 1, 0, 0, 1], 没有设置为0，只在设置为1
    end
    # else # 把这个else删掉后，就可以把两种情况放在一起用了
    # 说明原来的主节点还可以用，但信号出现异常，这个时候要么是MO，要么是FB
    # 我们只需要讨论一种情况，如果系统有一个MO出现
    cnt = 0
    idx = 0
    @inbounds for i = 1:NUM_NODE
        gN[i] == 3 && (cnt=cnt + 1, idx=i)
    end
    if cnt == 1
        master_flag[idx] = 1
        master_flag[master_node] = 0
        master_node = idx
    elseif cnt >= 2
        master_flag[master_node] = 0
        master_node = 0 # 这里随意设置，因为系统必然失效
    elseif cnt == 0 && (5 in gN)
        master_flag[master_node] = 0
        master_node = 0 # 这里随意设置，因为系统必然失效
    else
        # 也就是说现在信号是好的，不需要重选
    end
    # end
    nothing
end
# 我发现文本的解释好像有点问题，或者是我的理解有误

# > 有且仅有一个节点处于 gMO（注：按前文描述的系统工作机制，该节点必然担当主节点，虽然因为随机因素，过程可能曲折）

# 按理说，MO 节点应当是始终处于主模式的，无法退出到从模式，所以必然直接占据了主模式的信号，从而不需要重选，重选的情况应该只会发生在，主模式因故失效的情况下
function simulate!(gA, gB, gN, master_flag)
    master_node::Int8 = 1

    initialize!(gA, gB, gN, master_flag, master_node)
    Gsys::Int8 = 0  # 初始化后的系统可以正常工作
    state_system = false
    life_counter::Int32 = 0
    QPF::Int8 = QSO::Int8 = QDM::Int8 = QMO::Int8 = QDN::Int8 = QFB::Int8 = 0
    while !state_system && life_counter < LIFE_LIMIT
        # 更新各节点的状态
        update_node_state!(gA, gB, gN, master_flag, master_node)
        # 计算系统的状态，我是先写的这一部分，没有按着顺序从头写到尾
        QPF = QSO = QDM = QMO = QDN = QFB = 0
        @inbounds for elem in gN
            elem == 0 && (QPF += 1)
            elem == 1 && (QSO += 1)
            elem == 2 && (QDM += 1)
            elem == 3 && (QMO += 1)
            elem == 4 && (QDN += 1)
            elem == 5 && (QFB += 1)
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

function julia_main()
    # variables definition for each simulation
    gA = zeros(Int8, NUM_NODE)
    gB = zeros(Int8, NUM_NODE)
    gN = zeros(Int8, NUM_NODE)
    master_flag = zeros(Bool, NUM_NODE)
    system_life = zeros(Int32, NUM_SYSTEM)
    # system_life = SharedArray{Int32}(NUM_SYSTEM)
    # Gsys = 0 # valid only during the operation of current system
    # QPF = 0  # QN0 = QPF
    # QSO = 0  # QN1 = QSO
    # QDM = 0  # QN2 = QDM
    # QMO = 0  # QN3 = QMO
    # QDN = 0  # QN4 = QDN
    # QFB = 0  # QN5 = QFB
    # run simulation
    @time @inbounds for i = 1:NUM_SYSTEM
        @printf("\nsystem number:%8d\n", i)
        system_life[i] = simulate!(gA, gB, gN, master_flag)
        @printf("system life:%8d\n", system_life[i])
    end
    bugval = count(i -> (i == LIFE_LIMIT), system_life)
    writedlm("out/csv/$NUM_SYSTEM-systems-life-(limit-$LIFE_LIMIT-hours)-has-bugval.csv", system_life, ',')
    @printf("Invalid simulations: %8d (%5.2f%%)\n", bugval, bugval * 100.0 / NUM_SYSTEM)
    deleteat!(system_life, system_life .== LIFE_LIMIT)
    writedlm("out/csv/$NUM_SYSTEM-systems-life-(limit-$LIFE_LIMIT-hours)-no-bugval.csv", system_life, ',')
    # readdlm("out/csv/FileName.csv",',',Int32)
    @printf("Avg Life: %12.6f\n", mean(system_life))
    p = histogram(system_life, bins=min(NUM_SYSTEM, 128), legend=:topright, title="Monto-Carlo Simulation", label="$(NUM_SYSTEM-bugval) samples")
    display(p)
    savefig(p, "out/figure/$NUM_SYSTEM-systems-(limit-$LIFE_LIMIT-hours).png")
end

# if abspath(PROGRAM_FILE) == @__FILE__ # 会导致debuger失效
@time julia_main()
# end