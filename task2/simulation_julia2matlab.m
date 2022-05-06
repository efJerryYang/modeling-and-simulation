main();

function main()
    global NUM_SYSTEM TIME_STEP LIFE_LIMIT STATE_NUM_NODE k w
    global lambda_A PA0 PEA1 PEA2 PEA3
    global lambda_B PB0 PEB1 PEB2

    NUM_SYSTEM = 100000;
    TIME_STEP = 1;
    LIFE_LIMIT = 200000;
    STATE_NUM_NODE = 6;
    k = 3;
    w = 30000;
    lambda_A = 1/5.90e4;
    PA0 = exp(-lambda_A * TIME_STEP);
    PEA1 = 0.20 * (1 - PA0);
    PEA2 = 0.15 * (1 - PA0);
    PEA3 = 0.65 * (1 - PA0);
    lambda_B = 1/2.20e5;
    PB0 = exp(-lambda_B * TIME_STEP);
    PEB1 = 0.45 * (1 - PB0);
    PEB2 = 0.55 * (1 - PB0);

    avg_life_max = 0.0;
    avg_life_idx = 0;
    reliability_max = 0.0;
    reliability_idx = 0;

    for i = 3:20
        [avg_life_max, avg_life_idx, reliability_max, reliability_idx] = julia_main_varia(i, avg_life_max, avg_life_idx, reliability_max, reliability_idx);
    end

    fprintf("MTTF: %3d%16.4f\n", avg_life_idx, avg_life_max);
    fprintf("R(w): %3d%15.2f%%\n", reliability_idx, reliability_max * 100);
end

function [avg_life_max, avg_life_idx, reliability_max, reliability_idx] = julia_main_varia(NUM_NODE, avg_life_max, avg_life_idx, reliability_max, reliability_idx)
    system_life = zeros(100000, 1);
    % lifeA = zeros(NUM_NODE);
    % lifeB = zeros(NUM_NODE);
    reliability_counter = 0;
    w = 30000;

    for i = 1:100000
        system_life(i) = simulate_variable_timestep(NUM_NODE);

        if system_life(i) >= w
            reliability_counter = reliability_counter + 1;
        end

    end

    avg_life = mean(system_life);
    reliability = reliability_counter / 100000;
    fprintf("NUM_NODE:%3d\tMTTF: %12.6f\tReliability: %7.3f%%\n", NUM_NODE, avg_life, reliability * 100)
    avg_life_max = max(avg_life_max, avg_life);

    if avg_life_max == avg_life
        avg_life_idx = NUM_NODE;
    end
    
    reliability_max = max(reliability_max, reliability);

    if reliability_max == reliability
        reliability_idx = NUM_NODE;
    end
end

function [life_counter] = simulate_variable_timestep(NUM_NODE)
    % [gA, gB, gN] = initialize(gA, gB, gN);
    gA = zeros(NUM_NODE, 1);
    gB = zeros(NUM_NODE, 1);
    gN = zeros(NUM_NODE, 1);
    life_counter = 0;
    [gA, gB, lifeA, lifeB] = estimate_switch_state(NUM_NODE, gA, gB);

    for i = 1:2 * NUM_NODE
        minA = min(lifeA);
        idxA = find(lifeA == minA);
        minB = min(lifeB);
        idxB = find(lifeB == minB);
        min_life = min(minA, minB);

        if min_life > 200000
            life_counter = 200000;
            break;
        end

        switch_tag = true;

        if min_life == minA
            switch_tag = true;
            idx = idxA;
        else
            switch_tag = false;
            idx = idxB;
        end

        gN = estimate_node_performance_state(gA, gB, gN, lifeA, lifeB, switch_tag, idx);

        if switch_tag
            lifeA(idxA) = +inf;
        else
            lifeB(idxB) = +inf;
        end

        Gsys = estimate_system_state(gN);

        if Gsys == 2 || Gsys == 3
            life_counter = min_life;
        else
            life_counter = min_life;
            break;
        end

    end

end

function [gA, gB, gN] = initialize(gA, gB, gN)
    gA = gA * 0;
    gB = gB * 0;
    gN = gN * 0;
end

function [gA, gB, lifeA, lifeB] = estimate_switch_state(NUM_NODE, gA, gB)
    lambda_A = 1/5.90e4;
    TIME_STEP = 1;
    PA0 = exp(-lambda_A * TIME_STEP);
    PEA1 = 0.20 * (1 - PA0);
    PEA2 = 0.15 * (1 - PA0);
    PEA3 = 0.65 * (1 - PA0);
    lambda_B = 1/2.20e5;
    PB0 = exp(-lambda_B * TIME_STEP);
    PEB1 = 0.45 * (1 - PB0);
    PEB2 = 0.55 * (1 - PB0);

    lifeA = exprnd(1 / lambda_A, NUM_NODE, 1);
    lifeB = exprnd(1 / lambda_B, NUM_NODE, 1);

    for i = 1:NUM_NODE
        tolA = rand * (1 - PA0);

        if tolA < PEA1
            gA(i) = 1;
        elseif tolA < PEA1 + PEA2
            gA(i) = 2;
        else
            gA(i) = 3;
        end

        tolB = rand * (1 - PB0);

        if tolB < PEB1
            gB(i) = 1;
        else
            gB(i) = 2;
        end

    end

end

function [gN] = estimate_node_performance_state(gA, gB, gN, lifeA, lifeB, switch_tag, idx)

    if switch_tag

        if gA(idx) == 1

            if lifeB(idx) ~= +inf
                gN(idx) = 1;
            elseif gB(idx) == 1
                gN(idx) = 5;
            elseif gB(idx) == 2
                gN(idx) = 1;
            end

        elseif gA(idx) == 2

            if lifeB(idx) ~= +inf
                gN(idx) = 2;
            elseif gB(idx) == 1
                gN(idx) = 3;
            elseif gB(idx) == 2
                gN(idx) = 4;
            end

        elseif gA(idx) == 3

            if lifeB(idx) ~= +inf
                gN(idx) = 4;
            elseif gB(idx) == 1
                gN(idx) = 4;
            elseif gB(idx) == 2
                gN(idx) = 4;
            end

        end

    else

        if gB(idx) == 1

            if lifeA(idx) ~= +inf
                gN(idx) = 3;
            elseif gA(idx) == 1
                gN(idx) = 5;
            elseif gA(idx) == 2
                gN(idx) = 3;
            elseif gA(idx) == 3
                gN(idx) = 4;
            end

        elseif gB(idx) == 2

            if lifeA(idx) ~= +inf
                gN(idx) = 1;
            elseif gA(idx) == 1
                gN(idx) = 1;
            elseif gA(idx) == 2
                gN(idx) = 4;
            elseif gA(idx) == 3
                gN(idx) = 4;
            end

        end

    end

end

function [Gsys] = estimate_system_state(gN)
    QPF = sum(gN(:) == 0);
    QSO = sum(gN(:) == 1);
    QDM = sum(gN(:) == 2);
    QMO = sum(gN(:) == 3);
    QDN = sum(gN(:) == 4);
    QFB = sum(gN(:) == 5);
    k = 3;
    C1 = (QFB >= 1);
    C2 = (QMO >= 2);
    C3 = (QPF + QMO + QDM == 0);
    C4 = ((QPF + QSO + ((QMO + QDM) > 0)) < k);
    C5 = (QFB == 0);
    C6 = (QMO == 1 && QPF + QSO >= k - 1);
    C7 = ((QMO == 0 && QPF >= 1 && QPF + QSO >= k) || (QMO == 0 && QPF == 0 && QDM >= 1 && QSO >= k - 1));
    C8 = (QFB + QMO == 0);
    C9 = (QPF >= 1 && (QPF + QSO == k - 1) && QDM >= 1);

    if C1 || C2 || C3 || C4
        Gsys = 1;
    elseif C5 && (C6 || C7)
        Gsys = 2;
    elseif C8 && C9
        cond = QDM / (QDM + QPF);

        if rand < cond
            Gsys = 3;
        else
            Gsys = 4;
        end

    end

end
