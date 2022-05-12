main();

function [C] = initializeConstants()
    C.NUM_SYSTEM = 100000;
    C.TIME_STEP = 1;
    C.LIFE_LIMIT = 200000;
    C.STATE_NUM_NODE = 6;
    C.k = 3;
    C.w = 30000;
    C.Lambda_A = 1/5.90e4;
    C.PA0 = exp(-C.Lambda_A * C.TIME_STEP);
    C.PEA1 = 0.20 * (1 - C.PA0);
    C.PEA2 = 0.15 * (1 - C.PA0);
    C.PEA3 = 0.65 * (1 - C.PA0);
    C.Lambda_B = 1/2.20e5;
    C.PB0 = exp(-C.Lambda_B * C.TIME_STEP);
    C.PEB1 = 0.45 * (1 - C.PB0);
    C.PEB2 = 0.55 * (1 - C.PB0);
    C.NUM_NODE = 10;
end

function main()
    C = initializeConstants();
    result.averageLifeIndex = 0;
    result.averageLifeMaxVal = 0.0;
    result.reliabilityIndex = 0;
    result.reliabilityMaxVal = 0.0;

    for i = 3:20
        C.NUM_NODE = i;
        [result] = simulate(C, result);
    end

    fprintf("MTTF: %3d%16.4f\n", result.averageLifeIndex, result.averageLifeMaxVal);
    fprintf("R(w): %3d%15.2f%%\n", result.reliabilityIndex, result.reliabilityMaxVal * 100);
end

function [result] = simulate(C, result)
    systemLife = zeros(C.NUM_SYSTEM, 1);
    reliabilityCounter = 0;

    for i = 1:C.NUM_SYSTEM
        systemLife(i) = simulateVariableTimestep(C);

        if systemLife(i) >= C.w
            reliabilityCounter = reliabilityCounter + 1;
        end

    end

    averageLife = mean(systemLife);
    reliability = reliabilityCounter / C.NUM_SYSTEM;
    fprintf("NUM_NODE:%3d\tMTTF: %12.6f\tReliability: %7.3f%%\n", C.NUM_NODE, averageLife, reliability * 100)

    % update return result
    result = updateResult(C, result, averageLife, reliability);
end

function [result] = updateResult(C, result, averageLife, reliability)

    if max(result.averageLifeMaxVal, averageLife) == averageLife
        result.averageLifeMaxVal = averageLife;
        result.averageLifeIndex = C.NUM_NODE;
    end

    if max(result.reliabilityMaxVal, reliability) == reliability
        result.reliabilityMaxVal = reliability;
        result.reliabilityIndex = C.NUM_NODE;
    end

end

function [lifeCounter] = simulateVariableTimestep(C)
    gA = zeros(C.NUM_NODE, 1);
    gB = zeros(C.NUM_NODE, 1);
    gN = zeros(C.NUM_NODE, 1);
    gR = zeros(C.NUM_NODE, 1); % not used currently

    iMasterNode = 1;
    gR(iMasterNode) = 1;

    lifeCounter = 0;
    [gA, gB, lifeA, lifeB] = computeSwitchState(C, gA, gB);

    for i = 1:2 * C.NUM_NODE
        minA = min(lifeA);
        idxA = find(lifeA(:) == minA);
        minB = min(lifeB);
        idxB = find(lifeB(:) == minB);
        minLife = min(minA, minB);

        if minLife > C.LIFE_LIMIT
            lifeCounter = C.LIFE_LIMIT;
            break;
        end

        % switchTag = true;

        if minLife == minA
            switchTag = true;
            idx = idxA;
        else
            switchTag = false;
            idx = idxB;
        end

        gN = computeNodePerformanceState(gA, gB, gN, lifeA, lifeB, switchTag, idx);

        if switchTag
            lifeA(idxA) = +inf;
        else
            lifeB(idxB) = +inf;
        end

        [gR, iMasterNode] = computeNodeRoleState(C, gN, gR, iMasterNode);

        Gsys = computeSystemState(C, gN, iMasterNode);

        if Gsys == 2 || Gsys == 3
            lifeCounter = minLife;
        else
            lifeCounter = minLife;
            break;
        end

    end

end

function [gA, gB, lifeA, lifeB] = computeSwitchState(C, gA, gB)

    lifeA = exprnd(1 / C.Lambda_A, C.NUM_NODE, 1);
    lifeB = exprnd(1 / C.Lambda_B, C.NUM_NODE, 1);

    for i = 1:C.NUM_NODE
        tolA = rand * (1 - C.PA0);

        if tolA < C.PEA1
            gA(i) = 1;
        elseif tolA < C.PEA1 + C.PEA2
            gA(i) = 2;
        else
            gA(i) = 3;
        end

        tolB = rand * (1 - C.PB0);

        if tolB < C.PEB1
            gB(i) = 1;
        else
            gB(i) = 2;
        end

    end

end

function [gN] = computeNodePerformanceState(gA, gB, gN, lifeA, lifeB, switchTag, idx)

    if switchTag

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

function [isOk] = isOkForMaster(gNi)

    if gNi == 0
        isOk = true; % PF
    elseif gNi == 1
        isOk = false; % SO
    elseif gNi == 2
        isOk = true; % DM
    elseif gNi == 3
        isOk = true; % MO
    elseif gNi == 4
        isOk = false; % DN
    elseif gNi == 5
        isOk = false; % FB
    end

end

function [gR, iMasterNode] = computeNodeRoleState(C, gN, gR, iMasterNode)
    % if previous master node corrupts
    if ~isOkForMaster(gN(iMasterNode))

        alertMinTime = 1.0; % rand range (0, 1)
        minTimeIndex = 0; % assign invalid index by default
        alertCounter = rand(C.NUM_NODE, 1);

        for i = 1:C.NUM_NODE

            if isOkForMaster(gN(i))
                tmp = alertCounter(i);

                if min(alertMinTime, tmp) == tmp
                    alertMinTime = tmp;
                    minTimeIndex = i;
                end

            end

        end

        if minTimeIndex == 0
            return
        end

        iMasterNode = minTimeIndex;

    end

    % now, master node is available
    if sum(gN(:) == 3) == 1
        iMasterNode = find(gN == 3);
        return
    end

    % The other cases will corrupt the system later
    % when computing system state, so we only need to
    % handle the only case above.
end

function [Gsys] = computeSystemState(C, gN, iMasterNode)
    QPF = sum(gN(:) == 0);
    QSO = sum(gN(:) == 1);
    QDM = sum(gN(:) == 2);
    QMO = sum(gN(:) == 3);
    QDN = sum(gN(:) == 4); % not used
    QFB = sum(gN(:) == 5);

    C1 = (QFB >= 1);
    C2 = (QMO >= 2);
    C3 = (QPF + QMO + QDM == 0);
    C4 = ((QPF + QSO + ((QMO + QDM) > 0)) < C.k);
    C5 = (QFB == 0);
    C6 = (QMO == 1 && QPF + QSO >= C.k - 1);
    C7 = ((QMO == 0 && QPF >= 1 && QPF + QSO >= C.k) || (QMO == 0 && QPF == 0 && QDM >= 1 && QSO >= C.k - 1));
    C8 = (QFB + QMO == 0);
    C9 = (QPF >= 1 && (QPF + QSO == C.k - 1) && QDM >= 1);

    if C1 || C2 || C3 || C4
        Gsys = 1;
    elseif C5 && (C6 || C7)
        Gsys = 2;
    elseif C8 && C9
        % cond = QDM / (QDM + QPF);
        % if rand < cond
        if gN(iMasterNode) == 2
            Gsys = 3;
        elseif gN(iMasterNode) == 0
            Gsys = 4;
        end

    end

end
