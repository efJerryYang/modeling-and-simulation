C = initializeConstants();
PPF = C.PA0 * C.PB0;
PMO = C.PA0 * C.PB1 + C.PA2 * C.PB1;
PSO = C.PA0 * C.PB2 + C.PA1 * C.PB0 + C.PA1 * C.PB2;
PFB = C.PA1 * C.PB1;
PDM = C.PA2 * C.PB0;
PDN = C.PA2 * C.PB2 + C.PA3 * (C.PB0 + C.PB1 + C.PB2);
availabilityVector = zeros(20,2);
availabilityVector(:,1) = computeAvailability(C, PPF, PMO, PSO, PFB, PDM, PDN);
C.PB1 = C.PB1^2 * (1 - C.PB0)^2;
C.PB2 = 2 * C.PB2 * (1 - C.PB0) - C.PB2^2 * (1 - C.PB0)^2;
C.PB0 = 1 - C.PB1 -C.PB2;

PPF = C.PA0 * C.PB0;
PMO = C.PA0 * C.PB1 + C.PA2 * C.PB1;
PSO = C.PA0 * C.PB2 + C.PA1 * C.PB0 + C.PA1 * C.PB2;
PFB = C.PA1 * C.PB1;
PDM = C.PA2 * C.PB0;
PDN = C.PA2 * C.PB2 + C.PA3 * (C.PB0 + C.PB1 + C.PB2);

availabilityVector(:,2) = computeAvailability(C, PPF, PMO, PSO, PFB, PDM, PDN);
availabilityVector
function [C] = initializeConstants()
    C.NUM_SYSTEM = 100000;
    C.TIME_STEP = 30000;
    C.LIFE_LIMIT = 200000;
    C.STATE_NUM_NODE = 6;
    C.k = 3;
    C.w = 30000;
    C.Lambda_A = 1/5.90e4;
    C.PA0 = exp(-C.Lambda_A * C.TIME_STEP);
    C.PA1 = 0.20 * (1 - C.PA0);
    C.PA2 = 0.15 * (1 - C.PA0);
    C.PA3 = 0.65 * (1 - C.PA0);
    C.Lambda_B = 1/2.20e5;
    C.PB0 = exp(-C.Lambda_B * C.TIME_STEP);
    C.PB1 = 0.45 * (1 - C.PB0);
    C.PB2 = 0.55 * (1 - C.PB0);
    C.NUM_NODE = 10;
end

function [availabilityVector] = computeAvailability(C, PPF, PMO, PSO, PFB, PDM, PDN)
    availabilityVector = zeros(20, 1);

    for NUM_NODE = 3:20

        for QPF = 0:NUM_NODE

            for QMO = 0:(NUM_NODE - QPF)

                for QSO = 0:(NUM_NODE - QPF - QMO)

                    for QFB = 0:(NUM_NODE - QPF - QMO - QSO)

                        for QDM = 0:(NUM_NODE - QPF - QMO - QSO - QFB)
                            QDN = NUM_NODE - QPF - QMO - QSO - QFB - QDM;

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
                                continue;
                            end
                            
                            if C5 && (C6 || C7)
                                availabilityVector(NUM_NODE) = availabilityVector(NUM_NODE) + nchoosek(NUM_NODE, QPF) * nchoosek(NUM_NODE - QPF, QMO) * nchoosek(NUM_NODE - QPF - QMO, QSO) * nchoosek(NUM_NODE - QPF - QMO - QSO, QFB) * nchoosek(NUM_NODE - QPF - QMO - QSO - QFB, QDM) * nchoosek(NUM_NODE - QPF - QMO - QSO - QFB - QDM, QDN) * PPF^(QPF) * PMO^(QMO) * PSO^(QSO) * PFB^(QFB) * PDM^(QDM) * PDN^(QDN);
                            else
                                availabilityVector(NUM_NODE) = availabilityVector(NUM_NODE) + nchoosek(NUM_NODE, QPF) * nchoosek(NUM_NODE - QPF, QMO) * nchoosek(NUM_NODE - QPF - QMO, QSO) * nchoosek(NUM_NODE - QPF - QMO - QSO, QFB) * nchoosek(NUM_NODE - QPF - QMO - QSO - QFB, QDM) * nchoosek(NUM_NODE - QPF - QMO - QSO - QFB - QDM, QDN) * PPF^(QPF) * PMO^(QMO) * PSO^(QSO) * PFB^(QFB) * PDM^(QDM) * PDN^(QDN) * QDM / (QDM + QPF);
                            end

                        end

                    end

                end

            end

        end

    end

end
