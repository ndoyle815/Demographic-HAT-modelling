function  R0Contributions(Dir, IDStr, Data, Paras, fitted_para_names, ourCluster)

    % Projection vectors; 6 humans, 2 animals, 4 tsetse 
    P_H1=diag([ones(3,1); zeros(3,1); zeros(2,1); ones(4,1)]);
    P_H2=diag([zeros(3,1); ones(3,1); zeros(2,1); ones(4,1)]);
    P_A=diag([zeros(6,1); ones(2,1); ones(4,1)]);
    P_H=diag([ones(6,1); zeros(2,1); ones(4,1)]);
    P_H1A=diag([ones(3,1); zeros(3,1); ones(2,1); ones(4,1)]);
    P_H2A=diag([zeros(3,1); ones(3,1); ones(2,1); ones(4,1)]);


    load([Dir, 'Posterior', Data.FileStr, strtok(Data.LocStr, '('), IDStr, '.mat']);

    NumPosterior = size(Posterior, 1);
    if NumPosterior ~= 1  
         
        R0_All = zeros(NumPosterior, 1);
        R0_H1 = zeros(NumPosterior, 1);
        R0_H2 = zeros(NumPosterior, 1);
        R0_A = zeros(NumPosterior, 1);
        R0_H = zeros(NumPosterior, 1);
        R0_H1A = zeros(NumPosterior, 1);
        R0_H2A = zeros(NumPosterior, 1);

        pc = parcluster;   %('local');
        pc.NumWorkers = ourCluster.NumWorkers;
        pc.JobStorageLocation = ourCluster.JobStorageLocation;
        parpool(pc);            
        parfor P = 1 : NumPosterior
            Paraz(P) = Paras;
            for i = 1 : length(fitted_para_names)
                Paraz(P).(fitted_para_names{i}) = Posterior{P, i};
            end
            % update meff
            meff = 1;
            K = NGM(Data.N_H, meff, Paraz(P));
            meff = (Paraz(P).R0 / max(abs(eig(K))))^2 * meff;
            % R0 contributions
            K = NGM(Data.N_H, meff, Paraz(P));
            R0_All(P) = max(abs(eig(K)));
            R0_H1(P) = max(abs(eig(P_H1*K)));
            R0_H2(P) = max(abs(eig(P_H2*K)));
            R0_A(P) = max(abs(eig(P_A*K)));
            R0_H(P) = max(abs(eig(P_H*K)));
            R0_H1A(P) = max(abs(eig(P_H1A*K)));
            R0_H2A(P) = max(abs(eig(P_H2A*K)));
        end
        delete(gcp('nocreate'));
    
        R0 = [R0_All R0_H1 R0_H2 R0_A R0_H R0_H1A R0_H2A];
        VariableNames = {'R0_All' ,'R0_H1' ,'R0_H2' ,'R0_A' ,'R0_H' ,'R0_H1A' ,'R0_H2A'};
        %save([Dir, 'R0_', Model, '_', LocStr, '_', IDStr, '.mat'], 'R0');
        %save([Dir, 'R0', Data.FileStr, Data.LocStr, IDStr, '.mat'], 'R0', 'VariableNames');

        Persistence = sum(R0_All > 1);
        Human = sum(R0_All > 1 & R0_H > 1 & R0_A < 1) / Persistence;
        Animal = sum(R0_All > 1 & R0_H < 1 & R0_A > 1) / Persistence;
        Both = sum(R0_All > 1 & R0_H < 1 & R0_A < 1) / Persistence;
        Either = sum(R0_All > 1 & R0_H > 1 & R0_A > 1) / Persistence;
        Persistence = Persistence / NumPosterior;

        save([Dir, 'R0Contributions', Data.FileStr, Data.LocStr, IDStr, '.mat'],...
             'VariableNames', 'R0', 'NumPosterior',...
             'Persistence', 'Human', 'Animal', 'Both', 'Either');
    end
end

function K = NGM(N_H, meff, Paras)
    
    k4 = 1 - Paras.k1 - Paras.k2 - Paras.k3;
    Plow = Paras.k1 + Paras.k3;
    Phigh = Paras.k2 + k4;

    N_Hlow = N_H * Plow;
    N_Hhigh = N_H * Phigh;
    f_Hlow = Paras.f_H * Plow / (Plow + Paras.r * Phigh);
    f_Hhigh = Paras.f_H * Paras.r * Phigh / (Plow + Paras.r * Phigh);

    N_A = N_H * Paras.k_A;

    if N_Hhigh == 0
        N_Hhigh = 1;
        f_Hhigh = 0;
    end
    if N_A == 0
        N_A = 1;
        Paras.f_A = 0;
    end

    eta_H0 = Paras.eta_H * Paras.b_eta_H0;
    
    
    %%% Compute m_eff using NGM approach, where the order of matrix elements is
    %%% (E_H, I1_H, I2_H) for low risk, (E_H, I1_H, I2_H) for high risk, E_A, I1_A, E1_V, E2_V, E3_V, I_V
    T = zeros(12,12);
    S = zeros(12,12);
    S_Vstar = Paras.mu_V * N_H / (Paras.alpha + Paras.mu_V); % equilibrium S_V
    G_Vfrom0 = Paras.alpha * N_H / (Paras.alpha + Paras.mu_V); % equilibrium (no infection) G_V

    %meff = 1; % assign a value for computing the corresponding R_0

    %T is transmissions
    T(1,12) = Paras.alpha * meff * f_Hlow; % I_V infects E_Hlow
    T(4,12) = Paras.alpha * meff * f_Hhigh; % I_V infects E_Hhigh
    T(7,12) = Paras.alpha * meff * Paras.f_A; % I_V infects E_A
    T(9,2) = Paras.alpha * f_Hlow * Paras.p_V * (S_Vstar + Paras.epsilon * G_Vfrom0) / N_Hlow; % I1_Hlow infects E1_V
    T(9,3) = T(9,2); % I2_Hlow infects E1_V
    T(9,5) = Paras.alpha * f_Hhigh * Paras.p_V * (S_Vstar + Paras.epsilon * G_Vfrom0) / N_Hhigh; % I1_Hhigh infects E1_V
    T(9,6) = T(9,5); % I2_Hhigh infects E1_V
    T(9,8) = Paras.alpha * Paras.f_A * Paras.p_V * (S_Vstar + Paras.epsilon * G_Vfrom0) / N_A; % I1_A infects E1_V

    %S is transissions (including passive stage1 detection)
    S(1,1) = - Paras.sigma_H - Paras.mu_H; % leaves E_Hlow
    S(2,1) = Paras.sigma_H; % enters I1_Hlow from E_Hlow
    S(2,2) = - eta_H0 - Paras.phi_H - Paras.mu_H; % leaves I1_Hlow
    S(3,2) = Paras.phi_H; % enters I2_Hlow from I1_Hlow
    S(3,3) = - Paras.gamma_H0 - Paras.mu_H; % leaves I2_Hlow
    S(4,4) = - Paras.sigma_H - Paras.mu_H; % leaves E_Hhigh
    S(5,4) = Paras.sigma_H; % enters I1_Hhigh from E_Hhigh
    S(5,5) = - eta_H0 - Paras.phi_H - Paras.mu_H; % leaves I1_Hhigh
    S(6,5) = Paras.phi_H; % enters I2_Hhigh from I1_Hhigh
    S(6,6) = - Paras.gamma_H0 - Paras.mu_H; % leaves I2_Hhigh
    S(7,7) = - Paras.sigma_A - Paras.mu_A; % leaves E_A
    S(8,7) = Paras.sigma_A; % enters I_A from E_A
    S(8,8) = - Paras.mu_A; % leaves I1_A (no treatment and no staging) 
    S(9,9) = - 3 * Paras.sigma_V - Paras.mu_V; % leaves E1_V
    S(10,9) =  3 * Paras.sigma_V; % enters E2_V from E1_V
    S(10,10) = - 3 * Paras.sigma_V - Paras.mu_V; % leaves E2_V
    S(11,10) = 3 * Paras.sigma_V; % enters E3_V from E2_V
    S(11,11) = - 3 * Paras.sigma_V - Paras.mu_V; % leaves E3_V
    S(12,11) = 3 * Paras.sigma_V; % enters I_V from E3_V
    S(12,12) = - Paras.mu_V; % leaves I_V

    K = - T * inv(S);
    R0_current = max(abs(eig(K)));
end
