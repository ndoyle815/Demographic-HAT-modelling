
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                               %
%   This code deals with    %
%                                                                                                               %
%   Inputs:                                                                                                     %
%       Data - structure containing location-specific historical data                                           %
%       Paras - structure containing location-specific parameters (fixed, fitted and intervention parameters)   %
%       FittedPrior - structure containing prior distributions of fitted parameters                             %
%       CstrParas - structure containing constrained parameters                                                 %
%       ProjStrat - structure containing parameters associated with future strategy                             %
%       locFittedPar - a table containing location-specific fitted parameters for all models                    %
%       MachineSetting - structure containing parameters to control parcluster                                  %
%                                                                                                               %
%   Main output files:                                                                                          %
%       VC - a table containing ReductionPct and TargetDie when it's fitted                                     %
%       R0 - a table containing R0 contributions from different hosts                                           %
%       DIC - a table containing values of all DIC variants                                                     %
%       Evidence - a table model evidence                                                                       %
%       Ensemble - a file containing 2000 selected posteriors for each DIC variant                              %
%                                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PostProcessor(RunInfo, Data, Paras, FittedPrior, CstrParas, ProjStrat, locFittedPar, MachineSetting) % fitted_para_names is stored in Paras.FittedNames
    if strcmp(MachineSetting.Type, 'cluster') 
        pc = parcluster;  %('local');
        pc.NumWorkers = MachineSetting.NumWorkers;
        pc.JobStorageLocation = MachineSetting.JobStorageLocation;
    end
    
    if isfile(RunInfo.FitFilePath('Posterior'))
        input = load(RunInfo.FitFilePath('Posterior'));
        Posterior = input.Posterior;
        NumPosterior = size(Posterior, 1);

        %%% Compute VCreduction if TargetDie is a fitted parameter and VC table doesn't exist in Posterior.mat
        if sum(ismember(Paras.FittedNames, 'TargetDie')) == 1 && exist('VC') == 0
            TargetDie = Posterior.TargetDie;
            ReductionPct = GetVCReductionPct(Paras, Paras.VCwgt, TargetDie, MachineSetting);

            VC = array2table([round(ReductionPct, 2) TargetDie], 'VariableNames', {'ReductionPct', 'TargetDie'});
            save(RunInfo.FitFilePath('Posterior'), 'VC', '-append')                
        end
        

        %%% Compute R0 of different categories, the maintenance host and meff if R0 table doesn't exist in Posterior.mat
%         if exist('R0') == 0
%             if sum(Posterior.R0 < 1) == 0
%                 % Projection vectors; 6 humans (E, I1, I2 for low/high risk), 2 animals (E, I), 4 tsetse (E1, E2, E3, I) 
%                 P_H1=diag([ones(3,1); zeros(3,1); zeros(2,1); ones(4,1)]);
%                 P_H2=diag([zeros(3,1); ones(3,1); zeros(2,1); ones(4,1)]);
%                 P_A=diag([zeros(6,1); ones(2,1); ones(4,1)]);
%                 P_H=diag([ones(6,1); zeros(2,1); ones(4,1)]);
%                 P_H1A=diag([ones(3,1); zeros(3,1); ones(2,1); ones(4,1)]);
%                 P_H2A=diag([zeros(3,1); ones(3,1); ones(2,1); ones(4,1)]);
%                 
%                 R0_All = zeros(NumPosterior, 1);
%                 R0_H1 = zeros(NumPosterior, 1);
%                 R0_H2 = zeros(NumPosterior, 1);
%                 R0_A = zeros(NumPosterior, 1);
%                 R0_H = zeros(NumPosterior, 1);
%                 R0_H1A = zeros(NumPosterior, 1);
%                 R0_H2A = zeros(NumPosterior, 1);
%                 maintenance = zeros(NumPosterior, 4); % human, animal, both, either
%                 meff = zeros(NumPosterior, 1);
%                 code = [1 2 0 3]; % human, animal, both, either
%                 meff0 = 1;
% 
%                 if strcmp(MachineSetting.Type, 'cluster') 
%                     parpool(pc);
%                 end
%                 parfor p = 1 : NumPosterior
%                     Paraz(p) = Paras;
%                     for i = 1 : length(Paras.FittedNames)
%                         Paraz(p).(Paras.FittedNames{i}) = Posterior{p, i};
%                     end
% 
%                     if size(CstrParas, 1) == 1
%                         k = [Paraz(p).k1 Paraz(p).k2 Paraz(p).k3 Paraz(p).k4];
%                         Paraz(p).(CstrParas.Notation{:}) = 1 - sum(k(~strcmp(CstrParas.Notation, {'k1','k2','k3','k4'})));
%                     end
% 
%                     % Update meff
%                     K = NGM(Data.N_H, meff0, Paraz(p), Data);
%                     meff(p) = (Paraz(p).R0 / max(abs(eig(K))))^2 * meff0;
%                     
%                     % R0 contributions
%                     K = NGM(Data.N_H, meff(p), Paraz(p), Data);
%                     R0_All(p) = max(abs(eig(K)));
%                     R0_H1(p) = max(abs(eig(P_H1*K)));
%                     R0_H2(p) = max(abs(eig(P_H2*K)));
%                     R0_A(p) = max(abs(eig(P_A*K)));
%                     R0_H(p) = max(abs(eig(P_H*K)));
%                     R0_H1A(p) = max(abs(eig(P_H1A*K)));
%                     R0_H2A(p) = max(abs(eig(P_H2A*K)));
%                     
%                     % Determining maintenance hosts
%                     maintenance(p, :) = code == (R0_H(p) > 1) + 2 * (R0_A(p) > 1);
%                 end
%                 delete(gcp('nocreate'));
%                 clear Paraz pc;
%                 
%                 R0 = array2table([R0_All R0_H R0_A R0_H1 R0_H2 R0_H1A R0_H2A maintenance meff], ...
%                                  'VariableNames', {'R0_All', 'R0_H', 'R0_A', 'R0_H1', 'R0_H2', 'R0_H1A', 'R0_H2A', 'Huamn', 'Animal', 'Both', 'Either', 'meff'});
%                 save(RunInfo.FitFilePath('Posterior'), 'R0', '-append')
%             else
%                 warning(['R0 < 1 in some posteriors in ', strtok(Data.FileStr, '_')])
%             end
%         end
            
            
        %%% Compute DIC values if DIC table doesn't exist in Posterior.mat
        if exist('DIC') == 0
            % Mean of all fitted parameters
            ThetaMean = mean(Posterior{:, 1 : end - 1});
            MeanLL = mean(-Posterior.NegLogLikelihood);
            VarLL = var(Posterior.NegLogLikelihood);
            LLThetaMean = - Get_log_Prob(Data, Paras, Paras.FittedNames, ThetaMean, FittedPrior, CstrParas, ProjStrat);
                
            DIC1 = -2 * MeanLL + (-2 * MeanLL + 2 * LLThetaMean);
            DIC2 = -2 * MeanLL + 2 * VarLL;
            DIC1v2 = -2 * MeanLL + 2 * (-2 * MeanLL + 2 * LLThetaMean);
            DICalt = -2 * LLThetaMean + 4 * VarLL;
            DIC = table(DICalt, DIC1, DIC2, DIC1v2);
            save(RunInfo.FitFilePath('Posterior'), 'DIC', '-append')
        end
        
        %%% Compute model evidence if Evidence table doesn't exist in Posterior.mat
        if exist('Evidence') == 0
%         [LogModelEvidence, DF, PostCorr, PostSD, PostMean, VariableNames] = deal('From placeholder Evidence function');
%         save(RunInfo.FitFilePath('Posterior'),...
%         'VariableNames', 'Evidence', 'DF',...
%         'PostCorr', 'PostSD', 'PostMean');
        end
    end
end






function dPop = TsetseDyn(t, pop, parameter)
    if parameter.TargetDie==0
        f_T = 0;
    else
        f_T = parameter.TargetDie * (1 - 1/(1+exp(-25/365*(mod(t,365/parameter.TargetFreq)-0.35*365))));
    end
    
    dPop = zeros(3,1);
    P_V = pop(1);
    S_V = pop(2);
    G_V = pop(3);
    N_V = S_V + G_V;

    % Pupa
    dPop(1) = parameter.B_V * N_V - (parameter.xi_V + P_V/parameter.k_V) * P_V;

    % Teneral (feed twice as fast as other adults)
    dPop(2) = parameter.xi_V * parameter.p_survive * P_V - parameter.alpha * S_V - parameter.mu_V * S_V;

    % Non-teneral
    dPop(3)= parameter.alpha * (1 - f_T) * S_V - parameter.alpha * f_T * G_V - parameter.mu_V * G_V;
    
end

function K = NGM(N_H, meff, Paras, Data)
    
    % bite rates
    k=[Data.N_FY./Data.N_H Data.N_FW./Data.N_H Data.N_FP./Data.N_H Data.N_MY./Data.N_H Data.N_MW./Data.N_H Data.N_MP./Data.N_H Paras.k_A 1];
    
    s_A = 0;
    s_N = 0;
    bitepref=[1 Paras.rFW Paras.rFP Paras.rMY Paras.rMW Paras.rMP s_A s_N];            %biting preference on hosts (given no animal reservoir)  
    f=(bitepref.*k)/sum(bitepref.*k);
    f(7) = [];
    f(8) = [];

    eta_H0 = Paras.eta_H * Paras.b_eta_H0; % Fixed as 0 in Paras.mat


    %%% Compute m_eff using NGM approach, where the order of matrix elements is
    %%% (E_H (of all demog), I1_H (of all demog), I2_H (of all demog) demog order is FY FW FP MY MW MP, then   E1_V, E2_V, E3_V, I_V
    T = zeros(22,22);
    S = zeros(22,22);
    S_Vstar = Paras.mu_V * N_H / (Paras.alpha + Paras.mu_V); % equilibrium S_V
    G_Vfrom0 = Paras.alpha * N_H / (Paras.alpha + Paras.mu_V); % equilibrium (no infection) G_V

    %T is transmissions
    T(1,22) = Paras.alpha * meff * f(1);
    T(2,22) = Paras.alpha * meff * f(2);
    T(3,22) = Paras.alpha * meff * f(3);
    T(4,22) = Paras.alpha * meff * f(4);
    T(5,22) = Paras.alpha * meff * f(5);
    T(6,22) = Paras.alpha * meff * f(6);

    T(19,7) = Paras.alpha * f(1) * Paras.p_V * (S_Vstar + Paras.epsilon * G_Vfrom0) / Data.N_FY;
    T(19,8) = Paras.alpha * f(2) * Paras.p_V * (S_Vstar + Paras.epsilon * G_Vfrom0) / Data.N_FW;
    T(19,9) = Paras.alpha * f(3) * Paras.p_V * (S_Vstar + Paras.epsilon * G_Vfrom0) / Data.N_FP;
    T(19,10) =  Paras.alpha * f(4) * Paras.p_V * (S_Vstar + Paras.epsilon * G_Vfrom0) / Data.N_MY;
    T(19,11) =  Paras.alpha * f(5) * Paras.p_V * (S_Vstar + Paras.epsilon * G_Vfrom0) / Data.N_MW;
    T(19,12) =  Paras.alpha * f(6) * Paras.p_V * (S_Vstar + Paras.epsilon * G_Vfrom0) / Data.N_MP;
    T(19,13:18) = T(19,7:12);


    %S is transitions (including passive stage1 detection)
    % NATHAN
    % E_H
    S(1,1) = - (Paras.sigma_H + Paras.mu_H_FY + Paras.l_Y);
    S(2,2) = - (Paras.sigma_H + Paras.mu_H_FW + Paras.l_W);
    S(3,3) = - (Paras.sigma_H + Paras.mu_H_FP);
    S(4,4) = - (Paras.sigma_H + Paras.mu_H_MY + Paras.l_Y);
    S(5,5) = - (Paras.sigma_H + Paras.mu_H_MW + Paras.l_W);
    S(6,6) = - (Paras.sigma_H + Paras.mu_H_MP);
    S(2,1) = Paras.l_Y;
    S(3,2) = Paras.l_W;
    S(5,4) = S(2,1);
    S(6,5) = S(3,2);

    % I1_H
    S(7,7) = - (eta_H0 + Paras.phi_H + Paras.mu_H_FY + Paras.l_Y);
    S(8,8) = - (eta_H0 + Paras.phi_H + Paras.mu_H_FW + Paras.l_W);
    S(9,9) = - (eta_H0 + Paras.phi_H + Paras.mu_H_FP);
    S(10,10) = - (eta_H0 + Paras.phi_H + Paras.mu_H_MY + Paras.l_Y);
    S(11,11) = - (eta_H0 + Paras.phi_H + Paras.mu_H_MW + Paras.l_W);
    S(12,12) = - (eta_H0 + Paras.phi_H + Paras.mu_H_MP);
    S(8,7) = S(2,1);
    S(9,8) = S(3,2);
    S(11,10) = S(2,1);
    S(12,11) = S(3,2);
    S(7,1) = Paras.sigma_H;
    S(8,2) = S(7,1);
    S(9,3) = S(7,1);
    S(10,4) = S(7,1);
    S(11,5) = S(7,1);
    S(12,6) = S(7,1);

    % I2_H
    S(13,13) = - (Paras.gamma_H0 + Paras.mu_H_FY + Paras.l_Y);
    S(14,14) = - (Paras.gamma_H0 + Paras.mu_H_FW + Paras.l_W);
    S(15,15) = - (Paras.gamma_H0 + Paras.mu_H_FP);
    S(16,16) = - (Paras.gamma_H0 + Paras.mu_H_MY + Paras.l_Y);
    S(17,17) = - (Paras.gamma_H0 + Paras.mu_H_MW + Paras.l_W);
    S(18,18) = - (Paras.gamma_H0 + Paras.mu_H_MP);
    S(14,13) = S(2,1);
    S(15,14) = S(3,2);
    S(17,16) = S(2,1);
    S(18,17) = S(3,2);
    S(13,7) = Paras.phi_H;
    S(14,8) = S(13,7);
    S(15,9) = S(13,7);
    S(16,10) = S(13,7);
    S(17,11) = S(13,7);
    S(18,12) = S(13,7);
    
    % Tsetse (transitions unchanged)
    S(19,19) = - 3 * Paras.sigma_V - Paras.mu_V; % leaves E1_V
    S(20,19) =  3 * Paras.sigma_V; % enters E2_V from E1_V
    S(20,20) = - 3 * Paras.sigma_V - Paras.mu_V; % leaves E2_V
    S(21,20) = 3 * Paras.sigma_V; % enters E3_V from E2_V
    S(21,21) = - 3 * Paras.sigma_V - Paras.mu_V; % leaves E3_V
    S(22,21) = 3 * Paras.sigma_V; % enters I_V from E3_V
    S(22,22) = - Paras.mu_V; % leaves I_V


    % Compute NGM
    K = - T * inv(S);
end

