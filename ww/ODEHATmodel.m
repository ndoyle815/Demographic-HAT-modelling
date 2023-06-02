
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                               %
%   This code solves ordinary differential equations of the Warwick HAT model                                   %
%                                                                                                               %
%   Inputs:                                                                                                     %
%       meff - number denoting effective vector density, calculated by the relation Paras.R0^2 ~ meff           %
%       ICs - cell array containing endemic equilibrium of given set of parameters                              %
%       Data - structure containing location-specific historical data                                           %
%       Paras - structure containing location-specific parameters (fixed, fitted and intervention parameters)   %
%       ProjStrat - structure containing parameters associated with future strategy                             %
%                                                                                                               %
%   Outputs:                                                                                                    %
%       Classes - table containing time series of model outputs (e.g. susceptible humans, infectious vectors)   %
%       Aggregate - table containing yearly aggregated outputs (e.g. active/passive stage 1/2 cases, deathss)   %
%                                                                                                               %
%   Note: hosts are (1) FY                                                           %
%                   (2) FW                                                          %
%                   (3) FP                                                              %
%                   (4) MY                                                             %
%                   (5) MW                                                                       %
%                   (6) MP                                      %
%                                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Classes, Aggregate] = ODEHATmodel(meff, ICs, Data, Paras, ProjStrat)
    ODEoptions = odeset('NonNegative', 1:37); % M9
    [S_H, E_H, I1_H, I2_H, R_H, P_V, S_V, G_V, E1_V, E2_V, E3_V, I_V] = ICs{:}; % M9
    
    %k4 = 1 - Paras.k1 - Paras.k2 - Paras.k3;
    N_A = Data.N_H * Paras.k_A;
    % K_V = Data.N_H * Paras.k_V;
    K_V = Data.N_H * Paras.mu_V / (Paras.xi_V^2 * Paras.p_survive * (Paras.p_survive * Paras.B_V / Paras.mu_V - 1));
    
    k=[Data.N_FY./Data.N_H Data.N_FW./Data.N_H Data.N_FP./Data.N_H Data.N_MY./Data.N_H Data.N_MW./Data.N_H Data.N_MP./Data.N_H Paras.k_A 1];
    
    % ELLIOT question for kat pending an answer, currently have just set
    % these to 0
    s_A = 0;
    s_N = 0;
    bitepref=[1 Paras.rFW Paras.rFP Paras.rMY Paras.rMW Paras.rMP s_A s_N];            %biting preference on hosts (given no animal reservoir)  
    f=(bitepref.*k)/sum(bitepref.*k).*Paras.f_H;
    f(8) = [];
    f(7) = [];
    
    NumberScreening = length(Data.ModelScreeningTime);
    %NumberScreening
    

    % Active screening: will move random participating Infectious (I1_H and I2_H) to Recovery at the end of the year
    % Default AS strategy 'traditional' (screen Paras.k1 and Paras.k2 only)

    %Data.ModelPeopleScreened_FY
    TurnOut_FY = Data.ModelPeopleScreened_FY / Data.N_FY;
    TurnOut_FW = Data.ModelPeopleScreened_FW / Data.N_FW;
    TurnOut_FP = Data.ModelPeopleScreened_FP / Data.N_FP;
    TurnOut_MY = Data.ModelPeopleScreened_MY / Data.N_MY;
    TurnOut_MW = Data.ModelPeopleScreened_MW / Data.N_MW;
    TurnOut_MP = Data.ModelPeopleScreened_MP / Data.N_MP;

    TurnOut = [TurnOut_FY ; TurnOut_FW ; TurnOut_FP ; TurnOut_MY ; TurnOut_MW ; TurnOut_MP ];
    
% ELLIOT I think the below should be ignored? We aren't
% changing turnout dependent on how screening happens anymore - we're just
% assuming turnout is what is in our data (as me if you dont know what I
% mean, its a bit hard to explain)

    % Change AS strategy
    mode = split(ProjStrat.NewASstrat, '_');
    M = find(~strcmp(mode, 'traditional')' & ProjStrat.NewASyear <= Data.Years(end));
%     if ~isempty(M) % NewASstrat has something other than 'traditional' and simulation uses new screening mode
%         NewASyear = [ProjStrat.NewASyear Data.Years(end) + 1];
%         for m = M
%             Y = find(Data.ModelScreeningTime >= NewASyear(m) & Data.ModelScreeningTime < NewASyear(m + 1));
%             switch mode{m}
%                 case 'equal' % door-to-door (all populations get equal probability to be screened)
%                     TurnOut1(Y) = Data.ModelPeopleScreened(Y) / Data.N_H;        
%                       TurnOut2(Y) = TurnOut1(Y);
%                     TurnOut3(Y) = TurnOut1(Y);
%                     TurnOut4(Y) = TurnOut1(Y);
%                 case 'high' % work place screening (k2 group gets screened first and then k4 then k1)
%                     TurnOut4(Y) = 0;
%                     TurnOut2(Y) = min(Data.ModelPeopleScreened(Y) / Data.N_H / Paras.k2, 1);
%                     if TurnOut2(Y) == 1
%                         TurnOut4(Y) = min((Data.ModelPeopleScreened(Y) / Data.N_H / - Paras.k2) / Paras.k4, 1);
%                     end
%                     TurnOut1(Y) = max((Data.ModelPeopleScreened(Y) / Data.N_H - Paras.k2 - Paras.k4) / Paras.k1, 0);
%             end
%         end
%     end
    
% %     if ~strcmp(ProjStrat.NewASstrat, 'traditional') && ~strcmp(ProjStrat.NewASstrat, 'stop') && Data.ModelScreeningTime(1) < ProjStrat.NewASyear && ProjStrat.NewASyear < Data.ModelScreeningTime(end)
% %         Y = find(Data.ModelScreeningTime == ProjStrat.NewASyear);
% %         
% %         switch ProjStrat.NewASstrat
% %         case 'equal' % door-to-door (all populations get equal probability to be screened)
% %             TurnOut1(Y:end) = Data.ModelPeopleScreened(Y:end) / Data.N_H;
% %             TurnOut2(Y:end) = TurnOut1(Y:end);
% %             TurnOut3(Y:end) = TurnOut1(Y:end);
% %             TurnOut4(Y:end) = TurnOut1(Y:end);
% %         case 'high' % work place screening (k4 group gets screened first and then equally screened Paras.k1 and Paras.k2)
% %             TurnOut4(Y:end) = min(Data.ModelPeopleScreened(Y:end) / Data.N_H / k4, 1);
% %             TurnOut1(Y:end) = max((Data.ModelPeopleScreened(Y:end) / Data.N_H - k4) / (Paras.k1 + Paras.k2), 0);
% %             TurnOut2(Y:end) = TurnOut1(Y:end);
% %         end
% %     end

    % Specificity by screening
    ScreeningSpecificity = repmat(Paras.specificity, 1, NumberScreening);
    % Lower specificity can exist in some years (eg MSF interventions in
    % Ango and Ganga HZ of Orientale former province).
    altSpec = find(Data.ModelScreeningTime <= Paras.Last_year);
    ScreeningSpecificity(altSpec) = Paras.specificityMSF; %ScreeningSpecificity(altSpec) * Paras.b_specificity;
    % Also change sensitivity
    ScreeningSensitivity = repmat(Paras.Sensitivity, 1, NumberScreening);
    ScreeningSensitivity(altSpec) = Paras.SensitivityMSF;
    %ScreeningSensitivity
    S = NumberScreening + 1;
    if ProjStrat.TTyear ~= 0 && Data.ModelScreeningTime(end) >= ProjStrat.TTyear
        S = find(Data.ModelScreeningTime >= ProjStrat.TTyear, 1);
        ScreeningSensitivity(S : end) = ProjStrat.SensRDT;
        ScreeningSpecificity(S : end) = ProjStrat.SpecRDT;
        if contains(ProjStrat.LabTest, 'N')
            ProjStrat.SensLab = 1;
            ProjStrat.SpecLab = 0;
        end
    end
    
    % Passive screening: move Infected (I1 and I2) to Recovery continuously
    if Data.ModelScreeningTime(1) < ProjStrat.RDTyear && ProjStrat.RDTyear < Data.ModelScreeningTime(end)
        Y = find(Data.ModelScreeningTime == ProjStrat.RDTyear);
    else
        Y = length(Data.ModelScreeningTime) + 1;
    end
    yearlyeta_H = [(1 + Paras.eta_H_amp ./ (1 + exp(-Paras.d_steep * (double(Data.ModelScreeningTime(1:Y-1)) - (Paras.d_change+Paras.eta_H_lag))))) * Paras.eta_H,...
                   (1 + ProjStrat.RDTincrease) * (1 + Paras.eta_H_amp ./ (1 + exp(-Paras.d_steep * (double(Data.ModelScreeningTime(Y-1)) - (Paras.d_change+Paras.eta_H_lag))))) * Paras.eta_H * ones(1, NumberScreening - (Y-1))];
    yearlygamma_H = [(1 + Paras.gamma_H_amp ./ (1 + exp(-Paras.d_steep * (double(Data.ModelScreeningTime(1:Y-1)) - Paras.d_change)))) * Paras.gamma_H,...
                   (1 + ProjStrat.RDTincrease) * (1 + Paras.gamma_H_amp ./ (1 + exp(-Paras.d_steep * (double(Data.ModelScreeningTime(Y-1)) - Paras.d_change)))) * Paras.gamma_H * ones(1, NumberScreening - (Y-1))];
    
    if Data.ModelScreeningTime(1) == ProjStrat.RDTyear
        yearlyeta_H = (1 + ProjStrat.RDTincrease) * (1 + Paras.eta_H_amp ./ (1 + exp(-Paras.d_steep * (double(Data.ModelScreeningTime(1) - 1) - (Paras.d_change+Paras.eta_H_lag))))) * Paras.eta_H * ones(1, NumberScreening);
        yearlygamma_H = (1 + ProjStrat.RDTincrease) * (1 + Paras.gamma_H_amp ./ (1 + exp(-Paras.d_steep * (double(Data.ModelScreeningTime(1) - 1) - Paras.d_change)))) * Paras.gamma_H * ones(1, NumberScreening);
    end
%     Data.ModelScreeningTime
%     yearlyeta_H
%     yearlygamma_H
    
    % Update passive interruptions %%% Change to fitted interruption 
    yearlyeta_H(Data.ModelScreeningTime >= Paras.PPSstart & Data.ModelScreeningTime < Paras.PPSend) = 0;
    yearlygamma_H(Data.ModelScreeningTime >= Paras.PPSstart & Data.ModelScreeningTime < Paras.PPSend) = Paras.gamma_H0;
    yearlyeta_H(Data.ModelScreeningTime >= Paras.NoPSstart & Data.ModelScreeningTime < Paras.NoPSend) = 0;
    yearlygamma_H(Data.ModelScreeningTime >= Paras.NoPSstart & Data.ModelScreeningTime < Paras.NoPSend) = Paras.death;
    
%     yearlyeta_H
%     yearlygamma_H
    % calculate yearly uVector maintaining a constant death rate
    %Paras.death = (1-Paras.u) * Paras.gamma_H;
    %uVector = 1- Paras.death ./ yearlygamma_H;
    S2treatment = (yearlygamma_H - Paras.death) * Paras.effcy;
    yearlygamma_H = Paras.death + S2treatment;
    %yearlyeta_H = yearlyeta_H * Paras.effcy;
    
    % Vector control
    p_targetdie(1 : NumberScreening+1) = 0;
    TargetFreq(1 : NumberScreening+1) = 1;
    if Paras.VCstart ~= 0
        Y1 = find([Data.ModelScreeningTime Data.ModelScreeningTime(end) + 1] >= Paras.VCstart, 1);
        p_targetdie(Y1:end) = Paras.TargetDie;
        TargetFreq(Y1:end) = Paras.TargetFreq;
    end
    if Paras.VC_scaleback ~= 0
        Y2 = find([Data.ModelScreeningTime Data.ModelScreeningTime(end) + 1] >= Paras.VC_scaleback, 1);
        p_targetdie(Y2:end) = Paras.TargetDie_scaleback;
        TargetFreq(Y2:end) = Paras.TargetFreq_scaleback;
    end
    if ProjStrat.NewVCyear ~= 0
        Y3 = find([Data.ModelScreeningTime Data.ModelScreeningTime(end) + 1] >= ProjStrat.NewVCyear, 1);
        p_targetdie(Y3:end) = ProjStrat.NewTargetDie;
        TargetFreq(Y3:end) = ProjStrat.NewTargetFreq;
    end
    
    
    Pop = [S_H E_H I1_H I2_H R_H P_V S_V G_V E1_V E2_V E3_V I_V]; % M9
    T = mod(Data.ModelScreeningTime(1), 1) * 365;
    
%     Time = Data.ModelScreeningTime - 2000;
%     Freq = Data.ModelScreeningFreq;
%     Number = Data.ModelPeopleScreened;
%     Time
%     Freq
%     Number
    
    Data.ModelScreeningTime = [Data.ModelScreeningTime floor(Data.ModelScreeningTime(end))+1];

    [ActiveM1, ActiveM2, ActiveMFP,ActiveMAgeY, ActiveMAgeW, ActiveMAgeP,ActiveMGenderM, ActiveMGenderF, PassiveM1, PassiveM2,PassiveMAgeY, PassiveMAgeW, PassiveMAgeP, PassiveMGenderM, PassiveMGenderF, DeathsM, PersonYrsM1, PersonYrsM2,PersonYrsMAgeP, PersonYrsMAgeY, PersonYrsMAgeW,PersonYrsMGenderM, PersonYrsMGenderF, NewInfM, PerfectSpec, Transmission, RDTM1, RDTM2, RDTMFP,RDTMGenderF,RDTMGenderM,RDTMAgeY,RDTMAgeW,RDTMAgeP] = deal(zeros(1, length(Data.Years)));
    [Active1, Active2, ActiveFP,ActiveFP, ActiveAgeY, ActiveAgeW, ActiveAgeP,ActiveGenderM, ActiveGenderF,Passive1, Passive2,PassiveAgeY, PassiveAgeW, PassiveAgeP, PassiveGenderM, PassiveGenderF, Deaths, PersonYrs1, PersonYrs2,PersonYrsAgeP, PersonYrsAgeY, PersonYrsAgeW,PersonYrsGenderM, PersonYrsGenderF,PersonYrsAgeP_1, PersonYrsAgeY_1, PersonYrsAgeW_1,PersonYrsGenderM_1, PersonYrsGenderF_1,PersonYrsAgeP_2, PersonYrsAgeY_2, PersonYrsAgeW_2,PersonYrsGenderM_2, PersonYrsGenderF_2,  NewInf, Spec, Trans, RDT1, RDT2, RDTFP,RDTGenderF, RDTGenderM, RDTAgeY, RDTAgeW, RDTAgeP] = deal(zeros(1, NumberScreening));
%    mod(Data.ModelScreeningTime(1), 1) * 365 + [sum(Data.ModelScreeningFreq(1 : 0)) sum(Data.ModelScreeningFreq(1 : 1))]
%    Y = length(Data.Years);
    Y = length(Data.Years); %%%%%%%%%% ProjStrat
    for s = 1 : NumberScreening
        % No transmissons after EOT
        if s > 1 && Paras.alpha ~= 0 && Data.ModelScreeningTime(s) == round(Data.ModelScreeningTime(s)) && ...
           sum(NewInf(floor(Data.ModelScreeningTime)==floor(Data.ModelScreeningTime(s-1)))) < Paras.EOTthreshold * Data.PopGrowth .^ double(Data.PopSizeYear - floor(Data.ModelScreeningTime(s-1))) && ...
           sum(I1_H(end, :) + I2_H(end, :)) < Paras.EOTthreshold * Data.PopGrowth .^ double(Data.PopSizeYear - floor(Data.ModelScreeningTime(s-1))) % M9
            Paras.alpha = 0;
        end
        Trans(s) = Paras.alpha ~= 0;
        
       
        % Expected active PassiveAgeY(s) = (yearlyeta_H(s) * PersonYrsAgeY_1(s) + S2treatment(s) * PersonYrsAgeY_2(s))* 365;screening detections (S1/S2) each time interval

        TruePos = TurnOut(:,s) * ScreeningSensitivity(s);
        FalsePos = TurnOut(:,s) * (1 - ScreeningSpecificity(s));
        
        % If there are less than 0.5 cases per 10,000 above the expect false+ves and if it is before 2018	RDTMGenderF,RDTMGenderM,RDTMAgeY,RDTMAgeW,RDTMAgeP
        % Or beyond year_spec_100pct
        if (ProjStrat.TTyear == 0 || Data.ModelScreeningTime(s) < ProjStrat.TTyear) && ...
            ScreeningSpecificity(s) ~= 1 && Data.ModelPeopleScreened(s) > 0 &&  ((Paras.year_spec_100pct == 0 && (((I1_H(end,:) + I2_H(end,:)) * (TruePos - FalsePos))* 10000 / Data.N_H < 0.5 && Data.ModelScreeningTime(s) >= 3000)) || ... % M9 update here if skin FP is considered
                (Paras.year_spec_100pct ~= 0 && Data.ModelScreeningTime(s) >= Paras.year_spec_100pct))
            ScreeningSpecificity(s : S-1) = 1;
            FalsePos = TurnOut(:,s) * (1 - ScreeningSpecificity(s));
        end
         
        
         % Update case numbers from lab test
        if s < S
       
           Active1(s) = I1_H(end,:) * TruePos;
           Active2(s) = I2_H(end,:) * TruePos;
           ActiveFP(s) = S_H(end,:) * FalsePos;
           ActiveGenderF(s) = (I1_H(end,1:3) + I2_H(end,1:3))* TruePos(1:3);
           ActiveGenderM(s) = (I1_H(end,4:6) + I2_H(end,4:6))* TruePos(4:6);
           ActiveAgeY(s) = (I1_H(end,[1,4]) + I2_H(end,[1,4]))* TruePos([1,4]);
           ActiveAgeW(s) = (I1_H(end,[2,5]) + I2_H(end,[2,5]))* TruePos([2,5]);
           ActiveAgeP(s) = (I1_H(end,[3,6]) + I2_H(end,[3,6]))* TruePos([3,6]);
               

        else
           
           RDT1(s) = I1_H(end,:) * TruePos;
            RDT2(s) = I2_H(end,:) * TruePos;
            RDTFP(s) = S_H(end,:) * FalsePos;

            RDTGenderF(s) = (I1_H(end,1:3) + I2_H(end,1:3))* TruePos(1:3);
           RDTGenderM(s) = (I1_H(end,4:6) + I2_H(end,4:6))* TruePos(4:6);
           RDTAgeY(s) = (I1_H(end,[1,4]) + I2_H(end,[1,4]))* TruePos([1,4]);
           RDTAgeW(s) = (I1_H(end,[2,5]) + I2_H(end,[2,5]))* TruePos([2,5]);
           RDTAgeP(s) = (I1_H(end,[3,6]) + I2_H(end,[3,6]))* TruePos([3,6]);


            Active1(s) = RDT1(s) * ProjStrat.SensLab;
            Active2(s) = RDT2(s) * ProjStrat.SensLab;
            ActiveFP(s) = RDTFP(s) * (1 - ProjStrat.SpecLab);

            ActiveGenderF(s) = RDTGenderF(s) * ProjStrat.SensLab;
           ActiveGenderM(s) = RDTGenderM(s) * ProjStrat.SensLab;
           ActiveAgeY(s) = RDTAgeY(s) * ProjStrat.SensLab;
           ActiveAgeW(s) = RDTAgeW(s) * ProjStrat.SensLab;
           ActiveAgeP(s) = RDTAgeP(s) * ProjStrat.SensLab;
            ScreeningSpecificity(s) = ScreeningSpecificity(s) + (1 - ScreeningSpecificity(s)) * ProjStrat.SpecLab;

        
        end
       

        Spec(s) = ScreeningSpecificity(s) == 1;
%         Spec(s) = ScreeningSpecificity(s) == 1;
%          %changed FP code here!
%         Active1(s) = I1_H(end,:) * TruePos; %+ 0.5*S_H(end,:) * FalsePos;
%         Active2(s) = I2_H(end,:) * TruePos;%+ 0.5*S_H(end,:) * FalsePos; %assume false positives are detected as 50:50 S2 and S1
%         ActiveFP(s)= S_H(end,:) * FalsePos;
        
        
        % Dynamics
        %TurnOut(:,s)'
        DandT = TurnOut(:,s)' * ScreeningSensitivity(s) * Paras.Compliance; % proportional change in different HUMAN group
        %I1_H(end,:)
        ICs = [S_H(end,:) E_H(end,:) I1_H(end,:).*(1-DandT) I2_H(end,:).*(1-DandT) R_H(end,:)+(E_H(end,:)+I1_H(end,:)+I2_H(end,:)).* DandT Pop(end,31:37)]; % M9
        
        % Tsetse reintroduction after VC stops (change ICs)
        if Data.ModelScreeningTime(s) == ProjStrat.NewVCyear && ProjStrat.NewTargetDie == 0 && Paras.VCstart ~= 0 && ProjStrat.TsetseReintro ~= 0
            ICs(32) = 0.01 * ProjStrat.TsetseReintro * Data.N_H * Paras.mu_V / (Paras.mu_V + Paras.alpha0);
            ICs(33) = 0.01 * ProjStrat.TsetseReintro * Data.N_H * Paras.alpha0 / (Paras.mu_V + Paras.alpha0);
        end % M9
        
        parameter = Paras;
        parameter.f = f';
        parameter.meff = meff;
        parameter.mu_H = [Paras.mu_H_FY Paras.mu_H_FW Paras.mu_H_FP Paras.mu_H_MY Paras.mu_H_MW Paras.mu_H_MP]';
        parameter.sigma_H = [Paras.sigma_H*ones(1,6)]';
        parameter.phi_H = [Paras.phi_H*ones(1,6)]';
        parameter.omega_H = [Paras.omega_H*ones(1,6)]';
        parameter.K_V = K_V;
        
        parameter.gamma_H = [yearlygamma_H(s)*ones(1,6)]';
        parameter.eta_H = [yearlyeta_H(s)*ones(1,6)]';
        parameter.p_targetdie = p_targetdie(s);
        parameter.TargetFreq = TargetFreq(s);
            
        if (Data.ModelScreeningTime(s) < Paras.VCstart && Paras.VCstart < Data.ModelScreeningTime(s+1)) || (Data.ModelScreeningTime(s) < Paras.VC_scaleback && Paras.VC_scaleback < Data.ModelScreeningTime(s+1)) || (Data.ModelScreeningTime(s) < ProjStrat.NewVCyear && ProjStrat.NewVCyear < Data.ModelScreeningTime(s+1))
            dT = [Paras.VCstart Paras.VC_scaleback ProjStrat.NewVCyear] - double(Data.ModelScreeningTime(s));
            Tbreak = 365 * min(dT(dT > 0));
            [t1, pop1] = ode45(@diffHATmodel, mod(Data.ModelScreeningTime(1), 1) * 365 + [sum(Data.ModelScreeningFreq(1 : s-1)) sum(Data.ModelScreeningFreq(1 : s-1)) + Tbreak], ICs, ODEoptions, parameter);
            
            %update VC
            parameter.p_targetdie = p_targetdie(s+1);
            parameter.TargetFreq = TargetFreq(s+1);
            [t2, pop2] = ode45(@diffHATmodel, mod(Data.ModelScreeningTime(1), 1) * 365 + [sum(Data.ModelScreeningFreq(1 : s-1)) + Tbreak sum(Data.ModelScreeningFreq(1 : s))], pop1(end,:), ODEoptions, parameter);
            
            pop = [pop1; pop2(2:end,:)];
            t = [t1; t2(2:end,:)];
        else

%             mod(Data.ModelScreeningTime(1), 1) * 365 + [sum(Data.ModelScreeningFreq(1 : s-1)) sum(Data.ModelScreeningFreq(1 : s))]
%             ICs
%             for ls = 1:37
%                 if ICs(ls) < 0
%                     ls
%                     ICs(ls)
%                 end
%             end
            [t, pop] = ode45(@diffHATmodel, mod(Data.ModelScreeningTime(1), 1) * 365 + [sum(Data.ModelScreeningFreq(1 : s-1)) sum(Data.ModelScreeningFreq(1 : s))], ICs, ODEoptions, parameter);
        end
        Pop = [Pop; pop];
        T = [T; t];
       
        S_H = pop(:, 1:6);
        E_H = pop(:, 7:12);
        I1_H = pop(:, 13:18);
        I2_H = pop(:, 19:24);
        R_H = pop(:, 25:30);
        I_V = pop(:, 37);
        dNH = S_H + E_H + I1_H + I2_H + R_H;
        dNH(dNH == 0) = 1;
                
        % Person years infected in S1 and S2 %% 365 what if multiple
        % screening per year
        PersonYrs1(s) = trapz(t, sum(I1_H,2)) / 365;
        PersonYrs2(s) = trapz(t, sum(I2_H,2)) / 365;

        PersonYrsAgeY_1(s) = trapz(t, sum(I1_H(:,[1,4]),2)) / 365;
        PersonYrsAgeW_1(s) = trapz(t, sum(I1_H(:,[2,5]),2)) / 365;
        PersonYrsAgeP_1(s) = trapz(t, sum(I1_H(:,[3,6]),2)) / 365;
        PersonYrsGenderM_1(s) = trapz(t, sum(I1_H(:,1:3),2)) / 365;
        PersonYrsGenderF_1(s) = trapz(t, sum(I1_H(:,4:6),2)) / 365;

        PersonYrsAgeY_2(s) = trapz(t, sum(I2_H(:,[1,4]),2)) / 365;
        PersonYrsAgeW_2(s) = trapz(t, sum(I2_H(:,[2,5]),2)) / 365;
        PersonYrsAgeP_2(s) = trapz(t, sum(I2_H(:,[3,6]),2)) / 365;
        PersonYrsGenderM_2(s) = trapz(t, sum(I2_H(:,1:3),2)) / 365;
        PersonYrsGenderF_2(s) = trapz(t, sum(I2_H(:,4:6),2)) / 365;

        PersonYrsGenderF(s) = PersonYrsGenderF_1(s) + PersonYrsGenderF_2(s) ;
        PersonYrsGenderM(s) = PersonYrsGenderM_1(s) + PersonYrsGenderM_2(s) ;
        PersonYrsAgeY(s) = PersonYrsAgeY_1(s) + PersonYrsAgeY_2(s) ;
        PersonYrsAgeW(s) = PersonYrsAgeW_1(s) + PersonYrsAgeW_2(s) ;
        PersonYrsAgeP(s) = PersonYrsAgeP_1(s) + PersonYrsAgeP_2(s) ;
        
        
        % Passive detections (S1/S2) each time interval
        Passive1(s) = yearlyeta_H(s) * PersonYrs1(s) * 365;
        Passive2(s) = S2treatment(s) * PersonYrs2(s) * 365;

        PassiveAgeY(s) = (yearlyeta_H(s) * PersonYrsAgeY_1(s) + S2treatment(s) * PersonYrsAgeY_2(s))* 365;
        PassiveAgeW(s) = (yearlyeta_H(s) * PersonYrsAgeW_1(s) + S2treatment(s) * PersonYrsAgeW_2(s))* 365;
        PassiveAgeP(s) = (yearlyeta_H(s) * PersonYrsAgeP_1(s) + S2treatment(s) * PersonYrsAgeP_2(s))* 365;
        PassiveGenderF(s) = (yearlyeta_H(s) * PersonYrsGenderF_1(s) + S2treatment(s) * PersonYrsGenderF_2(s))* 365;
        PassiveGenderM(s) = (yearlyeta_H(s) * PersonYrsGenderM_1(s) + S2treatment(s) * PersonYrsGenderM_2(s))* 365;

        % Deaths
        Deaths(s) = Paras.death * PersonYrs2(s) * 365;
        
        % New infections (influx into I_1H) %%all human combined???
        FOI = Paras.alpha * meff * bsxfun(@times, f, S_H./dNH) .* I_V;
        NewInf(s) = sum(trapz(t, FOI));

        
    end
    
% Output
    % All timepoints)
    Classes = array2table([floor(Data.ModelScreeningTime(1))+T/365 Pop],... 
              'VariableNames', {'Time', 'S_H_FY', 'S_H_FW', 'S_H_FP', 'S_H_MY', 'S_H_MW', 'S_H_MP',...
                                        'E_H_FY', 'E_H_FW', 'E_H_FP', 'E_H_MY', 'E_H_MW', 'E_H_MP',...
                                    'I1_H_FY', 'I1_H_FW', 'I1_H_FP', 'I1_H_MY', 'I1_H_MW', 'I1_H_MP',...
                                    'I2_H_FY', 'I2_H_FW', 'I2_H_FP', 'I2_H_MY', 'I2_H_MW', 'I2_H_MP',...
                                    'R_H_FY', 'R_H_FW', 'R_H_FP', 'R_H_MY', 'R_H_MW', 'R_H_MP',...                          
                                'P_V', 'S_V', 'G_V', 'E1_V', 'E2_V', 'E3_V', 'I_V'});
%    T(1)
%                            Classes([1 2 end],:)
    %Y
    for y = 1 : Y%length(Data.Years) % Y-1
        s = find(floor(double(Data.ModelScreeningTime)) == Data.Years(y));
        ActiveM1(y) = sum(Active1(s));
        ActiveM2(y) = sum(Active2(s));
        ActiveMFP(y) = sum(ActiveFP(s));
        PassiveM1(y) = sum(Passive1(s));
        PassiveM2(y) = sum(Passive2(s));

        ActiveMAgeY(y) = sum(ActiveAgeY(s));
        ActiveMAgeW(y) = sum(ActiveAgeW(s));
        ActiveMAgeP(y) = sum(ActiveAgeP(s));
        ActiveMGenderF(y) = sum(ActiveGenderF(s));
        ActiveMGenderM(y) = sum(ActiveGenderM(s));

        PassiveMAgeY(y) = sum(PassiveAgeY(s));
        PassiveMAgeW(y) = sum(PassiveAgeW(s));
        PassiveMAgeP(y) = sum(PassiveAgeP(s));
        PassiveMGenderF(y) = sum(PassiveGenderF(s));
        PassiveMGenderM(y) = sum(PassiveGenderM(s));

        DeathsM(y) = sum(Deaths(s));

        PersonYrsMAgeY(y) = sum(PersonYrsAgeY(s));
        PersonYrsMAgeW(y) = sum(PersonYrsAgeW(s));
        PersonYrsMAgeP(y) = sum(PersonYrsAgeP(s));
        PersonYrsMGenderF(y) = sum(PersonYrsGenderF(s));
        PersonYrsMGenderM(y) = sum(PersonYrsGenderM(s));

        NewInfM(y) = sum(NewInf(s));
        PerfectSpec(y) = mean(Spec(s));
        Transmission(y) = mean(Trans(s));
        RDTM1(y) = sum(RDT1(s));
        RDTM2(y) = sum(RDT2(s));
        RDTMFP(y) = sum(RDTFP(s));

        RDTMAgeY(y) = sum(RDTAgeY(s));
        RDTMAgeW(y) = sum(RDTAgeW(s));
        RDTMAgeP(y) = sum(RDTAgeP(s));
        RDTMGenderF(y) = sum(RDTGenderF(s));
        RDTMGenderM(y) = sum(RDTGenderM(s));

    end
   
   % disp(ActiveM1)
   % disp(ActiveMGenderM)
    Aggregate = table(Data.Years', ActiveM1', ActiveM2', ActiveMFP',ActiveMGenderM',ActiveMGenderF',ActiveMAgeY',ActiveMAgeW',ActiveMAgeP', PassiveM1', PassiveM2',PassiveMAgeY',PassiveMAgeW',PassiveMAgeP',PassiveMGenderM',PassiveMGenderF', DeathsM', PersonYrsM1', PersonYrsM2', PersonYrsMGenderM',PersonYrsMGenderF',PersonYrsMAgeY',PersonYrsMAgeW',PersonYrsMAgeP', NewInfM', Transmission', PerfectSpec', RDTM1', RDTM2', RDTMFP',RDTMAgeY',RDTMAgeW',RDTMAgeP',RDTMGenderM',RDTMGenderF', ...
                'VariableNames', {'Year', 'ActiveM1', 'ActiveM2', 'ActiveMFP','ActiveMGenderM','ActiveMGenderF','ActiveMAgeY','ActiveMAgeW','ActiveMAgeP', 'PassiveM1', 'PassiveM2','PassiveMAgeY','PassiveMAgeW','PassiveMAgeP','PassiveMGenderM','PassiveMGenderF', 'DeathsM', 'PersonYrsM1', 'PersonYrsM2', 'PersonYrsMGenderM','PersonYrsMGenderF','PersonYrsMAgeY','PersonYrsMAgeW','PersonYrsMAgeP', 'NewInfM', 'Transmission', 'PerfectSpec', 'RDTM1', 'RDTM2', 'RDTMFP','RDTMAgeY','RDTMAgeW','RDTMAgeP','RDTMGenderM','RDTMGenderF'}); % M9
    %Aggregate = struct('Year', Data.Year, 'ActiveM1', ActiveM1, 'ActiveM2', ActiveM2, 'PassiveM1', PassiveM1, 'PassiveM2',PassiveM2, 'DeathsM', DeathsM,...
    %                   'PersonYrsM1', PersonYrsM1, 'PersonYrsM2', PersonYrsM2, 'NewInfM', NewInfM);
    
   
% Main ODE code
function dPop = diffHATmodel(t, pop, parameter) % M9

%Compute vector reduction function
%f_T = parameter.p_targetdie * (1 - sigmf(mod(t,365/parameter.TargetFreq),[25/365 0.35*365]));
if parameter.p_targetdie==0
    f_T = 0;
else
    f_T = parameter.p_targetdie * (1 - 1/(1+exp(-25/365*(mod(t - parameter.VC_t0, 365/parameter.TargetFreq)-0.35*365))));
end % Fix_f_t

%Get populations from inputs
S_H = pop(1:6); E_H = pop(7:12); I1_H = pop(13:18); I2_H = pop(19:24); R_H = pop(25:30);
P_V = pop(31); S_V = pop(32); G_V = pop(33); E1_V=pop(34); E2_V=pop(35); E3_V=pop(36); I_V=pop(37);

N_H = S_H + E_H + I1_H + I2_H + R_H;
N_V = S_V + G_V + E1_V + E2_V + E3_V + I_V;

dNH = N_H;
dNH(N_H==0) = 1;

%Human infection dynamics 
age_out_rate = [parameter.l_Y parameter.l_W 0 parameter.l_Y parameter.l_W 0]';
age_in_rate = [0 parameter.l_Y parameter.l_W 0 parameter.l_Y parameter.l_W]';

in_idx = [3,1,2,6,4,5];

dS_H = age_in_rate .* S_H(in_idx) + parameter.omega_H .* R_H  + parameter.omega_IB .* I1_H - I_V * parameter.alpha * parameter.meff .* parameter.f .* S_H ./ dNH - parameter.mu_H .* S_H - age_out_rate .* S_H;
dE_H = age_in_rate .* E_H(in_idx) + I_V * parameter.alpha * parameter.meff .* parameter.f .* S_H ./ dNH - (parameter.sigma_H + parameter.mu_H) .* E_H - age_out_rate .* E_H;
dI1_H = age_in_rate .* I1_H(in_idx) + parameter.sigma_H .* E_H * parameter.p_BS - (parameter.eta_H + parameter.phi_H + parameter.mu_H + parameter.omega_IB) .* I1_H  - age_out_rate .* I1_H;
dI2_H = age_in_rate .* I2_H(in_idx) + parameter.phi_H .* I1_H -  (parameter.gamma_H + parameter.mu_H) .* I2_H  - age_out_rate .* I2_H;
dR_H =  age_in_rate .* R_H(in_idx) + parameter.eta_H .* I1_H + parameter.gamma_H .* I2_H - (parameter.omega_H + parameter.mu_H) .* R_H - age_out_rate .* R_H;

dS_H(1) = dS_H(1) + sum(parameter.mu_H(1:3) .* N_H(1:3));
dS_H(4) = dS_H(4) + sum(parameter.mu_H(4:6) .* N_H(4:6));


%Tsetse Infection dynamics
%Pupa
%                                                  
dP_V = parameter.B_V * N_V - (parameter.xi_V + P_V/parameter.K_V) * P_V;
%Teneral
dS_V = parameter.xi_V * parameter.p_survive * P_V - parameter.alpha * S_V - parameter.mu_V * S_V;
%Non-teneral
dG_V = parameter.alpha * (1 - f_T) * (1 - sum(parameter.f .* (I1_H + I2_H) ./ dNH) * parameter.p_V) * S_V - parameter.alpha * ((1 - f_T) * parameter.epsilon * sum(parameter.f .* (I1_H + I2_H) ./ dNH) * parameter.p_V + f_T) * G_V - parameter.mu_V * G_V;
%Exposed
dE1_V = parameter.alpha * (1 - f_T) * sum(parameter.f .* (I1_H + I2_H) ./ dNH) * parameter.p_V * (S_V + parameter.epsilon * G_V) - 3 * parameter.sigma_V * E1_V - (parameter.mu_V + parameter.alpha * f_T) * E1_V;
dE2_V = 3 * parameter.sigma_V * E1_V - (3 * parameter.sigma_V + parameter.mu_V + parameter.alpha * f_T) * E2_V;
dE3_V = 3 * parameter.sigma_V * E2_V - (3 * parameter.sigma_V + parameter.mu_V + parameter.alpha * f_T) * E3_V;
%Infected
dI_V= 3 * parameter.sigma_V * E3_V - (parameter.mu_V + parameter.alpha * f_T) * I_V;

dPop = [dS_H; dE_H; dI1_H; dI2_H; dR_H; dP_V; dS_V; dG_V; dE1_V; dE2_V; dE3_V; dI_V];

function demographic_pop = get_demog_partition(this_N)
float_pop = k*this_N;
demographic_pop = floor(float_pop);
pop_remainder = float_pop - demographic_pop;
[~,idx] = sort(pop_remainder);
j = 6;
while sum(demographic_pop ~= this_N)
 demographic_pop(idx(j)) = demographic_pop(idx(j)) + 1;
 j = j - 1;
end


