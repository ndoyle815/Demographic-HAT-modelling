
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                               %
%   This code runs forcasting of all strategies considered in the Warwick HAT model                             %
%                                                                                                               %
%   Inputs:                                                                                                     %
%       Data - structure containing location-specific historical data                                           %
%       Paras - structure containing location-specific parameters (fixed, fitted and intervention parameters)   %
%       Strategy - table containing parameters of all strategies                                                %
%       samples - number denoting sample size from ODE                                                          %
%       ReactiveParameters - table containing parameters of all reactive strategies                             %
%                                                                                                               %
%   Outputs:                                                                                                    %
%       Outputs - structure containing yearly aggregated outputs, reactive information and ICs for forcasting   %
%                                                                                                               %
%   Functions required: GetEndemicEq & ODEHATmodel & Cbetabinornd                                               %
%                                                                                                               %
%                                                                                                               %
%                                                                                                               %
%   Notes from Ching-I: WHO classic algorithm (no 4th year screening) in RS isn't an option in this version     %
%                                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Outputs = Projection(RunInfo, Data, Paras, Strategy, samples, ReactiveParameters)
    NumStrat = size(Strategy, 1) - RunInfo.StratMin;
    NumReact = size(ReactiveParameters, 1) - 1;
    Years = Data.Years(1) : max([Data.Years(end); Strategy.SIMyear]); % all strategis must have the same SIMyear
    MtoAbsScaling = Data.PopGrowth .^ double(Years - Data.PopSizeYear);
    Pop = round(Data.N_H * MtoAbsScaling);
    NumYear = length(Years);
    NumYear0 = length(Data.Years);
    
    [Active1, Active2, Active1FP, Active2FP,ActiveAgeY, ActiveAgeW, ActiveAgeP, ActiveGenderM,ActiveGenderF, Passive1, Passive2,PassiveAgeY, PassiveAgeW, PassiveAgeP, PassiveGenderM,PassiveGenderF, Deaths, PersonYrs1, PersonYrs2,PersonYrsAgeY, PersonYrsAgeW, PersonYrsAgeP, PersonYrsGenderM,PersonYrsGenderF, NewInf, PerfectSpec, Transmission, RDT1, RDT2, RDTFP,RDTAgeY, RDTAgeW,RDTAgeP, RDTGenderF,RDTGenderM] = deal(zeros(1, NumYear, max(NumStrat, 1))); % M9
    [YEPHP, YEOT] = deal(zeros(1, max(NumStrat, 1)));
    [SampledActive1, SampledActive2, SampledActive1FP, SampledActive2FP, SampledActiveAgeY, SampledActiveAgeW, SampledActiveAgeP, SampledActiveGenderM,SampledActiveGenderF, SampledPassive1, SampledPassive2,SampledPassiveAgeY, SampledPassiveAgeW, SampledPassiveAgeP, SampledPassiveGenderM,SampledPassiveGenderF, SampledDeaths, SampledRDT1, SampledRDT2, SampledRDTFP, SampledRDTAgeY, SampledRDTAgeW, SampledRDTAgeP, SampledRDTGenderF, SampledRDTGenderM] = deal(zeros(samples, NumYear, max(NumStrat, 1))); % M9
    SampledYEPHP = zeros(samples, max(NumStrat, 1));

    %%% Fitted part
    % Get equilibrium ICs
    [meff, ICs] = GetEndemicEq(Data.N_H, Paras, Data); 

    % Run fitted dynamics
    ProjStrat = table2struct(Strategy('Strat0',:));
    Paras.death = (1-Paras.u) * Paras.gamma_H;
    Paras.effcy = 1;
    
    %Data.ScreeningSpecificity
    [Classes0, Aggregate0] = ODEHATmodel(meff, ICs, Data, Paras, ProjStrat);
    MtoAbsScaling0 = Data.PopGrowth .^ double(Data.Years - Data.PopSizeYear);
    
    if RunInfo.StratMin == 1
        % Sampling from fitted dynamics
        Pop0 = round(Data.N_H * MtoAbsScaling0);
        ScaledPeopleScreened0 = Data.PeopleScreened ./ MtoAbsScaling0;
        AS = Data.PeopleScreened;
    

        Disp_act = Paras.disp_act*(Aggregate0.ActiveM1' + Aggregate0.ActiveM2') ./ (Aggregate0.ActiveM1' + Aggregate0.ActiveM2' + Aggregate0.ActiveMFP'); %M9
        ActiveS = Cbetabinornd(repmat(Data.PeopleScreened, samples, 1), repmat((Aggregate0.ActiveM1' + Aggregate0.ActiveM2' + Aggregate0.ActiveMFP') ./ ScaledPeopleScreened0, samples, 1), repmat(Disp_act, samples, 1)); % M9
        ActiveS1 = Cbinornd(ActiveS, repmat((Aggregate0.ActiveM1' + Paras.S1givenFP * Aggregate0.ActiveMFP') ./ (Aggregate0.ActiveM1' + Aggregate0.ActiveM2'+ Aggregate0.ActiveMFP'), samples, 1)); % M9
        ActiveS1FP = Cbinornd(ActiveS1, repmat(Paras.S1givenFP * Aggregate0.ActiveMFP' ./ (Aggregate0.ActiveM1' + Paras.S1givenFP * Aggregate0.ActiveMFP'), samples, 1));
        ActiveS2 = Cbinornd(ActiveS - ActiveS1, repmat(Aggregate0.ActiveM2' + (1-Paras.S1givenFP) * Aggregate0.ActiveMFP', samples, 1) ./ (repmat(Aggregate0.ActiveM2' + (1-Paras.S1givenFP) * Aggregate0.ActiveMFP', samples, 1))); % M9
        ActiveS2FP = Cbinornd(ActiveS2, repmat((1 - Paras.S1givenFP) * Aggregate0.ActiveMFP' ./ (Aggregate0.ActiveM2' + (1 - Paras.S1givenFP) * Aggregate0.ActiveMFP'), samples, 1));     
 
       
        ActiveSAgeY = Cbinornd(ActiveS, repmat(Aggregate0.ActiveMAgeY' ./ (Aggregate0.ActiveMAgeY' + Aggregate0.ActiveMAgeW' +  Aggregate0.ActiveMAgeP'), samples, 1));
        ActiveSAgeW = Cbinornd(ActiveS, repmat(Aggregate0.ActiveMAgeW' ./ (Aggregate0.ActiveMAgeY' + Aggregate0.ActiveMAgeW' +  Aggregate0.ActiveMAgeP'), samples, 1));
      
 
        PassiveS = Cbetabinornd(repmat(Pop0, samples, 1), repmat((Aggregate0.PassiveM1' + Aggregate0.PassiveM2') / Data.N_H, samples, 1), repmat(Paras.disp_pass, samples, NumYear0));
        PassiveS1 = Cbinornd(PassiveS, repmat(Aggregate0.PassiveM1' ./ (Aggregate0.PassiveM1' + Aggregate0.PassiveM2'), samples, 1));
        DeathsS = Cbinornd(repmat(Pop0, samples, 1), repmat(Aggregate0.DeathsM' / Data.N_H, samples, 1));
        
        ActiveSAgeP = ActiveS - ActiveSAgeY - ActiveSAgeW;
        PassiveSAgeY = Cbinornd(PassiveS, repmat(Aggregate0.PassiveMAgeY' ./ (Aggregate0.PassiveMAgeY' + Aggregate0.PassiveMAgeW' +  Aggregate0.PassiveMAgeP'), samples, 1));
        PassiveSAgeW = Cbinornd(PassiveS, repmat(Aggregate0.PassiveMAgeW' ./ (Aggregate0.PassiveMAgeY' + Aggregate0.PassiveMAgeW' +  Aggregate0.PassiveMAgeP'), samples, 1));
        PassiveSAgeP = PassiveS - PassiveSAgeY - PassiveSAgeW;
    
        ActiveSGenderF = Cbinornd(ActiveS, repmat(Aggregate0.ActiveMGenderF' ./ (Aggregate0.ActiveMGenderF' + Aggregate0.ActiveMGenderM'), samples, 1));
        ActiveSGenderM = ActiveS - ActiveSGenderF;

        PassiveSGenderF = Cbinornd(PassiveS, repmat(Aggregate0.PassiveMGenderF' ./ (Aggregate0.PassiveMGenderF' + Aggregate0.PassiveMGenderM'), samples, 1));
        PassiveSGenderM = PassiveS - PassiveSGenderF;
       


        for s = 1 : max(NumStrat, 1)
            % fitted ODE
            Active1(:,1:NumYear0,s) = (Aggregate0.ActiveM1 + Paras.S1givenFP * Aggregate0.ActiveMFP)' .* MtoAbsScaling0;
            Active2(:,1:NumYear0,s) = (Aggregate0.ActiveM2 + (1 - Paras.S1givenFP) * Aggregate0.ActiveMFP)' .* MtoAbsScaling0;
            Active1FP(:,1:NumYear0,s) = Paras.S1givenFP * Aggregate0.ActiveMFP' .* MtoAbsScaling0;
            Active2FP(:,1:NumYear0,s) = (1 - Paras.S1givenFP) * Aggregate0.ActiveMFP' .* MtoAbsScaling0;
            
            ActiveAgeY(:,1:NumYear0,s) = Aggregate0.ActiveMAgeY' .* MtoAbsScaling0;
            ActiveAgeW(:,1:NumYear0,s) = Aggregate0.ActiveMAgeW' .* MtoAbsScaling0;
            ActiveAgeP(:,1:NumYear0,s) = Aggregate0.ActiveMAgeP' .* MtoAbsScaling0;
            ActiveGenderF(:,1:NumYear0,s) = Aggregate0.ActiveMGenderF' .* MtoAbsScaling0;
            ActiveGenderM(:,1:NumYear0,s) = Aggregate0.ActiveMGenderM' .* MtoAbsScaling0;

            PassiveAgeY(:,1:NumYear0,s) = Aggregate0.PassiveMAgeY' .* MtoAbsScaling0;
            PassiveAgeW(:,1:NumYear0,s) = Aggregate0.PassiveMAgeW' .* MtoAbsScaling0;
            PassiveAgeP(:,1:NumYear0,s) = Aggregate0.PassiveMAgeP' .* MtoAbsScaling0;
            PassiveGenderF(:,1:NumYear0,s) = Aggregate0.PassiveMGenderF' .* MtoAbsScaling0;
            PassiveGenderM(:,1:NumYear0,s) = Aggregate0.PassiveMGenderM' .* MtoAbsScaling0;

            Passive1(:,1:NumYear0,s) = Aggregate0.PassiveM1' .* MtoAbsScaling0;
            Passive2(:,1:NumYear0,s) = Aggregate0.PassiveM2' .* MtoAbsScaling0;
            Deaths(:,1:NumYear0,s) = Aggregate0.DeathsM' .* MtoAbsScaling0;


            PersonYrs1(:,1:NumYear0,s) = Aggregate0.PersonYrsM1' .* MtoAbsScaling0;
            PersonYrs2(:,1:NumYear0,s) = Aggregate0.PersonYrsM2' .* MtoAbsScaling0;
            PersonYrsAgeY(:,1:NumYear0,s) = Aggregate0.PersonYrsMAgeY' .* MtoAbsScaling0;
            PersonYrsAgeW(:,1:NumYear0,s) = Aggregate0.PersonYrsMAgeW' .* MtoAbsScaling0;
            PersonYrsAgeP(:,1:NumYear0,s) = Aggregate0.PersonYrsMAgeP' .* MtoAbsScaling0;
            PersonYrsGenderF(:,1:NumYear0,s) = Aggregate0.PersonYrsMGenderF' .* MtoAbsScaling0;
            PersonYrsGenderM(:,1:NumYear0,s) = Aggregate0.PersonYrsMGenderM' .* MtoAbsScaling0;

            NewInf(:,1:NumYear0,s) = Aggregate0.NewInfM' .* MtoAbsScaling0;
            PerfectSpec(:,1:NumYear0,s) = Aggregate0.PerfectSpec';
            Transmission(:,1:NumYear0,s) = Aggregate0.Transmission';
            
            % fitted sampling
            SampledActive1(:, 1:NumYear0, s) = ActiveS1;
            SampledActive2(:, 1:NumYear0, s) = ActiveS2;
            SampledActive1FP(:, 1:NumYear0, s) = ActiveS1FP;
            SampledActive2FP(:, 1:NumYear0, s) = ActiveS2FP;

            SampledActiveAgeY(:, 1:NumYear0, s) = ActiveSAgeY;
            SampledActiveAgeW(:, 1:NumYear0, s) = ActiveSAgeW;
            SampledActiveAgeP(:, 1:NumYear0, s) = ActiveSAgeP;
            SampledActiveGenderF(:, 1:NumYear0, s) = ActiveSGenderF;
            SampledActiveGenderM(:, 1:NumYear0, s) = ActiveSGenderM;

            SampledPassive1(:, 1:NumYear0, s) = PassiveS1;
            SampledPassive2(:, 1:NumYear0, s) = PassiveS - PassiveS1;

            SampledPassiveAgeY(:, 1:NumYear0, s) = PassiveSAgeY;
            SampledPassiveAgeW(:, 1:NumYear0, s) = PassiveSAgeW;
            SampledPassiveAgeP(:, 1:NumYear0, s) = PassiveSAgeP;
            SampledPassiveGenderF(:, 1:NumYear0, s) = PassiveSGenderF;
            SampledPassiveGenderM(:, 1:NumYear0, s) = PassiveSGenderM;

            SampledDeaths(:, 1:NumYear0, s) = DeathsS;
        end
    else
        fitted = load(RunInfo.ProjFilePath('Fitted'));
        SampledActive1(:, 1:NumYear0, :) = repmat(fitted.SampledActive1(Paras.PostRows, :), 1, 1, max(NumStrat, 1));
        SampledActive2(:, 1:NumYear0, :) = repmat(fitted.SampledActive2(Paras.PostRows, :), 1, 1, max(NumStrat, 1));

        SampledActiveAgeY(:, 1:NumYear0, :) = repmat(fitted.SampledActiveAgeY(Paras.PostRows, :), 1, 1, max(NumStrat, 1));
        SampledActiveAgeW(:, 1:NumYear0, :) = repmat(fitted.SampledActiveAgeW(Paras.PostRows, :), 1, 1, max(NumStrat, 1));
        SampledActiveAgeP(:, 1:NumYear0, :) = repmat(fitted.SampledActiveAgeP(Paras.PostRows, :), 1, 1, max(NumStrat, 1));
        SampledActiveGenderM(:, 1:NumYear0, :) = repmat(fitted.SampledActiveGenderM(Paras.PostRows, :), 1, 1, max(NumStrat, 1));
        SampledActiveGenderF(:, 1:NumYear0, :) = repmat(fitted.SampledActiveGenderF(Paras.PostRows, :), 1, 1, max(NumStrat, 1));

        SampledPassive1(:, 1:NumYear0, :) = repmat(fitted.SampledPassive1(Paras.PostRows, :), 1, 1, max(NumStrat, 1));
        SampledPassive2(:, 1:NumYear0, :) = repmat(fitted.SampledPassive2(Paras.PostRows, :), 1, 1, max(NumStrat, 1));

        SampledPassiveAgeY(:, 1:NumYear0, :) = repmat(fitted.SampledPassiveAgeY(Paras.PostRows, :), 1, 1, max(NumStrat, 1));
        SampledPassiveAgeW(:, 1:NumYear0, :) = repmat(fitted.SampledPassiveAgeW(Paras.PostRows, :), 1, 1, max(NumStrat, 1));
        SampledPassiveAgeP(:, 1:NumYear0, :) = repmat(fitted.SampledPassiveAgeP(Paras.PostRows, :), 1, 1, max(NumStrat, 1));
        SampledPassiveGenderM(:, 1:NumYear0, :) = repmat(fitted.SampledPassiveGenderM(Paras.PostRows, :), 1, 1, max(NumStrat, 1));
        SampledPassiveGenderF(:, 1:NumYear0, :) = repmat(fitted.SampledPassiveGenderF(Paras.PostRows, :), 1, 1, max(NumStrat, 1));

    end

    

    
    %%% Projection part
    % Get ICs
    pop = Classes0(end, :);
    ICs = {[pop.S_H_FY, pop.S_H_FW, pop.S_H_FP, pop.S_H_MY,  pop.S_H_MW,  pop.S_H_MP], [pop.E_H_FY, pop.E_H_FW, pop.E_H_FP, pop.E_H_MY,  pop.E_H_MW,  pop.E_H_MP],...
           [pop.I1_H_FY, pop.I1_H_FW, pop.I1_H_FP, pop.I1_H_MY,  pop.I1_H_MW,  pop.I1_H_MP], [pop.I2_H_FY, pop.I2_H_FW, pop.I2_H_FP, pop.I2_H_MY,  pop.I2_H_MW,  pop.I2_H_MP],[pop.R_H_FY, pop.R_H_FW, pop.R_H_FP, pop.R_H_MY,  pop.R_H_MW,  pop.R_H_MP],...
            pop.P_V, pop.S_V, pop.G_V, pop.E1_V, pop.E2_V, pop.E3_V, pop.I_V};
   
    % Update alpha
    if Aggregate0.NewInfM(end) * MtoAbsScaling0(end) < Paras.EOTthreshold && ...
       (pop.I1_H_FY+ pop.I1_H_FW+ pop.I1_H_FP+ pop.I1_H_MY+  pop.I1_H_MW+  pop.I1_H_MP+ pop.I2_H_FY+ pop.I2_H_FW+ pop.I2_H_FP+ pop.I2_H_MY+  pop.I2_H_MW+  pop.I2_H_MP) * MtoAbsScaling0(end) < Paras.EOTthreshold % M9
        alpha = 0;
    else
        alpha = Paras.alpha;
    end
       
    Data.Years = Data.Years(end) + 1 : max([Data.Years(end); Strategy.SIMyear]); % all strategis must have the same SIMyear
    MtoAbsScaling1 = Data.PopGrowth .^ double(Data.Years - Data.PopSizeYear);
    Pop1 = round(Data.N_H * MtoAbsScaling1);
    
    ReactInfo = [];
    for s = 1 : NumStrat % 1:3:7
        S = ['Strat', num2str(s + RunInfo.StratMin - 1)];
        % Strategy Parameters
        ProjStrat = table2struct(Strategy(S,:));
        
        Data.PeopleScreened = [];
        ScreeningCapacity = [];
        NewASnum = split(ProjStrat.NewASnum, '_');
        NewASstrat = split(ProjStrat.NewASstrat, '_');
        NewASyear = [ProjStrat.NewASyear Data.Years(end)+1];

        for i = 1 : length(NewASnum)
            yrs = NewASyear(i+1) - NewASyear(i);
            switch NewASnum{i}
                case 'mean'
                    PeopleScreened = Data.MeanPeopleScreened * ones(1, yrs);
                case 'max'
                    PeopleScreened = Data.MaxPeopleScreened * ones(1, yrs);
                case 'intensified'
                    PeopleScreened = Data.IntensifiedPeopleScreened * ones(1, yrs);
                case 'off'
                    PeopleScreened = zeros(1, yrs);
                case 'stop'
                    PeopleScreened = zeros(1, yrs);
                otherwise
                    if contains(NewASnum{i}, '%') % percentage 'x%'
                        pct = NewASnum{i};
                        PeopleScreened = round(0.01 * str2num(pct(1 : end - 1)) * Data.N_H * MtoAbsScaling0(end - 2 - Data.DeltaGapYear)) * ones(1, yrs);
                    else
                        PeopleScreened = str2num(NewASnum{i}) * ones(1, yrs);
                    end % Isangi
            end

            Data.PeopleScreened = [Data.PeopleScreened PeopleScreened];

            % ELLIOT: The bit after "switch" here I'm not really sure how we fix, i
% have a feeling it is skipped though. Suspect SC  = screening capacity
%             switch NewASstrat{i}
%                 case 'traditional'
%                     SC = Paras.k1 * ones(1, yrs);
%                 case 'equal'
            SC = ones(1, yrs);
%                 case 'high'
%                     SC = (Paras.k1 + Paras.k2 + Paras.k4) * ones(1, yrs);
%             end
            ScreeningCapacity = [ScreeningCapacity SC];
        end
        
        
%         switch ProjStrat.NewASnum %ProjStrat{2}
%             case 'mean'
%                 Data.PeopleScreened = Data.MeanPeopleScreened * ones(1, length(Data.Years));
%             case 'max'
%                 Data.PeopleScreened = Data.MaxPeopleScreened * ones(1, length(Data.Years));
%             case 'intensified'
%                 Data.PeopleScreened = Data.IntensifiedPeopleScreened * ones(1, length(Data.Years));
%             %case 'off'
%             %    Data.PeopleScreened = zeros(1, length(Data.Years));
%             otherwise % percentage 
%                 Data.PeopleScreened = round(0.01 * str2num(ProjStrat.NewASnum(1 : end - 1)) * Data.N_H * MtoAbsScaling0(end - 2 - Data.DeltaGapYear)) * ones(1, length(Data.Years));
%         end
        
%         % Modify PeopleScreened when screening stops forever from NewASyear 
%         if strcmp(ProjStrat.NewASstrat, 'stop')
%             Data.PeopleScreened(Data.Years >= ProjStrat.NewASyear) = 0;
%         end

        ScaledPeopleScreened1 = Data.PeopleScreened ./ MtoAbsScaling1;
        ModelScreeningFreq = [];
        ModelScreeningTime = [];
        ModelPeopleScreened = [];
        ExpectedScreeningTimes = ceil(ScaledPeopleScreened1 ./ ScreeningCapacity / Data.N_H);
        ExpectedScreeningTimes(ExpectedScreeningTimes == 0) = 1;
        for y = 1 : length(Data.Years)
            freq = 365 / ExpectedScreeningTimes(y);
            time = 1 / ExpectedScreeningTimes(y);
            number = ScaledPeopleScreened1(y) / ExpectedScreeningTimes(y);
            
            ModelScreeningFreq = [ModelScreeningFreq freq * ones(1, ExpectedScreeningTimes(y))];
            ModelScreeningTime = [ModelScreeningTime Data.Years(y) + (0 : ExpectedScreeningTimes(y) - 1) * time];
            ModelPeopleScreened = [ModelPeopleScreened number * ones(1, ExpectedScreeningTimes(y))];
            
%             if ScaledPeopleScreened1(y) < Data.N_H * Paras.ScreeningCapacity
%                 ModelScreeningFreq = [ModelScreeningFreq 365];
%                 ModelScreeningTime = [ModelScreeningTime Data.Years(y)];
%                 ModelPeopleScreened = [ModelPeopleScreened round(ScaledPeopleScreened1(y))];
%             else
%                 ModelScreeningFreq = [ModelScreeningFreq 365/2 365/2];
%                 ModelScreeningTime = [ModelScreeningTime Data.Years(y) Data.Years(y)+0.5];
%                 ModelPeopleScreened = [ModelPeopleScreened round(ScaledPeopleScreened1(y)/2) round(ScaledPeopleScreened1(y)/2)];
%             end
        end
        Data.ModelScreeningFreq = ModelScreeningFreq;
        Data.ModelScreeningTime = ModelScreeningTime;
        Data.ModelPeopleScreened = ModelPeopleScreened;
        
        % Update year_spec_100pct if year_spec_100pct == 0 and NewSPECyear ~= 0
        Paras.year_spec_100pct = Paras.year_spec_100pct + (Paras.year_spec_100pct == 0) * ProjStrat.NewSPECyear;

        % Update PPSstart, PPSend, NoPSstart and NoPSend
        switch ProjStrat.NewPS 
            case 'full'
                Paras.PPSstart = 0;
                Paras.PPSend = 0;
                Paras.NoPSstart = 0;
                Paras.NoPSend = 0;
            case 'partial'
                Paras.PPSstart = NewASyear(1);
                Paras.PPSend = NewASyear(end);
                Paras.NoPSstart = 0;
                Paras.NoPSend = 0;
            case 'no'
                Paras.PPSstart = 0;
                Paras.PPSend = 0;
                Paras.NoPSstart = NewASyear(1);
                Paras.NoPSend = NewASyear(end);
        end % Isangi

        Paras.alpha = alpha;
        VCstarted = Paras.VCstart > 0;
        %Paras.VC_t0 = mod(ProjStrat.NewVCyear, 1) * 365; % Fix_f_t
% %         ProjStrat.NewTargetFreq = VCstarted * (ProjStrat.NewTargetDie ~= 0) * Paras.TargetFreq + (1 - VCstarted) * ProjStrat.NewTargetFreq;
%         ProjStrat.NewVCyear = max([ProjStrat.NewVCyear ceil(Paras.VCstart + Paras.MinVC)]);
%         vc = (ProjStrat.NewTargetFreq ~= 0) + (ProjStrat.NewTargetDie ~= 0);
%         if vc == 1
%             %warning(['NewTargetFreq or NewTargetDie is set to 0 as its value conficts with the other in ', S])
%             ProjStrat.NewTargetFreq = 0;
%             ProjStrat.NewTargetDie = 0;
%         elseif vc == 2 && VCstarted
%             %warning('NewTargetFreq and NewTargetDie are replaced by its known values as VC started prior to projections')
%             ProjStrat.NewTargetFreq = Paras.TargetFreq;
%             ProjStrat.NewTargetDie = Paras.TargetDie;
%         end % ScalingUpVC

        % NATHAN: I suspect that our earlier bugs is that we have updated
        % the length of "Data.ModelPeopleScreened" here in Projection for
        % the 30 years of further simulation up to 2053, but in
        % ODEHATmodel.m we determine Turnout from
        % Data.ModelPeopleScreened_GC, which as of now is not updated here
        % in Projection. I believe different strategies given more time can
        % be designed to scale screening differently across populations,
        % but for now I will gague that this is uniform in proportions so I
        % can run basic strategies.

        % MY STRATEGY: If "TargetTest" in Strategy = 1, we conduct
        % screening on GC people based to their relative bite risk rGC:
        
        % begin as usual - uniform across demographics
        Data.ModelPeopleScreened_FY = ModelPeopleScreened * (Data.N_FY/Data.N_H);
        Data.ModelPeopleScreened_FW = ModelPeopleScreened * (Data.N_FW/Data.N_H);
        Data.ModelPeopleScreened_FP = ModelPeopleScreened * (Data.N_FP/Data.N_H);
        Data.ModelPeopleScreened_MY = ModelPeopleScreened * (Data.N_MY/Data.N_H);
        Data.ModelPeopleScreened_MW = ModelPeopleScreened * (Data.N_MW/Data.N_H);
        Data.ModelPeopleScreened_MP = ModelPeopleScreened * (Data.N_MP/Data.N_H);
        
        % if we want to target (who we believe) are the highest risk cohort
        if ProjStrat.TargetTest == 1
            % copy of ModelPeopleScreened
            MPS = [Data.ModelPeopleScreened_FY; Data.ModelPeopleScreened_FW; Data.ModelPeopleScreened_FP; ...
                   Data.ModelPeopleScreened_MY; Data.ModelPeopleScreened_MW; Data.ModelPeopleScreened_MP;];
            
            % identify the highest risk group
            [~,argmax] = max([Paras.rFY, Paras.rFW, Paras.rFP, Paras.rMY, Paras.rMW, Paras.rMP]);
            
            for idx = 1:6
                % take (Strategy.TestFactor)% of screening from other cohorts and allocate to
                % highest risk group
                % set default to 20% (see Stratfiles)
                if idx ~= argmax
                    MPS(argmax,:) = MPS(argmax,:) + 0.1.*MPS(idx,:);
                    MPS(idx,:) = MPS(idx,:) - 0.1.*MPS(idx,:);
                end
            end

            % reconstruct ModelPeopleScreened_GC
            Data.ModelPeopleScreened_FY = MPS(1,:);
            Data.ModelPeopleScreened_FW = MPS(2,:);
            Data.ModelPeopleScreened_FP = MPS(3,:);
            Data.ModelPeopleScreened_MY = MPS(4,:);
            Data.ModelPeopleScreened_MW = MPS(5,:);
            Data.ModelPeopleScreened_MP = MPS(6,:);
        end

        % Run projected dynamics
        [Classes1, Aggregate1] = ODEHATmodel(meff, ICs, Data, Paras, ProjStrat);
        Aggregate = [Aggregate0; Aggregate1];
        Classes = [Classes0; Classes1];
        
        % Projected ODE dynamics
%         Active1(:,NumYear0+1:end,s) = Aggregate1.ActiveM1' .* MtoAbsScaling1;
%         Active2(:,NumYear0+1:end,s) = Aggregate1.ActiveM2' .* MtoAbsScaling1;
        
        Active1(:,NumYear0+1:end,s) = (Aggregate1.ActiveM1 + Paras.S1givenFP * Aggregate1.ActiveMFP)' .* MtoAbsScaling1;
        Active2(:,NumYear0+1:end,s) = (Aggregate1.ActiveM2 + (1 - Paras.S1givenFP) * Aggregate1.ActiveMFP)' .* MtoAbsScaling1;
        Active1FP(:,NumYear0+1:end,s) = Paras.S1givenFP * Aggregate1.ActiveMFP' .* MtoAbsScaling1;
        Active2FP(:,NumYear0+1:end,s) = (1 - Paras.S1givenFP) * Aggregate1.ActiveMFP' .* MtoAbsScaling1;
        
        PassiveAgeY(:,NumYear0+1:end,s) = Aggregate1.PassiveMAgeY' .* MtoAbsScaling1;
        PassiveAgeW(:,NumYear0+1:end,s) = Aggregate1.PassiveMAgeW' .* MtoAbsScaling1;
        PassiveAgeP(:,NumYear0+1:end,s) = Aggregate1.PassiveMAgeP' .* MtoAbsScaling1;
        PassiveGenderF(:,NumYear0+1:end,s) = Aggregate1.PassiveMGenderF' .* MtoAbsScaling1;
        PassiveGenderM(:,NumYear0+1:end,s) = Aggregate1.PassiveMGenderM' .* MtoAbsScaling1;

        ActiveAgeY(:,NumYear0+1:end,s) = Aggregate1.ActiveMAgeY' .* MtoAbsScaling1;
        ActiveAgeW(:,NumYear0+1:end,s) = Aggregate1.ActiveMAgeW' .* MtoAbsScaling1;
        ActiveAgeP(:,NumYear0+1:end,s) = Aggregate1.ActiveMAgeP' .* MtoAbsScaling1;
        ActiveGenderF(:,NumYear0+1:end,s) = Aggregate1.ActiveMGenderF' .* MtoAbsScaling1;
        ActiveGenderM(:,NumYear0+1:end,s) = Aggregate1.ActiveMGenderM' .* MtoAbsScaling1;
        

        Passive1(:,NumYear0+1:end,s) = Aggregate1.PassiveM1' .* MtoAbsScaling1;
        Passive2(:,NumYear0+1:end,s) = Aggregate1.PassiveM2' .* MtoAbsScaling1;
        
        Deaths(:,NumYear0+1:end,s) = Aggregate1.DeathsM' .* MtoAbsScaling1;

        
        PersonYrs1(:,NumYear0+1:end,s) = Aggregate1.PersonYrsM1' .* MtoAbsScaling1;
        PersonYrs2(:,NumYear0+1:end,s) = Aggregate1.PersonYrsM2' .* MtoAbsScaling1;

        PersonYrsAgeY(:,NumYear0+1:end,s) = Aggregate1.PersonYrsMAgeY' .* MtoAbsScaling1;
        PersonYrsAgeW(:,NumYear0+1:end,s) = Aggregate1.PersonYrsMAgeW' .* MtoAbsScaling1;
        PersonYrsAgeP(:,NumYear0+1:end,s) = Aggregate1.PersonYrsMAgeP' .* MtoAbsScaling1;
        PersonYrsGenderM(:,NumYear0+1:end,s) = Aggregate1.PersonYrsMGenderM' .* MtoAbsScaling1;
        PersonYrsGenderF(:,NumYear0+1:end,s) = Aggregate1.PersonYrsMGenderF' .* MtoAbsScaling1;


        
        NewInf(:,NumYear0+1:end,s) = Aggregate1.NewInfM' .* MtoAbsScaling1;
        PerfectSpec(:,NumYear0+1:end,s) = Aggregate1.PerfectSpec';
        Transmission(:,NumYear0+1:end,s) = Aggregate1.Transmission';
        RDT1(:,NumYear0+1:end,s) = Aggregate1.RDTM1' .* MtoAbsScaling1;
        RDT2(:,NumYear0+1:end,s) = Aggregate1.RDTM2' .* MtoAbsScaling1;
        RDTFP(:,NumYear0+1:end,s) = Aggregate1.RDTMFP' .* MtoAbsScaling1;

        RDTGenderM(:,NumYear0+1:end,s) = Aggregate1.RDTMGenderM' .* MtoAbsScaling1;
        RDTGenderF(:,NumYear0+1:end,s) = Aggregate1.RDTMGenderF' .* MtoAbsScaling1;
        RDTAgeY(:,NumYear0+1:end,s) = Aggregate1.RDTMAgeY' .* MtoAbsScaling1;
        RDTAgeW(:,NumYear0+1:end,s) = Aggregate1.RDTMAgeW' .* MtoAbsScaling1;
        RDTAgeP(:,NumYear0+1:end,s) = Aggregate1.RDTMAgeP' .* MtoAbsScaling1;
       
        
        % Sampling from projected dynamics
        
        %Amended to include switch in overdispersion once all FP and
        %seperate FP output from model
%         Disp_act = Paras.disp_act * (Aggregate1.ActiveM1' + Aggregate1.ActiveM2') ./ (Aggregate1.ActiveM1' + Aggregate1.ActiveM2' + Aggregate1.ActiveMFP');
%         ActiveS = CCbetbinornd(repmat(Data.PeopleScreened, samples, 1), repmat((Aggregate1.ActiveM1' + Aggregate1.ActiveM2' + Aggregate1.ActiveMFP') ./ ScaledPeopleScreened1, samples, 1), repmat(Disp_act, samples, 1));
%         ActiveS1 = Cbinornd(ActiveS, repmat((Aggregate1.ActiveM1' + Paras.S1givenFP * Aggregate1.ActiveMFP') ./ (Aggregate1.ActiveM1' + Aggregate1.ActiveM2'+ Aggregate1.ActiveMFP'), samples, 1));
%         PassiveS = CCbetbinornd(repmat(Pop1, samples, 1), repmat((Aggregate1.PassiveM1' + Aggregate1.PassiveM2') / Data.N_H, samples, 1), repmat(Paras.disp_pass, samples, NumYear - NumYear0));
%         PassiveS1 = Cbinornd(PassiveS, repmat(Aggregate1.PassiveM1' ./ (Aggregate1.PassiveM1' + Aggregate1.PassiveM2'), samples, 1));
%         DeathsS = Cbinornd(repmat(Pop1, samples, 1), repmat(Aggregate1.DeathsM' / Data.N_H, samples, 1));
        Aggregate1.ActiveMS = zeros(size(Aggregate1.ActiveM1));

        Disp_act = Paras.disp_act * (Aggregate1.ActiveMS' + Aggregate1.ActiveM1' + Aggregate1.ActiveM2') ./ (Aggregate1.ActiveMS' + Aggregate1.ActiveM1' + Aggregate1.ActiveM2' + Aggregate1.ActiveMFP'); % M9
        ActiveS = Cbetabinornd(repmat(Data.PeopleScreened, samples, 1), repmat((Aggregate1.ActiveMS' + Aggregate1.ActiveM1' + Aggregate1.ActiveM2' + Aggregate1.ActiveMFP') ./ ScaledPeopleScreened1, samples, 1), repmat(Disp_act, samples, 1)); % M9
        ActiveS1 = Cbinornd(ActiveS, repmat((Aggregate1.ActiveM1' + Paras.S1givenFP * Aggregate1.ActiveMFP') ./ (Aggregate1.ActiveMS' + Aggregate1.ActiveM1' + Aggregate1.ActiveM2'+ Aggregate1.ActiveMFP'), samples, 1)); % M9
        ActiveS1FP = Cbinornd(ActiveS1, repmat(Paras.S1givenFP * Aggregate1.ActiveMFP' ./ (Aggregate1.ActiveM1' + Paras.S1givenFP * Aggregate1.ActiveMFP'), samples, 1));
        ActiveS2 = Cbinornd(ActiveS - ActiveS1, repmat(Aggregate1.ActiveM2' + (1-Paras.S1givenFP) * Aggregate1.ActiveMFP',samples,1) ./ (Aggregate1.ActiveMS' + Aggregate1.ActiveM2' + (1-Paras.S1givenFP) * Aggregate1.ActiveMFP')); % M9
        ActiveS2FP = Cbinornd(ActiveS2, repmat((1 - Paras.S1givenFP) * Aggregate1.ActiveMFP' ./ (Aggregate1.ActiveM2' + (1 - Paras.S1givenFP) * Aggregate1.ActiveMFP'), samples, 1));
        PassiveS = Cbetabinornd(repmat(Pop1, samples, 1), repmat((Aggregate1.PassiveM1' + Aggregate1.PassiveM2') / Data.N_H, samples, 1), repmat(Paras.disp_pass, samples, NumYear - NumYear0));
        PassiveS1 = Cbinornd(PassiveS, repmat(Aggregate1.PassiveM1' ./ (Aggregate1.PassiveM1' + Aggregate1.PassiveM2'), samples, 1));
        DeathsS = Cbinornd(repmat(Pop1, samples, 1), repmat(Aggregate1.DeathsM' / Data.N_H, samples, 1));
        
        ActiveSAgeY = Cbinornd(ActiveS, repmat(Aggregate1.ActiveMAgeY' ./ (Aggregate1.ActiveMAgeY' + Aggregate1.ActiveMAgeW' +  Aggregate1.ActiveMAgeP'), samples, 1));
        ActiveSAgeW = Cbinornd(ActiveS, repmat(Aggregate1.ActiveMAgeW' ./ (Aggregate1.ActiveMAgeY' + Aggregate1.ActiveMAgeW' +  Aggregate1.ActiveMAgeP'), samples, 1));
        ActiveSAgeP = ActiveS - ActiveSAgeY - ActiveSAgeW;
        PassiveSAgeY = Cbinornd(PassiveS, repmat(Aggregate1.PassiveMAgeY' ./ (Aggregate1.PassiveMAgeY' + Aggregate1.PassiveMAgeW' +  Aggregate1.PassiveMAgeP'), samples, 1));
        PassiveSAgeW = Cbinornd(PassiveS, repmat(Aggregate1.PassiveMAgeW' ./ (Aggregate1.PassiveMAgeY' + Aggregate1.PassiveMAgeW' +  Aggregate1.PassiveMAgeP'), samples, 1));
        PassiveSAgeP = PassiveS - PassiveSAgeY - PassiveSAgeW;
    
        ActiveSGenderF = Cbinornd(ActiveS, repmat(Aggregate1.ActiveMGenderF' ./ (Aggregate1.ActiveMGenderF' + Aggregate1.ActiveMGenderM'), samples, 1));
        ActiveSGenderM = ActiveS - ActiveSGenderF;

        PassiveSGenderF = Cbinornd(PassiveS, repmat(Aggregate1.PassiveMGenderF' ./ (Aggregate1.PassiveMGenderF' + Aggregate1.PassiveMGenderM'), samples, 1));
        PassiveSGenderM = PassiveS - PassiveSGenderF;
       


        [RDTS1, RDTS2, RDTSFP, RDTSAgeP, RDTSAgeY, RDTSAgeW, RDTSGenderM, RDTSGenderF] = deal(zeros(size(ActiveS))); % M9

        if ProjStrat.TTyear ~= 0
            Y = find(Data.Years >= ProjStrat.TTyear, 1);
            RDTS1(:, Y:end) = ActiveS1(:, Y:end) + Cbinornd(repmat(Data.PeopleScreened(Y:end), samples, 1), repmat((Aggregate1.RDTM1(Y:end)' - Aggregate1.ActiveM1(Y:end)') ./ ScaledPeopleScreened1(Y:end), samples, 1));
            RDTS2(:, Y:end) = ActiveS2(:, Y:end) + Cbinornd(repmat(Data.PeopleScreened(Y:end), samples, 1), repmat((Aggregate1.RDTM2(Y:end)' - Aggregate1.ActiveM2(Y:end)') ./ ScaledPeopleScreened1(Y:end), samples, 1));
            RDTSFP(:, Y:end) = Cbinornd(repmat(Data.PeopleScreened(Y:end), samples, 1), repmat(1 - ProjStrat.SpecRDT, samples, length(Data.Years) - Y + 1));
        
            RDTSAgeY(:, Y:end) = ActiveSAgeY(:, Y:end) + Cbinornd(repmat(Data.PeopleScreened(Y:end), samples, 1), repmat((Aggregate1.RDTMAgeY(Y:end)' - Aggregate1.ActiveMAgeY(Y:end)') ./ ScaledPeopleScreened1(Y:end), samples, 1));
            RDTSAgeW(:, Y:end) = ActiveSAgeW(:, Y:end) + Cbinornd(repmat(Data.PeopleScreened(Y:end), samples, 1), repmat((Aggregate1.RDTMAgeW(Y:end)' - Aggregate1.ActiveMAgeW(Y:end)') ./ ScaledPeopleScreened1(Y:end), samples, 1));
            RDTSAgeP(:, Y:end) = ActiveSAgeP(:, Y:end) + Cbinornd(repmat(Data.PeopleScreened(Y:end), samples, 1), repmat((Aggregate1.RDTMAgeP(Y:end)' - Aggregate1.ActiveMAgeP(Y:end)') ./ ScaledPeopleScreened1(Y:end), samples, 1));
            RDTSGenderM(:, Y:end) = ActiveSGenderM(:, Y:end) + Cbinornd(repmat(Data.PeopleScreened(Y:end), samples, 1), repmat((Aggregate1.RDTMGenderM(Y:end)' - Aggregate1.ActiveMGenderM(Y:end)') ./ ScaledPeopleScreened1(Y:end), samples, 1));
            RDTSGenderF(:, Y:end) = ActiveSGenderF(:, Y:end) + Cbinornd(repmat(Data.PeopleScreened(Y:end), samples, 1), repmat((Aggregate1.RDTMGenderF(Y:end)' - Aggregate1.ActiveMGenderF(Y:end)') ./ ScaledPeopleScreened1(Y:end), samples, 1));
            
        end
        
        SampledActive1(:, NumYear0+1:end, s) = ActiveS1;
        SampledActive2(:, NumYear0+1:end, s) = ActiveS2; % M9
        SampledActive1FP(:, NumYear0+1:end, s) = ActiveS1FP;
        SampledActive2FP(:, NumYear0+1:end, s) = ActiveS2FP;

        SampledActiveAgeY(:, NumYear0+1:end, s) = ActiveSAgeY;
        SampledActiveAgeW(:, NumYear0+1:end, s) = ActiveSAgeW;
        SampledActiveAgeP(:, NumYear0+1:end, s) = ActiveSAgeP;
        SampledActiveGenderM(:, NumYear0+1:end, s) = ActiveSGenderM;
        SampledActiveGenderF(:, NumYear0+1:end, s) = ActiveSGenderF;

        SampledPassiveAgeY(:, NumYear0+1:end, s) = PassiveSAgeY;
        SampledPassiveAgeW(:, NumYear0+1:end, s) = PassiveSAgeW;
        SampledPassiveAgeP(:, NumYear0+1:end, s) = PassiveSAgeP;
        SampledPassiveGenderM(:, NumYear0+1:end, s) = PassiveSGenderM;
        SampledPassiveGenderF(:, NumYear0+1:end, s) = PassiveSGenderF;

        SampledPassive1(:, NumYear0+1:end, s) = PassiveS1;
        SampledPassive2(:, NumYear0+1:end, s) = PassiveS - PassiveS1;
        SampledDeaths(:, NumYear0+1:end, s) = DeathsS;

        SampledRDT1(:, NumYear0+1:end, s) = RDTS1;
        SampledRDT2(:, NumYear0+1:end, s) = RDTS2;
        SampledRDTFP(:, NumYear0+1:end, s) = RDTSFP;
        
        SampledRDTAgeY(:, NumYear0+1:end, s) = RDTSAgeY;
        SampledRDTAgeW(:, NumYear0+1:end, s) = RDTSAgeW;
        SampledRDTAgeP(:, NumYear0+1:end, s) = RDTSAgeP;
        SampledRDTGenderM(:, NumYear0+1:end, s) = RDTSGenderM;
        SampledRDTGenderF(:, NumYear0+1:end, s) = RDTSGenderF;
        
        % Elimination years
        YEPHP(s) = max([find(sum(Aggregate{:, {'ActiveM1','ActiveM2', 'ActiveMFP', 'PassiveM1','PassiveM2'}}, 2) * 10000 >= Data.N_H, 1, 'last') 0]) + Years(1); % M9
        YEOT(s) = max([find(Aggregate.NewInfM' .* MtoAbsScaling >= Paras.EOTthreshold, 1, 'last') 0]) + Years(1); % no transmission threshold = 1
        SampledCases = SampledActive1(:,:,s) + SampledActive2(:,:,s) + SampledPassive1(:,:,s) + SampledPassive2(:,:,s); % M9
        SampledYEPHP(:, s) = table2array(rowfun(@(SampledCases)(max([find(SampledCases * 10000 >= Pop, 1, 'last') 0])), table(SampledCases))) + double(Years(1)) * ones(samples, 1);
        
        % Post-processing of NewInf (NewInf = 0 after EoT rather than after EoI)
        %NewInf(NewInf(:, :, s) < Paras.EOTthreshold) = 0;
        
        switch ProjStrat.CaseDef
            case 'Parasitol'
                Cases = SampledCases - SampledActive1FP(:, :, s) - SampledActive2FP(:, :, s);
            case 'Serol'
                Cases = SampledCases;
        end
% %         AS = [AS Data.PeopleScreened];
        
        % Starting info for reactive interactions
        for r = 1 : NumReact
            R = ['React', num2str(r)];
            % Reactive Parameters
            React = table2struct(ReactiveParameters(R,:));
            
% %             Y0 = find(Years == React.Ryear);
% %             if sum(AS(Y0 - 4 : Y0 - 1)) == 0 % Historical foci (No AS for 4 years before Ryear & one-off screening ~= 0)
% %                 % Check if VC will cease as well
% %                 RyearVC = max((ProjStrat.NewVCyear + (1 - VCstarted) * Paras.MinVC) * ((VCstarted + ProjStrat.NewTargetFreq) ~= 0), React.Ryear);
% %                 NumZerosVC = React.StoppingVC + 100 * (ProjStrat.NewTargetDie == 0) * (RVC == 0);
% %                 if RyearVC == React.Ryear
% %                     y_VC = max([find(Years == RyearVC) - NumZerosVC 1]);
% %                     Y_VC = table2array(rowfun(@(Cases)(min([strfind(Cases(:, y_VC:end), zeros(1, NumZerosVC)) NumYear])), table(Cases))) + (y_VC - 1) + NumZerosVC;
% %                 
% %                     VC = (Y_VC > Y0) * (React.StoppingVC == NumZerosVC);
% %                 else
% %                     VC = ones(samples, 1) * (React.StoppingVC == NumZerosVC);
% %                 end
% %                 YearIDs = Y0;
% %                 AS = 0;
% %                 RS = 1;
% %                 OS = 1;
% % 
% %                 for i = 1 : samples
% %                     ReactInfo = [ReactInfo; [s r i React.Ryear AS RS OS VC(i) Aggregate.PerfectSpec(YearIDs) Aggregate.Transmission(YearIDs) meff Classes{Classes.Time == Years(YearIDs), 2 : end}(1, :)]];
% %                 end 
% %             else
                switch React.ReactiveASnum
                    case 'mean'
                        RS = Data.MeanPeopleScreened ~= 0;
                    case 'max'
                        RS = Data.MaxPeopleScreened ~= 0;
                    case 'off'
                        RS = 0;
                    otherwise
                        if ischar(React.ReactiveASnum) % percentage 'x%'
                            RS = str2num(React.ReactiveASnum(1 : end - 1)) ~= 0;
                        else
                            RS = React.ReactiveASnum ~= 0;
                        end
                end

                switch React.OneoffASnum
                    case 'mean'
                        OS = Data.MeanPeopleScreened ~= 0;
                    case 'max'
                        OS = Data.MaxPeopleScreened ~= 0;
                    case 'off'
                        OS = 0;
                    otherwise
                        if ischar(React.OneoffASnum) % percentage 'x%'
                            OS = str2num(React.OneoffASnum(1 : end - 1)) ~= 0;
                        else
                            OS = React.OneoffASnum ~= 0;
                        end
                end

                if isequal(React.ReactiveTargetFreq, 'AsStrat')
                    ReactTargetFreq = ProjStrat.NewTargetFreq;
                elseif isequal(React.ReactiveTargetFreq, 'scaleback')
                    ReactTargetFreq = Paras.TargetFreq_scaleback;
                else
                    ReactTargetFreq = React.ReactiveTargetFreq;
                end
                TargetFreq = VCstarted * Paras.TargetFreq + (1 - VCstarted) * ProjStrat.NewTargetFreq;

%                 if isequal(React.ReactiveTargetDie, 'AsStrat')
%                     React.ReactiveTargetDie = ProjStrat.NewTargetDie;
%                 end
                switch React.ReactiveTargetCoverage
                    case 'AsStrat'
                        ReactTargetDie = ProjStrat.NewTargetDie;
                    case 'full'
                        ReactTargetDie = round(GetTargetDie(ProjStrat.NewVC, TargetFreq, ProjStrat.NewTrapCycle, Paras), 4);
                    case 'realistic'
                        ReactTargetDie = round(GetTargetDie(ProjStrat.NewVC * Paras.RealisticVCwgt, TargetFreq, ProjStrat.NewTrapCycle, Paras), 4);
                    case 'scaleback'
                        ReactTargetDie = round(GetTargetDie(ProjStrat.NewVC * Paras.VCwgt_scaleback, TargetFreq, ProjStrat.NewTrapCycle, Paras), 4);
                    otherwise
                        ReactTargetDie = round(GetTargetDie(ProjStrat.NewVC * React.ReactiveTargetCoverage, TargetFreq, ProjStrat.NewTrapCycle, Paras), 4);
                end % ScaleBackVC

                RVC = ReactTargetFreq * ReactTargetDie;

                RyearAS = React.Ryear;
% %             if Data.PeopleScreened(1) ~= 0 && ProjStrat.NewASyear < ProjStrat.SIMyear && ProjStrat.NewASyear > React.Ryear
% %                 RyearAS = ProjStrat.NewASyear; % + MinAS if considered
% %             end
%                 RyearVC = max((ProjStrat.NewVCyear + (1 - VCstarted) * Paras.MinVC) * ((VCstarted + ProjStrat.NewTargetFreq) ~= 0), React.Ryear);
                RyearVC = max(ceil(VCstarted * Paras.VCstart + (1 - VCstarted) * ProjStrat.NewVCyear * (ProjStrat.NewTargetDie ~= 0)) + Paras.MinVC, React.Ryear); %ScalingUpVC
            
                NumZerosAS = React.StoppingAS + 100 * (sum(Data.PeopleScreened) == 0) * (RS == 0);
                NumZerosVC = React.StoppingVC + 100 * (ProjStrat.NewTargetDie == 0) * (RVC == 0);
            
                y_AS = max([find(Years == RyearAS) - NumZerosAS 1]);
                y_VC = max([find(Years == RyearVC) - NumZerosVC 1]);
            
                Y_AS = table2array(rowfun(@(Cases)(min([strfind(Cases(:, y_AS:end), zeros(1, NumZerosAS)) NumYear])), table(Cases))) + (y_AS - 1) + NumZerosAS;
                Y_VC = table2array(rowfun(@(Cases)(min([strfind(Cases(:, y_VC:end), zeros(1, NumZerosVC)) NumYear])), table(Cases))) + (y_VC - 1) + NumZerosVC;
            
                Y = min(Y_AS, Y_VC);
            
            
% % %             NumZeros = min([React.StoppingAS + 100 * (sum(Data.PeopleScreened) == 0 | strcmp(ProjStrat.NewASstrat, 'stop')) React.StoppingVC + 100 * (ProjStrat.NewTargetDie == 0)]);
% % %             %y = max([find(Years == React.Ryear) - NumZeros 1]);
% % %             y = max([find(Years == max([React.Ryear ProjStrat.NewVCyear + (Paras.VCstart == 0 & ProjStrat.NewTargetDie ~= 0) * Paras.MinVC])) - NumZeros 1]);
% % % 
% % %             Y = table2array(rowfun(@(SampledCases)(min([strfind(SampledCases(:, y:end), zeros(1, NumZeros)) NumYear])), table(SampledCases))) + (y - 1) + NumZeros;
                SampleIDs = find(Y <= NumYear);
                YearIDs = Y(SampleIDs);
                AS = (YearIDs ~= Y_AS(SampleIDs)) * (React.StoppingAS == NumZerosAS); % status of AS; 0:stop (including both AS and RS are zero) 1:conti
                VC = (YearIDs ~= Y_VC(SampleIDs)) * (React.StoppingVC == NumZerosVC) * ((VCstarted + ProjStrat.NewTargetDie) ~= 0); % status of VC; 0:stop (including 1. both NewVC and RVC are zero 2. VCstart = 0 and NewVC = 0) 1:conti
            
                for i = 1 : length(SampleIDs)
                    if AS(i) == 0 && OS == 0
                        rs = 1;
                    else
                        rs = 0;
                    end
                    os = 0;
                    ReactInfo = [ReactInfo; [s+RunInfo.StratMin-1 r SampleIDs(i) double(YearIDs(i)+Years(1)-1) AS(i) VC(i) rs os Aggregate.PerfectSpec(YearIDs(i)) Aggregate.Transmission(YearIDs(i)) meff Classes{Classes.Time == Years(YearIDs(i)), 2 : end}(1, :)]];
                end
            %%end
        end
        
        %%% If we haven't reached EOT, simulate this many extra years (step size)
        extra_sim_years = 50;
        maxYEOT = 2100; % don't go beyond here (prevent an infinte while loop if, eg, gHAT increasing)
        dYearsBackup = Data.Years;
        YearsBackup = Years;
        
        % Update ProjStrat
        ProjStrat_extra = ProjStrat;
        ProjStrat_extra.NewASnum = NewASnum{end};
        ProjStrat_extra.NewASstrat = NewASstrat{end};
        ProjStrat_extra.NewASyear = NewASyear(end);
        
        while YEOT(s) == length(Years)+Years(1) && Years(end) < maxYEOT
            Data.Years = Data.Years(end) + 1 : Data.Years(end) + extra_sim_years; % all strategis must have the same SIMyear
            Data.Years = Data.Years(find(Data.Years <= maxYEOT));
            Years = [Years, Data.Years];
            MtoAbsScalingA = Data.PopGrowth .^ double(Years - Data.PopSizeYear);
            MtoAbsScalingB = Data.PopGrowth .^ double(Data.Years - Data.PopSizeYear);
            
            switch NewASnum{end}
                case 'mean'
                    Data.PeopleScreened = Data.MeanPeopleScreened * ones(1, length(Data.Years));
                case 'max'
                    Data.PeopleScreened = Data.MaxPeopleScreened * ones(1, length(Data.Years));
                case 'intensified'
                    Data.PeopleScreened = Data.IntensifiedPeopleScreened * ones(1, length(Data.Years));
                case 'off'
                    Data.PeopleScreened = zeros(1, length(Data.Years));
                case 'stop'
                    Data.PeopleScreened = zeros(1, length(Data.Years));
                otherwise
                    if contains(NewASnum{end}, '%') % percentage 'x%'
                        pct = NewASnum{end};
                        Data.PeopleScreened = round(0.01 * str2num(pct(1 : end - 1)) * Data.N_H * MtoAbsScaling0(end - 2 - Data.DeltaGapYear)) * ones(1, length(Data.Years));
                    else
                        Data.PeopleScreened = str2num(NewASnum{end}) * ones(1, length(Data.Years));
                    end % Isangi
            end
            ScaledPeopleScreened1 = Data.PeopleScreened ./ MtoAbsScalingB;
            ModelScreeningFreq = [];
            ModelScreeningTime = [];
            ModelPeopleScreened = [];
            ExpectedScreeningTimes = ceil(ScaledPeopleScreened1 / Data.N_H / ScreeningCapacity(end));
            ExpectedScreeningTimes(ExpectedScreeningTimes == 0) = 1;
            for y = 1 : length(Data.Years)
                freq = 365 / ExpectedScreeningTimes(y);
                time = 1 / ExpectedScreeningTimes(y);
                number = ScaledPeopleScreened1(y) / ExpectedScreeningTimes(y);
            
                ModelScreeningFreq = [ModelScreeningFreq freq * ones(1, ExpectedScreeningTimes(y))];
                ModelScreeningTime = [ModelScreeningTime Data.Years(y) + (0 : ExpectedScreeningTimes(y) - 1) * time];
                ModelPeopleScreened = [ModelPeopleScreened number * ones(1, ExpectedScreeningTimes(y))];
            end
            Data.ModelScreeningFreq = ModelScreeningFreq;
            Data.ModelScreeningTime = ModelScreeningTime;
            Data.ModelPeopleScreened = ModelPeopleScreened;
            
            % NATHAN: see lines 325-372 for relevant comments
            Data.ModelPeopleScreened_FY = ModelPeopleScreened * (Data.N_FY/Data.N_H);
            Data.ModelPeopleScreened_FW = ModelPeopleScreened * (Data.N_FW/Data.N_H);
            Data.ModelPeopleScreened_FP = ModelPeopleScreened * (Data.N_FP/Data.N_H);
            Data.ModelPeopleScreened_MY = ModelPeopleScreened * (Data.N_MY/Data.N_H);
            Data.ModelPeopleScreened_MW = ModelPeopleScreened * (Data.N_MW/Data.N_H);
            Data.ModelPeopleScreened_MP = ModelPeopleScreened * (Data.N_MP/Data.N_H);
        
            if ProjStrat.TargetTest == 1
                MPS = [Data.ModelPeopleScreened_FY; Data.ModelPeopleScreened_FW; Data.ModelPeopleScreened_FP; ...
                       Data.ModelPeopleScreened_MY; Data.ModelPeopleScreened_MW; Data.ModelPeopleScreened_MP;];

                [~,argmax] = max([Paras.rFY, Paras.rFW, Paras.rFP, Paras.rMY, Paras.rMW, Paras.rMP]);
            
                for idx = 1:6
                    if idx ~= argmax
                        MPS(argmax,:) = MPS(argmax,:) + ProjStrat.TestFactor.*MPS(idx,:);
                        MPS(idx,:) = MPS(idx,:) - ProjStrat.TestFactor.*MPS(idx,:);
                    end
                end

                Data.ModelPeopleScreened_FY = MPS(1,:);
                Data.ModelPeopleScreened_FW = MPS(2,:);
                Data.ModelPeopleScreened_FP = MPS(3,:);
                Data.ModelPeopleScreened_MY = MPS(4,:);
                Data.ModelPeopleScreened_MW = MPS(5,:);
                Data.ModelPeopleScreened_MP = MPS(6,:);
            end
            
            % Update PPSstart, PPSend, NoPSstart and NoPSend
            switch ProjStrat.NewPS 
                case 'partial'
                    Paras.PPSstart = NewASyear(end);
                    Paras.PPSend = Data.Years(end) + 1;
                    Paras.NoPSstart = 0;
                    Paras.NoPSend = 0;
                case 'no'
                    Paras.PPSstart = 0;
                    Paras.PPSend = 0;
                    Paras.NoPSstart = NewASyear(end);
                    Paras.NoPSend = Data.Years(end) + 1;
            end % Isangi

            % Get ICs
            pop = Classes1(end, :);
            yICs = {[pop.S_H_FY, pop.S_H_FW, pop.S_H_FP, pop.S_H_MY,  pop.S_H_MW,  pop.S_H_MP], [pop.E_H_FY, pop.E_H_FW, pop.E_H_FP, pop.E_H_MY,  pop.E_H_MW,  pop.E_H_MP],...
           [pop.I1_H_FY, pop.I1_H_FW, pop.I1_H_FP, pop.I1_H_MY,  pop.I1_H_MW,  pop.I1_H_MP], [pop.I2_H_FY, pop.I2_H_FW, pop.I2_H_FP, pop.I2_H_MY,  pop.I2_H_MW,  pop.I2_H_MP],[pop.R_H_FY, pop.R_H_FW, pop.R_H_FP, pop.R_H_MY,  pop.R_H_MW,  pop.R_H_MP],...
            pop.P_V, pop.S_V, pop.G_V, pop.E1_V, pop.E2_V, pop.E3_V, pop.I_V};

            % Update alpha
            if Aggregate1.NewInfM(end) * MtoAbsScaling1(end) < Paras.EOTthreshold && ...
            (pop.I1_H_FY+ pop.I1_H_FW+ pop.I1_H_FP+ pop.I1_H_MY+  pop.I1_H_MW+  pop.I1_H_MP+ pop.I2_H_FY+ pop.I2_H_FW+ pop.I2_H_FP+ pop.I2_H_MY+  pop.I2_H_MW+  pop.I2_H_MP) * MtoAbsScaling1(end) < Paras.EOTthreshold % M9
          Paras.alpha = 0;
            else
                Paras.alpha = alpha;
            end

            
            [Classes1, Aggregate1] = ODEHATmodel(meff, yICs, Data, Paras, ProjStrat_extra);
            Aggregate = [Aggregate; Aggregate1];            
            YEOT(s) = max([find(Aggregate.NewInfM' .* MtoAbsScalingA >= Paras.EOTthreshold, 1, 'last') 0]) + Years(1); % no transmission threshold = 1
        end
        Data.Years = dYearsBackup;
        Years = YearsBackup;
    end
    
    if size(ReactInfo, 1) == 0
        ReactInfo = zeros(1, 48); % M9
    end
   
    % output here
    Outputs = struct('Years', Years, 'ReactInfo', ReactInfo, 'Transmission', Transmission, 'PerfectSpec', PerfectSpec, ...
                      'Active1', Active1, 'Active2', Active2, 'Active1FP', Active1FP, 'Active2FP', Active2FP, 'Passive1', Passive1, 'Passive2', Passive2, ...
                      'ActiveAgeY', ActiveAgeY,'ActiveAgeW', ActiveAgeW,'ActiveAgeP', ActiveAgeP,'ActiveGenderM', ActiveGenderM,'ActiveGenderF', ActiveGenderF,...
                      'PassiveAgeY', PassiveAgeY,'PassiveAgeW', PassiveAgeW,'PassiveAgeP', PassiveAgeP,'PassiveGenderM', PassiveGenderM,'PassiveGenderF', PassiveGenderF,... 
                      'PersonYrsAgeY', PersonYrsAgeY,'PersonYrsAgeW', PersonYrsAgeW,'PersonYrsAgeP', PersonYrsAgeP,'PersonYrsGenderM', PersonYrsGenderM,'PersonYrsGenderF', PersonYrsGenderF,...
                     'Deaths', Deaths,  'PersonYrs1', PersonYrs1, 'PersonYrs2', PersonYrs2, 'NewInf', NewInf, ...
                     'YEPHP', YEPHP, 'YEOT', YEOT, 'SampledYEPHP', SampledYEPHP,...
                     ...
                      'SampledActive1', SampledActive1, 'SampledActive2', SampledActive2,'SampledActive1FP', SampledActive1FP, 'SampledActive2FP', SampledActive2FP, ...                     
                      'SampledActiveAgeY', SampledActiveAgeY,'SampledActiveAgeW', SampledActiveAgeW,'SampledActiveAgeP', SampledActiveAgeP,'SampledActiveGenderM', SampledActiveGenderM,'SampledActiveGenderF', SampledActiveGenderF,...
                     'SampledPassiveAgeY', SampledPassiveAgeY,'SampledPassiveAgeW', SampledPassiveAgeW,'SampledPassiveAgeP', SampledPassiveAgeP,'SampledPassiveGenderM', SampledPassiveGenderM,'SampledPassiveGenderF', SampledPassiveGenderF,...         
                     'SampledPassive1', SampledPassive1, 'SampledPassive2', SampledPassive2, ...
                     'SampledDeaths', SampledDeaths, 'ProjectionICs', Classes0{end, 2:end}, ...
                     'RDTAgeY', RDTAgeY,'RDTAgeW', RDTAgeW,'RDTAgeP', RDTAgeP,'RDTGenderM', RDTGenderM,'RDTGenderF', RDTGenderF,...
                     'SampledRDTAgeY', SampledRDTAgeY,'SampledRDTAgeW', SampledRDTAgeW,'SampledRDTAgeP', SampledRDTAgeP,'SampledRDTGenderM', SampledRDTGenderM,'SampledRDTGenderF', SampledRDTGenderF,...
                      'RDT1', RDT1, 'RDT2', RDT2, 'RDTFP', RDTFP, 'SampledRDT1', SampledRDT1, 'SampledRDT2', SampledRDT2, 'SampledRDTFP', SampledRDTFP);

    
    

