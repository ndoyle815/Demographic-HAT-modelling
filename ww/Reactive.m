
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                               %
%   This code runs reactive forcasting of a given strategy in the Warwick HAT model                             %
%                                                                                                               %
%   Inputs:                                                                                                     %
%       Data - structure containing location-specific historical data                                           %
%       Paras - structure containing location-specific parameters (fixed, fitted and intervention parameters)   %
%       Strategy - table containing parameters of all strategies                                                %
%       s - number denoting strategy ID                                                                         %
%       ReactiveParameters - table containing parameters of all reactive strategies                             %
%       ReactInfo - table containing key inputs (e.g. year and initial conditions) of reactive interventions    %
%       RunCFS - a cell array determining to run specific CounterFactual scenarios or not                       %
%                                                                                                               %
%   Outputs:                                                                                                    %
%       Outputs - structure containing yearly aggregated outputs, reactive information and ICs for forcasting   %
%                                                                                                               %
%   Functions required: ODEHATmodel & Cbetabinornd                                                              %
%                                                                                                               %
%                                                                                                               %
%                                                                                                               %
%   Notes from Ching-I: WHO classic algorithm (no 4th year screening) in RS isn't an option in this version     %
%                                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Outputs = Reactive(RunInfo, Data, Paras, CstrParas, Strategy, s, ReactiveParameters, ReactInfo, RunCFS)
    %%% All strategies
    load(RunInfo.FitFilePath('Posterior'));
    Fitted = load(RunInfo.ProjFilePath('Fitted'));
    Elim = load(RunInfo.ProjFilePathR('Elimination', 0));
%     alpha0 = Paras.alpha;
%    Spec0 = Paras.specificity;
    Y0 = Fitted.Years(1);
    
    % Counter-factual correction in Posterior
    if sum(strcmp(RunCFS, '2')) == 1
        if sum(ismember(Paras.FittedNames, {'eta_H_amp', 'gamma_H_amp'})) == 2
            Posterior.eta_H_amp(:) = 0;
            Posterior.gamma_H_amp(:) = 0;
        else
            Paras.eta_H_amp = 0;
            Paras.gamma_H_amp = 0;
            warning('eta_H_amp and gamma_H_amp are NOT fitted parameters')
        end
    end
%     Posterior.eta_H_amp = (sum(strcmp(RunCFS, '2')) == 0) * Posterior.eta_H_amp;
%     Posterior.gamma_H_amp = (sum(strcmp(RunCFS, '2')) == 0) * Posterior.gamma_H_amp;

    % Update the constrained parameter
    if size(CstrParas, 1) == 1
        Posterior.(CstrParas.Notation{:}) = 1 - sum(Posterior{:, ismember(Paras.FittedNames, {'k1','k2','k3','k4'})}, 2);
        Paras.FittedNames = [Paras.FittedNames CstrParas.Notation];
    end

    
    %%% Strategy sepcific info
    S = ['Strat', num2str(s)];
    ProjStrat = table2struct(Strategy(S,:));
    load(RunInfo.ProjFilePathSR('Projection', s, 0)) % Years is stored here
    Yrs = [Years Years(end) + 1];
    
    % AS parameters in the current strategy
    StratPeopleScreened = [];
    ASstrat = {};
    NewASnum = split(ProjStrat.NewASnum, '_');
    NewASstrat = split(ProjStrat.NewASstrat, '_');
    NewASyear = [ProjStrat.NewASyear Years(end)+1];
    for i = 1 : length(NewASnum)
        yrs = NewASyear(i+1) - NewASyear(i);
        switch NewASnum{i}
            case 'mean'
                PS = Data.MeanPeopleScreened * ones(1, yrs);
            case 'max'
                PS = Data.MaxPeopleScreened * ones(1, yrs);
            case 'intensified'
                PS = Data.IntensifiedPeopleScreened * ones(1, yrs);
            case 'off'
                PS = zeros(1, yrs);
            case 'stop'
                PS = zeros(1, yrs);
            otherwise 
                %PS = round(0.01 * str2num(ProjStrat.NewASnum(1 : end - 1)) * Data.N_H * Data.PopGrowth ^ (Fitted.Years(end - 2 - Data.DeltaGapYear) - Data.PopSizeYear)) * ones(1, yrs);
                if contains(NewASnum{i}, '%') % percentage 'x%'
                    pct = NewASnum{i};
                    PS = round(0.01 * str2num(pct(1 : end - 1)) * Data.N_H * Data.PopGrowth ^ (Fitted.Years(end - 2 - Data.DeltaGapYear) - Data.PopSizeYear)) * ones(1, yrs);
                else
                    PS = str2num(NewASnum{i}) * ones(1, yrs);
                end % Isangi
        end
        StratPeopleScreened = [StratPeopleScreened PS];
        
        ASstrat = [ASstrat repmat(NewASstrat(i), 1, yrs)];
    end


%     switch ProjStrat.NewASnum
%         case 'mean'
%             StratPeopleScreened = Data.MeanPeopleScreened;
%         case 'max'
%             StratPeopleScreened = Data.MaxPeopleScreened;
%         case 'intensified'
%             StratPeopleScreened = Data.IntensifiedPeopleScreened;
%         %case 'off'
%         %    StratPeopleScreened = 0;
%         otherwise % percentage 
%             StratPeopleScreened = round(0.01 * str2num(ProjStrat.NewASnum(1 : end - 1)) * Data.N_H * Paras.PopGrowth ^ (Data.Years(end - 2 - Data.DeltaGapYear) - Data.PopSizeYear));
%     end
% 
%     % Modify StratPeopleScreened when screening stops forever from NewASyear 
%     if strcmp(ProjStrat.NewASstrat, 'stop')
%         StratPeopleScreened = 0;
%     end
    
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
    
    % VC parameters in historical and current strategy
%     ProjStrat.NewVCyear = max([ProjStrat.NewVCyear ceil(Paras.VCstart + Paras.MinVC)]); % ScalingUpVC
    StratTargetFreq = ProjStrat.NewTargetFreq;
    HisTargetFreq = (Paras.VCstart ~= 0) * Paras.TargetFreq;
    StratTargetDie = ProjStrat.NewTargetDie;
    HisTargetDie = (Paras.VCstart ~= 0) * Paras.TargetDie;
    Y_NewVC = find(Years == ceil(ProjStrat.NewVCyear));
    VCstarted = Paras.VCstart ~= 0;
    
    %StratVCyear = (Paras.VCstart ~= 0) * (ProjStrat.NewTargetDie ~= 0) * Paras.VCstart + (Paras.VCstart == 0 | ProjStrat.NewTargetDie == 0) * ProjStrat.NewVCyear;
    %StratTargetFreq = (Paras.VCstart ~= 0) * (ProjStrat.NewTargetDie ~= 0) * Paras.TargetFreq + (Paras.VCstart == 0) * (ProjStrat.NewTargetDie > 0) * ProjStrat.NewTargetFreq;
    %ProjStrat.NewTargetDie = (Paras.VCstart ~= 0) * Paras.TargetDie + (Paras.VCstart == 0) * ProjStrat.NewTargetDie;
    
    scaling = Data.PopGrowth .^ double([Fitted.Years(1) : Years(end)] - Data.PopSizeYear);
    
    YEPHP = Elim.YEPHP(:, s);
    YEOT = Elim.YEOT(:, s);
    SampledYEPHP = Elim.SampledYEPHP(:, s);

    ParasitolCases = [Fitted.SampledActiveS + Fitted.SampledActive1 + Fitted.SampledActive2 + Fitted.SampledPassive1 + Fitted.SampledPassive2 - Fitted.SampledActiveSFP - Fitted.SampledActive1FP - Fitted.SampledActive2FP, ...
                      SampledActiveS + SampledActive1 + SampledActive2 + SampledPassive1 + SampledPassive2 - SampledActiveSFP - SampledActive1FP - SampledActive2FP]; % M9
    SerolCases = [Fitted.SampledActiveS + Fitted.SampledActive1 + Fitted.SampledActive2 + Fitted.SampledPassive1 + Fitted.SampledPassive2, SampledActiveS + SampledActive1 + SampledActive2 + SampledPassive1 + SampledPassive2]; % M9
    
    %%% Reactive part
    NumReact = size(ReactiveParameters, 1) - 1;
    MtoAbsScaling = Data.PopGrowth .^ double(Years - Data.PopSizeYear);
    Pop = round(Data.N_H * MtoAbsScaling);
    
    ActiveIS = ActiveS; % M9
    for r = 1 : NumReact
        R = ['React', num2str(r)];
        % Reactive Parameters
        React = table2struct(ReactiveParameters(R,:));
        Outputs.(R) = struct('PostID', PostID, 'SampPostID', SampPostID, 'Years', Years, 'AS', AS, 'VC', VC,  'Transmission', Transmission, 'PerfectSpec', PerfectSpec, ...
                             'ActiveS', ActiveIS, 'Active1', Active1, 'Active2', Active2, 'ActiveSFP', ActiveSFP, 'Active1FP', Active1FP, 'Active2FP', Active2FP, 'Passive1', Passive1, 'Passive2', Passive2, ...
                             'Deaths', Deaths, 'PersonYrsS', PersonYrsS, 'PersonYrs1', PersonYrs1, 'PersonYrs2', PersonYrs2, 'NewInf', NewInf, 'RDTS', RDTS, 'RDT1', RDT1, 'RDT2', RDT2, 'RDTFP', RDTFP, ...
                             'YEPHP', YEPHP, 'YEOT', YEOT, 'SampledYEPHP', SampledYEPHP,...
                             'SampledActiveS', SampledActiveS, 'SampledActive1', SampledActive1, 'SampledActive2', SampledActive2, 'SampledActiveSFP', SampledActiveSFP, 'SampledActive1FP', SampledActive1FP, 'SampledActive2FP', SampledActive2FP, ...
                             'SampledPassive1', SampledPassive1, 'SampledPassive2', SampledPassive2, 'SampledDeaths', SampledDeaths, ...
                             'SampledRDTS', SampledRDTS, 'SampledRDT1', SampledRDT1, 'SampledRDT2', SampledRDT2, 'SampledRDTFP', SampledRDTFP); % M9
%         switch React.ReactiveASnum
%             case 'AsStrat'
%                 ReactPeopleScreened = StratPeopleScreened;
%             otherwise % all others, value has been updated in Run.m 
%                 ReactPeopleScreened = React.ReactiveASnum;
%         end
        ReactPeopleScreened = React.ReactiveASnum;
        OneoffPeopleScreened = React.OneoffASnum;
        
        ReactStrat = React.ReactiveASstrat;
        OneoffStrat = React.OneoffASstrat;

        switch React.ReactiveTargetFreq
            case 'AsStrat'
                ReactTargetFreq = StratTargetFreq;
            otherwise
                ReactTargetFreq = React.ReactiveTargetFreq;
        end % ScalingUpVC002only
        
        if ReactTargetFreq > 0
            TargetFreq = VCstarted * HisTargetFreq + (1 - VCstarted) * StratTargetFreq;
            switch React.ReactiveTargetCoverage
                case 'AsStrat'
                    ReactTargetDie = StratTargetDie;
                case 'full'
                    ReactTargetDie = round(GetTargetDie(ProjStrat.NewVC, TargetFreq, ProjStrat.NewTrapCycle, Paras), 4);
                case 'realistic'
                    ReactTargetDie = round(GetTargetDie(ProjStrat.NewVC * Paras.RealisticVCwgt, TargetFreq, ProjStrat.NewTrapCycle, Paras), 4);
                case 'scaleback'
                    ReactTargetDie = round(GetTargetDie(ProjStrat.NewVC * Paras.VCwgt_scaleback, TargetFreq, ProjStrat.NewTrapCycle, Paras), 4);
                otherwise
                    ReactTargetDie = round(GetTargetDie(ProjStrat.NewVC * React.ReactiveTargetCoverage, TargetFreq, ProjStrat.NewTrapCycle, Paras), 4);
            end % ScaleBackVC
        else
            ReactTargetDie = 0;
        end % ScalingUpVC002only

        
        
%         RyearVC = max((ProjStrat.NewVCyear + (1 - VCstarted) * Paras.MinVC) * ((VCstarted + StratTargetFreq) ~= 0), React.Ryear); 
        RyearVC = max(ceil(VCstarted * Paras.VCstart + (1 - VCstarted) * ProjStrat.NewVCyear * (StratTargetDie ~= 0)) + Paras.MinVC, React.Ryear); %ScalingUpVC
        Y_RVC = find(Years == RyearVC);
 
        RunReact = ReactInfo(ReactInfo.Reactive == r, :);
        ReactStats = [];
        for k = 1 : size(RunReact, 1)
            %[s r k]
            % Replace values of fitted parameters in Paras by the values from Posterior
            p = RunReact.Posterior(k);
            for i = 1 : length(Paras.FittedNames)
                Paras.(Paras.FittedNames{i}) = Posterior.(Paras.FittedNames{i})(p);
                %Paras.(Paras.FittedNames{i}) = Posterior{p, i};
            end
            
            % Rows for replacement
            row = find(PostID == p);
            Srow = find(SampPostID == p, 1) + RunReact.SampleID(k) - 1;
            
            [ASstatus, VCstatus] = deal(ones(1, length(Yrs)));
            ASstatus(RunReact.Year(k) - Years(1) + 1) = RunReact.AS(k);
            VCstatus(RunReact.Year(k) - Years(1) + 1) = RunReact.VC(k);
            RS = RunReact.RS(k);
            OS = RunReact.OS(k);
            SwitchToRS = Years(end) + 1;
            YearsOS = [];
            YearsRS = [];
            YearsRVC = [];

%             if React.ReactiveAS == 0
%                 ASstatus(RunReact.Year(k) - Years(1) + 2 : end) = 0;
%             end
%             if React.ReactiveVC == 0
%                 VCstatus(RunReact.Year(k) - Years(1) + 2 : end) = 0;
%             end
            
            %Paras.alpha = RunReact.alpha(k);
            if RunReact.PerfectSpec(k) == 1
                Paras.specificity = 1;
            end
            Paras.alpha = Paras.alpha0 * (RunReact.Transmission(k) == 1);
            
            Paras.death = (1-Paras.u) * Paras.gamma_H;
            Paras.effcy = 1;
            
            meff = RunReact.meff(k);
            ICs = {[RunReact.S_H1(k), RunReact.S_H2(k), RunReact.S_H3(k), RunReact.S_H4(k)], [RunReact.E_H1(k), RunReact.E_H2(k), RunReact.E_H3(k), RunReact.E_H4(k)], ...
                   [RunReact.IS_H1(k), RunReact.IS_H2(k), RunReact.IS_H3(k), RunReact.IS_H4(k)], ... % M9
                   [RunReact.I1_H1(k), RunReact.I1_H2(k), RunReact.I1_H3(k), RunReact.I1_H4(k)], [RunReact.I2_H1(k), RunReact.I2_H2(k), RunReact.I2_H3(k), RunReact.I2_H4(k)], ...
                   [RunReact.R_H1(k), RunReact.R_H2(k), RunReact.R_H3(k), RunReact.R_H4(k)], RunReact.S_A(k), RunReact.E_A(k), RunReact.I1_A(k), ...
                    RunReact.P_V(k), RunReact.S_V(k), RunReact.G_V(k), RunReact.E1_V(k), RunReact.E2_V(k), RunReact.E3_V(k), RunReact.I_V(k)};
            
            %Cases = [ones(1, abs(React.StoppingAS - React.StoppingVC)) zeros(1, min([React.StoppingAS React.StoppingVC]))];
            
            switch ProjStrat.CaseDef
                case 'Parasitol'
                    Cases = ParasitolCases(Srow, RunReact.Year(k) - Y0 - max([React.StoppingAS React.StoppingVC]) + 1 : RunReact.Year(k) - Y0);
                case 'Serol'
                    Cases = SerolCases(Srow, RunReact.Year(k) - Y0 - max([React.StoppingAS React.StoppingVC]) + 1 : RunReact.Year(k) - Y0);
            end
                        
            % run ODE year by year
            ReactYears = RunReact.Year(k) : ProjStrat.SIMyear;
            for y = 1 : length(ReactYears)
                %ReactYears(y)
                Data.Years = ReactYears(y);
                Y = Data.Years - Years(1) + 1; % column for replacement
                
%                 % No transmissons after EOT 
%                 if Paras.alpha ~= 0 && Outputs.(R).NewInf(Srow, Y-1) < 1
%                     Paras.alpha = 0;
%                 end
                
                % Update AS
                ProjStrat.NewASyear = Data.Years; % Multi-period AS
                if ASstatus(Y) == 0 % no AS this year
                    Data.PeopleScreened = 0;
                    ProjStrat.NewASstrat = 'traditional'; % Multi-period AS
                elseif RS == 1
                    if OS == 1 % one-off screening
                        Data.PeopleScreened = OneoffPeopleScreened;
                        ProjStrat.NewASstrat = OneoffStrat;
                        %ProjStrat.NewASyear = Data.Years;
                        if OneoffPeopleScreened ~= 0
                            YearsOS = [YearsOS Data.Years];
                        end
                    else % RS
                        Data.PeopleScreened = ReactPeopleScreened;
                        ProjStrat.NewASstrat = ReactStrat;
                        %ProjStrat.NewASyear = Data.Years;
                        if ReactPeopleScreened ~= 0
                            YearsRS = [YearsRS Data.Years];
                        end
                    end
                else % RS = 0; either 1st-3rd regular AS or the 5th year (assuming the same as AS, not one-off)
                    Data.PeopleScreened = StratPeopleScreened(Y);
                    ProjStrat.NewASstrat = ASstrat{Y}; % Multi-period AS
                    %ProjStrat.NewASyear = Data.Years; % Multi-period AS
                end

% %                 % Update AS
% %                 if ASstatus(Y) == 0 % no AS this year
% %                     Data.PeopleScreened = 0;
% %                     ActiveS = 0;
% %                     ActiveS1 = 0;
% %                     ActiveS1FP = 0;
% %                     ActiveS2FP = 0;
% %                     RDTS1 = 0;
% %                     RDTS2 = 0;
% %                     RDTSFP = 0;
% %                 elseif sum(ASstatus(1 : Y - 1)) < Y - 1 % AS has stopped before
% % %                     AsStrat = strcmp('AsStrat', React.ReactiveASnum);
% % %                     Data.PeopleScreened = AsStrat * StratPeopleScreened + (1 - AsStrat) * max([str2num(React.ReactiveASnum) 0]);
% %                     Data.PeopleScreened = ReactPeopleScreened;
% %                 else % AS has never stopped
% %                     Data.PeopleScreened = StratPeopleScreened;
% %                 end


                switch ProjStrat.NewASstrat
                    case 'traditional'
                        ScreeningCapacity = Paras.k1;
                    case 'equal'
                        ScreeningCapacity = 1;
                    case 'high'
                        ScreeningCapacity = Paras.k1 + Paras.k2 + Paras.k4;
                end
                ScaledPeopleScreened = Data.PeopleScreened / MtoAbsScaling(Y);
                ExpectedScreeningTimes = max([ceil(ScaledPeopleScreened / Data.N_H / ScreeningCapacity) 1]);
            
                Data.ModelScreeningFreq = 365 / ExpectedScreeningTimes * ones(1, ExpectedScreeningTimes);
                Data.ModelScreeningTime = Data.Years + (0 : ExpectedScreeningTimes - 1) / ExpectedScreeningTimes;
                Data.ModelPeopleScreened = ScaledPeopleScreened / ExpectedScreeningTimes * ones(1, ExpectedScreeningTimes);
            
%                 if ScaledPeopleScreened < Data.N_H * Paras.ScreeningCapacity
%                 	Data.ModelScreeningFreq = 365;
%                 	Data.ModelScreeningTime = Data.Years;
%                 	Data.ModelPeopleScreened = round(ScaledPeopleScreened);
%                 else
%                     Data.ModelScreeningFreq = [365/2 365/2];
%                     Data.ModelScreeningTime = [Data.Years Data.Years+0.5];
%                     Data.ModelPeopleScreened = [round(ScaledPeopleScreened/2) round(ScaledPeopleScreened/2)];
%                 end
            
                % Update VC
                if VCstatus(Y) == 0 % no VC this year
                    ProjStrat.NewTargetFreq = 0;
                    ProjStrat.NewTargetDie = 0;
                elseif sum(VCstatus(1 : Y - 1)) < Y - 1 % VC has stopped before
                    ProjStrat.NewTargetFreq = ReactTargetFreq;
                    ProjStrat.NewTargetDie = ReactTargetDie;
                    if ReactTargetFreq ~= 0
                        YearsRVC = [YearsRVC Data.Years];
                    end
                else % VC has never stopped (cound be before NewVCyear but ODE will use historical values)
                    ProjStrat.NewTargetFreq = StratTargetFreq;
                    ProjStrat.NewTargetDie = StratTargetDie;
                end
                
                Outputs.(R).AS(Srow, Y) = Data.PeopleScreened;
                Outputs.(R).VC(Srow, Y) = (Y < Y_NewVC) * HisTargetFreq + (Y >= Y_NewVC) * ProjStrat.NewTargetFreq;
                
                [Classes, Aggregate] = ODEHATmodel(meff, ICs, Data, Paras, ProjStrat);
                %Aggregate
                
                % Store projected dynamics
                Outputs.(R).ActiveS(Srow, Y) = Aggregate.ActiveMS * MtoAbsScaling(Y); % M9
                Outputs.(R).Active1(Srow, Y) = (Aggregate.ActiveM1 + Paras.S1givenFP * Aggregate.ActiveMFP) * MtoAbsScaling(Y);
                Outputs.(R).Active2(Srow, Y) = (Aggregate.ActiveM2 + (1 - Paras.S1givenFP) * Aggregate.ActiveMFP) * MtoAbsScaling(Y);
        	    %Outputs.(R).ActiveSFP(Srow, Y); % M9
                Outputs.(R).Active1FP(Srow, Y) = Paras.S1givenFP * Aggregate.ActiveMFP * MtoAbsScaling(Y);
                Outputs.(R).Active2FP(Srow, Y) = (1 - Paras.S1givenFP) * Aggregate.ActiveMFP * MtoAbsScaling(Y);
                Outputs.(R).Passive1(Srow, Y) = Aggregate.PassiveM1 * MtoAbsScaling(Y);
                Outputs.(R).Passive2(Srow, Y) = Aggregate.PassiveM2 * MtoAbsScaling(Y);
                Outputs.(R).Deaths(Srow, Y) = Aggregate.DeathsM * MtoAbsScaling(Y);
                Outputs.(R).PersonYrsS(Srow, Y) = Aggregate.PersonYrsMS * MtoAbsScaling(Y); % M9
                Outputs.(R).PersonYrs1(Srow, Y) = Aggregate.PersonYrsM1 * MtoAbsScaling(Y);
                Outputs.(R).PersonYrs2(Srow, Y) = Aggregate.PersonYrsM2 * MtoAbsScaling(Y);
                Outputs.(R).NewInf(Srow, Y) = Aggregate.NewInfM * MtoAbsScaling(Y);
                Outputs.(R).RDTS(Srow, Y) = Aggregate.RDTMS * MtoAbsScaling(Y); % M9
                Outputs.(R).RDT1(Srow, Y) = Aggregate.RDTM1 * MtoAbsScaling(Y);
                Outputs.(R).RDT2(Srow, Y) = Aggregate.RDTM2 * MtoAbsScaling(Y);
                Outputs.(R).RDTFP(Srow, Y) = Aggregate.RDTMFP * MtoAbsScaling(Y);
                 
                % Sampling from projected dynamics
                if Data.PeopleScreened ~= 0
                    Disp_act = Paras.disp_act * (Aggregate.ActiveMS + Aggregate.ActiveM1 + Aggregate.ActiveM2) ./ (Aggregate.ActiveMS + Aggregate.ActiveM1 + Aggregate.ActiveM2 + Aggregate.ActiveMFP); % M9
                    ActiveS = Cbetabinornd(Data.PeopleScreened, (Aggregate.ActiveMS + Aggregate.ActiveM1 + Aggregate.ActiveM2 + Aggregate.ActiveMFP) / ScaledPeopleScreened, Disp_act); % M9
                    ActiveS1 = Cbinornd(ActiveS, (Aggregate.ActiveM1 + Paras.S1givenFP * Aggregate.ActiveMFP) / (Aggregate.ActiveMS + Aggregate.ActiveM1 + Aggregate.ActiveM2+ Aggregate.ActiveMFP)); % M9
                    ActiveS1FP = Cbinornd(ActiveS1, Paras.S1givenFP * Aggregate.ActiveMFP / (Aggregate.ActiveM1 + Paras.S1givenFP * Aggregate.ActiveMFP));
                    ActiveS2 = Cbinornd(ActiveS - ActiveS1, (Aggregate.ActiveM2 + (1 - Paras.S1givenFP) * Aggregate.ActiveMFP) / (Aggregate.ActiveMS + Aggregate.ActiveM2 + (1 - Paras.S1givenFP) * Aggregate.ActiveMFP)); % M9
                    ActiveS2FP = Cbinornd(ActiveS2, (1 - Paras.S1givenFP) * Aggregate.ActiveMFP / (Aggregate.ActiveM2 + (1 - Paras.S1givenFP) * Aggregate.ActiveMFP));
                    ActiveSS = ActiveS - ActiveS1 - ActiveS2; % M9
                    %ActiveSSFP = zeros(size(ActiveSS)); % M9

                    if ProjStrat.TTyear ~= 0 && Data.Years >= ProjStrat.TTyear
                        RDTSS = ActiveSS + Cbinornd(Data.PeopleScreened, (Aggregate.RDTMS - Aggregate.ActiveMS) / ScaledPeopleScreened); % M9
                        RDTS1 = ActiveS1 + Cbinornd(Data.PeopleScreened, (Aggregate.RDTM1 - Aggregate.ActiveM1) / ScaledPeopleScreened);
                        RDTS2 = ActiveS2 + Cbinornd(Data.PeopleScreened, (Aggregate.RDTM2 - Aggregate.ActiveM2) / ScaledPeopleScreened);
                        RDTSFP = Cbinornd(Data.PeopleScreened, 1 - ProjStrat.SpecRDT);
                    else
                        RDTSS = 0; % M9
                        RDTS1 = 0;
                        RDTS2 = 0;
                        RDTSFP = 0;
                    end
                else
                    ActiveS = 0;
                    ActiveS1 = 0;
                    ActiveS1FP = 0;
                    ActiveS2 = 0; % M9
                    ActiveS2FP = 0;
                    ActiveSS = 0; % M9
                    %ActiveSSFP = 0; % M9
                    RDTSS = 0; % M9
                    RDTS1 = 0;
                    RDTS2 = 0;
                    RDTSFP = 0;
                end
                PassiveS = Cbetabinornd(Pop(Y), (Aggregate.PassiveM1 + Aggregate.PassiveM2) / Data.N_H, Paras.disp_pass);
                PassiveS1 = Cbinornd(PassiveS, Aggregate.PassiveM1 / (Aggregate.PassiveM1 + Aggregate.PassiveM2));
                DeathsS = Cbinornd(Pop(Y), Aggregate.DeathsM / Data.N_H);
                
                
                % Store sampled dynamics
                Outputs.(R).SampledActiveS(Srow, Y) = ActiveSS; % M9
                Outputs.(R).SampledActive1(Srow, Y) = ActiveS1;
                Outputs.(R).SampledActive2(Srow, Y) = ActiveS2; % M9
                %Outputs.(R).SampledActiveSFP(Srow, Y) = ActiveSSFP; % M9
                Outputs.(R).SampledActive1FP(Srow, Y) = ActiveS1FP;
                Outputs.(R).SampledActive2FP(Srow, Y) = ActiveS2FP;
                Outputs.(R).SampledPassive1(Srow, Y) = PassiveS1;
                Outputs.(R).SampledPassive2(Srow, Y) = PassiveS - PassiveS1;
                Outputs.(R).SampledDeaths(Srow, Y) = DeathsS;
                Outputs.(R).SampledRDTS(Srow, Y) = RDTSS; % M9
                Outputs.(R).SampledRDT1(Srow, Y) = RDTS1;
                Outputs.(R).SampledRDT2(Srow, Y) = RDTS2;
                Outputs.(R).SampledRDTFP(Srow, Y) = RDTSFP;

                % Check if cases = 0 and update AS and VC
                %SampledCases = [SampledCases ActiveS + PassiveS];
                switch ProjStrat.CaseDef
                    case 'Parasitol'
                        Cases = [Cases ActiveS - ActiveS1FP - ActiveS2FP + PassiveS];
                    case 'Serol'
                        Cases = [Cases ActiveS + PassiveS];
                end
                
                if RS == 1 % Reactive screening algorithm started already
                    if Cases(end) ~= 0 && ASstatus(Y) == 0 && OneoffPeopleScreened ~= 0 % No AS but cases from PS 
                        ASstatus(Y + 1) = 1;
                        OS = 1;
                    elseif Cases(end) ~= 0 && (ASstatus(Y) == 1 || OneoffPeopleScreened == 0) % Had AS (could be RS or OS) and cases from either AS or PS % need to check again
                        ASstatus(Y + 1) = 1;
                        OS = 0; % will use RS value next year
                    else % No cases
                        ASstatus(Y + 1) = (1 - OS) * (sum(Cases(end - React.ReactiveAS + 1 : end)) > 0);
                        %OS = ASstatus(Y + 1) == 0;
                    end
                else % Not achieved RS criterion yet (4th or 5th year of the WHO classic algorithm or found cases in the 4th or 5th year and try to get a new set of zero cases)
                    if ASstatus(Y) == 0 % 4th year
                        ASstatus(Y + 1) = 1;
                        OS = Cases(end) == 0;                        
                    elseif OS == 1 % 5th year
                        if Cases(end) == 0
                            ASstatus(Y + 1) = 0;
                            %OS = 0;
                            RS = 1;
                            SwitchToRS = Yrs(Y+1);
                        else
                            ASstatus(Y + 1) = 1;
                            OS = 0;
                        end
                    else % 1st-3rd year
                        ASstatus(Y + 1) = sum(Cases(end - React.StoppingAS + 1 : end)) > 0;
                        OS = ASstatus(Y + 1) == 0; % 
                    end
                end

                if sum(VCstatus(1 : Y)) ~= Y % VC has stopped before (including 1. both VC and RVC are zero 2. VCstart = 0 and NewVC = 0), 
                    if RS == 1
                        VCstatus(Y + 1) = (React.ReactiveVC ~= 0) * (ASstatus(Y + 1) == 1) * (OS == 0); % only when RS happens next year
                    else % use Reactive.VC to check if no VC next year (might happen when VC ceases before AS)
                        VCstatus(Y + 1) = sum(Cases(end - React.ReactiveVC + 1 : end)) > 0;
                    end
                elseif Y + 1 >= Y_RVC % VC hasn't stopped before & next year > NewVCyear + MinVC, use Stopping.VC to check if no VC next year
                    VCstatus(Y + 1) = sum(Cases(end - React.StoppingVC + 1 : end)) > 0;
                end

% % %                 %if ActiveS + PassiveS == 0
% % %                     % No cases this year then need to check if no interventions next year
% % %                     if sum(ASstatus(1 : Y)) ~= Y % AS has stopped before (including both AS and RS are zero), use Reactive.AS to check if no AS next year
% % %                         ASstatus(Y + 1) = sum(Cases(end - React.ReactiveAS + 1 : end)) > 0;
% % %                     else % AS hasn't stopped before, use Stopping.AS to check if no AS next year
% % %                         ASstatus(Y + 1) = sum(Cases(end - React.StoppingAS + 1 : end)) > 0;
% % %                     end
% % % 
% % %                     if sum(VCstatus(1 : Y)) ~= Y % VC has stopped before (including both VC and RVC are zero), use Reactive.VC to check if no VC next year
% % %                         VCstatus(Y + 1) = (sum(Cases(end - React.ReactiveVC + 1 : end)) > 0);
% % %                     elseif Y + 1 >= Y_RVC% VC hasn't stopped before, use Stopping.VC to check if no VC next year (need to be after NewVCyear + MinVC)
% % %                         VCstatus(Y + 1) = sum(Cases(end - React.StoppingVC + 1 : end)) > 0;
% % %                     end
% % %                 %else
% % %                     % Cases are found and interventions are needed next year
% % %                 %    ASstatus(Y + 1) = 1 * ();
% % %                 %    VCstatus(Y + 1) = 1;
% % %                 %end
                
                
                
                % Update ICs
                pop = Classes(end, :);
                ICs = {[pop.S_H1, pop.S_H2, pop.S_H3, pop.S_H4], [pop.E_H1, pop.E_H2, pop.E_H3, pop.E_H4], [pop.IS_H1, pop.IS_H2, pop.IS_H3, pop.IS_H4],... % M9
                       [pop.I1_H1, pop.I1_H2, pop.I1_H3, pop.I1_H4], [pop.I2_H1, pop.I2_H2, pop.I2_H3, pop.I2_H4], [pop.R_H1, pop.R_H2, pop.R_H3, pop.R_H4],...
                        pop.S_A, pop.E_A, pop.I1_A, pop.P_V, pop.S_V, pop.G_V, pop.E1_V, pop.E2_V, pop.E3_V, pop.I_V};
                    
                % Update PerfectSpec and Transmission
                Outputs.(R).Transmission(Srow, Y) = Aggregate.Transmission;
                Outputs.(R).PerfectSpec(Srow, Y) = Aggregate.PerfectSpec;
                if Paras.alpha ~= 0 && Outputs.(R).NewInf(Srow, Y) < Paras.EOTthreshold && ...
                   (pop.IS_H1 + pop.IS_H2 + pop.IS_H3 + pop.IS_H4 + pop.I1_H1 + pop.I1_H2 + pop.I1_H3 + pop.I1_H4 + pop.I2_H1 + pop.I2_H2 + pop.I2_H3 + pop.I2_H4) * MtoAbsScaling(Y) < Paras.EOTthreshold % M9
                    Paras.alpha = 0;
                end
                if Paras.specificity ~= 1 && Aggregate.PerfectSpec ~= 0
                    Paras.specificity = 1;
                end    
            end

            % Recalculate YEPHP, YEOT and SampledYEPHP
            FittedDynamics = Fitted.ActiveS(row, :) + Fitted.Active1(row, :) + Fitted.Active2(row, :) + Fitted.Passive1(row, :) + Fitted.Passive2(row, :); % M9
            FittedSamplings = Fitted.SampledActiveS(Srow, :) + Fitted.SampledActive1(Srow, :) + Fitted.SampledActive2(Srow, :) + Fitted.SampledPassive1(Srow, :) + Fitted.SampledPassive2(Srow, :); % M9
            ProjDyanmics = Outputs.(R).ActiveS(Srow, :) + Outputs.(R).Active1(Srow, :) + Outputs.(R).Active2(Srow, :) + Outputs.(R).Passive1(Srow, :) + Outputs.(R).Passive2(Srow, :); % M9
            ProjSamplings = Outputs.(R).SampledActiveS(Srow, :) + Outputs.(R).SampledActive1(Srow, :) + Outputs.(R).SampledActive2(Srow, :) + Outputs.(R).SampledPassive1(Srow, :) + Outputs.(R).SampledPassive2(Srow, :); % M9
            
            WholeDynamics = [FittedDynamics ProjDyanmics];
            WholeSamplings = [FittedSamplings ProjSamplings];
            WholeNewInf = [Fitted.NewInf(row, :) Outputs.(R).NewInf(Srow, :)];
            
            Outputs.(R).YEPHP(Srow) = max([find(WholeDynamics * 10000 >= Data.N_H * scaling, 1, 'last') 0]) + Fitted.Years(1);
            Outputs.(R).YEOT(Srow) = max([find(WholeNewInf >= Paras.EOTthreshold, 1, 'last') 0]) + Fitted.Years(1);
            Outputs.(R).SampledYEPHP(Srow) = max([find(WholeSamplings * 10000 >= Data.N_H * scaling, 1, 'last') 0]) + Fitted.Years(1);
            
            
            % ReactStats here
            %SwitchRS = (StratPeopleScreened ~= 0) * min([Years(find(Outputs.(R).AS(Srow, :) == 0, 1)) Years(end) + 1]); % 0 means 
            SwitchToRVC = min([Years(find(VCstatus(1:end-1) == 0, 1)) Years(end) + 1]); % Year that VC ceases
            
            % Fisrt time that OS, RS, RVC happen
            FirstOS = min([YearsOS Years(end) + 1]);
            FirstRS = min([YearsRS Years(end) + 1]);
            FirstRVC = min([YearsRVC Years(end) + 1]);
%             FirstRS = min([Years(YidRS) Years(end) + 1]);
%             FirstRVC = min([Years(YidRVC) Years(end) + 1]);
            
            % Times that OS happens
            TimesOS = length(YearsOS);

            % Times that RS happens
            if isempty(YearsRS) % No RS
                TimesRS = 0;
            else % RS happens
                TimesRS = sum(YearsRS - [0 YearsRS(1 : end-1)] > 1);
            end
            
            % Times that RVC happens
            if isempty(YearsRVC) % No RVC
                TimesRVC = 0;
            else % RVC happens
                TimesRVC = sum(YearsRVC - [0 YearsRVC(1 : end-1)] > 1);
            end


%             FirstRS = (React.ReactiveAS ~= 0 & ReactPeopleScreened ~= 0) * min([Years(find(ReactAS ~= 0, 1) + Y - 1) Years(end) + 1]);
%             FirstRVC = (React.ReactiveVC ~= 0 & ReactTargetFreq ~= 0) * min([Years(find(ReactVC ~= 0, 1) + Y - 1) Years(end) + 1]);
%             TimesRS = sum(ReactAS(2:end)-ReactAS(1:end-1) > 0); % times of restart AS
%             TimesRVC = sum(ReactVC(2:end)-ReactVC(1:end-1) > 0); % times of restart VC
            
%             ASrestarts = sum(AS(2:end)-AS(1:end-1) > 0); % times of restart AS
%             VCrestarts = sum(VC(2:end)-VC(1:end-1) > 0); % times of restart VC
%             YearsFirstStopAS = max([find(AS ~= 0, 1) - 1, 0]); % years of first stopping in AS
%             YearsFirstStopVC = max([find(VC ~= 0, 1) - 1, 0]); % years of first stopping in VC

            ReactStats = [ReactStats; [RunReact.Posterior(k) Srow RunReact.Year(k) SwitchToRS SwitchToRVC FirstOS FirstRS FirstRVC TimesOS TimesRS TimesRVC]];
        end

        % Post-processing of NewInf (NewInf = 0 after EoT rather than after EoI)
        %Outputs.(R).NewInf(Outputs.(R).NewInf < Paras.EOTthreshold) = 0;
        
        % Recalculate PEPHP and PEOT
        for y = 1 : length(Elim.ElimByYears)
            Y = Elim.ElimByYears(y);
            Outputs.(R).PEPHP(Y - Elim.ElimByYears(1) + 1) = mean(Outputs.(R).YEPHP <= Y);
            Outputs.(R).PEOT(Y - Elim.ElimByYears(1) + 1) = mean(Outputs.(R).YEOT <= Y);
            Outputs.(R).SampledPEPHP(Y - Elim.ElimByYears(1) + 1) = mean(Outputs.(R).SampledYEPHP <= Y);
        end
        
%         Outputs.(R).PostID = PostID;
%         Outputs.(R).SampPostID = SampPostID;
%         Outputs.(R).Years = Years;
        Outputs.(R).ElimByYears = Elim.ElimByYears;
        %Outputs.(R).ReactStats = ReactStats;
        if size(ReactStats, 1) == 0
            ReactStats = zeros(1, 11);
        end
        Outputs.(R).ReactStats = array2table(ReactStats, 'VariableNames', {'PostID', 'RowID', 'RunReactYear', 'SwitchToRS', 'SwitchToRVC', 'FirstOS', 'FirstRS', 'FirstRVC', 'TimesOS', 'TimesRS', 'TimesRVC'});
    end
    
  
    
    

