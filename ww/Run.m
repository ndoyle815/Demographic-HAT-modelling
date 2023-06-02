
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                     %
%   This code processes the raw data, reorganises the parameters and runs the simulations of the Warwick HAT model.   %
%                                                                                                                     %
%   Inputs:                                                                                                           %
%       Cloc - a string in the format of ALPHA-3 country codes                                                        %
%       Lv1loc - an integer, provine/district/focus index for DRC/CIV,GIN,TCD/UGA                                     %
%       Lv2loc - an integer, health zone/county index for DRC/UGA                                                     %
%       Lv3loc - an integer, health area index for DRC                                                                %
%       Model - a cell array containing the model indices                                                             %
%       EndCaseStr
%       EndScrStr
%       ParaStr - a string of 3-digits ID related to parameter settings                                               %
%       StratDefStr
%       MinimumData - an integer indicating the minimum requirement in Data                                           %
%       RunMCMC - an integer determining to run MCMC (1 - from scratch, 2 - re-adapt, or                              %
%                 3 - additional sampling from posterior) or no MCMC (0)                                              %
%       RunEvidence - integer number of importance samples to use, or 0(skip)                                         %
%       RunEnsProjection
%       RunEnsFromExisting
%       RunProjection - an integer denoting the number of realizations used in Projection                             %
%       StratMin - an integer denoting the first strategy simulated in (Ensemble) Projection                          %
%       RunSamples - an integer denoting the number of samples from ODE                                               %
%       RunReactive - a boolean determining to run Reactive or not                                                    %
%       RunCFS - a cell array determining to run specific CounterFactual scenarios or not                             %
%       MachineSetting - structure containing parameters to control parcluster                                        %
%                                                                                                                     %
%   Main output files:                                                                                                %
%       Output_MCMC - a file containing MCMC settings and outcomes                                                    %
%       Posterior - a table containing the potseriors for fitted parameters and corresponding likelihood              %
%       Fitted - matrices containing key outputs (e.g. cases, new infections, deaths) in data period                  %
%       ProjectionICs - a table containing initial conditions for forcasting                                          %
%       Projection - matrices containing key outputs (e.g. cases, new infections, deaths) in forcasting period        %
%       Elimination - matrices containing years and probabilities of EPHP and EOT under different strategies          %
%       ReactInfo - a table containing key inputs (e.g. year and initial conditions) of reactive interventions        %
%                                                                                                                     %
%   Note: hosts are (1) low-risk, random participants                                                                 %
%                   (2) high-risk, random participants                                                                %
%                   (3) low-risk, non-participants                                                                    %
%                   (4) high-risk, non-participants                                                                   %
%                   (5) reservoir animals                                                                             %
%                   (6) non-reservoir animals, no dynamics and is ignored                                             %
%                                                                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Run(Cloc, Lv1loc, Lv2loc, Lv3loc, Model, EndCaseStr, EndScrStr, ParaStr,...
    StratDefStr, MinimumData, RunMCMC, MCMCOptions, RunEvidence, RunEnsProjection,...
    RunEnsFromExisting, RunProjection, StratMin, RunSamples, RunDisaster,...
    RunReactive, RunCFS, MachineSetting)
    
    RunInfo = RunInformation(Cloc, Lv1loc, Lv2loc, Lv3loc,...
        Model, EndCaseStr, EndScrStr, ParaStr, StratDefStr,...
        RunMCMC, MCMCOptions, RunProjection, RunEnsProjection, RunEnsFromExisting,...
        RunSamples, RunReactive, RunCFS, StratMin);
    
    load(RunInfo.DataPath);
    load(RunInfo.ParasFile);
    copyfile(RunInfo.ParasFile, RunInfo.OutputParasPath);
    
    %%% Parameters for current location
    % FittedParameters for current location...
    locFittedPar = cell2table(cell(0,width(FittedParameters)),'VariableNames',FittedParameters.Properties.VariableNames);
    % work from most-local to least-local area name
    for i = 1:length(RunInfo.names)
        % Keep line if 1) Location matches current level, and
        %              2) parameter hasn't been matched at a more local level.
        choose = (strcmp(FittedParameters.Location, RunInfo.names{i}) .* ...
                  ~ismember(FittedParameters.Notation, char(locFittedPar.Notation))) == 1;
        locFittedPar = [locFittedPar; FittedParameters(choose,:)];
    end
    locFittedPar = locFittedPar(:,~strcmp(locFittedPar.Properties.VariableNames,'Location')); % drop location
    locFittedPar = [locFittedPar(~startsWith(locFittedPar.Notation,'active_neg_'),:);...                       % if fitting neg active tests,
                    sortrows(locFittedPar(startsWith(locFittedPar.Notation,'active_neg_'),:), 'Notation')];    % put that chunk, sorted, last.
                
    % FixedParameters for current location...
    x = 0;
    a = 0;
    while x == 0
        a = a + 1;
        x = strcmp(FixedParameters.Location, RunInfo.names{a});
    end
    locFixedPar = FixedParameters(strcmp(FixedParameters.Location, RunInfo.names{a}), ...
                                  ~strcmp(FixedParameters.Properties.VariableNames,'Location'));
    
    % InterventionParameters for current location
    x = 0;
    a = 0;
    while x == 0
        a = a + 1;
        x = contains(InterventionParameters.Location, RunInfo.names{a});
    end
    locIntPar = InterventionParameters(strcmp(InterventionParameters.Location, RunInfo.names{a}), ...
                                       ~strcmp(InterventionParameters.Properties.VariableNames,'Location'));
    IntrpParas.PPSstart = locIntPar.PPSstart;
    IntrpParas.PPSend = locIntPar.PPSend;
    IntrpParas.NoPSstart = locIntPar.NoPSstart;
    IntrpParas.NoPSend = locIntPar.NoPSend;

    %%% Calculate autual VC reduction = VC effectiveness (VC) * Coverage (VCwgt) for loaction with existing VC
    FullTargetDie = locIntPar.TargetDie;
    if sum(strcmp(InterventionParameters.Properties.VariableNames,'B_V')) ~= 0
        lB_V = locIntPar.B_V;
    else
        lB_V = locFixedPar.B_V;
    end
    TsetseParas = struct('xi_V', locFixedPar.xi_V, 'p_survive', locFixedPar.p_survive, 'alpha', locFixedPar.alpha, 'mu_V', locFixedPar.mu_V, 'B_V', lB_V);
    if locIntPar.VCstart ~= 0 && locIntPar.VCwgt ~= 1 % B_V
        WgtTargetDie = round(GetTargetDie(locIntPar.VC * locIntPar.VCwgt, locIntPar.TargetFreq, locIntPar.TrapCycle, TsetseParas), 4);
        locIntPar.TargetDie = WgtTargetDie;
    end% ScalingUpVC % B_V
    
    
    %%% Data
    location = RunInfo.locationIdx;
    % Overwrite SCREENED by NaN when its prevalence is too high	
    % ***this next line might be causing issues when the number of cases is
    %    low but the weighted mean AS is high***
    %SCREENED(location, (ACTIVE1(location, :)+ ACTIVE2(location, :)+ ACTIVENa(location, :)) ./ SCREENED(location, :)> 0.05) = NaN;
                              
    % Are there NaN in SCREENING ===> imputation of number of neg AS tests	
    unk_scr = isnan(SCREENED(location,:));	
    n_unk_scr = sum(unk_scr);	
    if n_unk_scr > 0	
        unk_scr_years = YEAR(unk_scr);	
        for i=1:n_unk_scr	
            to_add.Notation = ['active_neg_', num2str(unk_scr_years(i))];	
            to_add.Model = 'All';	
            to_add.Initial = 5000;	
            to_add.Lower = 0;	
            to_add.Upper = Inf;	
            to_add.Distribution = 'Geometric';
            scrYears = ~isnan(SCREENED(location,:)) & SCREENED(location,:) ~= 0;
            %WtMeanAS = SCREENED(location,scrYears) * (1./abs(YEAR(scrYears)-unk_scr_years(i))') / sum(1./abs(YEAR(scrYears)-unk_scr_years(i)));
            WtMeanAS = SCREENED(location,scrYears) * exp(-abs(YEAR(scrYears)-unk_scr_years(i))') / sum(exp(-abs(YEAR(scrYears)-unk_scr_years(i))));
            %WtMeanAS = mean(SCREENED(location,scrYears));
            lambda=1.0/(1.0+WtMeanAS);
            to_add.Parameters = [lambda];	
            to_add.Initial_sigma = 0;	
            to_add.Last_year = 0;	
            locFittedPar = [locFittedPar; struct2table(to_add, 'AsArray', true)];	
        end	
    end
                              
    % Find parameters for negative active tests...
    work = startsWith(locFittedPar.Notation,'active_neg_') & ~strcmp(locFittedPar.Distribution, 'Delta');
    ActiveNeg = struct('Fitting', sum(work));
    if ActiveNeg.Fitting > 0
        ActiveNeg.Notation = table2cell(locFittedPar(work,'Notation'));
        % ...convert name to year of interest...
        ActiveNeg.Year     = str2num(cell2mat(erase(table2cell(locFittedPar(work,'Notation')), 'active_neg_')));
        % ...and index of year of interest
        [work, ActiveNeg.YearIdx] = ismember(ActiveNeg.Year, YEAR);
        % ...ensure any present-but-unacceptable number screened is not counted as a true record	
        SCREENED(location,ActiveNeg.YearIdx) = 0;
    end
    clear work;
    
    
    %%% Strategy and Reactive
    NumStrat = 1;
    NumReact = 0;
    % Only need the StratDef*.mat contents if going beyond the fitting stage	
    if RunProjection ~= 0 || RunEnsProjection ~= 0 ||...	
       RunEnsFromExisting ~= 0 || RunReactive ~= 0 	
   	
       load(RunInfo.StratDefFile);	
       save(RunInfo.OutputStratDefPath, 'ReactiveParameters', 'Strategy', 'VCwgt'); % save a copy to Result directory	
       
        % Strategy for current location
        locStrategy = cell2table(cell(0,width(Strategy)),'VariableNames',Strategy.Properties.VariableNames);
        SplitRowNames = extractBefore(strcat(Strategy.Properties.RowNames,'_'),'_');
    
        % work from most-local to least-local area name
        for i = 1:length(RunInfo.names)
            % Keep line if 1) Location matches current level, and
            %              2) parameter hasn't been matched at a more local level.
            choose = (strcmp(Strategy.Location, RunInfo.names{i}) .* ...        
                      ~ismember(SplitRowNames, locStrategy.Properties.RowNames)) == 1;
            StratToKeep = Strategy(choose,:);
            StratToKeep.Properties.RowNames = SplitRowNames(choose);
            locStrategy = [locStrategy; StratToKeep];
        end % MultiLocInStrat
        locStrategy = locStrategy(:,~strcmp(locStrategy.Properties.VariableNames,'Location')); % drop location
        locStrategy = sortrows(locStrategy,'RowNames');
        NumStrat = size(locStrategy,1) - StratMin;
        NewVCyear0 = locStrategy.NewVCyear;
        
        row = find(locStrategy.NewVC > 0 & locStrategy.NewTargetFreq == 0);
        if ~isempty(row)
            error(['Strategy.NewTargetFreq = 0 but Strategy.NewVC > 0 in ', strjoin(locStrategy.Properties.RowNames(row)), '.'])
        end
        locStrategy.NewTargetDie(locStrategy.NewVC == 0) = 0;
        
        if locIntPar.VCstart ~= 0
            if floor(locIntPar.VC_scaleback) > max([YEAR GapYear])
                error('Scaling back VC has to be in the past!')
            end % ScaleBackVC

            if locIntPar.VC_scaleback > 0 && (locIntPar.VCwgt_scaleback >= locIntPar.VCwgt) + (locIntPar.TargetFreq_scaleback >= locIntPar.TargetFreq) == 2
                error('Either InterventionParameters.VCwgt_scaleback or locIntPar.TargetFreq_scaleback should be smaller than the existing values.')
            end % ScaleBackVC

            if mod(locIntPar.VCstart, 1) > 0 && mod(locIntPar.VCstart, 1) * locIntPar.TargetFreq < 1
                if locIntPar.VC_scaleback > 0 %&& mod(locIntPar.VCstart, 1) ~= mod(locIntPar.VC_scaleback, 1)
                    locIntPar.VC_scaleback = floor(locIntPar.VC_scaleback) + mod(locIntPar.VCstart, 1);
                    warning('Scale-back VC started from the same time of the year as VC.')
                else
                    locStrategy.NewVCyear = locStrategy.NewVCyear + mod(locIntPar.VCstart, 1);
                    warning('Strategy.NewVCyear will start from the same time of the year as Intervention.VCstart')
                end    
            end % CompleteVCcycle + ScaleBackVC

            if locIntPar.VC_scaleback > 0 && locIntPar.VCwgt_scaleback == 0 
                locIntPar.TargetFreq_scaleback = 0;
                if locIntPar.VC_scaleback < locIntPar.VCstart + locIntPar.MinVC
                    error('VC scaled back too soon! VC has to last MinVC years.')
                end
            end % ScaleBackVC

            row = find(locStrategy.NewVCyear > 0 & locStrategy.NewVC == 0 & locStrategy.NewVCyear < locIntPar.VCstart + locIntPar.MinVC);
            if ~isempty(row)
                locStrategy.NewVCyear(row) = locIntPar.VCstart + locIntPar.MinVC;
                warning(['NewVCyear is delayed because of the requirement of MinVC in ', strjoin(locStrategy.Properties.RowNames(row)), '.'])
            end

            row = find(locStrategy.NewVC > 0 & locStrategy.NewVC ~= locIntPar.VC);
            if ~isempty(row)
                locStrategy.NewVC(row) = locIntPar.VC;
                warning(['Strategy.NewVC is set as InterventionParameter.VC in ', strjoin(locStrategy.Properties.RowNames(row)), '.'])
            end

            row = find(locStrategy.NewVC > 0); % & strcmp(locStrategy.NewTargetCoverage, 'full') & locStrategy.NewTargetDie ~= FullTargetDie);
            if ~isempty(row)
                locStrategy.NewTargetDie(row) = FullTargetDie;
                %warning(['Strategy.NewTargetDie is different from InterventionParameter.TargetDie in ', strjoin(locStrategy.Properties.RowNames(row)), '.'])
            end
            
            row = find(locStrategy.NewVC > 0 & locStrategy.NewTargetFreq ~= locIntPar.TargetFreq);
            if ~isempty(row)
                %locStrategy.NewTargetFreq(row) = locIntPar.TargetFreq;
                warning(['Strategy.NewTargetFreq is different from InterventionParameter.TargetFreq in ', strjoin(locStrategy.Properties.RowNames(row)), '.'])
            end

            row = find(locStrategy.NewVC > 0 & locStrategy.NewTrapCycle ~= locIntPar.TrapCycle);
            if ~isempty(row)
                locStrategy.NewTrapCycle(row) = locIntPar.TrapCycle;
                warning(['Strategy.NewTrapCycle is set as InterventionParameter.TrapCycle in ', strjoin(locStrategy.Properties.RowNames(row)), '.'])
            end

            if VCwgt.Realistic(VCwgt.LocStr == RunInfo.LocStr) ~= locIntPar.VCwgt
                %error('Different VC weights defined in Paras (InterventionParameter.VCwgt) and StratDef (VCwgt.Realistic).')
                VCwgt.Realistic(VCwgt.LocStr == RunInfo.LocStr) = locIntPar.VCwgt;
                warning('Replaced VCwgt.Realistic by InterventionParameter.VCwgt.')
            end % ScalingUpVC % If warning shows up, let Sam know about it
            
            row1 = find(strcmp(locStrategy.NewTargetCoverage, 'existing'));
            row2 = find(strcmp(locStrategy.NewTargetCoverage, 'realistic'));
            if ~isempty(row1) && ~isempty(row2)
                locStrategy(row2,:) = [];
                NumStrat = NumStrat - sum(row2 > StratMin);
            end % Reduce identical strategies
        else
            locIntPar.VC_scaleback = 0;
            locIntPar.VCwgt_scaleback = 0;
            locIntPar.TargetFreq_scaleback = 0;

            row1 = find(strcmp(locStrategy.NewTargetCoverage, 'existing'));
            row2 = find(strcmp(locStrategy.NewTargetCoverage, 'realistic'));
            if ~isempty(row1) && ~isempty(row2) && VCwgt.Realistic(VCwgt.LocStr == RunInfo.LocStr) == 0
                locStrategy(row2,:) = [];
                NumStrat = NumStrat - sum(row2 > StratMin);
            end % Reduce identical strategies
        end % ScaleBackVC
        
        for s = 1 : NumStrat
            % adjusting AS
            if length(locStrategy.NewASyear{s + StratMin}) == 0
                locStrategy.NewASyear{s + StratMin} = max([YEAR GapYear]) + 1;
            elseif locStrategy.NewASyear{s + StratMin}(1) > max([YEAR GapYear]) + 1
                if length(split(locStrategy.NewASnum{s + StratMin}, '_')) == length(locStrategy.NewASyear{s + StratMin})
                    locStrategy.NewASnum{s + StratMin} = append('mean_', locStrategy.NewASnum{s + StratMin}); 
                    locStrategy.NewASstrat{s + StratMin} = append('traditional_', locStrategy.NewASstrat{s + StratMin});
                    warning('Traditional MeanAS is assumed from the beginning of the projections.')
                end
                locStrategy.NewASyear{s + StratMin} = [max([YEAR GapYear]) + 1 locStrategy.NewASyear{s + StratMin}];
            else
                if length(split(locStrategy.NewASnum{s + StratMin}, '_')) ~= length(locStrategy.NewASyear{s + StratMin})
                    error('Missing AS num and strat from the beginning of the projections.')
                end
                locStrategy.NewASyear{s + StratMin}(1) = max([YEAR GapYear]) + 1;
            end

            % adjusting VC
            locStrategy{s + StratMin, 'NewVCyear'} = max([YEAR+1 GapYear+1 locStrategy{s + StratMin, 'NewVCyear'}]);
            TargetFreq = (locIntPar.VCstart > 0)* locIntPar.TargetFreq + (locIntPar.VCstart == 0) * locStrategy.NewTargetFreq(s + StratMin);
            if locStrategy.NewVC(s + StratMin) > 0 && strcmp(locStrategy.NewTargetCoverage{s + StratMin}, 'realistic') && VCwgt.Realistic(VCwgt.LocStr == RunInfo.LocStr) ~= 1 % B_V
%                 if locStrategy.NewVC(s + StratMin) ~= locIntPar.VC
%                     error(['Strategy.NewVC is different from InterventionParameter.VC in ', locStrategy.Properties.RowNames(s + StratMin), '.'])
%                 end
                WgtTargetDie = round(GetTargetDie(locStrategy.NewVC(s + StratMin) * VCwgt.Realistic(VCwgt.LocStr == RunInfo.LocStr), TargetFreq, locStrategy.NewTrapCycle(s + StratMin), TsetseParas), 4);
                locStrategy.NewTargetDie(s + StratMin) = WgtTargetDie;
            end
            if locStrategy.NewVC(s + StratMin) > 0 && strcmp(locStrategy.NewTargetCoverage{s + StratMin}, 'scaleback')
                WgtTargetDie = round(GetTargetDie(locStrategy.NewVC(s + StratMin) * locIntPar.VCwgt_scaleback, TargetFreq, locStrategy.NewTrapCycle(s + StratMin), TsetseParas), 4);
                locStrategy.NewTargetDie(s + StratMin) = WgtTargetDie;
            end
            if locStrategy.NewVC(s + StratMin) > 0 && strcmp(locStrategy.NewTargetCoverage{s + StratMin}, 'existing')
                WgtTargetDie = round(GetTargetDie(locStrategy.NewVC(s + StratMin) * locIntPar.VCwgt, TargetFreq, locStrategy.NewTrapCycle(s + StratMin), TsetseParas), 4);
                locStrategy.NewTargetDie(s + StratMin) = WgtTargetDie;
            end % ExistingVC
        end

        locStrategy.NewTargetFreq(locStrategy.NewTargetDie == 0) = 0; % including realistic coverage = 0

    
        % ReactiveParameters for current location
    	locReactivePar = cell2table(cell(0,width(ReactiveParameters)),'VariableNames',ReactiveParameters.Properties.VariableNames);
        SplitRowNames = extractBefore(strcat(ReactiveParameters.Properties.RowNames,'_'),'_');

        % work from most-local to least-local area name
        for i = 1:length(RunInfo.names)
            % Keep line if 1) Location matches current level, and
            %              2) parameter hasn't been matched at a more local level.
            choose = (strcmp(ReactiveParameters.Location, RunInfo.names{i}) .* ...        
                      ~ismember(SplitRowNames, locReactivePar.Properties.RowNames)) == 1;
            StratToKeep = ReactiveParameters(choose,:);
            StratToKeep.Properties.RowNames = SplitRowNames(choose);
            locReactivePar = [locReactivePar; StratToKeep];
        end % MultiLocInStrat

        locReactivePar = locReactivePar(:,~strcmp(locReactivePar.Properties.VariableNames,'Location')); % drop location
        locReactivePar = sortrows(locReactivePar,'RowNames');

        NumReact = size(locReactivePar, 1) - 1;	
        for r = 1 : NumReact	
            locReactivePar{r + 1, 'Ryear'} = max([YEAR+1 GapYear+1 locReactivePar{r + 1, 'Ryear'}]);	
        end
    end
    
    
    if not(isfolder(RunInfo.Dir))
        mkdir(RunInfo.Dir);
    end
    do_not_run = 0;
    no_run_msg = {'No Transmission', ...
                  'No Data', ...
                  'No Detections', ...
                  ['Data < ', num2str(MinimumData)], ...
                  'Errors in Data'};
    %%% Skip simulating locations if
    % (1) there's no local transmission (based on information from national programme)
    if do_not_run == 0 && Transmission(location) == 0
        do_not_run = 1;
    end
    
    % (2) there's no data at all (no Active Screening and no Passive Detection)
    if do_not_run == 0 && Present(location) == 0
        do_not_run = 2;
    end
    
    % (3) there's no detection at all (no Active and Passive Detections in the past)
    if do_not_run == 0 && sum(ACTIVE1(location,:)+ACTIVE2(location,:)+ACTIVENa(location,:)+PASSIVE1(location,:)+PASSIVE2(location,:)+PASSIVENa(location,:)) == 0
        do_not_run = 3;
    end
    
    % (4) data quality is poor (the total number of data in Active Scrrening and Passive Detection is less than minimum requirement)
    if do_not_run == 0 && sum(SCREENED(location, :) >= 20) + sum(PASSIVE1(location,:)+PASSIVE2(location,:)+PASSIVENa(location,:) ~= 0) < MinimumData
        do_not_run = 4;
    end
    
    % (5) there're errors in raw data (more cases than tested population)
    if ActiveNeg.Fitting == 0
        if do_not_run == 0 && sum(ACTIVE1(location,:)+ACTIVE2(location,:)+ACTIVENa(location,:) > SCREENED(location, :)) > 0
            do_not_run = 5;
        end
    else
        OnlyCols = ~ismember(1:size(SCREENED, 2), ActiveNeg.YearIdx);
        if do_not_run == 0 && sum(ACTIVE1(location,OnlyCols)+ACTIVE2(location,OnlyCols)+ACTIVENa(location,OnlyCols) > SCREENED(location, OnlyCols)) > 0
            do_not_run = 5;
        end
    end

    % if we aren't going to run the analysis for any of the above reasons:
    if do_not_run ~= 0
        for m = 1 : RunInfo.NumModel
            [Posterior, YEPHP, YEOT, PEPHP, PEOT, SampledYEPHP, SampledPEPHP, Active1, Active2, Active1FP, Active2FP,Passive1, Passive2, Deaths,  PersonYrs1, PersonYrs2, NewInf, SampledActive1, SampledActive2, SampledActive1FP, SampledActive2FP, SampledPassive1, SampledPassive2, SampledDeaths, RDT1, RDT2, RDTFP, SampledRDT1, SampledRDT2, SampledRDTFP, ReactInfo, LogModelEvidence,...
               ActiveGenderF, ActiveGenderM, ActiveAgeY, ActiveAgeW, ActiveAgeP, PassiveGenderF, PassiveGenderM, PassiveAgeY, PassiveAgeW, PassiveAgeP, PersonYrsGenderF, PersonYrsGenderM, PersonYrsAgeY, PersonYrsAgeW, PersonYrsAgeP,RDTGenderF, RDTGenderM, RDTAgeY, RDTAgeW, RDTAgeP,...
               SampledActiveGenderF, SampledActiveGenderM, SampledActiveAgeY, SampledActiveAgeW, SampledActiveAgeP, SampledPassiveGenderF, SampledPassiveGenderM, SampledPassiveAgeY,SampledPassiveAgeW, SampledPassiveAgeP, SampledRDTGenderF, SampledRDTGenderM, SampledRDTAgeY, SampledRDTAgeW, SampledRDTAgeP...
                ] = deal(no_run_msg{do_not_run});
            save(RunInfo.FitFilePath('Posterior'), 'Posterior');
            save(RunInfo.ProjFilePath('Fitted'), 'Active1', 'Active2', 'Active1FP', 'Active2FP', 'Passive1', 'Passive2', 'PersonYrs1', 'PersonYrs2', 'NewInf', ...
                'ActiveGenderF', 'ActiveGenderM', 'ActiveAgeY', 'ActiveAgeW', 'ActiveAgeP', 'PassiveGenderF', 'PassiveGenderM', 'PassiveAgeY', 'PassiveAgeW', 'PassiveAgeP', 'PersonYrsGenderF', 'PersonYrsGenderM', 'PersonYrsAgeY', 'PersonYrsAgeW', 'PersonYrsAgeP',...
                'SampledActiveGenderF', 'SampledActiveGenderM', 'SampledActiveAgeY', 'SampledActiveAgeW', 'SampledActiveAgeP', 'SampledPassiveGenderF', 'SampledPassiveGenderM', 'SampledPassiveAgeY','SampledPassiveAgeW', 'SampledPassiveAgeP',...
                                                                  'SampledActive1', 'SampledActive2', 'SampledActive1FP', 'SampledActive2FP', 'SampledPassive1', 'SampledPassive2', 'SampledDeaths');
            save(RunInfo.ProjFilePath('ReactInfo'), 'ReactInfo');
            
            for r = 0 : NumReact
                save(RunInfo.ProjFilePathR('Elimination', r), 'YEPHP', 'YEOT', 'SampledYEPHP', 'PEPHP', 'PEOT', 'SampledPEPHP');
                for s = 1 : NumStrat
                    save(RunInfo.ProjFilePathSR('Projection', s + StratMin - 1, r),  'Active1', 'Active2', 'Active1FP', 'Active2FP', 'Passive1', 'Passive2', 'Deaths', 'PersonYrs1', 'PersonYrs2', 'NewInf', ...
                                                                                  'SampledActive1', 'SampledActive2', 'SampledActive1FP', 'SampledActive2FP', 'SampledPassive1', 'SampledPassive2', 'SampledDeaths', ...
                                                                                   'ActiveGenderF', 'ActiveGenderM', 'ActiveAgeY', 'ActiveAgeW', 'ActiveAgeP', 'PassiveGenderF', 'PassiveGenderM', 'PassiveAgeY', 'PassiveAgeW', 'PassiveAgeP', 'PersonYrsGenderF', 'PersonYrsGenderM', 'PersonYrsAgeY', 'PersonYrsAgeW', 'PersonYrsAgeP',...
                                                                                     'SampledActiveGenderF', 'SampledActiveGenderM', 'SampledActiveAgeY', 'SampledActiveAgeW', 'SampledActiveAgeP', 'SampledPassiveGenderF', 'SampledPassiveGenderM', 'SampledPassiveAgeY','SampledPassiveAgeW', 'SampledPassiveAgeP', 'Deaths', 'RDT1', 'RDT2', 'RDTFP','SampledRDT1', 'SampledRDT2', 'SampledRDTFP',...
                                                                                    'RDTAgeY','RDTAgeW','RDTAgeP','RDTGenderM','RDTGenderF','SampledRDTAgeY','SampledRDTAgeW','SampledRDTAgeP','SampledRDTGenderM','SampledRDTGenderF',...
                                                                                  'RDT1', 'RDT2', 'RDTFP', 'SampledRDT1', 'SampledRDT2', 'SampledRDTFP');
                end
            end
        end
        return;
    end
    
    
          
    %%% Data processing
    N_H = PopSize(location);
    PopGrowth = PopGrowthRate(location);
    MeanPeopleScreened = round(mean(SCREENED(location, end-4 : end))); % mean of last 5 years
    MaxPeopleScreened = max(SCREENED(location,:)); % max of all period
    IntensifiedPeopleScreened = max([round(0.3 * N_H * PopGrowth ^ (YEAR(end - 2) - PopSizeYear)) MaxPeopleScreened]); % Ron
    %ModelYears = double(YEAR(1)- InterventionParameters.ScreeningBeforeData(Irow) - 1) : double(YEAR(1)-1); % simulation years
    if locIntPar.ScreeningBeforeData >= 0
        ModelYears = double(YEAR(1)- locIntPar.ScreeningBeforeData) : double(YEAR(1)-1); % pre-data simulation years
        for i = 1:length(YEAR)
            %ModelYears(i+InterventionParameters.ScreeningBeforeData(Irow)+1) =  YEAR(i);
            ModelYears(i+locIntPar.ScreeningBeforeData) =  YEAR(i);
        end
    else
        ModelYears = YEAR((1-locIntPar.ScreeningBeforeData):length(YEAR));
    end % Isangi
    
    %ScaledPeopleScreened = [0 SCREENED(location,1) * ones(1, InterventionParameters.ScreeningBeforeData(Irow)) SCREENED(location,:)] .* (InterventionParameters.PopGrowth(Irow) .^ double(PopSizeYear - ModelYears));
    ActivePosTests = ACTIVE1(location,:) + ACTIVE2(location,:) + ACTIVENa(location,:);
    ActiveNegTests = max(SCREENED(location,:) - ActivePosTests, 0);
    if ActiveNeg.Fitting ~= 0
        ActiveNeg.ObservedCases = ActivePosTests(ActiveNeg.YearIdx);
        ActiveNegTests(ActiveNeg.YearIdx) = table2array(locFittedPar(startsWith(locFittedPar.Notation,'active_neg_') & ~strcmp(locFittedPar.Distribution, 'Delta'),'Initial')); % Initial values for fitted neg test numbers.
    end
    ScreeningInfo = GetScaledScreening(ActivePosTests + ActiveNegTests, ModelYears,...
                                       locIntPar.ScreeningBeforeData,...
                                       locIntPar.ScreeningCapacity,...
                                       PopGrowth, N_H, PopSizeYear, IntrpParas);


    screening_demog_total = SCREENED_FY(location,end) + SCREENED_FW(location,end) + SCREENED_FP(location,end) + SCREENED_MY(location,end) + SCREENED_MW(location,end) + SCREENED_MP(location,end);
    float_ModelPeopleScreened_GC = ScreeningInfo.ModelPeopleScreened .* [SCREENED_FY(location,end) SCREENED_FW(location,end) SCREENED_FP(location,end) SCREENED_MY(location,end) SCREENED_MW(location,end) SCREENED_MP(location,end)]' /screening_demog_total ;

    % Makes sure we're screening integer people
    % Actually I'm not sure it even cares about this being integer
    ModelPeopleScreened_GC = floor(float_ModelPeopleScreened_GC);
    pop_remainder = float_ModelPeopleScreened_GC - ModelPeopleScreened_GC ;
    for index = 1:size(ModelPeopleScreened_GC,2)
        [~,idx] = sort(pop_remainder(:,index));
        j = 6;
        while sum(ModelPeopleScreened_GC(:,index)) - ScreeningInfo.ModelPeopleScreened(index) < 0
            ModelPeopleScreened_GC(idx(j),index) = ModelPeopleScreened_GC(idx(j),index) + 1;
            j = j - 1;
        end
    end
    
     Data = struct('N_FY', PopSize_FY(location), 'N_FW', PopSize_FW(location),'N_FP', PopSize_FP(location),'N_MY', PopSize_MY(location), 'N_MW', PopSize_MW(location),'N_MP', PopSize_MP(location), ...
                  'Years', YEAR, 'N_H', N_H, 'PopSizeYear', PopSizeYear, 'PopGrowth', PopGrowth, 'GapYear', GapYear, 'ModelYears', ModelYears,...
                  'ModelScreeningTime', ScreeningInfo.ModelScreeningTime, 'ModelScreeningFreq', ScreeningInfo.ModelScreeningFreq, 'ModelPeopleScreened', ScreeningInfo.ModelPeopleScreened,...
                  'ModelPeopleScreened_FY', ModelPeopleScreened_GC(1,:),'ModelPeopleScreened_FW', ModelPeopleScreened_GC(2,:),'ModelPeopleScreened_FP', ModelPeopleScreened_GC(3,:),'ModelPeopleScreened_MY', ModelPeopleScreened_GC(4,:),'ModelPeopleScreened_MW', ModelPeopleScreened_GC(5,:),'ModelPeopleScreened_MP', ModelPeopleScreened_GC(6,:),...
                  'PeopleScreened_FY', SCREENED_FY(location,:),'PeopleScreened_FW', SCREENED_FW(location,:),'PeopleScreened_FP', SCREENED_FP(location,:),'PeopleScreened_MY', SCREENED_MY(location,:),'PeopleScreened_MW', SCREENED_MW(location,:),'PeopleScreened_MP', SCREENED_MP(location,:),...
                  'PeopleScreened', SCREENED(location,:), 'MeanPeopleScreened', MeanPeopleScreened, 'MaxPeopleScreened', MaxPeopleScreened, 'IntensifiedPeopleScreened', IntensifiedPeopleScreened, ...
                  'ActiveDY', ACTIVEC1(location,:), 'ActiveDW', ACTIVEC2(location,:),'ActiveDP', ACTIVEC3(location,:),'ActiveDCNa', ACTIVECNa(location,:), 'ActiveDGF', ACTIVEGF(location,:),'ActiveDGM', ACTIVEGM(location,:),'ActiveDGNa', ACTIVEGNa(location,:),...
                  'ActiveD1',ACTIVE1(location,:),'ActiveD2', ACTIVE2(location,:),'ActiveDNa', ACTIVENa(location,:),...
                  'PassiveD1',PASSIVE1(location,:), 'PassiveD2',PASSIVE2(location,:),'PassiveDNa', PASSIVENa(location,:),...
                  'PassiveDY', PASSIVEC1(location,:), 'PassiveDW', PASSIVEC2(location,:),'PassiveDP', PASSIVEC3(location,:),'PassiveDCNa', PASSIVECNa(location,:), 'PassiveDGF', PASSIVEGF(location,:),'PassiveDGM', PASSIVEGM(location,:),'PassiveDGNa', PASSIVEGNa(location,:),...
                  'ActivePos', ActivePosTests, 'ActiveNeg', ActiveNegTests); % input data for main functions      

              
%     Wrapper = struct('Cloc', Cloc, 'Lv1loc', Lv1loc, 'Lv2loc', Lv2loc, 'Lv3loc', Lv3loc, 'DataStr', DataStr, 'ParaStr', ParaStr, 'MinimumData', MinimumData, ...
%                      'RunEnsProjection', RunEnsProjection, 'StratMin', StratMin, 'RunSamples', RunSamples, 'RunReactive', RunReactive);
%     Wrapper.RunCFS = RunCFS;

    if RunEnsProjection ~= 0
        EnsembleInfo = EnsembleSetup(Data, Model, RunEnsProjection, 'A');
    end
   
    DataYears = YEAR;
    for m = 1 : RunInfo.NumModel
        RunInfo.currentModelIdx = m;
        YEAR = DataYears;
        M = Model;
        %M = Model{m};
        if startsWith(M, 'Ens')
            M0 = 'Ens';
        else
            M0 = 'ens'; % contains is case sensitive
        end
        % FittedParameters for current model only
        %%% Parameter reorganisation - all parameters    
        FittedAll = locFittedPar(contains(locFittedPar.Model, {'All', M, M0}), {'Notation','Initial','Distribution','Lower','Upper'});
        
        
        % check TargetDie
        if ismember('TargetDie', FittedAll.Notation) % has TargetDie in FittedParameters for current location and current model
            if locIntPar.VCstart ~= 0
                WgtTargetDie = round(GetTargetDie(locIntPar.VC * locIntPar.VCwgt, locIntPar.TargetFreq, locIntPar.TrapCycle, TsetseParas), 4);
                FittedAll.Initial(strcmp(FittedAll.Notation, 'TargetDie')) = WgtTargetDie;
            else
                error('No historical VC therefore can not fit TargetDie.')
            end
            %if strcmp(FittedAll.Distribution(strcmp(FittedAll.Notation, 'TargetDie')), 'Delta') == 1 && locIntPar.TargetDie ~= FittedAll.Initial(strcmp(FittedAll.Notation, 'TargetDie'))
            %    error('Different TargetDie values assigned from FittedParameters and InterventionParameters')
            %elseif ismember('TargetDie', locIntPar.Properties.VariableNames) % might be removed already when running multiple models 
                MlocIntPar = removevars(locIntPar, 'TargetDie');
            %end
        else
            MlocIntPar = locIntPar;
        end % ScalingUpVC
            
        Paras = table2struct([cell2table(num2cell(FittedAll.Initial)', 'VariableNames', FittedAll.Notation'), locFixedPar, MlocIntPar]); % input parameters for main functions
        Paras.mu_H_FP = Paras.l_W*Data.N_FW/Data.N_FP;
        Paras.mu_H_MP = Paras.l_W*Data.N_MW/Data.N_MP;
        Paras.mu_H_FW = (Paras.l_Y * Data.N_FY - Paras.l_W*Data.N_FW)/Data.N_FW;
        Paras.mu_H_MW = (Paras.l_Y * Data.N_MY - Paras.l_W*Data.N_MW)/Data.N_MW;
        
        % No change in eta_H or gamma_H if d_change == 0
        Paras.eta_H_amp = Paras.eta_H_amp * (Paras.d_change ~=  0);     % set eta_H_amp   = 0 if d_change = 0
        Paras.gamma_H_amp = Paras.gamma_H_amp * (Paras.d_change ~=  0); % set gamma_H_amp = 0 if d_change = 0
        Paras.alpha0 = Paras.alpha;
        % Final year in which the active screening specifity can be lower
        % (eg for MSF screenings)
        %Paras.Last_year = locFittedPar.Last_year(strcmp(locFittedPar.Notation, 'b_specificity'));
        Paras.Last_year = locFittedPar.Last_year(strcmp(locFittedPar.Notation, 'specificityMSF'));
        Paras.ActiveNeg = ActiveNeg;
        %Paras.StratMin = StratMin;
        if exist('VCwgt') == 1
            % Save realistic target coverage to Paras
            Paras.RealisticVCwgt = VCwgt.Realistic(VCwgt.LocStr == RunInfo.LocStr); %ScalingUpVC002only
        end
        Paras.VC_t0 = mod(Paras.VCstart, 1) * 365; % Fix_f_t
%         Paras.TargetDie_full
%         Paras.TargetDie_realistic
        Paras.TargetDie_scaleback = round(GetTargetDie(Paras.VC * Paras.VCwgt_scaleback, Paras.TargetFreq_scaleback, Paras.TrapCycle, TsetseParas), 4); % ScaleBackVC
        
        %%% Parameter reorganisation - fitted parameters, extra info for MCMC
        FittedParas = locFittedPar(contains(locFittedPar.Model, {'All', M, M0}) & ~contains(locFittedPar.Distribution, {'Delta', 'Constraint'}), {'Notation', 'Initial', 'Lower', 'Upper', 'Distribution', 'Parameters', 'Initial_sigma'});
        fitted_para_names = FittedParas.Notation';
        Paras.FittedNames = fitted_para_names;
        FittedInitial = FittedParas.Initial';
        for i=1:length(fitted_para_names)
            FittedPrior.(fitted_para_names{i})={[FittedParas.Lower(i) FittedParas.Upper(i)], FittedParas.Distribution{i}, FittedParas.Parameters{i}};  %, FittedParas.Last_year{i}};
        end 
        Initial_sigma = FittedParas.Initial_sigma;
      
        % Constrained parameters
        CstrParas = FittedAll(contains(FittedAll.Distribution, 'Constraint'), {'Notation', 'Lower','Upper'});

        %%% Skip simulating locations if
        % (6) initial likelihood is insane, usually implies something is wrong in data 
        d = length(fitted_para_names(~startsWith(fitted_para_names, 'active_neg_')));	
        ProjStrat = table2struct(Strategy('Strat0',:)); % use no projection strategy 'Strat0' to run MCMC	
        if Paras.ActiveNeg.Fitting ~= 0		
            Data.chain = 0;	
            [FittedInitial, Data, AP] = sample_active_neg(d, Data, Paras, fitted_para_names, FittedInitial, FittedPrior, CstrParas, ProjStrat);	
        end
        prob = Get_log_Prob(Data, Paras, fitted_para_names, FittedInitial, FittedPrior, CstrParas, ProjStrat);
        %prob
        if prob == Inf || isnan(prob)
           % do_not_run
            [Posterior, YEPHP, YEOT, PEPHP, PEOT, SampledYEPHP, SampledPEPHP, Active1, Active2, Active1FP, Active2FP,Passive1, Passive2, Deaths,  PersonYrs1, PersonYrs2, NewInf, SampledActive1, SampledActive2, SampledActive1FP, SampledActive2FP, SampledPassive1, SampledPassive2, SampledDeaths, RDT1, RDT2, RDTFP, SampledRDT1, SampledRDT2, SampledRDTFP, ReactInfo, LogModelEvidence,...
               ActiveGenderF, ActiveGenderM, ActiveAgeY, ActiveAgeW, ActiveAgeP, PassiveGenderF, PassiveGenderM, PassiveAgeY, PassiveAgeW, PassiveAgeP, PersonYrsGenderF, PersonYrsGenderM, PersonYrsAgeY, PersonYrsAgeW, PersonYrsAgeP,RDTGenderF, RDTGenderM, RDTAgeY, RDTAgeW, RDTAgeP,...
               SampledActiveGenderF, SampledActiveGenderM, SampledActiveAgeY, SampledActiveAgeW, SampledActiveAgeP, SampledPassiveGenderF, SampledPassiveGenderM, SampledPassiveAgeY,SampledPassiveAgeW, SampledPassiveAgeP, SampledRDTGenderF, SampledRDTGenderM, SampledRDTAgeY, SampledRDTAgeW, SampledRDTAgeP] = deal(no_run_msg{do_not_run});
            save(RunInfo.FitFilePath('Posterior'), 'Posterior');
            save(RunInfo.ProjFilePath('Fitted'), 'Active1', 'Active2', 'Active1FP', 'Active2FP', 'Passive1', 'Passive2', 'PersonYrs1', 'PersonYrs2', 'NewInf', ...
                'ActiveGenderF', 'ActiveGenderM', 'ActiveAgeY', 'ActiveAgeW', 'ActiveAgeP', 'PassiveGenderF', 'PassiveGenderM', 'PassiveAgeY', 'PassiveAgeW', 'PassiveAgeP', 'PersonYrsGenderF', 'PersonYrsGenderM', 'PersonYrsAgeY', 'PersonYrsAgeW', 'PersonYrsAgeP',...
                'SampledActiveGenderF', 'SampledActiveGenderM', 'SampledActiveAgeY', 'SampledActiveAgeW', 'SampledActiveAgeP', 'SampledPassiveGenderF', 'SampledPassiveGenderM', 'SampledPassiveAgeY','SampledPassiveAgeW', 'SampledPassiveAgeP',...
                                                                  'SampledActive1', 'SampledActive2', 'SampledActive1FP', 'SampledActive2FP', 'SampledPassive1', 'SampledPassive2', 'SampledDeaths');
            save(RunInfo.ProjFilePath('ReactInfo'), 'ReactInfo');       
            for r = 0 : NumReact
                save(RunInfo.ProjFilePathR('Elimination', r), 'YEPHP', 'YEOT', 'SampledYEPHP', 'PEPHP', 'PEOT', 'SampledPEPHP');
                for s = 1 : NumStrat
                    save(RunInfo.ProjFilePathSR('Projection', s + StratMin - 1, r), 'Active1', 'Active2','Active1FP', 'Active2FP', 'Passive1', 'Passive2', 'Deaths', 'PersonYrs1', 'PersonYrs2', 'NewInf', ...
                                                                                   'SampledActive1', 'SampledActive2', 'SampledActive1FP', 'SampledActive2FP', 'SampledPassive1', 'SampledPassive2', 'SampledDeaths', ...
                                                                                   'ActiveGenderF', 'ActiveGenderM', 'ActiveAgeY', 'ActiveAgeW', 'ActiveAgeP', 'PassiveGenderF', 'PassiveGenderM', 'PassiveAgeY', 'PassiveAgeW', 'PassiveAgeP', 'PersonYrsGenderF', 'PersonYrsGenderM', 'PersonYrsAgeY', 'PersonYrsAgeW', 'PersonYrsAgeP',...
                                                                                     'SampledActiveGenderF', 'SampledActiveGenderM', 'SampledActiveAgeY', 'SampledActiveAgeW', 'SampledActiveAgeP', 'SampledPassiveGenderF', 'SampledPassiveGenderM', 'SampledPassiveAgeY','SampledPassiveAgeW', 'SampledPassiveAgeP', 'Deaths', 'RDT1', 'RDT2', 'RDTFP','SampledRDT1', 'SampledRDT2', 'SampledRDTFP',...
                                                                                    'RDTAgeY','RDTAgeW','RDTAgeP','RDTGenderM','RDTGenderF','SampledRDTAgeY','SampledRDTAgeW','SampledRDTAgeP','SampledRDTGenderM','SampledRDTGenderF',...
                                                                                  'RDT1', 'RDT2', 'RDTFP',  'SampledRDT1', 'SampledRDT2', 'SampledRDTFP');
                end
            end
            return;
        end
%         Data.FileStr = ['_', M, '_']; % for MCMC, Fitted dynamics and ReactInfo
%         FileStr = ['_', M, '_React0_'];
        
        %%% Run simulation - MCMC part
        if RunMCMC ~= 0
            tic
            HAT_MCMC_Wrapper(RunInfo, Data, Paras, fitted_para_names, FittedInitial, FittedPrior, Initial_sigma, CstrParas, ProjStrat, MachineSetting);
            time=toc;
            
            fileID = fopen('../Result/MCMCtimes','a');
            fprintf(fileID,'%s %d %d %d %s %s %g\n', Cloc, Lv1loc, Lv2loc, Lv3loc, M, ParaStr, time);
            fclose(fileID);
            
            PostProcessor(RunInfo, Data, Paras, FittedPrior, CstrParas, ProjStrat, locFittedPar, MachineSetting)
        end
        %PostProcessor(Data, Paras, FittedPrior, CstrParas, ProjStrat, locFittedPar, Wrapper, MachineSetting)
        if RunEnsProjection ~= 0
            PostProcessor(RunInfo, Data, Paras, FittedPrior, CstrParas, ProjStrat, locFittedPar, MachineSetting)
        end
        
        %%% Run importance sampling to obtain model evidence
        if RunEvidence ~= 0
            Evidence(RunInfo, RunEvidence, Data, Paras, CstrParas, FittedPrior, ProjStrat, MachineSetting);
            % for models with an animal reservoir, calculate the
            % contributions of humans and animals to R0
            %if Paras.f_A ~= 0
            %    R0Contributions(Dir, IDStr, Data, Paras, fitted_para_names, MachineSetting);
            %end
        end
        
        % Get PostID for counter-factual scenario(CFS); using the same PostID as non-CSF simulations
        clear ProjectionICs
        ProjectionICs.PostID = 0;
        if sum(strcmp(RunCFS, '0')) == 0
            load(RunInfo.ProjFilePath('ProjectionICs'));
            if sum(strcmp(RunCFS, '1')) == 1 && Paras.VCstart ~= 0
                Paras.VCstart = 0;
                Paras.TargetDie = 0;
                locStrategy.NewVCyear = NewVCyear0;
%                 if VCwgt(location) >= 0.1
%                     for i = 2 : length(locStrategy.NewVC)
%                         if locStrategy.NewVC(i) > 0
%                             WgtTargetDie = round(GetTargetDie(locStrategy.NewVC(i) * VCwgt(location), locStrategy.NewTargetFreq(i), locStrategy.NewTrapCycle(i)), 4);
%                             locStrategy.NewTargetDie(i) = WgtTargetDie;
%                         end
%                     end
%                 end
            else
                error('Can not run CFS_NoVC without historical VC')
            end
        end % ScalingUpVC


        % Merge GapYear and GapScreening into YEAR, Data.Years, Data.PeopleScreened
        ScaledPeopleScreened = [];
        GapPeopleScreened = [];
        DeltaGapYear = 0;
            
        if length(GapYear) ~= 0
            %extras = sum(strcmp(GapScreening(location, :), 'mean') + strcmp(GapScreening(location, :), 'max') == 0);
            gy = 1;
            while gy <= length(GapYear) && ~ischar(GapScreening{location, gy})
                gy = gy + 1;
            end
            extras = gy - 1;
            DeltaGapYear = length(GapYear) - extras;

            MeanPeopleScreened = round(mean([SCREENED(location, end - 4 + extras : end) GapScreening{location, 1 : extras}]));
            MaxPeopleScreened = max([SCREENED(location, :) GapScreening{location, 1 : extras}]);
            IntensifiedPeopleScreened = max([round(0.3 * N_H * PopGrowth ^ (YEAR(end - 2) + extras - PopSizeYear)) MaxPeopleScreened]);
            Data.MeanPeopleScreened = MeanPeopleScreened;
            Data.MaxPeopleScreened = MaxPeopleScreened;
            Data.IntensifiedPeopleScreened = IntensifiedPeopleScreened;
            
            YEAR = [YEAR GapYear];
            Data.Years = YEAR;
            
            for y = 1 : length(GapYear)
                switch GapScreening{location, y}
                    case 'mean'
                        GapPeopleScreened = [GapPeopleScreened MeanPeopleScreened];
                    case 'max'
                        GapPeopleScreened = [GapPeopleScreened MaxPeopleScreened];
                    otherwise % real screening numbers
                        GapPeopleScreened = [GapPeopleScreened GapScreening{location, y}];
                end
            end

            ScaledPeopleScreened = GapPeopleScreened .* PopGrowth .^ double(PopSizeYear - GapYear);
        end
        Data.DeltaGapYear = DeltaGapYear;
        

        NumPosterior = RunProjection; % number of realizations used in Projection
        samples = RunSamples; %10000 / NumPosterior; % 10,000 samples in total by drawing the same number of samples for each realization, fixed value for all strategies
        %%% Run simulation - Projection part
        if RunProjection ~= 0
            if sum(strcmp(locReactivePar.ReactiveASnum, 'AsStrat') + strcmp(locReactivePar.ReactiveASstrat, 'AsStrat')) ~= 0
                error('AsStrat is no longer an acceptable value for ReactiveASnum and ReactiveASstrat')
            end
            
            err = 0;
            for row = StratMin + 1 : size(locStrategy, 1)
                SizeASnum = length(split(locStrategy.NewASnum{row}, '_'));
                SizeASstrat = length(split(locStrategy.NewASstrat{row}, '_'));
                SizeASyear = length(locStrategy.NewASyear{row});
                if SizeASnum ~= SizeASyear || SizeASstrat ~= SizeASyear
                    warning(['Check the settings of NewASnum, NewASstrat, and NewASyear in ', locStrategy.Properties.RowNames{row}])
                    err = err + 1;
                end
            end
            if err
                return
            end
            
            if strcmp(MachineSetting.Type, 'cluster') 
                pc = parcluster;  %('local');
                pc.NumWorkers = MachineSetting.NumWorkers;
                pc.JobStorageLocation = MachineSetting.JobStorageLocation;
            end

            %NATHAN: ORIGINALLY 2000
            Pstr = datasample(1:2000, NumPosterior, 'Replace', false); % randomly select realizations from Posteior
            %Pstr = [1973 1830 102 847 724 523 400 184 440 1029];
            %ImportID = readtable('SampPostID.csv');
            %Pstr = ImportID{:,1}';
            PostID = (sum(strcmp(RunCFS, '0')) == 1) * Pstr' + (sum(strcmp(RunCFS, '0')) == 0) * ProjectionICs.PostID;
            SampPostID = reshape(repmat(PostID', samples, 1), [], 1);
            

            
            
            input = load(RunInfo.FitFilePath('Posterior'));
            Posterior = input.Posterior;
            NumYear = length(YEAR(1) : locStrategy{'Strat1', 'SIMyear'}); % fixed value for all strategies
            NumYear0 = length(YEAR);
            ElimByYears = YEAR(1) : locStrategy{'Strat1', 'SIMyear'};
            Years0 = Data.Years(1) - Paras.ScreeningBeforeData - 1 : Data.Years(end) - length(GapYear);

    
            [Active1All, Active2All, Active1FPAll, Active2FPAll,Passive1All, Passive2All, DeathsAll, PersonYrs1All, PersonYrs2All, NewInfAll, SampledActive1All, SampledActive2All, SampledActive1FPAll, SampledActive2FPAll, SampledPassive1All, SampledPassive2All, SampledDeathsAll, RDT1All, RDT2All, RDTFPAll, SampledRDT1All, SampledRDT2All, SampledRDTFPAll,...
            NewInfAll, TransmissionAll, PerfectSpecAll,    ActiveGenderFAll, ActiveGenderMAll, ActiveAgeYAll, ActiveAgeWAll, ActiveAgePAll, PassiveGenderFAll, PassiveGenderMAll, PassiveAgeYAll, PassiveAgeWAll, PassiveAgePAll, PersonYrsGenderFAll, PersonYrsGenderMAll, PersonYrsAgeYAll, PersonYrsAgeWAll, PersonYrsAgePAll,RDTGenderFAll, RDTGenderMAll, RDTAgeYAll, RDTAgeWAll, RDTAgePAll,...
            SampledActiveGenderFAll, SampledActiveGenderMAll, SampledActiveAgeYAll,SampledActiveAgeWAll, SampledActiveAgePAll, SampledPassiveGenderFAll,...
            SampledPassiveGenderMAll, SampledPassiveAgeYAll,SampledPassiveAgeWAll, SampledPassiveAgePAll, SampledRDTGenderFAll, SampledRDTGenderMAll, SampledRDTAgeYAll, SampledRDTAgeWAll, SampledRDTAgePAll]  = deal(zeros(NumPosterior * samples, NumYear, NumStrat));       
            [YEPHP, YEOT, SampledYEPHP] = deal(zeros(NumPosterior * samples, NumStrat + StratMin - 1));
            ReactInfo = [];
            if StratMin ~= 1
                load(RunInfo.ProjFilePathR('Elimination', 0));
                Pstr = PostID';
                Yzeros = zeros(NumPosterior * samples, NumStrat);
                YEPHP = [YEPHP(:, 1:StratMin-1) Yzeros];
                YEOT = [YEOT(:, 1:StratMin-1) Yzeros];
                SampledYEPHP = [SampledYEPHP(:, 1:StratMin-1) Yzeros];
                ElimByYears = ElimByYears';

                load(RunInfo.ProjFilePath('ReactInfo'));
                ReactInfo = ReactInfo{ReactInfo.Strategy < StratMin, :};

                Paras.PostRows = 0; 
            end
            [PEPHP, PEOT, SampledPEPHP] = deal(zeros(length(ElimByYears), NumStrat + StratMin - 1));

            ProjectionICs = zeros(NumPosterior, 38); % PostID + 37 state variables ('S_H'*6, 'E_H'*6 'I1_H'*6, 'I2_H'*6, 'R_H'*6,'P_V', 'S_V', 'G_V', 'E1_V', 'E2_V', 'E3_V', 'I_V') % M9
            [ProjAS, ProjPS1, ProjPS2, ProjTargetDie, ProjTargetFreq] = deal(zeros(NumPosterior, NumYear-NumYear0, NumStrat));
            ResrgYearAll = zeros(NumPosterior, length(Years0), NumStrat);
            
            if sum(strcmp(RunCFS, '2')) == 1
                if sum(ismember(fitted_para_names, {'eta_H_amp', 'gamma_H_amp'})) == 2
                    Posterior.eta_H_amp(:) = 0;
                    Posterior.gamma_H_amp(:) = 0;
                else
                    Paras.eta_H_amp = 0;
                    Paras.gamma_H_amp = 0;
                    warning('eta_H_amp and gamma_H_amp are NOT fitted parameters')
                end
            end
%             Posterior.eta_H_amp = (sum(strcmp(RunCFS, '2')) == 1) * Posterior.eta_H_amp;
%             Posterior.gamma_H_amp = (sum(strcmp(RunCFS, '2')) == 1) * Posterior.gamma_H_amp;
            
            if locStrategy{'Strat1', 'SIMyear'} <= Data.Years(end)
                locStrategy = locStrategy('Strat0', :);
                NumStrat = 0;
            end

            if strcmp(MachineSetting.Type, 'cluster') 
                parpool(pc);
            end
            clear Paraz Dataz
            for p = 1 : NumPosterior
                       
                ['Parameter', num2str(PostID(p))]; 
                Paraz(p) = Paras;
                Dataz(p) = Data;

                if StratMin ~= 1
                    Paraz(p).PostRows = samples * (p - 1) + 1 : samples * p;
                end


                % Replace values of fitted parameters in Paras by the values from Posterior
                for i = 1 : length(fitted_para_names)
                    
                    Paraz(p).(fitted_para_names{i}) = Posterior.(fitted_para_names{i})(PostID(p));
%                     Paraz(p).(fitted_para_names{i}) = Posterior{PostID(p), i};
                end

                % Update screeening here....
                if Paraz(p).ActiveNeg.Fitting ~= 0
                                       % Dataz(p).ModelScreeningTime
                    for i = 1 : Paraz(p).ActiveNeg.Fitting
                        Dataz(p).ActiveNeg(Paraz(p).ActiveNeg.YearIdx(i)) = Paraz(p).(Paraz(p).ActiveNeg.Notation{i});
                    end
                    scr = GetScaledScreening(Dataz(p).ActivePos+Dataz(p).ActiveNeg, Dataz(p).ModelYears, Paraz(p).ScreeningBeforeData,...
                                             Paraz(p).ScreeningCapacity, PopGrowth, N_H, Dataz(p).PopSizeYear, Paraz(p));

                    Dataz(p).ModelScreeningTime  = scr.ModelScreeningTime;
                    Dataz(p).ModelScreeningFreq  = scr.ModelScreeningFreq;
                    Dataz(p).ModelPeopleScreened = scr.ModelPeopleScreened;
                                       % Dataz(p).ModelScreeningTime
                    Dataz(p).PeopleScreened(Paraz(p).ActiveNeg.YearIdx) = Dataz(p).ActivePos(Paraz(p).ActiveNeg.YearIdx) +...
                                                                          Dataz(p).ActiveNeg(Paraz(p).ActiveNeg.YearIdx);
                end
                
                % Calculate Data.ModelScreeningTime,
                % Data.ModelScreeningFreq, Data.ModelPeopleScreened based on k1 in Gap years
                ExpectedScreeningTimes = ceil(ScaledPeopleScreened / N_H / Paraz(p).k1);
                ExpectedScreeningTimes(ExpectedScreeningTimes == 0) = 1;
                
                for y = 1 : length(Dataz(p).GapYear)
                    freq = 365 / ExpectedScreeningTimes(y);
                    time = 1 / ExpectedScreeningTimes(y);
                    number = ScaledPeopleScreened(y) / ExpectedScreeningTimes(y);

                    Dataz(p).ModelScreeningFreq = [Dataz(p).ModelScreeningFreq freq * ones(1, ExpectedScreeningTimes(y))];
                    Dataz(p).ModelScreeningTime = [Dataz(p).ModelScreeningTime Dataz(p).GapYear(y) + (0 : ExpectedScreeningTimes(y) - 1) * time];
                    Dataz(p).ModelPeopleScreened = [Dataz(p).ModelPeopleScreened number * ones(1, ExpectedScreeningTimes(y))];
                end
                Dataz(p).PeopleScreened = [Dataz(p).PeopleScreened GapPeopleScreened];
                
                screening_demog_total = SCREENED_FY(location,end) + SCREENED_FW(location,end) + SCREENED_FP(location,end) + SCREENED_MY(location,end) + SCREENED_MW(location,end) + SCREENED_MP(location,end);
                float_ModelPeopleScreened_GC = Dataz(p).ModelPeopleScreened .* [SCREENED_FY(location,end) SCREENED_FW(location,end) SCREENED_FP(location,end) SCREENED_MY(location,end) SCREENED_MW(location,end) SCREENED_MP(location,end)]' /screening_demog_total ;

                % Makes sure we're screening integer people
                % Actually I'm not sure it even cares about this being integer
                ModelPeopleScreened_GC = floor(float_ModelPeopleScreened_GC);
                pop_remainder = float_ModelPeopleScreened_GC - ModelPeopleScreened_GC ;
                for index = 1:size(ModelPeopleScreened_GC,2)
                    [~,idx] = sort(pop_remainder(:,index));
                    j = 6;
                    while sum(ModelPeopleScreened_GC(:,index)) - Dataz(p).ModelPeopleScreened(index) < 0
                        ModelPeopleScreened_GC(idx(j),index) = ModelPeopleScreened_GC(idx(j),index) + 1;
                        j = j - 1;
                    end
                end
                
                Dataz(p).ModelPeopleScreened_FY = ModelPeopleScreened_GC(1,:);
                Dataz(p).ModelPeopleScreened_FW = ModelPeopleScreened_GC(2,:);
                Dataz(p).ModelPeopleScreened_FP = ModelPeopleScreened_GC(3,:);
                Dataz(p).ModelPeopleScreened_MY = ModelPeopleScreened_GC(4,:);
                Dataz(p).ModelPeopleScreened_MW = ModelPeopleScreened_GC(5,:);
                Dataz(p).ModelPeopleScreened_MP = ModelPeopleScreened_GC(6,:);

                if RunDisaster == 0
                    ProjectionOutputs(p) = Projection(RunInfo, Dataz(p), Paraz(p), locStrategy, samples, locReactivePar); % single realization and all strategies
                else
                    ProjectionOutputs(p) = Disaster(RunInfo, Dataz(p), Paraz(p), locStrategy, samples, locReactivePar);
                end
            end
            delete(gcp('nocreate'));	
            clear Paraz Dataz pc;
                
            for p = 1 : NumPosterior
                Outputs = ProjectionOutputs(p);
                % Post-processing of NewInf (NewInf = 0 after EoT rather than after EoI)
% % %                 for s = 1 : NumStrat
% % %                     NewInf = Outputs.NewInf(:, :, s);
% % %                     idx = max([find(NewInf >= Paras.EOTthreshold, 1, 'last') 0]) + 1;
% % %                     Outputs.NewInf(:, idx : end, s) = 0;
% % %                 end
                rsamples = max([samples 1]);
                range = (p-1) * rsamples + 1 : p * rsamples;

                Active1All(range,:,:) = repmat(Outputs.Active1, rsamples, 1);
                Active2All(range,:,:) = repmat(Outputs.Active2, rsamples, 1);
                Active1FPAll(range,:,:) = repmat(Outputs.Active1FP, rsamples, 1);
                Active2FPAll(range,:,:) = repmat(Outputs.Active2FP, rsamples, 1);
                Passive1All(range,:,:) = repmat(Outputs.Passive1, rsamples, 1);
                Passive2All(range,:,:) = repmat(Outputs.Passive2, rsamples, 1);
                DeathsAll(range,:,:) = repmat(Outputs.Deaths, rsamples, 1);
                PersonYrs1All(range,:,:) = repmat(Outputs.PersonYrs1, rsamples, 1);
                PersonYrs2All(range,:,:) = repmat(Outputs.PersonYrs2, rsamples, 1);
                NewInfAll(range,:,:) = repmat(Outputs.NewInf, rsamples, 1);
                RDT1All(range,:,:) = repmat(Outputs.RDT1, rsamples, 1);
                RDT2All(range,:,:) = repmat(Outputs.RDT2, rsamples, 1);
                RDTFPAll(range,:,:) = repmat(Outputs.RDTFP, rsamples, 1);
                TransmissionAll(range,:,:) = repmat(Outputs.Transmission, rsamples, 1);
                PerfectSpecAll(range,:,:) = repmat(Outputs.PerfectSpec, rsamples, 1);
                YEPHP(range, StratMin:end) = repmat(Outputs.YEPHP, rsamples, 1);
                YEOT(range, StratMin:end) = repmat(Outputs.YEOT, rsamples, 1);

                ActiveGenderFAll(range,:,:) = repmat(Outputs.ActiveGenderF, rsamples, 1);
                ActiveGenderMAll(range,:,:) = repmat(Outputs.ActiveGenderM, rsamples, 1);
                ActiveAgeYAll(range,:,:) = repmat(Outputs.ActiveAgeY, rsamples, 1);
                ActiveAgeWAll(range,:,:) = repmat(Outputs.ActiveAgeW, rsamples, 1);
                ActiveAgePAll(range,:,:) = repmat(Outputs.ActiveAgeP, rsamples, 1);
                PassiveGenderFAll(range,:,:) = repmat(Outputs.PassiveGenderF, rsamples, 1);
                PassiveGenderMAll(range,:,:) = repmat(Outputs.PassiveGenderM, rsamples, 1);
                PassiveAgeYAll(range,:,:) = repmat(Outputs.PassiveAgeY, rsamples, 1);
                PassiveAgeWAll(range,:,:) = repmat(Outputs.PassiveAgeW, rsamples, 1);
                PassiveAgePAll(range,:,:) = repmat(Outputs.PassiveAgeP, rsamples, 1);
                PersonYrsGenderFAll(range,:,:) = repmat(Outputs.PersonYrsGenderF, rsamples, 1);
                PersonYrsGenderMAll(range,:,:) = repmat(Outputs.PersonYrsGenderM, rsamples, 1);
                PersonYrsAgeYAll(range,:,:) = repmat(Outputs.PersonYrsAgeY, rsamples, 1);
                PersonYrsAgeWAll(range,:,:) = repmat(Outputs.PersonYrsAgeW, rsamples, 1);
                PersonYrsAgePAll(range,:,:) = repmat(Outputs.PersonYrsAgeP, rsamples, 1);
                RDTGenderFAll(range,:,:) = repmat(Outputs.RDTGenderF, rsamples, 1);
                RDTGenderMAll(range,:,:) = repmat(Outputs.RDTGenderM, rsamples, 1);
                RDTAgeYAll(range,:,:) = repmat(Outputs.RDTAgeY, rsamples, 1);
                RDTAgeWAll(range,:,:) = repmat(Outputs.RDTAgeW, rsamples, 1);
                RDTAgePAll(range,:,:) = repmat(Outputs.RDTAgeP, rsamples, 1);


                range = (p-1) * samples + 1 : p * samples;
                SampledActive1All(range,:,:) = Outputs.SampledActive1;
                SampledActive2All(range,:,:) = Outputs.SampledActive2;
                SampledActive1FPAll(range,:,:) = Outputs.SampledActive1FP;
                SampledActive2FPAll(range,:,:) = Outputs.SampledActive2FP;
                SampledPassive1All(range,:,:) = Outputs.SampledPassive1;
                SampledPassive2All(range,:,:) = Outputs.SampledPassive2;
                SampledDeathsAll(range,:,:) = Outputs.SampledDeaths;
                SampledRDT1All(range,:,:) = Outputs.SampledRDT1;
                SampledRDT2All(range,:,:) = Outputs.SampledRDT2;
                SampledRDTFPAll(range,:,:) = Outputs.SampledRDTFP;
                SampledYEPHP(range, StratMin:end) = Outputs.SampledYEPHP;

                SampledActiveGenderFAll(range,:,:) = Outputs.SampledActiveGenderF;
                SampledActiveGenderMAll(range,:,:) = Outputs.SampledActiveGenderM;
                SampledActiveAgeYAll(range,:,:) = Outputs.SampledActiveAgeY;
                SampledActiveAgeWAll(range,:,:) = Outputs.SampledActiveAgeW;
                SampledActiveAgePAll(range,:,:) = Outputs.SampledActiveAgeP;
                SampledPassiveGenderFAll(range,:,:) = Outputs.SampledPassiveGenderF;
                SampledPassiveGenderMAll(range,:,:) = Outputs.SampledPassiveGenderM;
                SampledPassiveAgeYAll(range,:,:) = Outputs.SampledPassiveAgeY;
                SampledPassiveAgeWAll(range,:,:) = Outputs.SampledPassiveAgeW;
                SampledPassiveAgePAll(range,:,:) = Outputs.SampledPassiveAgeP;
                SampledRDTGenderFAll(range,:,:) = Outputs.SampledRDTGenderF;
                SampledRDTGenderMAll(range,:,:) = Outputs.SampledRDTGenderM;
                SampledRDTAgeYAll(range,:,:) = Outputs.SampledRDTAgeY;
                SampledRDTAgeWAll(range,:,:) = Outputs.SampledRDTAgeW;
                SampledRDTAgePAll(range,:,:) = Outputs.SampledRDTAgeP;
 
         
                
                ProjectionICs(p, :) = [PostID(p) Outputs.ProjectionICs];
                ReactInfo = [ReactInfo; [repmat(PostID(p), size(Outputs.ReactInfo, 1), 1) Outputs.ReactInfo]];
                
                if RunDisaster == 1
                    ProjAS(p,:,:) = Outputs.AS;
                    ProjPS1(p,:,:) = Outputs.PS1;
                    ProjPS2(p,:,:) = Outputs.PS2;
                    ProjTargetDie(p,:,:) = Outputs.VCTargetDie;
                    ProjTargetFreq(p,:,:) = Outputs.VCTargetFreq;
                    ResrgYearAll(p,:,:) = Outputs.ResrgYear;
                end
            end
            
            % Post-processing of NewInf (NewInf = 0 after EoT rather than after EoI)
%             NewInfAll(NewInfAll < Paras.EOTthreshold) = 0; 
% %             for s = 1 : NumStrat
% %                 NewInf = NewInfAll(:, :, s);
% %                 NoNewInfidx = table2array(rowfun(@(NewInf)(max([find(NewInf >= Paras.EOTthreshold, 1, 'last') 0])), table(NewInf))) + 1;
% %                 for row = 1 : NumPosterior * samples
% %                     NewInfAll(row, NoNewInfidx(row) : end, s) = 0;
% %                 end
% %             end
            max(locStrategy.SIMyear)
            Data.Years(end)
            if max(locStrategy.SIMyear) > Data.Years(end)
                % Calculate elimination probabilities
                for Y = ElimByYears
                    PEPHP(Y-ElimByYears(1)+1, :) = mean(YEPHP <= Y);
                    PEOT(Y-ElimByYears(1)+1, :) = mean(YEOT <= Y);
                    SampledPEPHP(Y-ElimByYears(1)+1, :) = mean(SampledYEPHP <= Y);
                end
            ElimByYears = ElimByYears';
            save(RunInfo.ProjFilePathR('Elimination', 0), 'PostID', 'SampPostID', 'ElimByYears', 'YEPHP', 'YEOT', 'SampledYEPHP', 'PEPHP', 'PEOT', 'SampledPEPHP');
       
            % Reactive infomation
            ReactInfo = array2table(ReactInfo, 'VariableNames', {'Posterior', 'Strategy', 'Reactive', 'SampleID', 'Year', 'AS', 'VC', 'RS', 'OS', 'PerfectSpec', 'Transmission', 'meff', ...
                                                              'S_H_FY', 'S_H_FW', 'S_H_FP', 'S_H_MY', 'S_H_MW', 'S_H_MP',...
                                        'E_H_FY', 'E_H_FW', 'E_H_FP', 'E_H_MY', 'E_H_MW', 'E_H_MP',...
                                    'I1_H_FY', 'I1_H_FW', 'I1_H_FP', 'I1_H_MY', 'I1_H_MW', 'I1_H_MP',...
                                    'I2_H_FY', 'I2_H_FW', 'I2_H_FP', 'I2_H_MY', 'I2_H_MW', 'I2_H_MP',...
                                    'R_H_FY', 'R_H_FW', 'R_H_FP', 'R_H_MY', 'R_H_MW', 'R_H_MP',...                          
                                'P_V', 'S_V', 'G_V', 'E1_V', 'E2_V', 'E3_V', 'I_V'});
            save(RunInfo.ProjFilePath('ReactInfo'), 'ReactInfo')
            end
             
            % Dynamics - fitted part (strategy 1)
            if StratMin == 1
                Years = YEAR;
                rsamples = max([samples 1]);
                Active1 = Active1All(1:rsamples:end, 1:NumYear0, 1);
                Active2 = Active2All(1:rsamples:end, 1:NumYear0, 1);
                Active1FP = Active1FPAll(1:rsamples:end, 1:NumYear0, 1);
                Active2FP = Active2FPAll(1:rsamples:end, 1:NumYear0, 1);
                Passive1 = Passive1All(1:rsamples:end, 1:NumYear0, 1);
                Passive2 = Passive2All(1:rsamples:end, 1:NumYear0, 1);
                Deaths = DeathsAll(1:rsamples:end, 1 : NumYear0, 1);
                PersonYrs1 = PersonYrs1All(1:rsamples:end, 1:NumYear0, 1);
                PersonYrs2 = PersonYrs2All(1:rsamples:end, 1:NumYear0, 1);
                NewInf = NewInfAll(1:rsamples:end, 1:NumYear0, 1);
                Transmission = TransmissionAll(1:rsamples:end, 1:NumYear0, 1);
                PerfectSpec = PerfectSpecAll(1:rsamples:end, 1:NumYear0, 1);

                ActiveAgeY = ActiveAgeYAll(1:rsamples:end, 1:NumYear0, 1);
                ActiveAgeW = ActiveAgeWAll(1:rsamples:end, 1:NumYear0, 1);
                ActiveAgeP = ActiveAgePAll(1:rsamples:end, 1:NumYear0, 1);
                ActiveGenderF = ActiveGenderFAll(1:rsamples:end, 1:NumYear0, 1);
                ActiveGenderM = ActiveGenderMAll(1:rsamples:end, 1:NumYear0, 1);
                PassiveAgeY = PassiveAgeYAll(1:rsamples:end, 1:NumYear0, 1);
                PassiveAgeW = PassiveAgeWAll(1:rsamples:end, 1:NumYear0, 1);
                PassiveAgeP = PassiveAgePAll(1:rsamples:end, 1:NumYear0, 1);
                PassiveGenderF = PassiveGenderFAll(1:rsamples:end, 1:NumYear0, 1);
                PassiveGenderM = PassiveGenderMAll(1:rsamples:end, 1:NumYear0, 1);
                PersonYrsAgeY = PersonYrsAgeYAll(1:rsamples:end, 1:NumYear0, 1);
                PersonYrsAgeW = PersonYrsAgeWAll(1:rsamples:end, 1:NumYear0, 1);
                PersonYrsAgeP = PersonYrsAgePAll(1:rsamples:end, 1:NumYear0, 1);
                PersonYrsGenderF = PersonYrsGenderFAll(1:rsamples:end, 1:NumYear0, 1);
                PersonYrsGenderM = PersonYrsGenderMAll(1:rsamples:end, 1:NumYear0, 1);
                
                
                SampledActive1 = SampledActive1All(:, 1:NumYear0, 1);
                SampledActive2 = SampledActive2All(:, 1:NumYear0, 1);
                SampledActive1FP = SampledActive1FPAll(:, 1:NumYear0, 1);
                SampledActive2FP = SampledActive2FPAll(:, 1:NumYear0, 1);
                SampledPassive1 = SampledPassive1All(:, 1:NumYear0, 1);
                SampledPassive2 = SampledPassive2All(:, 1:NumYear0, 1);
                SampledDeaths = SampledDeathsAll(:, 1:NumYear0, 1);

                SampledActiveAgeY = SampledActiveAgeYAll(:, 1:NumYear0, 1);
                SampledActiveAgeW = SampledActiveAgeWAll(:, 1:NumYear0, 1);
                SampledActiveAgeP = SampledActiveAgePAll(:, 1:NumYear0, 1);
                SampledActiveGenderF = SampledActiveGenderFAll(:, 1:NumYear0, 1);
                SampledActiveGenderM = SampledActiveGenderMAll(:, 1:NumYear0, 1);
                SampledPassiveAgeY = SampledPassiveAgeYAll(:, 1:NumYear0, 1);
                SampledPassiveAgeW = SampledPassiveAgeWAll(:, 1:NumYear0, 1);
                SampledPassiveAgeP = SampledPassiveAgePAll(:, 1:NumYear0, 1);
                SampledPassiveGenderF = SampledPassiveGenderFAll(:, 1:NumYear0, 1);
                SampledPassiveGenderM = SampledPassiveGenderMAll(:, 1:NumYear0, 1);
            
                AS = [Data.PeopleScreened GapPeopleScreened];
                VC = zeros(1, NumYear0);
                if Paras.VCstart > 0
                    VC(YEAR >= floor(Paras.VCstart)) = Paras.TargetFreq;
                    if Paras.VC_t0 > 0
                        VC(YEAR == floor(Paras.VCstart)) = floor((1 - Paras.VC_t0/365) * Paras.TargetFreq);
                    end
                    if Paras.VC_scaleback > 0 
                        VC(YEAR > Paras.VC_scaleback) = Paras.TargetFreq_scaleback;
                        VC(YEAR == floor(Paras.VC_scaleback)) = floor(mod(Paras.VC_scaleback, 1) *  Paras.TargetFreq + (1 - mod(Paras.VC_scaleback, 1)) * Paras.TargetFreq_scaleback);
                    end %ScaleBackVC
                end
                
                % Generate required samples of screening numbers	
                if Paras.ActiveNeg.Fitting > 0	
                    SampledScreeningYear = Paras.ActiveNeg.Year';	
                    yrIdx = Paras.ActiveNeg.YearIdx;	
                    dd = length(fitted_para_names(~startsWith(fitted_para_names, 'active_neg_')));	
                    SampledScreening = SampledActive1(:,yrIdx) + ...	
                        SampledActive2(:,yrIdx) + SampledActive1FP(:,yrIdx) + ...	
                        SampledActive2FP(:,yrIdx) + ...	
                        table2array(Posterior(SampPostID,(dd+1):(dd+Paras.ActiveNeg.Fitting)));	
                else	
                    [SampledScreening, SampledScreeningYear] = deal('No missing screening number sampling');	
                end
                save(RunInfo.ProjFilePath('Fitted'), 'PostID', 'SampPostID', 'Years', 'AS', 'VC', 'Transmission', 'PerfectSpec', ... % use Data.FileStr because of no reactive intervention in the past
                    'Active1', 'Active2', 'Active1FP', 'Active2FP', 'Passive1', 'Passive2', 'PersonYrs1', 'PersonYrs2', 'NewInf', ...
                'ActiveGenderF', 'ActiveGenderM', 'ActiveAgeY', 'ActiveAgeW', 'ActiveAgeP', 'PassiveGenderF', 'PassiveGenderM', 'PassiveAgeY', 'PassiveAgeW', 'PassiveAgeP', 'PersonYrsGenderF', 'PersonYrsGenderM', 'PersonYrsAgeY', 'PersonYrsAgeW', 'PersonYrsAgeP',...
                'SampledActiveGenderF', 'SampledActiveGenderM', 'SampledActiveAgeY', 'SampledActiveAgeW', 'SampledActiveAgeP', 'SampledPassiveGenderF', 'SampledPassiveGenderM', 'SampledPassiveAgeY','SampledPassiveAgeW', 'SampledPassiveAgeP',...
                                                                  'SampledActive1', 'SampledActive2', 'SampledActive1FP', 'SampledActive2FP', 'SampledPassive1', 'SampledPassive2', 'SampledDeaths','SampledScreening', 'SampledScreeningYear');
             
                        
                % Initial condition for projections
                ProjectionICs = array2table(ProjectionICs, 'VariableNames', {'PostID', 'S_H_FY', 'S_H_FW', 'S_H_FP', 'S_H_MY', 'S_H_MW', 'S_H_MP',...
                                        'E_H_FY', 'E_H_FW', 'E_H_FP', 'E_H_MY', 'E_H_MW', 'E_H_MP',...
                                    'I1_H_FY', 'I1_H_FW', 'I1_H_FP', 'I1_H_MY', 'I1_H_MW', 'I1_H_MP',...
                                    'I2_H_FY', 'I2_H_FW', 'I2_H_FP', 'I2_H_MY', 'I2_H_MW', 'I2_H_MP',...
                                    'R_H_FY', 'R_H_FW', 'R_H_FP', 'R_H_MY', 'R_H_MW', 'R_H_MP',...                          
                                'P_V', 'S_V', 'G_V', 'E1_V', 'E2_V', 'E3_V', 'I_V'});
                save(RunInfo.ProjFilePath('ProjectionICs'), 'ProjectionICs')
            end
             
            % Dynamics - projection part
            Years = YEAR(end) + 1 : max(locStrategy.SIMyear);
            for s = 1 : NumStrat
                S = ['Strat', num2str(s + StratMin - 1)];
                
                
                Active1 = Active1All(:, NumYear0+1:end, s);
                Active2 = Active2All(:, NumYear0+1:end, s);               
                Active1FP = Active1FPAll(:, NumYear0+1:end, s);
                Active2FP = Active2FPAll(:, NumYear0+1:end, s);
                Passive1 = Passive1All(:, NumYear0+1:end, s);
                Passive2 = Passive2All(:, NumYear0+1:end, s);
                Deaths = DeathsAll(:, NumYear0+1:end, s);
                PersonYrs1 = PersonYrs1All(:, NumYear0+1:end, s);
                PersonYrs2 = PersonYrs2All(:, NumYear0+1:end, s);
                NewInf = NewInfAll(:, NumYear0+1:end, s);
                RDT1 = RDT1All(:, NumYear0+1:end, s);
                RDT2 = RDT2All(:, NumYear0+1:end, s);
                RDTFP = RDTFPAll(:, NumYear0+1:end, s);
                Transmission = TransmissionAll(:, NumYear0+1:end, s);
                PerfectSpec = PerfectSpecAll(:, NumYear0+1:end, s);

                ActiveAgeY = ActiveAgeYAll(:, NumYear0+1:end, s);
                ActiveAgeW = ActiveAgeWAll(:, NumYear0+1:end, s);
                ActiveAgeP = ActiveAgePAll(:, NumYear0+1:end, s);
                ActiveGenderF = ActiveGenderFAll(:, NumYear0+1:end, s);
                ActiveGenderM = ActiveGenderMAll(:, NumYear0+1:end, s);
                PassiveAgeY = PassiveAgeYAll(:, NumYear0+1:end, s);
                PassiveAgeW = PassiveAgeWAll(:, NumYear0+1:end, s);
                PassiveAgeP = PassiveAgePAll(:, NumYear0+1:end, s);
                PassiveGenderF = PassiveGenderFAll(:, NumYear0+1:end, s);
                PassiveGenderM = PassiveGenderMAll(:, NumYear0+1:end, s);
                PersonYrsAgeY = PersonYrsAgeYAll(:, NumYear0+1:end, s);
                PersonYrsAgeW = PersonYrsAgeWAll(:, NumYear0+1:end, s);
                PersonYrsAgeP = PersonYrsAgePAll(:, NumYear0+1:end, s);
                PersonYrsGenderF = PersonYrsGenderFAll(:, NumYear0+1:end, s);
                PersonYrsGenderM = PersonYrsGenderMAll(:, NumYear0+1:end, s);
                RDTAgeY = RDTAgeYAll(:, NumYear0+1:end, s);
                RDTAgeW = RDTAgeWAll(:, NumYear0+1:end, s);
                RDTAgeP = RDTAgePAll(:, NumYear0+1:end, s);
                RDTGenderF = RDTGenderFAll(:, NumYear0+1:end, s);
                RDTGenderM = RDTGenderMAll(:, NumYear0+1:end, s);
                
                SampledActive1 = SampledActive1All(:, NumYear0+1:end, s);
                SampledActive2 = SampledActive2All(:, NumYear0+1:end, s);
                 SampledActive1FP = SampledActive1FPAll(:, NumYear0+1:end, s);
                SampledActive2FP = SampledActive2FPAll(:, NumYear0+1:end, s);
                SampledPassive1 = SampledPassive1All(:, NumYear0+1:end, s);
                SampledPassive2 = SampledPassive2All(:, NumYear0+1:end, s);
                 SampledRDT1 = SampledRDT1All(:, NumYear0+1:end, s);
                SampledRDT2 = SampledRDT2All(:, NumYear0+1:end, s);
                SampledRDTFP = SampledRDTFPAll(:, NumYear0+1:end, s);
                SampledDeaths = SampledDeathsAll(:, NumYear0+1:end, s);

                SampledActiveAgeY = SampledActiveAgeYAll(:, NumYear0+1:end, s);
                SampledActiveAgeW = SampledActiveAgeWAll(:, NumYear0+1:end, s);
                SampledActiveAgeP = SampledActiveAgePAll(:, NumYear0+1:end, s);
                SampledActiveGenderF = SampledActiveGenderFAll(:, NumYear0+1:end, s);
                SampledActiveGenderM = SampledActiveGenderMAll(:, NumYear0+1:end, s);
                SampledPassiveAgeY = SampledPassiveAgeYAll(:, NumYear0+1:end, s);
                SampledPassiveAgeW = SampledPassiveAgeWAll(:, NumYear0+1:end, s);
                SampledPassiveAgeP = SampledPassiveAgePAll(:, NumYear0+1:end, s);
                SampledPassiveGenderF = SampledPassiveGenderFAll(:, NumYear0+1:end, s);
                SampledPassiveGenderM = SampledPassiveGenderMAll(:, NumYear0+1:end, s);                
                SampledRDTAgeY = SampledRDTAgeYAll(:, NumYear0+1:end, s);
                SampledRDTAgeW = SampledRDTAgeWAll(:, NumYear0+1:end, s);
                SampledRDTAgeP = SampledRDTAgePAll(:, NumYear0+1:end, s);
                SampledRDTGenderF = SampledRDTGenderFAll(:, NumYear0+1:end, s);
                SampledRDTGenderM = SampledRDTGenderMAll(:, NumYear0+1:end, s);
                
                AS = zeros(1, length(Years));
                NewASnum = split(locStrategy.NewASnum{S}, '_');
                NewASyear = locStrategy.NewASyear{S};
                for y = 1 : length(NewASyear)
                    Y = find(Years == NewASyear(y));
                    switch NewASnum{y}
                        case 'mean'
                            AS(Y:end) = Data.MeanPeopleScreened;
                        case 'max'
                            AS(Y:end) = Data.MaxPeopleScreened;
                        case 'intensified'
                            AS(Y:end) = Data.IntensifiedPeopleScreened;
                        case 'off'
                            AS(Y:end) = 0;
                        case 'stop'
                            AS(Y:end) = 0;
                        otherwise 
                            if contains(NewASnum{y}, '%') % percentage 'x%'
                                pct = NewASnum{y};
                                AS(Y:end) = round(0.01 * str2num(pct(1 : end - 1)) * N_H * PopGrowth ^ (Data.Years(end - 2 - Data.DeltaGapYear) - Data.PopSizeYear));
                            else
                                AS(Y:end) = str2num(NewASnum{y});
                            end % Isangi
                    end
                end
                AS = repmat(AS, length(PostID)*max([samples 1]), 1);

%                 if strcmp(locStrategy.NewASstrat{S}, 'stop')
%                     AS(:, Years >= locStrategy.NewASyear(S)) = 0;
%                 end
                            
                VC = zeros(1, length(Years));
                VC(Years >= locStrategy{S, 'NewVCyear'}) = locStrategy{S, 'NewTargetFreq'};
                if Paras.VCstart > 0
                    if Paras.VC_scaleback == 0
                        VC(Years < locStrategy{S, 'NewVCyear'}) = Paras.TargetFreq;
                    elseif Paras.VCwgt_scaleback > 0
                        VC(Years < locStrategy{S, 'NewVCyear'}) = Paras.TargetFreq_scaleback;
                    end
                    t0 = mod(locStrategy{S, 'NewVCyear'}, 1);
                    if t0 > 0
                        VC(Years == floor(locStrategy{S, 'NewVCyear'})) = floor(t0 * Paras.TargetFreq + (1 - t0) * locStrategy{S, 'NewTargetFreq'});
                    end
                end % CompleteCycle & ScaleBackVC
                VC = repmat(VC, length(PostID)*max([samples 1]), 1);        
                save(RunInfo.ProjFilePathSR('Projection', s + StratMin - 1, 0), 'PostID', 'SampPostID', 'Years', 'AS', 'VC', 'Transmission', 'PerfectSpec', ...
                    'Active1', 'Active2','Active1FP', 'Active2FP', 'Passive1', 'Passive2', 'Deaths', 'PersonYrs1', 'PersonYrs2', 'NewInf', ...
                                                                                   'SampledActive1', 'SampledActive2',  'SampledActive1FP', 'SampledActive2FP', 'SampledPassive1', 'SampledPassive2', 'SampledDeaths', ...
                                                                                   'ActiveGenderF', 'ActiveGenderM', 'ActiveAgeY', 'ActiveAgeW', 'ActiveAgeP', 'PassiveGenderF', 'PassiveGenderM', 'PassiveAgeY', 'PassiveAgeW', 'PassiveAgeP', 'PersonYrsGenderF', 'PersonYrsGenderM', 'PersonYrsAgeY', 'PersonYrsAgeW', 'PersonYrsAgeP',...
                                                                                     'SampledActiveGenderF', 'SampledActiveGenderM', 'SampledActiveAgeY', 'SampledActiveAgeW', 'SampledActiveAgeP', 'SampledPassiveGenderF', 'SampledPassiveGenderM', 'SampledPassiveAgeY','SampledPassiveAgeW', 'SampledPassiveAgeP', 'Deaths', 'RDT1', 'RDT2', 'RDTFP','SampledRDT1', 'SampledRDT2', 'SampledRDTFP',...
                                                                                    'RDTAgeY','RDTAgeW','RDTAgeP','RDTGenderM','RDTGenderF','SampledRDTAgeY','SampledRDTAgeW','SampledRDTAgeP','SampledRDTGenderM','SampledRDTGenderF',...
                                                                                  'RDT1', 'RDT2', 'RDTFP','SampledRDT1', 'SampledRDT2', 'SampledRDTFP');
                
                % modify AS and VC 
                if RunDisaster == 1 % Need to check StratMin
%                     ProjAS(:, :, s)
%                     ProjPS(:, :, s)
%                     ProjTargetDie(:, :, s)
%                     ProjTargetFreq(:, :, s)
                    if isequal(ProjAS(1:end-1, :, s), ProjAS(2:end, :, s)) * isequal(ProjTargetDie(1:end-1, :, s), ProjTargetDie(2:end, :, s)) * isequal(ProjTargetFreq(1:end-1, :, s), ProjTargetFreq(2:end, :, s)) == 1
                        InterventionCapacity = [Years; ProjAS(1, :, s); ProjPS1(1, :, s); ProjPS2(1, :, s); ProjTargetDie(1, :, s); ProjTargetFreq(1, :, s)];
                        
                        AS = ProjAS(1, :, s);
                        VC = ProjTargetFreq(1, :, s);
                        save(RunInfo.ProjFilePathSR('Projection', s + StratMin - 1, 0), 'AS', 'VC', 'InterventionCapacity', '-append')
                        if locStrategy.CheckResrgYears(S) ~= 0 && locStrategy.DisasterEnd(S) + locStrategy.PostYears{S} - locStrategy.SIMyear(S) > 0 % simulating resurgence
                            ResrgYear = ResrgYearAll(:, :, s);
                            save(RunInfo.ProjFilePathSR('Projection', s + StratMin - 1, 0), 'ResrgYear', '-append')
                        end
                    else
                        error(['Different suspend patterens in Strategy', num2str(s + StratMin - 1)])
                    end
                end
            end
            
            
            
            
            
%             if samples == 0
%                 PlotODE(Data, 'React0')
%             end
        end
               
% % % %         %%% Run simulation - Reacive part
% % % %         if RunReactive ~= 0 && RunEnsProjection == 0
% % % %             if strcmp(MachineSetting.Type, 'cluster') 
% % % %                 pc = parcluster;  %('local');
% % % %                 pc.NumWorkers = MachineSetting.NumWorkers;
% % % %                 pc.JobStorageLocation = MachineSetting.JobStorageLocation;
% % % %             end
% % % % 
% % % %             input = load(RunInfo.ProjFilePath('ReactInfo'));
% % % %             ReactInfo = input.ReactInfo;
% % % %             %InTabs = load(['Paras_', ParaStr, '.mat']);
% % % %             for r  = 1 : NumReact
% % % %                 R = ['React', num2str(r)];
% % % %                 switch locReactivePar.ReactiveASnum{R}
% % % %                     case 'mean'
% % % %                         locReactivePar.ReactiveASnum{R} = Data.MeanPeopleScreened;
% % % %                     case 'max'
% % % %                         locReactivePar.ReactiveASnum{R} = Data.MaxPeopleScreened;
% % % %                     case 'off'
% % % %                         locReactivePar.ReactiveASnum{R} = 0;
% % % %                     otherwise
% % % %                         if ischar(locReactivePar.ReactiveASnum{R}) % percentage 'x%'
% % % %                             locReactivePar.ReactiveASnum{R} = round(0.01 * str2num(locReactivePar.ReactiveASnum{R}(1 : end - 1)) * N_H * PopGrowth ^ (Data.Years(end - 2 - Data.DeltaGapYear) - Data.PopSizeYear));
% % % %                         end
% % % %                 end
% % % % 
% % % %                 switch locReactivePar.OneoffASnum{R}
% % % %                     case 'mean'
% % % %                         locReactivePar.OneoffASnum{R} = Data.MeanPeopleScreened;
% % % %                     case 'max'
% % % %                         locReactivePar.OneoffASnum{R} = Data.MaxPeopleScreened;
% % % %                     case 'off'
% % % %                         locReactivePar.OneoffASnum{R} = 0;
% % % %                     otherwise
% % % %                         if ischar(locReactivePar.OneoffASnum{R}) % percentage 'x%'
% % % %                             locReactivePar.OneoffASnum{R} = round(0.01 * str2num(locReactivePar.OneoffASnum{R}(1 : end - 1)) * N_H * PopGrowth ^ (Data.Years(end - 2 - Data.DeltaGapYear) - Data.PopSizeYear));
% % % %                         end
% % % %                 end
% % % % 
% % % % %                 AsStrat = strcmp('AsStrat', locReactivePar{R, 'ReactiveASnum'});
% % % % %                 if AsStrat == 0
% % % % %                     locReactivePar.ReactiveASnum{R} = num2str(round(str2num(locReactivePar{R, 'ReactiveASnum'}{:}(1:2)) * N_H * Paras.PopGrowth ^ double(YEAR(end-2) - PopSizeYear) / 100));
% % % % %                 end
% % % %             end
% % % % 
% % % % %            % Save realistic target coverage to Paras
% % % % %            Paras.RealisticVCwgt = VCwgt.Realistic(VCwgt.LocStr == RunInfo.LocStr); %ScalingUpVC002only
% % % % 
% % % % 
% % % %             % Keep VC efficiency as it was if VC has started already
% % % %             % apply for all stratgies
% % % % %             VCstarted = Paras.VCstart > 0;
% % % % %             locStrategy.NewTargetDie = (VCstarted * Paras.TargetDie * (locStrategy.NewTargetDie ~= 0) .* ones(NumStrat + StratMin, 1) + (1 - VCstarted) * locStrategy.NewTargetDie) .* (locStrategy.NewTargetFreq ~= 0);
% % % % %             locStrategy.NewTargetFreq = (VCstarted * Paras.TargetFreq .* ones(NumStrat + StratMin, 1) + (1 - VCstarted) * locStrategy.NewTargetFreq) .* (locStrategy.NewTargetDie ~= 0);
% % % %             % ScalingUpVC
% % % %             %locReactivePar;
% % % % 
% % % %             %ReactOutputs(NumStrat, NumReact) = struct('PostID', [], 'SampPostID', [], 'Years', [], 'ElimByYears', [], ...
% % % %             %                                          'Active1', [], 'Active2', [], 'Passive1', [], 'Passive2', [], 'Deaths', [], 'PersonYrs1', [], 'PersonYrs2', [], 'NewInf', [], ...
% % % %             %                                          'SampledActive1', [], 'SampledActive2', [], 'SampledPassive1', [], 'SampledPassive2', [], 'SampledDeaths', [], ...
% % % %             %                                          'YEPHP', [], 'YEOT', [], 'SampledYEPHP', [], 'PEPHP', [], 'PEOT', [], 'SampledPEPHP', []);
% % % % 
% % % %             if strcmp(MachineSetting.Type, 'cluster') 
% % % %                 parpool(pc);
% % % %             end
% % % %             parfor s = 1 : NumStrat
% % % %                 %Rinfo = ReactInfo(ReactInfo.Strategy == s, :);
% % % %                 %['Strat', num2str(s)]
% % % %                 ReactOutputs(s) = Reactive(RunInfo, Data, Paras, CstrParas, locStrategy, s + StratMin - 1, locReactivePar, ReactInfo(ReactInfo.Strategy == s + StratMin - 1, :), RunCFS);            
% % % %             end
% % % %             delete(gcp('nocreate'));	
% % % %             clear pc;
% % % % 
% % % %             ElimByYears = ReactOutputs(1).(R).ElimByYears;
% % % %             PostID = ReactOutputs(1).(R).PostID;
% % % %             SampPostID = ReactOutputs(1).(R).SampPostID;
% % % %             for r = 1 : NumReact
% % % %                 R = ['React', num2str(r)];
% % % %                 [YEPHP, YEOT, SampledYEPHP] = deal(zeros(length(SampPostID), NumStrat + StratMin - 1));
% % % %                 [PEPHP, PEOT, SampledPEPHP] = deal(zeros(length(ElimByYears), NumStrat + StratMin - 1));
% % % % 
% % % %                 if StratMin ~= 1
% % % %                     input = load(RunInfo.ProjFilePathR('Elimination', r));
% % % %                     YEPHP(:, 1:StratMin) = input.YEPHP(:, 1:StratMin);
% % % %                     YEOT(:, 1:StratMin) = input.YEOT(:, 1:StratMin);
% % % %                     SampledYEPHP(:, 1:StratMin) = input.SampledYEPHP(:, 1:StratMin);
% % % %                     PEPHP(:, 1:StratMin) = input.PEPHP(:, 1:StratMin);
% % % %                     PEOT(:, 1:StratMin) = input.PEOT(:, 1:StratMin);
% % % %                     SampledPEPHP(:, 1:StratMin) = input.SampledPEPHP(:, 1:StratMin);
% % % %                 end
% % % % 
% % % %                 for s = 1 : NumStrat
% % % %                     Outputs = ReactOutputs(s).(R);
% % % % 
% % % %                     % Post-processing of NewInf (NewInf = 0 after EoT rather than after EoI)
% % % % %                     Outputs.NewInf(Outputs.NewInf < Paras.EOTthreshold) = 0;
% % % %                     NewInf = Outputs.NewInf;
% % % %                     NoNewInfidx = table2array(rowfun(@(NewInf)(max([find(NewInf >= Paras.EOTthreshold, 1, 'last') 0])), table(NewInf))) + 1;
% % % %                     for row = 1 : NumPosterior * samples
% % % %                         Outputs.NewInf(row, NoNewInfidx(row) : end) = 0;
% % % %                     end
% % % % 
% % % %                     save(RunInfo.ProjFilePathSR('Projection', s+StratMin-1, r), '-struct', 'Outputs', 'PostID', 'SampPostID', 'Years', 'AS', 'VC', 'Transmission', 'PerfectSpec', 'ReactStats', ...
% % % %                          'ActiveS', 'Active1', 'Active2', 'ActiveSFP', 'Active1FP', 'Active2FP', 'Passive1', 'Passive2', 'Deaths', 'PersonYrsS', 'PersonYrs1', 'PersonYrs2', 'NewInf', ...
% % % %                          'SampledActiveS', 'SampledActive1', 'SampledActive2', 'SampledActiveSFP', 'SampledActive1FP', 'SampledActive2FP', 'SampledPassive1', 'SampledPassive2', 'SampledDeaths', ...
% % % %                          'RDTS', 'RDT1', 'RDT2', 'RDTFP', 'SampledRDTS', 'SampledRDT1', 'SampledRDT2', 'SampledRDTFP') % M9
% % % % 
% % % %                     YEPHP(:, s+StratMin-1) = ReactOutputs(s).(R).YEPHP;
% % % %                     YEOT(:, s+StratMin-1) = ReactOutputs(s).(R).YEOT;
% % % %                     SampledYEPHP(:, s+StratMin-1) = ReactOutputs(s).(R).SampledYEPHP;
% % % %                     PEPHP(:, s+StratMin-1) = ReactOutputs(s).(R).PEPHP;
% % % %                     PEOT(:, s+StratMin-1) = ReactOutputs(s).(R).PEOT;
% % % %                     SampledPEPHP(:, s+StratMin-1) = ReactOutputs(s).(R).SampledPEPHP;
% % % %                 end
% % % %                 save(RunInfo.ProjFilePathR('Elimination', r), 'PostID', 'SampPostID', 'ElimByYears', 'YEPHP', 'YEOT', 'SampledYEPHP', 'PEPHP', 'PEOT', 'SampledPEPHP');
% % % %             end
% % % %         end
        
        %%% Plot
%         if RunPlot ~= 0
%             for r = 0 : NumReact
%                 R = ['React', num2str(r)];
%                 Plot(Data, R); 
%             end
%         end
    end
end

