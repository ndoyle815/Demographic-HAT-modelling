%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NB: This has been adapted from the original Get_log_Prob.m function to                                                %
% incorporate a composite likelihood function including the 2000-2019                                                   %
% long-term data as well as the 2020-2022 demographic data for Guinea                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                       %
%   This code computes the (negative) Log posterior probability for UpdatedParas selected in MCMC                       %
%                                                                                                                       %
%   Inputs:                                                                                                             %
%       Data - structure containing location-specific historical data                                                   %
%       Paras - structure containing location-specific parameters (fixed, fitted and intervention parameters)           %
%       fitted_para_names - cell array containing names of fitted parameters in order they appear in FittedParameters   %
%       UpdatedParas - vector containing values of fitted parameters selected in MCMC                                   %
%       FittedPrior - structure containing prior distributions of fitted parameters                                     %
%       CstrParas - structure containing constrained parameters                                                         %
%       ProjStrat - structure containing parameters associated with future strategy                                     %
%                                                                                                                       %
%   Outputs:                                                                                                            %
%       neg_log_prob - negative log posterior probability                                                               %
%                                                                                                                       %
%   Functions required: GetEndemicEq & ODEHATmodel & log_betabinopdf & log_binopdf & GetScaledScreening                 %
%                                                                                                                       %
%   Note: future strategy is always Strat0 for fitting, which skips forcasting part                                     %
%                                                                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function neg_log_Prob = Get_log_Prob(Data, Paras, fitted_para_names, UpdatedParas, FittedPrior, CstrParas, ProjStrat)
    neg_log_Prob=0;
    NumFittedParas = length(fitted_para_names);
    
    % NATHAN: CHANGE FOR NEW ki's
    %KS = [Paras.k1WY, Paras.k1WW, Paras.k1WP, Paras.k1MY, Paras.k1MW, Paras.k1MP; ...
    %     Paras.k2WY, Paras.k2WW, Paras.k2WP, Paras.k2MY, Paras.k2MW, Paras.k2MP; ...
    %      Paras.k3WY, Paras.k3WW, Paras.k3WP, Paras.k3MY, Paras.k3MW, Paras.k3MP; ...
    %      Paras.k4WY, Paras.k4WW, Paras.k4WP, Paras.k4MY, Paras.k4MW, Paras.k4MP];

    % Replace values of fitted parameters in Paras with proposed values from sampling
    for i = 1 : NumFittedParas
        Paras.(fitted_para_names{i}) = UpdatedParas(i);
    end

    % Check if proposed fitted parameters should be immediately rejected
    for i=1 : NumFittedParas
        if strcmp(FittedPrior.(fitted_para_names{i}){2}, 'Beta') || strcmp(FittedPrior.(fitted_para_names{i}){2}, 'Beta_shifted')
            if Paras.(fitted_para_names{i})>=FittedPrior.(fitted_para_names{i}){1}(2) || Paras.(fitted_para_names{i})<=FittedPrior.(fitted_para_names{i}){1}(1)%parameters <0
                neg_log_Prob=Inf;
                break
            end
        elseif Paras.(fitted_para_names{i})>FittedPrior.(fitted_para_names{i}){1}(2) || Paras.(fitted_para_names{i})<FittedPrior.(fitted_para_names{i}){1}(1)%parameters <0
            neg_log_Prob=Inf;
            break
        end
    end

    % ELLIOT mabye add a rejection criteria up there instead
    
    % Change the value of constrained parameter and reject if they sit outside of given ranges
    % NATHAN: ONLY APPLIES TO k4 - MIGHT WANT TO ADJUST FOR ki^GC EXTENSION
    %if size(CstrParas, 1) == 1
    %    k = [Paras.k1 Paras.k2 Paras.k3 Paras.k4];
    %    Paras.(CstrParas.Notation{:}) = 1 - sum(k(~strcmp(CstrParas.Notation, {'k1','k2','k3','k4'})));
    %    if Paras.(CstrParas.Notation{:}) > CstrParas.Upper || Paras.(CstrParas.Notation{:}) < CstrParas.Lower
    %        neg_log_Prob=Inf;
    %    end
    %end
    
    % Update screening info in data
    if Paras.ActiveNeg.Fitting ~= 0
        for i = 1 : Paras.ActiveNeg.Fitting
            Data.ActiveNeg(Paras.ActiveNeg.YearIdx(i)) = Paras.(Paras.ActiveNeg.Notation{i});
        end
        scr = GetScaledScreening(Data.ActivePos+Data.ActiveNeg, Data.ModelYears, Paras.ScreeningBeforeData,...
                                 Paras.ScreeningCapacity, Data.PopGrowth, Data.N_H, Data.PopSizeYear, Paras);
        Data.ModelScreeningTime  = scr.ModelScreeningTime;
        Data.ModelScreeningFreq  = scr.ModelScreeningFreq;
        Data.ModelPeopleScreened = scr.ModelPeopleScreened;
    end
    
    % Other constraints
    % NATHAN: ADJUST TO NOT INCLUDE ki^GC
    Z = max(Data.ModelPeopleScreened) / Data.N_H; %%%%%%%%%% NEED TO CHNAGE IF DOOR-TO-DOOR HAPPENED

    %Data.ModelPeopleScreened/Data.N_H



    % ELLIOT here onwards
    % the bad call was happening from trying to put too many people
    % screened into ODEHATmodel so i have moved it up here
           % [peopleScreened(5),Data.N_H * Data.PopGrowth .^ double(Data.Years(5) - Data.PopSizeYear)]

           % [peopleScreened(10), Data.N_H * Data.PopGrowth .^ double(Data.Years(10) - Data.PopSizeYear)]
            
            % ELLIOT:
            % 1. population equilibrium was bad
            % 2. it would choose bad random pairs of numberscreened and
            % actual infections, this issue was mentioned in Run.m as
            % something that happened when prevalence is low, but they
            % never fixed it
            % 3. testing Z > 1 is a red herring in our case, Z is
            % deterministic for us I think. we wanted peopleScreened not
            % Data.ModelPeopleScreened

        % augment PeopleScreened with missing values
        peopleScreened = Data.PeopleScreened;
        
        if Paras.ActiveNeg.Fitting ~= 0
            peopleScreened(Paras.ActiveNeg.YearIdx) = Data.ActivePos(Paras.ActiveNeg.YearIdx) + Data.ActiveNeg(Paras.ActiveNeg.YearIdx);
        end

    if Z > 1 ...    %screening exceeds the participating population
       || Paras.f_H + Paras.f_A > 1 ...         %or biting on reservoir hosts exceeds total biting minus biting on humans
       || Paras.gamma_H0 < (1 - Paras.u) * Paras.gamma_H ...     %or gamma_H0 cannot explain the number of deaths
       || Paras.gamma_H0 > Paras.gamma_H ...    %or health system got worse around 2000
       || Paras.gamma_H0 < Paras.mu_H ...       %or we get cases hanging around.
       || Paras.specificityMSF > Paras.specificity...   % MSF spec > local programme spec
        || length(Data.ModelPeopleScreened) ~= length(Data.Years) + 3
        neg_log_Prob=Inf; 
    end

    
    
    if neg_log_Prob ~= Inf
	    % Get ICs	
        [meff, ICs] = GetEndemicEq(Data.N_H, Paras, Data); % NATHAN: NEW GetEndemicEq() 	
        
        % Run Model from IC	
        Paras.death = (1-Paras.u) * Paras.gamma_H;	
        Paras.effcy = 1;	
        [Classes, Aggregate] = ODEHATmodel(meff, ICs, Data, Paras, ProjStrat); % NATHAN: NEW ODEHATmodel()


        % ELLIOT my check
        these_cases = Aggregate.ActiveM1' + Aggregate.ActiveM2' + Aggregate.ActiveMFP';        
        if any(these_cases(Paras.ActiveNeg.YearIdx) - peopleScreened(Paras.ActiveNeg.YearIdx) > 0)  
            % note I don't think it should ever actually hit this case but
            % it would be annoying to check
            neg_log_prob = Inf;
        else


        % Compute log prob
        % Weights from prior (_negative_ log prior probability)
        for i = 1 : NumFittedParas
            ParaInfo = FittedPrior.(fitted_para_names{i}); % {[Lower, Upper], Distribution, Parameters}
            if strcmp(ParaInfo{2}, 'Gamma_shifted')
                k = ParaInfo{3}(1);
                theta = ParaInfo{3}(2);
                shift = ParaInfo{1}(1);
                neg_log_Prob_gamma = -(k-1) * log(Paras.(fitted_para_names{i}) - shift) + (Paras.(fitted_para_names{i}) - shift) / theta + k * log(theta) + gammaln(k);
                neg_log_Prob = neg_log_Prob + neg_log_Prob_gamma;
                
            elseif strcmp(ParaInfo{2}, 'Beta')
                a=ParaInfo{3}(1);
                b=ParaInfo{3}(2);
                neg_log_Prob_beta = -(a-1) * log(Paras.(fitted_para_names{i})) - (b-1) * log(1 - Paras.(fitted_para_names{i})) + betaln(a,b);
                neg_log_Prob = neg_log_Prob + neg_log_Prob_beta;
            
            elseif strcmp(ParaInfo{2}, 'Beta_shifted') %shifted and scaled prior (eg for specificity)
                a=ParaInfo{3}(1);
                b=ParaInfo{3}(2);
                lb = ParaInfo{1}(1);
                ub = ParaInfo{1}(2);
                para01 = (Paras.(fitted_para_names{i}) - lb) / (ub - lb);
                neg_log_Prob_beta = -(a-1) * log(para01) - (b-1) * log(1 - para01) + betaln(a,b) - log(ub - lb);
                neg_log_Prob = neg_log_Prob + neg_log_Prob_beta;

            elseif strcmp(ParaInfo{2}, 'Exp_shifted') %shifted exponential prior (eg for R0, starting at 1)
                lambda = ParaInfo{3}(1);
                lb = ParaInfo{1}(1);
                neg_log_Prob_exp =  -log(lambda) + (lambda * (Paras.(fitted_para_names{i}) - lb));
                neg_log_Prob = neg_log_Prob + neg_log_Prob_exp;
                
            elseif strcmp(ParaInfo{2}, 'NegBinomial') %negative binomial prior (eg for number of negative active tests when screening number is absent)
                r = ParaInfo{3}(1);
                p = ParaInfo{3}(2);
                f = Paras.(fitted_para_names{i});
                neg_log_Prob_nbin = -gammaln(r + f) + gammaln(r) + gammaln(f + 1) - (r * log(p)) - (f * log(1 - p));
                neg_log_Prob = neg_log_Prob + neg_log_Prob_nbin;


            elseif strcmp(ParaInfo{2}, 'Geometric') %Geometric prior (eg for number of negative active tests when screening number is absent)                lambda = ParaInfo{3}(1);
                f = Paras.(fitted_para_names{i});
                lambda = ParaInfo{3}(1);
                neg_log_Prob_geom = -f*log(1-lambda)-log(lambda);
                neg_log_Prob = neg_log_Prob + neg_log_Prob_geom;
                % ELLIOT this line is the issue, its the one commented on
                % in Run.m about WtMeanAS but I dont really know what to do
                % about it
                % just make small f a rejection criteria i guess. f =
                % 200ish

            %elseif strcmp(ParaInfo{2}, 'Normal')
            %    mu = ParaInfo{3}(1);
            %    sigma = ParaInfo{3}(2);
            %    neg_log_Prob_normal = (Paras.(fitted_para_names{i}) - mu)^2 / (2 * sigma^2) + log(sigma) + log(2*pi) / 2;
            %    neg_log_Prob = neg_log_Prob + neg_log_Prob_normal;
            end
        end

        

        % Weight from data
        PassiveD = Data.PassiveD1 + Data.PassiveD2 + Data.PassiveDNa;
        if Paras.Last_year ~= 0
            % Using Last_year as a flag for the MSF active HZ in Isangi
            % in which we assume the stage 1 passive detections are FP
            PassiveD_S1FP = Data.PassiveD2 + Data.PassiveDNa;
        end
        ActiveD =  Data.ActiveD1 + Data.ActiveD2 + Data.ActiveDNa;
        AbsPop = Data.N_H * Data.PopGrowth .^ double(Data.Years - Data.PopSizeYear);
        PassiveProb = (Aggregate.PassiveM1' + Aggregate.PassiveM2') / Data.N_H;
        %(Aggregate.ActiveM1' + Aggregate.ActiveM2' + Aggregate.ActiveMFP')
        %peopleScreened
        ActiveProb = (Aggregate.ActiveM1' + Aggregate.ActiveM2' + Aggregate.ActiveMFP') .* Data.PopGrowth .^ double(Data.Years - Data.PopSizeYear) ./ peopleScreened;        
        PassiveS1Prob = Aggregate.PassiveM1' ./ (Aggregate.PassiveM1' + Aggregate.PassiveM2');
        %ActiveS1Prob = Aggregate.ActiveM1' ./ (Aggregate.ActiveM1' + Aggregate.ActiveM2');
        
        %Tweek for FP
        ProbS1givenTP = Aggregate.ActiveM1' ./ (Aggregate.ActiveM1' + Aggregate.ActiveM2');
        ProbS1givenFP = Paras.S1givenFP;
        ProbTP = (Aggregate.ActiveM1' + Aggregate.ActiveM2') ./ (Aggregate.ActiveM1' + Aggregate.ActiveM2' + Aggregate.ActiveMFP');
        ProbTP(ProbTP < 10^(-5)) = 0;
        
        ActiveS1Prob = ProbS1givenTP.*ProbTP + ProbS1givenFP.*(1-ProbTP);
        Disp_act = Paras.disp_act*ProbTP ; %need to discuss disp_act is currenlty 1e-3; 
        
        firstscreening = find([peopleScreened 1] >= 20, 1); % changed from ~= 0, taking 20 as cut-off for "real" screening
        firstactive = find([ActiveD 1] ~= 0, 1); 
        firstpassive = find([PassiveD 1] ~= 0, 1);
          
        for y = min([firstscreening, firstactive, firstpassive]) : length(Data.Years)
            % Passive               
            if Paras.Last_year ~= 0
                passive_ll = log_betabinopdf(PassiveD_S1FP(y), AbsPop(y), PassiveProb(y), Paras.disp_pass); 
            else
                %disp(['passive d ', num2str(PassiveD(y))])
                %disp(['scaled nh ', num2str(AbsPop(y))])
                %disp(['passive p ', num2str(PassiveProb(y))])
                passive_ll = log_betabinopdf(PassiveD(y), AbsPop(y), PassiveProb(y), Paras.disp_pass)... 
                           + log_binopdf(Data.PassiveD1(y), Data.PassiveD1(y) + Data.PassiveD2(y), PassiveS1Prob(y));
            end

            % Active
            %Tweek for FP in overdispersion
            active_ll = log_betabinopdf(ActiveD(y), peopleScreened(y), ActiveProb(y), Disp_act(y)) ...
                      + log_binopdf(Data.ActiveD1(y), Data.ActiveD1(y) + Data.ActiveD2(y), ActiveS1Prob(y));

            neg_log_Prob = neg_log_Prob - active_ll - passive_ll;
                        
        end
        
        % NATHAN: NEW COMPUTATIONS
        if Aggregate.PassiveMGenderM' + Aggregate.PassiveMGenderF' > 0
            PassiveGMProb = Aggregate.PassiveMGenderM' ./ (Aggregate.PassiveMGenderM' + Aggregate.PassiveMGenderF');
        else
            PassiveGMProb = zeros(1,23);
        end
        if Aggregate.PassiveMAgeY' + Aggregate.PassiveMAgeW' + Aggregate.PassiveMAgeP' > 0
            PassiveC1Prob = Aggregate.PassiveMAgeY' ./ (Aggregate.PassiveMAgeY' + Aggregate.PassiveMAgeW' + Aggregate.PassiveMAgeP');
            PassiveC2Prob = Aggregate.PassiveMAgeW' ./ (Aggregate.PassiveMAgeY' + Aggregate.PassiveMAgeW' + Aggregate.PassiveMAgeP');
        else
            PassiveC1Prob = zeros(1,23);
            PassiveC2Prob = zeros(1,23);
        end

        % NB: ASSUMING "ActiveMFP" IS DISTINCT FROM Active1,2 AND HENCE IS
        % DISTINCT FROM ActiveG/C
        if Aggregate.ActiveMGenderM' + Aggregate.ActiveMGenderF' == 0
            ActiveGMProb = zeros(1,23);
        else
            ProbGMgivenTP = Aggregate.ActiveMGenderM' ./ (Aggregate.ActiveMGenderM' + Aggregate.ActiveMGenderF');
            ProbGMgivenFP = 0.5;
            ProbTP = (Aggregate.ActiveMGenderM' + Aggregate.ActiveMGenderF') ./ (Aggregate.ActiveMGenderM' + Aggregate.ActiveMGenderF' + Aggregate.ActiveMFP');
            ProbTP(ProbTP < 10^(-5)) = 0;
            ActiveGMProb = ProbGMgivenTP.*ProbTP + ProbGMgivenFP.*(1-ProbTP);
        end

        if Aggregate.ActiveMAgeY' + Aggregate.ActiveMAgeW' + Aggregate.ActiveMAgeP' == 0
            ActiveC1Prob = zeros(1,23);
            ActiveC2Prob = zeros(1,23);
        else    
            ProbC1givenTP = Aggregate.ActiveMAgeY' ./ (Aggregate.ActiveMAgeY' + Aggregate.ActiveMAgeW' + Aggregate.ActiveMAgeP');
            ProbC1givenFP = 0.216; % CRUDE APPROXIMATIONS BASED OFF CASEDEMO
            ProbTP = (Aggregate.ActiveMAgeY' + Aggregate.ActiveMAgeW' + Aggregate.ActiveMAgeP') ./ (Aggregate.ActiveMAgeY' + Aggregate.ActiveMAgeW' + Aggregate.ActiveMAgeP' + Aggregate.ActiveMFP');
            ProbTP(ProbTP < 10^(-5)) = 0;
            ActiveC1Prob = ProbC1givenTP.*ProbTP + ProbC1givenFP.*(1-ProbTP);

            ProbC2givenTP = Aggregate.ActiveMAgeW' ./ (Aggregate.ActiveMAgeY' + Aggregate.ActiveMAgeW' + Aggregate.ActiveMAgeP');
            ProbC2givenFP = 0.763; % CRUDE APPROXIMATIONS BASED OFF CASEDEMO
            ActiveC2Prob = ProbC2givenTP.*ProbTP + ProbC2givenFP.*(1-ProbTP);
        end

        % NATHAN: HERE IS WHERE WE APPEND THE LOG-LIKELIHOOD FOR AGE AND
        % GENDER DISTRIBUTION CHECKS
       % neg_log_Prob
        for y = 21:23
%             y
%             Data.PassiveDGM
%             Data.PassiveDGF
%             PassiveGMProb
            % PASSIVE
            passive_ll = log_binopdf(Data.PassiveDGM(y), Data.PassiveDGM(y) + Data.PassiveDGF(y), PassiveGMProb(y)) ...
                       + log_binopdf(Data.PassiveDY(y), Data.PassiveDY(y) + Data.PassiveDW(y) + Data.PassiveDP(y), PassiveC1Prob(y)) ...
                       + log_binopdf(Data.PassiveDW(y), Data.PassiveDY(y) + Data.PassiveDW(y) + Data.PassiveDP(y), PassiveC2Prob(y));
        %    passive_ll
            % ACTIVE
            active_ll = log_binopdf(Data.ActiveDGM(y), Data.ActiveDGM(y) + Data.ActiveDGF(y), ActiveGMProb(y)) ...
                      + log_binopdf(Data.ActiveDY(y), Data.ActiveDY(y) + Data.ActiveDW(y) + Data.ActiveDP(y), ActiveC1Prob(y)) ...
                      + log_binopdf(Data.ActiveDW(y), Data.ActiveDY(y) + Data.ActiveDW(y) + Data.ActiveDP(y), ActiveC2Prob(y));
        %    active_ll
            neg_log_Prob = neg_log_Prob - active_ll - passive_ll;

        end
       % neg_log_Prob
        if neg_log_Prob < 0
            errfile = ['Get_log_Prob_input-dt',datestr(datetime,'yyyymmdd_HHMMSS'),'-',num2str(random('poisson',100)),'.mat'];
            save(errfile,'neg_log_Prob','Data', 'Paras', 'fitted_para_names', 'UpdatedParas', 'FittedPrior', 'CstrParas', 'ProjStrat');
            error(['Get_log_Prob: neg_log_Prob < 0. See ',errfile])
        end
        end% ELLIOT this is my end
    end
