%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                        %
%   This code runs generates samples of the number of negative active
%   screening results for imputation of missing numbers screened as
%   part of the Warwick HAT model,                                                                       %                                                                                %               
%                                                                                                        %
%   Inputs:                                                                                              %
%       d - the number of parameters being fitted; excluding the active_neg_* parameters                 %                              %
%       Data - structure containing historical data of the location (updated within function)            %
%       Paras - structure containing all parameters (fixed and fitted)                                   %
%       fitted_para_names - cell array containing the names of the parameters being fitted               %
%       NewP - vector of proposed parameter values (updated within function)                             %
%       FittedPrior - structure containing prior distributions for fitted parameters                     %
%       CstrParas - structure containing constrained parameters                                                 %
%       ProjStrat - structure containing parameters associated with future strategy                      %
%                                                                                                        %
%   Outputs:                                                                                             %
%       updated Data and NewP                                                                            %
%       AP - Probability of a postive active screening test (Inf indicates invalid parameter set         %
%                                                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [NewP, Data, AP] = sample_active_neg(d, Data, Paras, fitted_para_names, NewP, FittedPrior, CstrParas, ProjStrat)
    %fn = sprintf('%s_%d%s','sample_active_neg', Data.chain, '.mat');
    %save(fn, 'd', 'Data', 'Paras', 'fitted_para_names', 'NewP', 'FittedPrior', 'ProjStrat');
    
    % NATHAN: CHANGE FOR NEW ki's
    %KS = [Paras.k1WY, Paras.k1WW, Paras.k1WP, Paras.k1MY, Paras.k1MW, Paras.k1MP; ...
    %      Paras.k2WY, Paras.k2WW, Paras.k2WP, Paras.k2MY, Paras.k2MW, Paras.k2MP; ...
    %      Paras.k3WY, Paras.k3WW, Paras.k3WP, Paras.k3MY, Paras.k3MW, Paras.k3MP; ...
    %      Paras.k4WY, Paras.k4WW, Paras.k4WP, Paras.k4MY, Paras.k4MW, Paras.k4MP];

    NumFittedParas = length(fitted_para_names);
    % Replace values of fitted parameters in Paras with proposed values from sampling
    for i = 1 : NumFittedParas
        Paras.(fitted_para_names{i}) = NewP(i);
    end    
    valid_parameters = 1;
    % Check if proposed fitted parameters should be immediately rejected
    for i=1 : NumFittedParas
        if strcmp(FittedPrior.(fitted_para_names{i}){2}, 'Beta') || strcmp(FittedPrior.(fitted_para_names{i}){2}, 'Beta_shifted')
            if Paras.(fitted_para_names{i})>=FittedPrior.(fitted_para_names{i}){1}(2) || Paras.(fitted_para_names{i})<=FittedPrior.(fitted_para_names{i}){1}(1)%parameters <0
                valid_parameters = 0;
                break
            end
        elseif Paras.(fitted_para_names{i})>FittedPrior.(fitted_para_names{i}){1}(2) || Paras.(fitted_para_names{i})<FittedPrior.(fitted_para_names{i}){1}(1)%parameters <0
            valid_parameters = 0;
            break
        end
    end
    
    % Change the value of constrained parameter and reject if they sit outside of given ranges
    % NATHAN: CHANGE FOR NEW ki's
    %if size(CstrParas, 1) == 6
    %    for j = 1:6
    %        %k = [Paras.k1 Paras.k2 Paras.k3 Paras.k4];
    %        k = KS(:,j);
    %        %Paras.(CstrParas.Notation{:}) = 1 - sum(k(~strcmp(CstrParas.Notation, {'k1','k2','k3','k4'})));
    %        %if Paras.(CstrParas.Notation{:}) > CstrParas.Upper || Paras.(CstrParas.Notation{:}) < CstrParas.Lower
    %        if k(4) > 1 || k(4) < 0
    %            valid_parameters = 0;
    %        end
    %    end
    %end

    % Other constraints
    % NATHAN: CHANGE FOR NEW ki's
    %Z = max(Data.ModelPeopleScreened) / (Data.N_H * (K1 + K2)); %%%%%%%%%% NEED TO CHNAGE IF DOOR-TO-DOOR HAPPENED
    Z = max(Data.ModelPeopleScreened) / Data.N_H;
    if Z > 1   ... %screening exceeds the participating population
       || Paras.f_H + Paras.f_A > 1 ...         %or biting on reservoir hosts exceeds total biting minus biting on humans
       || Paras.gamma_H0 < (1 - Paras.u) * Paras.gamma_H ...     %or gamma_H0 cannot explain the number of deaths
       || Paras.gamma_H0 > Paras.gamma_H ...  %or we get cases hanging around.
       || Paras.gamma_H < Paras.mu_H
        valid_parameters = 0;        
    end

    if valid_parameters == 1
        NewAN = 10;
        for j=1:Paras.ActiveNeg.Fitting
%            if isnan(NewAN)
%               save('NewAN_is_nan.mat','j','p','AP','ODA','Paras','FittedPrior','Data','NumFittedParas','fitted_para_names','NewP','ProjStrat');
%            end
            [AP, ODA] = Get_Active_Prob(j, Data, Paras, ProjStrat);
            if isinf(AP)
                break
            end                
            %New proposal based on prevalence, overdispersion for current year, 
            %and expected number of people screened based on average of other years.
            p = betarnd(AP*(1/ ODA - 1), (1-AP)*(1/ ODA- 1));
            NewAN = nbinrnd(Paras.ActiveNeg.ObservedCases(j)+1,...
                            1-(1-p)*(1-FittedPrior.(Paras.ActiveNeg.Notation{j}){3}(1)));
            NewP(d+j) = NewAN;
            %NATHAN: DON'T LET OUR ESTIMATE GET TOO BIG
            if NewAN < 0.5*Data.N_H 
                Data.ActiveNeg(Paras.ActiveNeg.YearIdx(j)) = NewAN; 
            end
        end
    else
        AP = Inf;
    end        
end
