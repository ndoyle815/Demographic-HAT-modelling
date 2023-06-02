%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                        %
%   This code is a wrapper for the MCMC for the Warwick HAT model, it calls HAT_Auto_MCMC                %
%                                                                                                        %
%   Inputs:                                                                                              %
%       RunInfo - includes RunMCMC: Fresh run (1) or re-adapt (2) or additional samples (3)                                %
%       Data - structure containing historical data of the location                                      %
%       Model - Name of the model being fitted                                                           %
%       fitted_para_names - cell array containing the names of the parameters being fitted               %
%       Para - structure containing all parameters (fixed and fitted)                                    %
%       Intervention - structure containing parameters associated with intervention                      %
%       ProjStrat - structure containing parameters associated with future strategy                      %
%       FittedInitial - vector containing initial values for fitted parameters                           %
%       FittedPrior - structure containing prior distributions for fitted parameters                     %
%       Initial_sigma - structure containing variances used for proposal during initial/learning period  %
%       ParaStr - character string ID for the parameter input                                            %
%                                                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function HAT_MCMC_Wrapper(RunInfo, Data, Paras, fitted_para_names, FittedInitial, FittedPrior, Initial_sigma, CstrParas, ProjStrat, MachineSetting)
% Take the MCMC settings from RunInfo
MCMCSettings = RunInfo.MCMCSettings;
%  ... and add Initial_sigma
MCMCSettings.Initial_sigma = Initial_sigma;

% filename for the MCMC output
filename = RunInfo.FitFilePath('Output');

Best=zeros(MCMCSettings.n_chains,length(fitted_para_names));
if RunInfo.RunMCMC == 1 % from scratch "full" run
    %randomly perturb the initial values of the fitted parameters and check that they result in a finite log-posterior
    log_prob=Inf(1,MCMCSettings.n_chains);
    d=length(fitted_para_names(~startsWith(fitted_para_names, 'active_neg_')));
    AP = 1;
    % Find closest constraint to supplied initial value
    nearest_constraint = zeros(d,1);
    for i=1:d
        [~, idx] = min(abs(FittedPrior.(fitted_para_names{i}){1}-FittedInitial(i)));
        nearest_constraint(i) = FittedPrior.(fitted_para_names{i}){1}(idx);
    end
    n_tries = 0;
    while (sum(isinf(log_prob)+isnan(log_prob))>0)
        n_tries = n_tries+1;
        if (n_tries > 1000)
            error('HATMEPP error: failed to obtain valid initial values for log_prob')
        end
        u=exp(0.2*rand(MCMCSettings.n_chains,length(fitted_para_names))-0.1);
        for j=1:d    %length(fitted_para_names)
            Best(:,j)=nearest_constraint(j)+((FittedInitial(j)-nearest_constraint(j))*u(:,j)); %[u(j);1/u(j)]);
        end
        for i=1:MCMCSettings.n_chains        
            if Paras.ActiveNeg.Fitting ~= 0
                Data.chain = i;
                [Best(i,:), Data, AP] = sample_active_neg(d, Data, Paras, fitted_para_names, Best(i,:), FittedPrior, CstrParas, ProjStrat);
            end
            if ~isinf(AP)
                log_prob(i) = Get_log_Prob(Data, Paras, fitted_para_names, Best(i,:), FittedPrior, CstrParas, ProjStrat);
            else
                log_prob(i) = Inf;
            end
        end
    end
else  % restarting from last posterior sample from each chain to date
    load(filename, 'MCMCAdapt');
    d=length(fitted_para_names(~startsWith(fitted_para_names, 'active_neg_')));
    for i=1:MCMCSettings.n_chains
        Best(i,1:d)=MCMCAdapt.Mu;
    end
    if Paras.ActiveNeg.Fitting ~= 0
        Best(:,(d+1):(d+Paras.ActiveNeg.Fitting)) = ones(MCMCSettings.n_chains,Paras.ActiveNeg.Fitting)*10000;
    end
end

% run the MCMC
[burnin, thin_factor, final_ess, final_conv] = HAT_Auto_MCMC(Best, MCMCSettings, Data, Paras, fitted_para_names, FittedPrior, CstrParas, ProjStrat, filename, MachineSetting);

% save the posterior in its own file
%%%%%%%% changes were made here
load(filename, 'posterior', 'neg_log_likelihood', 'Chains');
Posterior = [array2table(posterior, 'VariableNames', fitted_para_names) array2table(neg_log_likelihood, 'VariableNames', {'NegLogLikelihood'})];
save(RunInfo.FitFilePath('Posterior'), 'Posterior', 'Chains');

% output MCMC diagnostics
fileID = fopen(RunInfo.FitFilePath('Diagnostics', ''), 'w');
fprintf(fileID,'burnin = %g\n',burnin);
fprintf(fileID,'%thin factor = %g\n',thin_factor);
fprintf(fileID,'effective sample size = %g\n',final_ess);
fprintf(fileID,'convergence diagnostic = %g\n',final_conv);
fclose(fileID);

% create a report if the final effective sample size is below a quarter of the total sample size and the final convergence diagnostic is above the threshold
if (final_ess<MCMCSettings.sample_size/2 || final_conv>MCMCSettings.convergence_threshold)
    RepDir = '../Result/MCMC_Reports/';
    if not(isfolder(RepDir))
        mkdir(RepDir);
    end
    fileID = fopen([RepDir RunInfo.FitFileName('Report', '.txt')],'w');
    fprintf(fileID,'effective sample size = %g\n',final_ess);
    fprintf(fileID,'convergence diagnostic = %g\n',final_conv);
    fprintf(fileID,'CHECK MCMC OUTPUT\n');
end

end
