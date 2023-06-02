%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                        %
%   This code automates the burnin and run length for the MCMC for the Warwick HAT model,                % 
%   it calls HAT_Adaptive_MCMC                                                                           %               
%                                                                                                        %
%   Inputs:                                                                                              %
%       Best - two row array with the initial values of the fitted parameters for each chain             %
%       MCMCSettings - structure containing all the settings for the MCMC                                %
%       Data - structure containing historical data of the location                                      %
%       Paras - structure containing all parameters (fixed and fitted)                                   %
%       Intervention - structure containing parameters associated with intervention                      %
%       fitted_para_names - cell array containing the names of the parameters being fitted               %
%       FittedPrior - structure containing prior distributions for fitted parameters                     %
%       ProjStrat - structure containing parameters associated with future strategy                      %
%       filename - character string of the filename for saving output                                    %
%       MachineSetting - object of class parcluster usig in parallelisation                                  %
%                                                                                                        %
%   Outputs (all single values):                                                                         %
%       burnin - the length of the burnin                                                                %
%       thin_factor - the amount the chains have been thinned                                            %
%       final_ess - the effective sample size of the posterior sample                                    %
%       final_conv - the (between-chain) convergence diagnostic of the posterior sample                  %
%                                                                                                        %
%   MCMCSettings includes run_type, which can take values:                                               %
%       1 - run MCMC from scratch                                                                        %
%            - single parameter updates, adaptation and sampling                                         %
%       2 - restart MCMC                                                                                 %
%            - adaptation and sampling                                                                   %
%            - starting with a merged chain covariance matrix                                            %
%            - increased max_thin                                                                        %
%       3 - restart MCMC                                                                                 %
%            - additional sampling with increased max_thin                                               %
%                                                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [burnin,thin_factor,final_ess,final_conv]=HAT_Auto_MCMC(Best,MCMCSettings,Data,Paras,fitted_para_names,FittedPrior,CstrParas,ProjStrat,filename,MachineSetting)

% d is the number of parameters being fitted *adaptively*
d=length(fitted_para_names(~startsWith(fitted_para_names,'active_neg_')));

%initialise the mean and variance of the chain (these are updated iteratively)
MCMCAdapt = repmat(struct('Sigma',zeros(d),'Mu',zeros(1,d),'lambda',1,...
    'BlockUpdates',0,'Acceptance',0),1,MCMCSettings.n_chains);
if MCMCSettings.run_type == 1
    for j=1:MCMCSettings.n_chains
        MCMCAdapt(j).Sigma=diag(MCMCSettings.Initial_sigma(1:d));
        MCMCAdapt(j).Mu=Best(j,1:d);
        MCMCAdapt(j).lambda=1;
        MCMCAdapt(j).BlockUpdates=0;
        MCMCAdapt(j).Acceptance=0;
    end
    thin_factor=0;
elseif MCMCSettings.run_type == 2
    input = load(filename,'MCMCAdapt','parameters','log_post');
    MCMCAdapt = input.MCMCAdapt;
    parameters = input.parameters;
    log_post = input.log_post;
    sizep = size(parameters(MCMCSettings.initial_period+1:end,:,1));
    workp = zeros(sizep(1)*MCMCSettings.n_chains,sizep(2));
    for j = 1:MCMCSettings.n_chains
        workp((((j-1)*sizep(1))+1):(j*sizep(1)),:) = parameters(MCMCSettings.initial_period+1:end,:,j);
    end
    across_chain_sigma = cov(workp);
    clear workp sizep
    for j=1:MCMCSettings.n_chains
        MCMCAdapt(j).Sigma=across_chain_sigma(1:d,1:d);
        MCMCAdapt(j).Mu=Best(j,1:d);
        MCMCAdapt(j).lambda=1;
        MCMCAdapt(j).BlockUpdates=0;
        MCMCAdapt(j).Acceptance=0;
    end
    iters1=MCMCSettings.initial_period+MCMCSettings.learning_period+MCMCSettings.sample_size+round(MCMCSettings.sample_size/2);
    parameters=parameters([1:iters1,end],:,:);
    log_post=log_post([1:iters1,end],:,:);
    thin_factor=0;
    MCMCSettings.max_thin=MCMCSettings.max_thin*1.5;
    MCMCSettings.n0=MCMCSettings.n0*5; % slow down movement of cov matrix to empirical
    MCMCSettings.m0=MCMCSettings.m0*2; % slow down convergence of lambda
elseif MCMCSettings.run_type == 3
    input = load(filename,'MCMCAdapt','parameters','log_post','thin_factor','burnin');
    MCMCAdapt = input.MCMCAdapt;
    parameters = input.parameters;
    log_post = input.log_post;
    thin_factor = input.thin_factor;
    burnin = input.burnin;
    MCMCSettings.max_thin=max(MCMCSettings.max_thin,thin_factor)*2.5;
end


%first part is to run the chains for the initial period, then the learning
%period and then long enough to start checking the burn-in
fprintf('MCMC starts running at %s\n', datetime('now','Format','d-MMM-y HH:mm:ss'))

iterations=MCMCSettings.initial_period+MCMCSettings.learning_period+MCMCSettings.sample_size+round(MCMCSettings.sample_size/2);

if strcmp(MachineSetting.Type, 'cluster') 
    pc = parcluster;  %('local');
    pc.NumWorkers = MachineSetting.NumWorkers;
    pc.JobStorageLocation = MachineSetting.JobStorageLocation;
end
if MCMCSettings.run_type == 1
    fprintf('MCMC initial and learning period starts running at %s\n', datetime('now','Format','d-MMM-y HH:mm:ss'))
    from_iteration=0;
    if strcmp(MachineSetting.Type, 'cluster') 
        parpool(pc);
    end
    Datap = Data;
    Datap.chain = 0;
    Dataz = repmat(Datap,1,MCMCSettings.n_chains);
    for j=1:MCMCSettings.n_chains
        %Dataz(j) = Datap;
        Dataz(j).chain=j;
    end
    clear Datap
    parfor j=1:MCMCSettings.n_chains
        rng(sum(j*clock));
        [parameters(:,:,j),log_post(:,j),MCMCAdapt(j)]=HAT_Adaptive_MCMC(from_iteration,iterations,Best(j,:),[],MCMCSettings,MCMCAdapt(j),Dataz(j),Paras,fitted_para_names,FittedPrior,CstrParas,ProjStrat);
    end
    delete(gcp('nocreate'));
    save(filename,'parameters','log_post','MCMCAdapt','MCMCSettings','fitted_para_names');
    fprintf([num2str(iterations), ' MCMC iterations have been done in the initial and the learning period \n'])
end

if MCMCSettings.run_type <= 2
    %second part is to check the burnin, and then run an additional 100
    %iterations while it does not pass the diagnostics test
    fprintf('MCMC diagnostic test starts running at %s\n', datetime('now','Format','d-MMM-y HH:mm:ss'))
    burnin=MCMCSettings.initial_period+MCMCSettings.learning_period+round(MCMCSettings.sample_size/2)-100;
    burned_in=0;
    count = 0;
    if strcmp(MachineSetting.Type, 'cluster')
        parpool(pc);
    end
    while (burned_in==0 && burnin<MCMCSettings.max_burnin) 
    
        burnin=burnin+100;
        y=zeros(d,MCMCSettings.n_chains);
        for i=1:d
            for j=1:MCMCSettings.n_chains
                %y(i,j) is the within-chain convergence diagnostic for parameter
                %i and chain j
                %y(i,j)=convergence_i(parameters(burnin-round(MCMCSettings.sample_size/2)+1:burnin-round(MCMCSettings.sample_size/2)+MCMCSettings.sample_size,i,j),parameters(burnin+1:burnin+MCMCSettings.sample_size,i,j));
                y(i,j)=convergence_i([parameters(burnin-round(MCMCSettings.sample_size/2)+1:burnin+round(MCMCSettings.sample_size/4),i,j),...
                                      parameters(burnin+round(MCMCSettings.sample_size/4)+1:burnin+MCMCSettings.sample_size,i,j)]);
            end
        end
        yy=zeros(1,d);
        samples_i = zeros(MCMCSettings.sample_size,MCMCSettings.n_chains);
        for i=1:d
            %yy(i) is the between-chain convergence diagnostic for parameter i
            for j = 1:MCMCSettings.n_chains
                samples_i(:,j) = parameters(burnin+1:burnin+MCMCSettings.sample_size,i,j);
            end
            yy(i)=convergence_i(samples_i);
        end
    
        burnin_convergence=max(max(y(:,1)),max(y(:,2)));
        cross_burnin_convergence=max(yy);
        if (burnin_convergence<=MCMCSettings.burnin_threshold && cross_burnin_convergence<=MCMCSettings.cross_burnin_threshold )
            burned_in=1;
            break
        end
    
    
        from_iteration=burnin+MCMCSettings.sample_size;
        iterations=100;
    
        Datap = Data;
        Datap.chain = 0;
        %Dataz(MCMCSettings.n_chains) = Datap;
        Dataz = repmat(Datap,1,MCMCSettings.n_chains);
        for j=1:MCMCSettings.n_chains
            %Dataz(j) = Datap;
            Dataz(j).chain=j;
        end
        clear Datap
        parfor j=1:MCMCSettings.n_chains
            rng(sum(j*clock));
            [parameters_new(:,:,j),new_log_post(:,j),MCMCAdapt(j)]=HAT_Adaptive_MCMC(from_iteration,iterations,parameters(end,:,j),log_post(end,j),MCMCSettings,MCMCAdapt(j),Dataz(j),Paras,fitted_para_names,FittedPrior,CstrParas,ProjStrat);
        end
        parameters=[parameters;parameters_new];
        log_post=[log_post;new_log_post];
        clear parameters_new log_post_new
    
        save(filename,'parameters','log_post','MCMCAdapt','MCMCSettings','fitted_para_names');
        count = count + 1;
    end
    delete(gcp('nocreate'));

    fprintf([num2str(count * iterations), ' MCMC iterations have been done in the diagnostic test \n'])
end

%third part is to check how long the chains are run for post-burnin before
%we have our sample. While it does not pass the diagnostic tests we
%increase the thinning factor by 1 and run the chains for an additional
%length of the sample size
fprintf('MCMC sampling starts running at %s\n', datetime('now','Format','d-MMM-y HH:mm:ss'))
end_chains=0;
count = 0;
if strcmp(MachineSetting.Type, 'cluster')
parpool(pc);
end
while (end_chains==0 && thin_factor<MCMCSettings.max_thin) 
    
    thin_factor=thin_factor+1;
    yy=zeros(1,d);
    xx=zeros(1,d);
    samples_i = zeros(MCMCSettings.sample_size,MCMCSettings.n_chains);
    for i=1:d
        for j = 1:MCMCSettings.n_chains
            samples_i(:,j) = parameters(burnin+1:thin_factor:burnin+MCMCSettings.sample_size*thin_factor,i,j);
        end
        yy(i)=convergence_i(samples_i);
        xx(i)=eff_samp_size(samples_i);
        % yy(i)=convergence_i(parameters(burnin+1:thin_factor:burnin+MCMCSettings.sample_size*thin_factor,i,1),parameters(burnin+1:thin_factor:burnin+MCMCSettings.sample_size*thin_factor,i,2));
    end
    
    convergence=max(yy);
    ess=min(xx);
    convergence_diagnostic(thin_factor)=convergence;
    effective_size(thin_factor)=ess;
    if (convergence<=MCMCSettings.convergence_threshold && ess>=MCMCSettings.min_ess )
        end_chains=1;
        break
    end
    
    
    from_iteration=burnin+MCMCSettings.sample_size*thin_factor;
    iterations=MCMCSettings.sample_size;
    
    parameters_new=[];
    Datap = Data;
    Datap.chain = 0;
    %Dataz(MCMCSettings.n_chains) = Datap;
    Dataz = repmat(Datap,1,MCMCSettings.n_chains);
    for j=1:MCMCSettings.n_chains
        %Dataz(j) = Datap;
        Dataz(j).chain=j;
    end
    clear Datap
    parfor j=1:MCMCSettings.n_chains
        rng(sum(j*clock));
        [parameters_new(:,:,j),log_post_new(:,j),MCMCAdapt(j)]=HAT_Adaptive_MCMC(from_iteration,iterations,parameters(end,:,j),log_post(end,j),MCMCSettings,MCMCAdapt(j),Dataz(j),Paras,fitted_para_names,FittedPrior,CstrParas,ProjStrat);
    end
    parameters=[parameters;parameters_new];
    log_post=[log_post;log_post_new];
    clear parameters_new log_post_new
    
    save(filename,'parameters','log_post','effective_size','convergence_diagnostic','burnin','MCMCAdapt','MCMCSettings','fitted_para_names');
    count = count + 1;
end
delete(gcp('nocreate'));
clear pc;

fprintf([num2str(count * iterations), ' MCMC iterations have been done in sampling part (including tuning the thinning factor) \n'])

final_ess=effective_size(end);
final_conv=convergence_diagnostic(end);

posterior = zeros(MCMCSettings.sample_size*MCMCSettings.n_chains,size(parameters,2));
neg_log_likelihood = zeros(MCMCSettings.sample_size*MCMCSettings.n_chains,1);
Chains = zeros(MCMCSettings.sample_size*MCMCSettings.n_chains,1);
for j = 1:MCMCSettings.n_chains
    posterior(((j-1)*MCMCSettings.sample_size)+1:((j-1)*MCMCSettings.sample_size)+MCMCSettings.sample_size,:) = ...
        parameters(burnin+1:thin_factor:burnin+thin_factor*MCMCSettings.sample_size,:,j);
    neg_log_likelihood(((j-1)*MCMCSettings.sample_size)+1:((j-1)*MCMCSettings.sample_size)+MCMCSettings.sample_size) = ...
        log_post(burnin+1:thin_factor:burnin+thin_factor*MCMCSettings.sample_size,j);
    Chains(((j-1)*MCMCSettings.sample_size)+1:((j-1)*MCMCSettings.sample_size)+MCMCSettings.sample_size) = j;
end
%posterior=[parameters(burnin+1:thin_factor:burnin+thin_factor*MCMCSettings.sample_size,:,1);parameters(burnin+1:thin_factor:burnin+thin_factor*MCMCSettings.sample_size,:,2)];
%neg_log_likelihood = [log_post(burnin+1:thin_factor:burnin+thin_factor*MCMCSettings.sample_size,1);log_post(burnin+1:thin_factor:burnin+thin_factor*MCMCSettings.sample_size,2)];
save(filename,'parameters','posterior','log_post','Chains', 'neg_log_likelihood', 'effective_size','convergence_diagnostic','burnin','thin_factor','MCMCAdapt','MCMCSettings','fitted_para_names');

end

function ess=eff_samp_size_j(samp)
%this function estimates effective sample size using autocorrelation
samp_size=max(size(samp));
max_lag=samp_size-1;
A=zeros(max_lag,1);
for i=1:max_lag
    A(i)=autocorrelation(samp,i);
    if (i>1)
        if ((A(i)<2/sqrt(samp_size)) || (A(i)-A(i-1)>0))
            A(i)=0;
            break
        end
    end
end

ess=samp_size/(1+2*sum(A));

end

function ess = eff_samp_size(smpls)
n = size(smpls,2);
ess = 0;
for j = 1:n
    ess = ess + eff_samp_size_j(smpls(:,j));
end
end

function autocorr=autocorrelation(samp,k)
% this function calculates the autocorrelation of samp at lag k
    autocorr=(samp(1:end-k)-mean(samp))'*(samp(1+k:end)-mean(samp))/(max(size(samp))*var(samp));

end

% function R=convergence_i(sample1,sample2)
% %this function calculates the convergence diagnostic
% n=max(size(sample1));
% 
% mean1=mean(sample1);
% mean2=mean(sample2);
% mean_both=(mean1+mean2)/2;
% B=n*((mean1-mean_both)^2+(mean2-mean_both)^2);
% 
% var1=var(sample1);
% var2=var(sample2);
% W=(var1+var2)/2;
% 
% est_var=(n-1)/n*W+B/n;
% R=est_var/W;
% 
% end

function R=convergence_i(samples)
%this function calculates the Gelman-Rubin convergence diagnostic
n=max(size(samples));

means=mean(samples);
mean_all=mean(means);
B=n*sum((means-mean_all).^2);

vars=var(samples);
W=mean(vars);

est_var=W*(n-1)/n+B/n;
R=est_var/W;

end
