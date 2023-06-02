%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                        %
%   This code runs the adaptive MCMC for the Warwick HAT model,                                          % 
%   it calls Get_log_Prob                                                                                %               
%                                                                                                        %
%   Inputs:                                                                                              %
%       from_iteration - the number of iterations that have been run so far                              %
%       iterations - the number of iterations to be run                                                  %
%       Initial_Param - vector of the current (i.e. initial) values of the fitted parameter              %
%       Initial_log_post - the log-posterior for the current fitted parameters                           %
%       MCMCSettings - structure containing all the settings for the MCMC                                %
%       Data - structure containing historical data of the location                                      %
%       Paras - structure containing all parameters (fixed and fitted)                                   %
%       Intervention - structure containing parameters associated with intervention                      %
%       fitted_para_names - cell array containing the names of the parameters being fitted               %
%       FittedPrior - structure containing prior distributions for fitted parameters                     %
%       ProjStrat - structure containing parameters associated with future strategy                      %
%       MCMCAdapt - structure containing the current mean/covariance of the chain                        %
%                                                                                                        %
%   Outputs:                                                                                             %
%       parameters - array of the sampled fitted parameters                                              %
%       log_post - array of the log-posterior values corresponding to parameters                         %
%       MCMCAdapt - structure containing the current mean/covariance of the chain                        %
%                                                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [parameters,log_post,MCMCAdapt]=HAT_Adaptive_MCMC(from_iteration,iterations,Initial_Param,Initial_log_post,MCMCSettings,MCMCAdapt,Data,Paras,fitted_para_names,FittedPrior,CstrParas,ProjStrat)

% n0 controls the rate at which the covariance matrix drifts towards the empirical covariance matrix
n0=MCMCSettings.n0;
% m0 controls the rate at which MCMCAdapt.lambda tends to one.
m0=MCMCSettings.m0;

%d is number of parameters being fitted *adaptively*
d=length(fitted_para_names(~startsWith(fitted_para_names, 'active_neg_')));

%SD vector used for single parameter proposal during initial period and learning period
Sigma_initial=MCMCSettings.Initial_sigma(1:d);
% Initialise block update covariance matrix
%MCMCAdapt.Sigma = diag(MCMCSettings.Initial_sigma);

Param=Initial_Param;

%calculate log-posterior if we do not already have it
BL=Initial_log_post;
if (from_iteration==0)
    BL=Get_log_Prob(Data, Paras, fitted_para_names, Param, FittedPrior, CstrParas, ProjStrat);
end
parameters=zeros(iterations,length(fitted_para_names));
log_post=zeros(iterations,1);
Accept=0;
BlockUpdates = 0;
AP = 1;

MCMCAdapt.Sigma = symmetricPD(MCMCAdapt.Sigma);   %MCMCAdapt.Sigma + MCMCAdapt.Sigma')/2;

for runnumber=from_iteration+1:from_iteration+iterations
    if mod(runnumber,500)==0
        runnumber
    end
    if runnumber > MCMCSettings.initial_period
        BlockUpdates = BlockUpdates + 1;
        %%BLOCK UPDATE (after initial period)
        %if sum(diag(MCMCAdapt.Sigma) <= 0) ~= 0 || sum(eig(MCMCAdapt.Sigma) <= 0) ~= 0
        %   save([Data.DirStr, Data.LocStr, Data.IDStr, '_neg.mat'],'Param','MCMCAdapt','from_iteration','runnumber','MCMCSettings','parameters','log_post');
        %end
        [Sigma1, SD] = corrcov(MCMCAdapt.Sigma);
        NewP = Param;  %mvnrnd(zeros(d, 1), Sigma1) .* SD' + Param(1:d);
        if mod(runnumber,2)==0
            NewP(1:d) = mvnrnd(zeros(d, 1), Sigma1) .* (MCMCAdapt.lambda*(2.38/sqrt(d))*SD') + Param(1:d);
        else
            NewP(1:d) = mvnrnd(zeros(d, 1), Sigma1) .* ((2.38/sqrt(d))*SD') + Param(1:d);
        end
        % Update Negative active tests
        if Paras.ActiveNeg.Fitting ~= 0
            [NewP, Data, AP] = sample_active_neg(d, Data, Paras, fitted_para_names, NewP, FittedPrior, CstrParas, ProjStrat);
            if ~isinf(AP)
                %Calculate log-posterior for proposed values and then accept/reject
                L=Get_log_Prob(Data, Paras, fitted_para_names, NewP, FittedPrior, CstrParas, ProjStrat);
            else
                L = Inf;
            end
        else
            %Calculate log-posterior for proposed values and then accept/reject
            L=Get_log_Prob(Data, Paras, fitted_para_names, NewP, FittedPrior, CstrParas, ProjStrat);
        end

        jj=runnumber-MCMCSettings.initial_period;
        logmh = log(rand(1,1));
        %fprintf('BL-L %20.10g logmh %20.10g BL %20.10g L %20.10g\n', BL-L, logmh, BL, L);
        if BL-L>logmh     %log(rand(1,1))
            Param=NewP; BL=L;
            Accept=Accept+1;
            if mod(runnumber,2)==0
                MCMCAdapt.lambda=MCMCAdapt.lambda*(1+m0/(m0+jj));
            end
        elseif mod(runnumber,2)==0
            MCMCAdapt.lambda=MCMCAdapt.lambda*(1+m0/(m0+jj))^(-0.305483);
        end
        
        %calculate the mean and covariance of the chain after the initial
        %period, this is updated iteratively
        newMu=(jj/(jj+1))*MCMCAdapt.Mu + Param(1:d)/(jj+1);
        % Contributions (sums of squares) of mean, new mean and current parameters to covariance matrix
        SS_Mu    = MCMCAdapt.Mu'*MCMCAdapt.Mu;
        SS_newMu = newMu'*newMu;
        SS_Param = Param(1:d)'*Param(1:d);
        % Updated formula reducing the weight of initial Sigma
        %MCMCAdapt.Sigma= ((jj-1+n0)*MCMCAdapt.Sigma+jj*(MCMCAdapt.Mu'*MCMCAdapt.Mu) - (jj+1)*newMu'*newMu +Param(1:d)'*Param(1:d))/(jj+n0);
        MCMCAdapt.Sigma= ((jj-1+n0)*MCMCAdapt.Sigma + jj*SS_Mu - (jj+1)*SS_newMu + SS_Param)/(jj+n0);
        % ensure that matrix is symmetric
        %MCMCAdapt.Sigma = (triu(MCMCAdapt.Sigma,1))' + triu(MCMCAdapt.Sigma);
        MCMCAdapt.Sigma = symmetricPD(MCMCAdapt.Sigma);
%         MCMCAdapt.Sigma = (MCMCAdapt.Sigma + MCMCAdapt.Sigma')/2;
%         if sum(eig(MCMCAdapt.Sigma) <= 0) ~= 0
%             add_to_diag = abs(min(diag(MCMCAdapt.Sigma)));
%             n_adjust = 0;
%             while sum(eig(MCMCAdapt.Sigma) <= 0) ~= 0
%                 n_adjust = n_adjust+1;
%                 if (n_adjust > 1000)
%                    error(['HAT_Adaptive_MCMC eigenvalue check: MCMCAdapt.Sigma not +ve semidefinite after ', ...
%                          num2str(n_adjust), ' tries'])
%                 end
%                 MCMCAdapt.Sigma = MCMCAdapt.Sigma + add_to_diag*eye(d);
%             end
%             warning(['HAT_Adaptive_MCMC eigenvalue check: diagonals of MCMCAdapt.Sigma increased by ', num2str(add_to_diag)])
%         end
        MCMCAdapt.Mu=newMu;
    end
    
    if runnumber<=MCMCSettings.initial_period + MCMCSettings.learning_period
        %%SINGLE SITE UPDATES (only during initial & learning periods)
        for j=1:d
            NewP=Param;
            NewP(j)=normrnd(Param(j),sqrt(Sigma_initial(j)));
            % New set of parameters --> new AP --> new active_neg_* proposal too
            if Paras.ActiveNeg.Fitting ~= 0
                [NewP, Data, AP] = sample_active_neg(d, Data, Paras, fitted_para_names, NewP, FittedPrior, CstrParas, ProjStrat);
                if ~isinf(AP)              
                    %Calculate log-posterior for proposed values and then accept/reject
                    L=Get_log_Prob(Data, Paras, fitted_para_names, NewP, FittedPrior, CstrParas, ProjStrat);
                else
                    L = Inf;
                end
            else
                %Calculate log-posterior for proposed values and then accept/reject
                L=Get_log_Prob(Data, Paras, fitted_para_names, NewP, FittedPrior, CstrParas, ProjStrat);
            end
            if BL-L>log(rand(1,1))
                Param=NewP;
                BL=L;
                Sigma_initial(j)=Sigma_initial(j)*1.4;
            else 
                Sigma_initial(j)=Sigma_initial(j)*1.4^(-0.7857143);
            end
        end
        if runnumber == MCMCSettings.initial_period
             % Initial values for block update
             MCMCAdapt.Mu = Param(1:d);
             %MCMCAdapt.Sigma = diag(Sigma_initial);
        end
    end
    
    parameters(runnumber-from_iteration,:)=Param;
    log_post(runnumber-from_iteration)=BL;        

end
if BlockUpdates ~= 0
    MCMCAdapt.BlockUpdates = BlockUpdates;
    MCMCAdapt.Acceptance = Accept / BlockUpdates;
    %fprintf('From = %d, To = %d, BlockUpdates = %d, Acceptance rate = %10.5g, lambda = %10.5g, AP = %10.5g\n', from_iteration+1, from_iteration+iterations, BlockUpdates, Accept / BlockUpdates, MCMCAdapt.lambda, AP);
else
    MCMCAdapt.BlockUpdates = 0;
    MCMCAdapt.Acceptance = 0;
end

end

function Gstar = symmetricPD(G)
%PDforce test for negative eigenvalues of G.
%  If any are found, transform G to ensure PD
%
% Input: G a covariance matrix
% Output: G or a PD transformation of G
G = (G + G')/2;  % Ensure that G is symmetric
[U,D] = eig(G);  % Create eigenvector decomposition
if sum(diag(D)<0) ~= 0 % if G is not positive semi-deifnite, transform G to PD
    Gstar = LRS14(U,D);
else                   % return G
    Gstar = G;
end
end

function Gstar = LRS14(U,D)
%LRS14 transform a covariance matrix G to have +ve eigenvalues
%   Using Larry Schaeffer's method to deal with
%   singularities in a covariance matrix.
%   It is a simple, unweighted procedure that
%   does not "move" the matrix too much from it's
%   input value.
%   Negative eigenvalues are replaced by a function
%   of the smallest positive eigenvalue and a new
%   covariance matrix is generated.
%
%   http://animalbiosciences.uoguelph.ca/~lrs/ELARES/PDforce.pdf
%   Steps 1 to 3 are as in PDforce.pdf
%
% Inputs: U eigenvectors of G
%         D diagonal matrix containing eigenvalues of G
% Step 1
s = sum(D(D<=0))*2;
t = (s*s*100)+1;
% Step 2
d = diag(D);
p = min(d(d>0));     % smallest positive eigenvalue of G
n = D(D<=0);         % negative eigenvalues of G
nStar = p.*(s-n).*(s-n)/t;   % "new" eigenvalues for G^* (Gstar)
D(D<=0) = nStar;             % D^*
% Step 3
Gstar = U*D*U';

% Possible diagnostic outputs
G = U*diag(d)*U'; % recreate input G
dfile = ['LRS14-dt',datestr(datetime,'yyyymmdd_HHMMSS'),'-',num2str(random('poisson',100)),'.mat'];
save(dfile,'U','D','G','Gstar');
warning(['HAT_Adaptive_MCMC:LRS14 -Covariance matrix transformed to PD, see ',dfile])

end
