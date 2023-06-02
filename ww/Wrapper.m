function Wrapper(a,b,c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                       %
%   This is the wrapper to run Warwick HAT model on your local machine                                  %
%                                                                                                       %
%   Inputs:                                                                                             %
%       Cloc - a string in the format of ALPHA-3 country codes                                          %
%       Lv1loc - an integer, provine/district/focus index for DRC/CIV,GIN,TCD/UGA                       %
%       Lv2loc - an integer, health zone/county index for DRC/UGA                                       %
%       Lv3loc - an integer, health area index for DRC                                                  %
%       Model - a cell array containing the model indices                                               %
%       EndCaseStr - 
%       EndScrStr -                               %
%       ParaStr - a string of 3-digits ID related to parameter settings                                 %
%       StratDefStr - 
%       MinimumData - an integer indicating the minimum requirement in Data                             %
%       RunMCMC - an integer determining to run MCMC (1 - from scratch, 2 - re-adapt, or                %
%                 3 - additional sampling from posterior) or no MCMC (0)                                %
%       RunEvidence -
%       RunEnsProjection - an integer denoting the number of realizations used in Ensemble Projection   %
%       RunEnsFromExisting - 
%       RunProjection - an integer denoting the number of realizations used in Projection               %
%       StratMin - an integer denoting the first strategy simulated in (Ensemble) Projection            %
%       RunSamples - an integer denoting the number of samples from ODE                                 %
% .     RunDisaster - a boolean determining to use Disaster instead of Projection or not                %
%       RunReactive - a boolean determining to run Reactive or not                                      %
%       RunCFS - a cell array determining to run specific CounterFactual scenarios or not               %
%       MachineSetting - structure containing parameters to control parcluster                          %
%                                                                                                       %
%   Note: (1) Data.mat stored in ../Data/$Cloc$EndCaseStr[$EndScrStr]/ or subdirectory is required      %
%         (2) Paras_$Cloc$EndCaseStr[$EndScrStr]_$ParaStr.mat in current directory is required          %
%         (3) StratDef_$Cloc$EndCaseStr[$EndScrStr]_$StratDefStr.mat in current directory is required   %
%                                                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cloc = 'GIN'; % options: 'CIV', 'DRC', 'GIN', 'TCD', 'UGA'
Lv1loc = a; % set as 0 if a country level run is desired
Lv2loc = b; % set as 0 if a province level run is desired
Lv3loc = c; % set as 0 if a heath zone level run is desired
Model = {'M5'}; % use {'M4', 'M7'} for multiple models
EndCaseStr = '22';
EndScrStr = '23';
ParaStr = '001';
StratDefStr = '001';
MinimumData = 10; % minimum data needed to run MCMC and Projection
RunMCMC = 0; % options: 1(run MCMC) and 0
MCMCOptions = struct('n_chains',2); % see RunInformation.m for all options that can be set here
RunEvidence = 0;
RunEnsProjection = 0; % options: any postive integer(number of realizations for Ensemble Projection) and 0(skip)
RunEnsFromExisting = 0; % options: any positive integer(total number of samples from existing model projections) and 0(skip)
RunProjection = 500; % options: any postive integer(number of realizations for Projection) and 0(skip)
StratMin = 4; % options: any positive integer (first strategy to run in the expanded Strategy table; 1 for new simulation)
RunSamples = 100; % options: any positive integer (number of samples per realisation, usually 10) or 0 (Projection only)
RunDisaster = 0; % options: 1(use Disaster.m to run projections) and 0 (use Projection.m to run projecitons)
RunReactive = 0; % options: 1(run all reactive strategies) and 0
RunCFS = {'0'}; % options: {'0'} for actual, {'1'} for no vector control, {'2'} for no passive improvement, {'1', '2'} for no vector control and no passive improvement
MachineSetting.Type = 'local'; % 'local' or 'cluster'
MachineSetting.NumWorkers = 2; % cluster only
MachineSetting.JobStorageLocation = 'tmp'; % cluster only

fprintf([Cloc, '(', num2str(a), ', ', num2str(b), ', ', num2str(c), ') simulation starts at ', datestr(datetime('now','Format','d-MMM-y HH:mm:ss')), '\n'])
Run(Cloc, Lv1loc, Lv2loc, Lv3loc, Model, EndCaseStr, EndScrStr, ParaStr, ...
    StratDefStr, MinimumData, RunMCMC, MCMCOptions, RunEvidence, RunEnsProjection, ...
    RunEnsFromExisting, RunProjection, StratMin, RunSamples, RunDisaster, ...
    RunReactive, RunCFS, MachineSetting)
fprintf([Cloc, '(', num2str(a), ', ', num2str(b), ', ', num2str(c), ') simulation ended successfully at ', datestr(datetime('now','Format','d-MMM-y HH:mm:ss')), '\n'])
