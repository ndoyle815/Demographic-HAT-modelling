classdef RunInformation
    % Primarily designed for use in applying revised file naming
    % conventions:
    %   - Store components of filenames and filepaths
    %   - Provide methods to create filenames and filepaths    
    properties
        LocStr
        Ctry
        Model
        NumModel
        currentModelIdx {mustBeInteger, mustBePositive}
        FitDataStr
        FitDataID
        DataStr
        DataID
        FitMethod
        ParaStr
        ParaID
        ProjMethod
        CFS
        StratDefStr
        StratID
        Lvloc
        DataPath
        names
        locationIdx
        Dir
        ParasFile
        OutputParasPath
        StratDefFile
        OutputStratDefPath
        RunMCMC
        MCMCSettings
        RunProjection
        RunEnsProjection
        RunEnsFromExisting
        RunSamples
        RunReactive
        RunCFS
        StratMin
    end
    
    methods
        function obj = RunInformation(Cloc, Lv1loc, Lv2loc, Lv3loc,...
                Model, EndCaseStr, EndScrStr, ParaStr, StratDefStr,...
                RunMCMC, MCMCOptions, RunProjection, RunEnsProjection, RunEnsFromExisting,...
                RunSamples, RunReactive, RunCFS, StratMin)
            %RunInformation Construct an instance of this class.
            %   Stores a lot of run-related information used in
            %   constructing filenames/filepaths
            
            % Fitting and Projection methods hard-coded for now
            obj.FitMethod = 'MCMC';
            obj.ProjMethod = 'DetProj';
            % hard-coded results directory
            ResDir = '../Result/';
            if exist(ResDir,'dir') == 0
                mkdir(ResDir);
            end

            CFS = {'_NoVC', '_NoPSImp'};
            cfs = '';
            if sum(strcmp(RunCFS, '0')) == 0
                cfs = '(CFS';
                if sum(strcmp(RunCFS, '1')) == 1
                    cfs = strcat(cfs, CFS{1});
                end
                if sum(strcmp(RunCFS, '2')) == 1
                    cfs = strcat(cfs, CFS{2});
                end
                cfs = strcat(cfs, ')');
            end
            obj.CFS = cfs;

            obj.Model = Model;
            obj.NumModel = length(Model);
            obj.currentModelIdx = 1;
            obj.Ctry = Cloc;
            EndCaseStr = strtrim(EndCaseStr);
            obj.FitDataStr = EndCaseStr;
            obj.FitDataID  = ['_Data', Cloc, obj.FitDataStr];
            EndScrStr = strtrim(EndScrStr);
            if isempty(EndScrStr) || strcmp(EndScrStr, EndCaseStr)
                DataStr = strtrim(EndCaseStr);
            else
                DataStr = [strtrim(EndCaseStr) '[' EndScrStr ']'];
            end
            obj.DataStr = DataStr;
            obj.DataID  = ['_Data', Cloc, DataStr];
            obj.ParaStr = strtrim(ParaStr);
            obj.ParaID  = ['_Paras', obj.ParaStr];
            obj.StratDefStr = strtrim(StratDefStr);
            obj.StratID = ['_StratDef', obj.StratDefStr];
            obj.Lvloc = [Lv1loc, Lv2loc, Lv3loc];
            switch Cloc
                case {'GIN', 'TCD'}
                    if Lv1loc ~= 0 % focus level simulation
                        DataPath = ['../Data/', Cloc, DataStr, '/Data.mat'];
                        load(DataPath, 'CCLOC', 'LOCSTR');
                        f.name = strcat('F', num2str(Lv1loc), '_', CCLOC{Lv1loc});
                        names = {f.name, Cloc};
                        location = Lv1loc;
                        Dir = [ResDir, Cloc, '/', f.name, '/'];
                    end
                case 'CIV'
                    if Lv2loc ~= 0 % sub-prefecture level simulation
                        d = dir(['../Data/', Cloc, DataStr, '/D', num2str(Lv1loc), '_*']);
                        DataPath = ['../Data/', Cloc, DataStr, '/', d.name, '/Data.mat'];
                        load(DataPath, 'CCLOC', 'LOCSTR');
                        s.name = strcat('S', num2str(Lv2loc), '_', CCLOC{Lv2loc});
                        names = {s.name, d.name, Cloc};
                        location = Lv2loc;
                        Dir = [ResDir, Cloc, '/', d.name, '/', s.name, '/'];

                    elseif Lv1loc ~= 0 % health district level simulation
                        DataPath = ['../Data/', Cloc, DataStr, '/Data.mat'];
                        load(DataPath, 'CCLOC', 'LOCSTR');
                        d.name = strcat('D', num2str(Lv1loc), '_', CCLOC{Lv1loc});
                        names = {d.name, Cloc};
                        location = Lv1loc;
                        Dir = [ResDir, Cloc, '/', d.name, '/'];

                    end
                case 'DRC'
                    if Lv3loc ~= 0 % health Area level simulation
                        c = dir(['../Data/', Cloc, DataStr, '/C', num2str(Lv1loc), '_*']);
                        z = dir(['../Data/', Cloc, DataStr, '/', c.name, '/Z', num2str(Lv2loc), '_*']);
                        DataPath = ['../Data/', Cloc, DataStr, '/', c.name, '/', z.name, '/Data.mat'];
                        load(DataPath, 'CCLOC', 'LOCSTR');
                        a.name = strcat('A', num2str(Lv3loc), '_', CCLOC{Lv3loc});
                        names = {a.name, z.name, c.name, Cloc};
                        location = Lv3loc;
                        Dir = [ResDir, Cloc, '/', c.name, '/', z.name, '/', a.name, '/'];
                    elseif Lv2loc ~= 0 % health Zone level simulation
                        c = dir(['../Data/', Cloc, DataStr, '/C', num2str(Lv1loc), '_*']);
                        DataPath = ['../Data/', Cloc, DataStr, '/', c.name, '/Data.mat'];
                        load(DataPath, 'CCLOC', 'LOCSTR');
                        z.name = strcat('Z', num2str(Lv2loc), '_', CCLOC{Lv2loc});
                        names = {z.name, c.name, Cloc};
                        location = Lv2loc;
                        Dir = [ResDir, Cloc, '/', c.name, '/', z.name, '/'];
                    elseif Lv1loc ~= 0 % coordination level simulation
                        DataPath = ['../Data/', Cloc, DataStr, '/Data.mat'];
                        load(DataPath, 'CCLOC', 'LOCSTR');
                        c.name = strcat('P', num2str(Lv1loc), '_', CCLOC{Lv1loc});
                        names = {c.name, Cloc};
                        location = Lv1loc;
                        Dir = [ResDir, Cloc, '/', c.name, '/'];
                    end
                case 'UGA'
                    if Lv2loc ~= 0 % county level simulation
                        d = dir(['../Data/', Cloc, DataStr, '/D', num2str(Lv1loc), '_*']);
                        DataPath = ['../Data/', Cloc, DataStr, '/', d.name, '/Data.mat'];
                        load(DataPath, 'CCLOC', 'LOCSTR');
                        c.name = strcat('C', num2str(Lv2loc), '_', CCLOC{Lv2loc});
                        names = {c.name, d.name, Cloc};
                        location = Lv2loc;
                        Dir = [ResDir, Cloc, '/', d.name, '/', c.name, '/'];
                    elseif Lv1loc ~= 0 % district level simulation
                        DataPath = ['../Data/', Cloc, DataStr, '/Data.mat'];
                        load(DataPath, 'CCLOC', 'LOCSTR');
                        d.name = strcat('D', num2str(Lv1loc), '_', CCLOC{Lv1loc});
                        names = {d.name, Cloc};
                        location = Lv1loc;
                        Dir = [ResDir, Cloc, '/', d.name, '/'];
                    end
                otherwise
                    error(['ERROR unknown country code: ', Cloc])
            end
            obj.LocStr = LOCSTR{location};
            obj.DataPath = DataPath;
            obj.names = names;
            obj.locationIdx = location;
            obj.Dir = Dir;
            obj.ParasFile = ['Paras_', obj.Ctry, obj.DataStr, '_', obj.ParaStr, '.mat'];
            obj.OutputParasPath = [ResDir, obj.ParasFile];
            if not(isempty(obj.StratDefStr))
                obj.StratDefFile =['StratDef_', obj.Ctry, obj.DataStr, '_', obj.StratDefStr, '.mat'];
                obj.OutputStratDefPath = [ResDir, obj.StratDefFile];
            end
            obj.RunMCMC = RunMCMC;
            % *****   MCMC options   *****
            %  - set default values for:
            %    n_chains - number of MCMC chains
            %    sample_size - number of samples from each chain, total sample size=n_chains*sample_size
            %    initial_period - number of iterations with non-adaptive updates at the beginning, before we start learning the covariance matrix
            %    learning_period - number of iterations with non-adaptive updates after the initial_period, during which we start learning the covariance matrix
            %    max_burnin - the maximum length of burnin
            %    thresholds for convergence diagnostic, the aim is to drop below these thresholds
            %      burnin_threshold - threshold for the convergence diagnostic when applied within chain during the burnin
            %      cross_burnin_threshold: threshold for the convergence diagnostic when applied between chain during the burnin
            %      convergence_threshold: threshold for the convergence diagnostic when applied to the post-burnin, post-thinned posterior sample
            %      min_ess: minimum effective sample size, set to half of the total sample size. MCMC can be ended when min_ess has been reached
            %    max_thin: the maximum factor for thinning
            %    Adaptive MCMC controls
            %      n0 - controls the rate at which the covariance matrix drifts towards the empirical covariance matrix
            %      m0 - controls the rate at which MCMCAdapt.lambda tends to one.
            n_chains = 2;       %number of chains
            sample_size = 1000; %number of samples from each chain
            obj.MCMCSettings=struct(...
                'run_type',  RunMCMC,...
                'n_chains',  n_chains,...
                'sample_size',  sample_size,...
                'initial_period',  1000,...
                'learning_period',  0,...
                'burnin_threshold',  1.2,...
                'cross_burnin_threshold',  1.5, ...
                'convergence_threshold',  1.2,...
                'max_burnin',  200000,...
                'min_ess',  round(0.5*n_chains*sample_size),...
                'max_thin',  200,...
                'n0',  10, ...
                'm0',  50);
            % - overwrite with supplied values
            option_names = fieldnames(MCMCOptions);
            d = length(option_names);
            if d ~= 0
                for i = 1:d
                    obj.MCMCSettings.(option_names{i}) = MCMCOptions.(option_names{i});
                end
                if ~isfield(MCMCOptions,'min_ess')
                    obj.MCMCSettings.min_ess = round(0.5*obj.MCMCSettings.n_chains*obj.MCMCSettings.sample_size);
                end
            end            
            % *****   end of MCMC options   *****
            obj.RunProjection = RunProjection;
            obj.RunEnsProjection = RunEnsProjection;
            obj.RunEnsFromExisting = RunEnsFromExisting;
            obj.RunSamples = RunSamples;
            obj.RunReactive = RunReactive;
            obj.RunCFS = RunCFS;
            obj.StratMin = StratMin;
        end
        
        function fn = FitFileName(obj, descriptor, varargin)
            %FitFileName file *name* generator for fitting of model
            %   Descriptor : eg 'Posterior', 'Fitted'
            %   optional file suffix: defaults to ".mat"
            if nargin == 2
                sfx = '.mat';
            else
                sfx = varargin{1};
            end
            MethMod = ['_', obj.FitMethod, '_', obj.Model{obj.currentModelIdx}, '_'];
            fn = [descriptor, MethMod, obj.LocStr, obj.FitDataID, obj.ParaID, sfx];
        end
        
        function fn = ProjFileName(obj, descriptor)
            %FileNameR file *name* generator including Reactive stategy
            %   Descriptor : eg 'Fitted'
            %   ReactJ : integer reactive strategy identifier
            MethMod = ['_', obj.FitMethod, '_', obj.ProjMethod, obj.CFS, '_',...
                obj.Model{obj.currentModelIdx}, '_'];
            fn = [descriptor, MethMod, obj.LocStr, obj.DataID,...
                  obj.ParaID, obj.StratID, '.mat'];
        end
        
        function fn = ProjFileNameR(obj, descriptor, ReactJ)
            %FileNameR file *name* generator including Reactive stategy
            %    - across Strategies, within Reactive strategy eg Elimination*.mat
            %   Descriptor : eg 'Elimination'
            %   ReactJ : integer reactive strategy identifier
            MethMod = ['_', obj.FitMethod, '_', obj.ProjMethod, obj.CFS, '_', obj.Model{obj.currentModelIdx},...
                '_React', num2str(ReactJ), '_'];
            fn = [descriptor, MethMod, obj.LocStr, obj.DataID,...
                  obj.ParaID, obj.StratID, '.mat'];
        end
        
        function fn = ProjFileNameSR(obj, descriptor, StratI, ReactJ)
            %ProjFileName file *name* generator including Strategy and Reactive strategy
            %    - within Strategy and within Reactive strategy eg Projection outputs
            %   Descriptor : eg 'Projection'
            %   StratI : integer strategy identifier
            %   ReactJ : integer reactive strategy identifier
            MethMod = ['_', obj.FitMethod, '_', obj.ProjMethod, obj.CFS, '_', obj.Model{obj.currentModelIdx},...
                '_Strat', num2str(StratI), '_React', num2str(ReactJ), '_'];
            fn = [descriptor, MethMod, obj.LocStr, obj.DataID,...
                  obj.ParaID, obj.StratID, '.mat'];
        end

        function fp = FitFilePath(obj, descriptor, varargin)
            %FitFilePath file *path* generator for fitting of model
            %   (path equivalent to FitFileName)
            %   Descriptor : eg 'Posterior', 'Fitted'
            %   optional file suffix: defaults to ".mat"
            if nargin == 2
                sfx = '.mat';
            else
                sfx = varargin{1};
            end
            fn = obj.FitFileName(descriptor, sfx);
            fp = [obj.Dir fn];
        end

        function fp = ProjFilePath(obj, descriptor)
            %ProjFilePathR file *path* generator for fitting of model
            %   (path equivalent to ProjFileNameR)
            %   Descriptor : eg 'Fitted'
            fn = obj.ProjFileName(descriptor);
            fp = [obj.Dir fn];
        end

        function fp = ProjFilePathR(obj, descriptor, ReactJ)
            %ProjFilePathR file *path* generator for fitting of model
            %   (path equivalent to ProjFileNameR)
            %   Descriptor : eg 'Posterior', 'Fitted'
            %   ReactJ : integer reactive strategy identifier
            fn = obj.ProjFileNameR(descriptor, ReactJ);
            fp = [obj.Dir fn];
        end

        function fp = ProjFilePathSR(obj, descriptor, StratI, ReactJ)
            %ProjFilePathSR file *path* generator for fitting of model
            %   (path equivalent to ProjFileNameSR)
            %   Descriptor : eg 'Posterior', 'Fitted'
            %   StratI : integer strategy identifier
            %   ReactJ : integer reactive strategy identifier
            fn = obj.ProjFileNameSR(descriptor, StratI, ReactJ);
            fp = [obj.Dir fn];
        end
        
    end
end

