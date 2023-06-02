
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                %
%   This code generates the scaled screening numbers used in ODE (for fitted part only)                          %
%                                                                                                                %
%   Inputs:                                                                                                      %
%       Screened - vector containing current (i.e. including sampled missing negative values) numbers screened   %
%       modelYears - vector containing years expanded as required to cover pre-data period                       %
%       ScreeningBeforeData - number denoting number of years of AS before data (IP)                             %
%       ScreeningCapacity - number denoting how much of popn can they screen (IP)                                %
%       PopGrowth - number denoting annual growth rate (IP)                                                      %
%       N_H - number denoting human population size used in ODE (Data)                                           %
%       PopSizeYear - number denoting the year in which recorded popn size occured (Data)                        %
%                                                                                                                %
%   Outputs:                                                                                                     %
%       output - struture containing screening information needed in ODE                                         %
%                                                                                                                %
%   Note: only use this code for fitted part because screening capacity is ScreeningCapacity                     %
%                                                                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function output = GetScaledScreening(Screened, modelYears, ScreeningBeforeData, ScreeningCapacity, PopGrowth, N_H, PopSizeYear, IntrpParas)
    if ScreeningBeforeData >= 0
        Scr = [Screened(1) * ones(1, ScreeningBeforeData) Screened];
    else
        Scr = Screened((1-ScreeningBeforeData):length(Screened));
    end % Isangi
    %ScaledPeopleScreened = [Screened(1) * ones(1, ScreeningBeforeData) Screened]...
    %                        .* (PopGrowth .^ double(PopSizeYear - modelYears));
    ScaledPeopleScreened = Scr .* (PopGrowth .^ double(PopSizeYear - modelYears));
    
    ExpectedScreeningTimes = ceil(ScaledPeopleScreened / N_H / ScreeningCapacity);
    ExpectedScreeningTimes(ExpectedScreeningTimes == 0) = 1;
    maxMT = sum(ExpectedScreeningTimes);
    
    ModelScreeningFreq = zeros(1, maxMT);
    ModelScreeningTime = zeros(1, maxMT);
    ModelPeopleScreened = zeros(1, maxMT);
    
    nModelTimes = 0;
    for y = 1 : length(modelYears)
        if ScaledPeopleScreened(y) < N_H * ScreeningCapacity
            nModelTimes = nModelTimes + 1;
            ModelScreeningFreq(nModelTimes)  = 365;
            ModelScreeningTime(nModelTimes)  = modelYears(y);
            ModelPeopleScreened(nModelTimes) = ScaledPeopleScreened(y);
        else
            nModelTimes = nModelTimes + 2;
            ModelScreeningFreq(nModelTimes-1:nModelTimes)  = [365/2 365/2];
            ModelScreeningTime(nModelTimes-1:nModelTimes)  = [modelYears(y) modelYears(y) + 0.5];
            ModelPeopleScreened(nModelTimes-1:nModelTimes) = [ScaledPeopleScreened(y)/2 ScaledPeopleScreened(y)/2];
        end
    end

    
    if (IntrpParas.PPSstart + IntrpParas.PPSend) > 0 && (IntrpParas.PPSstart * IntrpParas.PPSend) == 0
        error('Check PPS related parameters in InterventionParameters')
    end
    if (IntrpParas.NoPSstart + IntrpParas.NoPSend) > 0 && (IntrpParas.NoPSstart * IntrpParas.NoPSend) == 0
        error('Check NoPS related parameters in InterventionParameters')
    end
    
    if IntrpParas.PPSstart > 0 && IntrpParas.PPSstart < ModelScreeningTime(1)
        IntrpParas.PPSstart = ModelScreeningTime(1);
    end
    if IntrpParas.NoPSstart > 0 && IntrpParas.NoPSstart < ModelScreeningTime(1)
        IntrpParas.NoPSstart = ModelScreeningTime(1);
    end
    if IntrpParas.PPSend > 0 && IntrpParas.PPSend > floor(ModelScreeningTime(nModelTimes) + 1)
        IntrpParas.PPSend = floor(ModelScreeningTime(nModelTimes) + 1);
    end
    if IntrpParas.NoPSend > 0 && IntrpParas.NoPSend > floor(ModelScreeningTime(nModelTimes) + 1)
        IntrpParas.NoPSend = floor(ModelScreeningTime(nModelTimes) + 1);
    end
    
    if IntrpParas.PPSstart ~= 0
        if isempty(find(ModelScreeningTime == IntrpParas.PPSstart))
            i = find(ModelScreeningTime > IntrpParas.PPSstart, 1);
            ModelScreeningFreq = [ModelScreeningFreq(1 : i - 2),...
                                  (IntrpParas.PPSstart - ModelScreeningTime(i - 1)) * 365,...
                                  (ModelScreeningTime(i) - IntrpParas.PPSstart) * 365, ...
                                  ModelScreeningFreq(i : end)];
            ModelScreeningTime = [ModelScreeningTime(1 : i - 1) IntrpParas.PPSstart ModelScreeningTime(i : end)];
            ModelPeopleScreened = [ModelPeopleScreened(1 : i - 1) 0 ModelPeopleScreened(i : end)];
            nModelTimes =  nModelTimes + 1;
        end
        
        if isempty(find(ModelScreeningTime == IntrpParas.PPSend)) && IntrpParas.PPSend ~= ceil(ModelScreeningTime(end))
            j = find(ModelScreeningTime(1 : nModelTimes) < IntrpParas.PPSend, 1, 'last');
            ModelScreeningFreq = [ModelScreeningFreq(1 : j - 1),...
                                  (IntrpParas.PPSend - ModelScreeningTime(j)) * 365,...
                                  (ModelScreeningTime(j + 1) - IntrpParas.PPSend) * 365, ...
                                  ModelScreeningFreq(j + 1 : end)];
            ModelScreeningTime = [ModelScreeningTime(1 : j) IntrpParas.PPSend ModelScreeningTime(j + 1 : end)];
            ModelPeopleScreened = [ModelPeopleScreened(1 : j) 0 ModelPeopleScreened(j + 1 : end)];
            nModelTimes =  nModelTimes + 1;
        end
    end
    
    if IntrpParas.NoPSstart ~= 0
        if isempty(find(ModelScreeningTime == IntrpParas.NoPSstart))
            i = find(ModelScreeningTime > IntrpParas.NoPSstart, 1);
            ModelScreeningFreq = [ModelScreeningFreq(1 : i - 2),...
                                  (IntrpParas.NoPSstart - ModelScreeningTime(i - 1)) * 365,...
                                  (ModelScreeningTime(i) - IntrpParas.NoPSstart) * 365, ...
                                  ModelScreeningFreq(i : end)];
            ModelScreeningTime = [ModelScreeningTime(1 : i - 1) IntrpParas.NoPSstart ModelScreeningTime(i : end)];
            ModelPeopleScreened = [ModelPeopleScreened(1 : i - 1) 0 ModelPeopleScreened(i : end)];
            nModelTimes =  nModelTimes + 1;
        end
        
        if isempty(find(ModelScreeningTime == IntrpParas.NoPSend)) && IntrpParas.NoPSend ~= floor(ModelScreeningTime(nModelTimes) + 1)
            j = find(ModelScreeningTime(1 : nModelTimes) < IntrpParas.NoPSend, 1, 'last');
            if isempty(find(ModelScreeningTime(1 : nModelTimes) > IntrpParas.NoPSend, 1))
                ModelScreeningTime(j+1) = ceil(IntrpParas.NoPSend);
            end
            ModelScreeningFreq = [ModelScreeningFreq(1 : j - 1),...
                                  (IntrpParas.NoPSend - ModelScreeningTime(j)) * 365,...
                                  (ModelScreeningTime(j + 1) - IntrpParas.NoPSend) * 365, ...
                                  ModelScreeningFreq(j + 1 : end)];
            ModelScreeningTime = [ModelScreeningTime(1 : j) IntrpParas.NoPSend ModelScreeningTime(j + 1 : end)];
            ModelPeopleScreened = [ModelPeopleScreened(1 : j) 0 ModelPeopleScreened(j + 1 : end)];
            nModelTimes =  nModelTimes + 1;
        end
    end
    
%     ModelScreeningFreq(1:nModelTimes)
%     ModelScreeningTime(1:nModelTimes)
%     ModelPeopleScreened(1:nModelTimes)

    
    output = struct('ScaledPeopleScreened', ScaledPeopleScreened,...
                    'ModelScreeningFreq',   ModelScreeningFreq(1:nModelTimes),...
                    'ModelScreeningTime',   ModelScreeningTime(1:nModelTimes),...
                    'ModelPeopleScreened',  ModelPeopleScreened(1:nModelTimes));
end
