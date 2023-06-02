%
% Calculate the probability of finding a case with active screening
% for a given year of interest (YoI)
function [active_prob, Disp_act] = Get_Active_Prob(YoIidx, Data, Paras, ProjStrat)
    %k12gc = Paras.k1WY + Paras.k1WW + Paras.k1WP + Paras.k1MY + Paras.k1MW + Paras.k1MP + Paras.k2WY + Paras.k2WW + Paras.k2WP + Paras.k2MY + Paras.k2MW + Paras.k2MP;
    k12gc = 1;

    YoI=Paras.ActiveNeg.Year(YoIidx);
    aneg = Data.N_H * (Data.PopGrowth^double(Data.PopSizeYear - YoI)) * k12gc;
    Paras.(Paras.ActiveNeg.Notation{YoIidx}) = aneg * 0.25;
    Data.ActiveNeg(Paras.ActiveNeg.YearIdx(YoIidx)) = aneg * 0.25;
    scr = GetScaledScreening(Data.ActivePos+Data.ActiveNeg, Data.ModelYears, Paras.ScreeningBeforeData,...
                             Paras.ScreeningCapacity, Data.PopGrowth, Data.N_H, Data.PopSizeYear, Paras);
    % limit Data to year of interest or earlier
    in_MST = floor(scr.ModelScreeningTime) <= YoI;
    in_Y  = floor(Data.Years) <= YoI; 
    Data.ModelScreeningTime  = scr.ModelScreeningTime(in_MST);
    Data.ModelScreeningFreq  = scr.ModelScreeningFreq(in_MST);
    Data.ModelPeopleScreened = scr.ModelPeopleScreened(in_MST);
    
    Data.ModelYears = Data.ModelYears(floor(Data.ModelYears) <= YoI);
    Data.ActiveD1 = Data.ActiveD1(in_Y);
    Data.ActiveD2 = Data.ActiveD2(in_Y);
    Data.ActiveDNa = Data.ActiveDNa(in_Y);
    Data.ActivePos = Data.ActivePos(in_Y);
    Data.ActiveNeg = Data.ActiveNeg(in_Y);
    Data.PassiveD1 = Data.PassiveD1(in_Y);
    Data.PassiveD2 = Data.PassiveD2(in_Y);
    Data.PassiveDNa = Data.PassiveDNa(in_Y);
    Data.PeopleScreened = Data.ActivePos(in_Y)+Data.ActiveNeg(in_Y);
    Data.Years = Data.Years(Data.Years <= YoI);
    
    % NATHAN: ADD NEW DATA CATEGORIES
    Data.ActiveDGM = Data.ActiveDGM(in_Y);
    Data.ActiveDGF = Data.ActiveDGF(in_Y);
    Data.ActiveDGNa = Data.ActiveDGNa(in_Y);
    Data.ActiveDY = Data.ActiveDY(in_Y);
    Data.ActiveDW = Data.ActiveDW(in_Y);
    Data.ActiveDP = Data.ActiveDP(in_Y);
    Data.ActiveDCNa = Data.ActiveDCNa(in_Y);
    Data.PassiveDGM = Data.PassiveDGM(in_Y);
    Data.PassiveDGF = Data.PassiveDGF(in_Y);
    Data.PassiveDGNa = Data.PassiveDGNa(in_Y);
    Data.PassiveDY = Data.PassiveDY(in_Y);
    Data.PassiveDW = Data.PassiveDW(in_Y);
    Data.PassiveDP = Data.PassiveDP(in_Y);
    Data.PassiveDCNa = Data.PassiveDCNa(in_Y);
    
    % 2 lines added REC 20201030
    Paras.death = (1-Paras.u) * Paras.gamma_H;
    Paras.effcy = 1;
    %
    % Get initial conditions
    [meff, ICs] = GetEndemicEq(Data.N_H, Paras, Data); % NATHAN: NEW GetEndemicEq()
    %
    % Run model forward to YoI from ICs
    try
        [Classes, Aggregate] = ODEHATmodel(meff, ICs, Data, Paras, ProjStrat); % NATHAN: NEW ODEHATmodel()
        %
        active_prob = (Aggregate.ActiveM1(end) + Aggregate.ActiveM2(end) + Aggregate.ActiveMFP(end)) * (Data.PopGrowth ^ double(Data.Years(end) - Data.PopSizeYear)) / Data.PeopleScreened(end);
        %
        %Deals with changing FP proportion over time and therefore impacts overdispersion
        ProbTP = (Aggregate.ActiveM1(end) + Aggregate.ActiveM2(end)) /...
            (Aggregate.ActiveM1(end) + Aggregate.ActiveM2(end) + Aggregate.ActiveMFP(end));
        %if ProbTP < 10^(-4)
        %    ProbTP = 0;
        %end        
        Disp_act = Paras.disp_act*ProbTP ;
           % if Disp_act == 0
           %    save('Disp_act_zero.mat','Aggregate','Paras','ProbTP','Data');
           % end
    catch ME
        active_prob = Inf;
        Disp_act = Inf;
        %save(['err', Data.FileStr, Data.LocStr, Data.IDStr, '.mat'], 'meff', 'ICs', 'Data', 'Paras', 'ProjStrat');
        %rethrow(ME)
    end
 
end
