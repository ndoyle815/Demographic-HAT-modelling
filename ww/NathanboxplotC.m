
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                     %
%   This code plots fitting results of the age adapted Warwick HAT model                              %
%                                                                                                     %
%   Inputs:                                                                                           %
%       Cloc - string containing ALPHA-3 country codes                                                %
%       Ploc - number denoting the provine index                                                      %
%       Zloc - number denoting the health zone index                                                  %
%       Aloc - number denoting the health area index                                                  %
%       Model - cell array containing the model indices                                               %
%       ParaStr - string containing ALPHA-3 country code and 3-digits related to parameter settings   %
%                                                                                                     %
%   File required: Data.mat & Fitted_*.mat                                                            %
%                                                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Pick/change options here to plot different locations/data/colours/language

Cloc = 'GIN'; % options: 'CIV', 'DRC', 'GIN', 'TCD', 'UGA'
Lv1loc = 4; %options: focus for GIN, TCD, district for UGA, health district for CIV, province for DRC
Lv2loc = 0; %options: health zone for DRC, subprefecture for CIV
Lv3loc = 0; %options: health area for DRC
Model = 'M5'; %options: M1, M2,....,M8
ModelFitType = 'MCMC'; %options: MCMC or pMCMC
ModelProjType = 'DetProj'; % options: DetProj or StochProj
LastDataYear = 22; %enter last year for data set e.g. 22 for 2022
LastGapYear = 23; %enter last year for "gap year" e.g. 23 for 2023
ParaStr = '001';
StratStr = '001';

Lang = 'EN'; %options: EN or FR to make french or english figures
ShowGapYear = 1; %options: 1 to show gap year projection or 0 to only show fitted years
AddNewInf = 0; %options: 1 for 4th line with new infections or 0 for seperate new infection plot

%HAT MEPP colours
MyBlue = [67 147 195]/255;
MyGreen = [27 120 55]/255;
MyRed = [178 24 42]/255;
MyPurple = [118,42,131]/255;
MyGrey = [240 100 10]/256; %actually orange - grey is [0.6 0.6 0.6];
MyBlack = [0 0 0];

%Choose box colour
MyColor = MyBlue;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates values/positions for plotting boxes

%Number of years with data (adds on 1 due to starting at 2000)
Y1 = LastDataYear+1; 
%Number of gap years simulated 
GapY = LastGapYear-LastDataYear;
%Number of years to show in plots
if ShowGapYear==1
    ShowYears = Y1+GapY;
else
    ShowYears = Y1;
end

%Gives first year of data
firstyear = 2000;

%Defines the number of years to show
xmax = Y1 + GapY;
%Defines the positions for boxes for each year (all offset by -0.5)
xall = 0.5 : xmax - 0.5;
% Box width
width = 0.2;

%Creates a vector of years to show
yearall = firstyear : firstyear + xmax - 1;
%Defines only xticks for every 5 years
xpos = xall(mod(yearall, 5) == 0);
%Prints tick labels for the xticks for every 5 years
xlabels = num2cell(yearall(mod(yearall, 5) == 0));
%Prints minor ticks for every year
xminorpos = (0 : Y1 + GapY) + 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gets data and model outputs

DataStr = strcat(num2str(LastDataYear),'[',num2str(LastGapYear),']'); 
DataOnlyStr = num2str(LastDataYear); 

switch Cloc
    case {'CIV', 'GIN', 'TCD'}
        if Lv1loc ~= 0 % focus level simulation
            load(['../Data/', Cloc, DataStr, '/Data.mat']);
            f.name = strcat('F', num2str(Lv1loc), '_', CCLOC{Lv1loc});
            names = {f.name, Cloc};
            location = Lv1loc;
            Dir = ['../Result/', Cloc, '/', f.name, '/'];
        end
    case 'DRC'
        if Lv3loc ~= 0 % health Area level simulation
            p = dir(['../Data/', Cloc, DataStr, '/P', num2str(Lv1loc), '_*']);
            z = dir(['../Data/', Cloc, DataStr, '/', p.name, '/Z', num2str(Lv2loc), '_*']);
            load(['../Data/', Cloc, DataStr, '/', p.name, '/', z.name, '/Data.mat']);
            a.name = strcat('A', num2str(Lv3loc), '_', CCLOC{Lv3loc});
            names = {a.name, z.name, p.name, Cloc};
            location = Lv3loc;
            Dir = ['../Result/', Cloc, '/', p.name, '/', z.name, '/', a.name, '/'];
        
        elseif Lv2loc ~= 0 % health Zone level simulation
            p = dir(['../Data/', Cloc, DataStr, '/P', num2str(Lv1loc), '_*']);
            load(['../Data/', Cloc, DataStr, '/', p.name, '/Data.mat']);
            z.name = strcat('Z', num2str(Lv2loc), '_', CCLOC{Lv2loc});
            names = {z.name, p.name, Cloc};
            location = Lv2loc;
            Dir = ['../Result/', Cloc, '/', p.name, '/', z.name, '/'];
        
        elseif Lv1loc ~= 0 % province level simulation
            load(['../Data/', Cloc, DataStr, '/Data.mat']);
            p.name = strcat('P', num2str(Lv1loc), '_', CCLOC{Lv1loc});
            names = {p.name, Cloc};
            location = Lv1loc;
            Dir = ['../Result/', Cloc, '/', p.name, '/'];
        end
    case 'UGA'
        if Lv2loc ~= 0 % county level simulation
            d = dir(['../Data/', Cloc, DataStr, '/D', num2str(Lv1loc), '_*']);
            load(['../Data/', Cloc, DataStr, '/', d.name, '/Data.mat']);
            c.name = strcat('C', num2str(Lv2loc), '_', CCLOC{Lv2loc});
            names = {c.name, d.name, Cloc};
            location = Lv2loc;
            Dir = ['../Result/', Cloc, '/', d.name, '/', c.name, '/'];
    
        elseif Lv1loc ~= 0 % district level simulation
            load(['../Data/', Cloc, DataStr, '/Data.mat']);
            d.name = strcat('D', num2str(Lv1loc), '_', CCLOC{Lv1loc});
            names = {d.name, Cloc};
            location = Lv1loc;
            Dir = ['../Result/', Cloc, '/', d.name, '/'];
        end
end
    
LocStr = LOCSTR{location}; % location info string ending with the name of smallest scale
        
% Data
%Active = [ACTIVE1(location, :) + ACTIVE2(location, :) + ACTIVENa(location, :)];% GapActive(location, :)];
%Passive = [PASSIVE1(location, :) + PASSIVE2(location, :) + PASSIVENa(location, :)];% GapPassive(location, :)];
%load('../orig_screened.mat')
%Screened = [SCREENED(location, :)]; %modify here for extra GapYears: GapScreening{location, 1:GapY}];
ActiveY = ACTIVEC1(location, :);
ActiveW = ACTIVEC2(location, :);
ActiveP = ACTIVEC3(location, :);
PassiveY = PASSIVEC1(location, :);
PassiveW = PASSIVEC2(location, :);
PassiveP = PASSIVEC3(location, :);

% Fitted part
load([Dir, 'Fitted_', ModelFitType, '_', ModelProjType, '_', Model, '_', LocStr, '_Data', Cloc, DataStr, '_Paras', ParaStr,'_StratDef',StratStr, '.mat']);
%ActiveQ1 = quantile(SampledActive1 + SampledActive2, [0.025 0.25 0.5 0.5 0.75 0.975]);
%PassiveQ1 = quantile(SampledPassive1 + SampledPassive2, [0.025 0.25 0.5 0.5 0.75 0.975]);
%NewInfQ1 = quantile(NewInf, [0.025 0.25 0.5 0.5 0.75 0.975]);
ActiveYQ1 = quantile(SampledActiveAgeY, [0.025 0.25 0.5 0.5 0.75 0.975]);
ActiveWQ1 = quantile(SampledActiveAgeW, [0.025 0.25 0.5 0.5 0.75 0.975]);
ActivePQ1 = quantile(SampledActiveAgeP, [0.025 0.25 0.5 0.5 0.75 0.975]);
PassiveYQ1 = quantile(SampledPassiveAgeY, [0.025 0.25 0.5 0.5 0.75 0.975]);
PassiveWQ1 = quantile(SampledPassiveAgeW, [0.025 0.25 0.5 0.5 0.75 0.975]);
PassivePQ1 = quantile(SampledPassiveAgeP, [0.025 0.25 0.5 0.5 0.75 0.975]);


% Posterior of missing screening data
%load([Dir, 'Posterior_', ModelFitType, '_', Model, '_', LocStr, '_Data', Cloc, DataOnlyStr, '_Paras', ParaStr, '.mat']);
%ix=find(isnan(Screened));
%clear ScreenQ1
%for k=1:length(ix)
%    a=['Posterior.active_neg_',num2str(1999+ix(k))];
%    Z = eval(a);
%ScreenQ1(:,k) = quantile(Z, [0.025 0.25 0.5 0.5 0.75 0.975]);
%end
%PosScreen=ix-0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% Figure setting
figure('Name', ['Fitting_', LocStr], 'NumberTitle', 'off')

% ===== Screening =====
%subplot(2+AddNewInf, 1, 1)
%stairs(0 : length(Screened), [Screened Screened(end)], 'LineWidth', 1, 'Color', 'k')

%Draw on missing screening coverage from posterior
%if ~isempty(ix)
%    hold on
%    Katboxplot(ScreenQ1,Screened,ShowYears,MyColor,PosScreen,xpos,xlabels,xminorpos,0.5)
%else
%    plt = gca;
%    set(plt, 'children', flipud(get(gca, 'children')))
%    xticks(xpos)
%    xticklabels(xlabels)
%    plt.XRuler.TickLabelGapOffset = -2;
%    plt.XAxis.TickDirection = 'out';
%    plt.XAxis.TickLength = [0.018 1];
%    plt.YAxis.TickLength = [0.012 1];
%    plt.XAxis.MinorTick = 'on';
%    plt.XAxis.MinorTickValues = xminorpos;
%    box off
        
%    xlim([0, ShowYears])
%    ylim([0, Inf*1.15])
    %ylim([0, 1.15 * max(Screening)])
%end


%if strcmp(Lang,'EN')==1    
%    ylabel('# screened')
%else
%    ylabel('# testée')
%end
        
% ===== Active =====
subplot(2+AddNewInf, 1, 1)
s1 = scatter([20:22]+0.5-width, ActiveY(20:22), 80, MyGreen, 'Filled', 'o', 'MarkerEdgeColor', 'k');
hold on
s2 = scatter([20:22]+0.5, ActiveW(20:22), 80, MyPurple, 'Filled', 'o', 'MarkerEdgeColor', 'k');
hold on
s3 = scatter([20:22]+0.5+width, ActiveP(20:22), 80, MyGrey, 'Filled', 'o', 'MarkerEdgeColor', 'k');
hold on
%Calls Kat's boxplot function
k1 = Katboxplot(ActiveYQ1, ActiveY, ShowYears, MyGreen, xall-width, xpos, xlabels, xminorpos,width);
hold on
k2 = Katboxplot(ActiveWQ1, ActiveW, ShowYears, MyPurple, xall, xpos, xlabels, xminorpos,width);
hold on
k3 = Katboxplot(ActivePQ1, ActiveP, ShowYears, MyGrey, xall+width, xpos, xlabels, xminorpos,width);

legend([k1 k2 k3 s1 s2 s3],{'Model (Y)', 'Model (W)', 'Model (P)', 'Data (Y)', 'Data (W)', 'Data (P)'},'FontSize',14)
xlim([-0.5 23])
ylim([0, 1.15 * max([max(ActiveYQ1) max(ActiveWQ1) max(ActivePQ1)])])

if strcmp(Lang,'EN')==1 
    ylabel('Active cases')
else
    ylabel('Cas actives')
end
title(CCLOC{location}, 'FontSize', 12)

% ===== Passive =====
subplot(2+AddNewInf, 1, 2)
s1 = scatter([20:22]+0.5-width, PassiveY(20:22), 80, MyGreen, 'Filled', 'o', 'MarkerEdgeColor', 'k');
hold on
s2 = scatter([20:22]+0.5, PassiveW(20:22), 80, MyPurple, 'Filled', 'o', 'MarkerEdgeColor', 'k');
hold on
s3 = scatter([20:22]+0.5+width, PassiveP(20:22), 80, MyGrey, 'Filled', 'o', 'MarkerEdgeColor', 'k');
hold on
%Calls Kat's boxplot function
k1 = Katboxplot(PassiveYQ1, PassiveY, ShowYears, MyGreen, xall-width, xpos, xlabels, xminorpos, width);
hold on
k2 = Katboxplot(PassiveWQ1, PassiveW, ShowYears, MyPurple, xall, xpos, xlabels, xminorpos, width);
hold on
k3 = Katboxplot(PassivePQ1, PassiveP, ShowYears, MyGrey, xall+width, xpos, xlabels, xminorpos,width);

legend([k1 k2 k3 s1 s2 s3],{'Model (Y)', 'Model (W)', 'Model (P)', 'Data (Y)', 'Data (W)', 'Data (P)'},'FontSize',14)
xlim([-0.5 23])
ylim([0, 1.15 * max([max(PassiveYQ1) max(PassiveWQ1) max(PassivePQ1)])])

if strcmp(Lang,'EN')==1 
    if AddNewInf==0
        xlabel('Year')
    end
    ylabel('Passive cases')
else
    if AddNewInf==0
        xlabel('Année')
    end
    ylabel('Cas passives')
end

% ===== New infections =====
%if AddNewInf==1
%    subplot(3+AddNewInf, 1, 4)
%else
%    figure('Name', ['NewInf_', LocStr], 'NumberTitle', 'off')
%end

%Calls Kat's boxplot function
%Katboxplot(NewInfQ1, 0, ShowYears, MyColor, xall, xpos, xlabels, xminorpos);

%if strcmp(Lang,'EN')==1 
%    xlabel('Year')
%    ylabel('New infections')
%else
%    xlabel('Année')
%    ylabel('Nouvelle infections')
%end

%Changes size of figure to be longer if including 4 rows
%if AddNewInf==1
%   set(gcf,'position',[1069         775         560         562])
%else
%    set(gcf,'position',[508   908   560   420])
%end


%To add in
%(1) legend underneath fitted/observed
%(2) gap year screening dashed

     