% Code to plot and compare YEOT and YEPHP for Guinea regions as output from
% projection strategies arising from fitted posterior samples
clear all

% Plotting preferences
set(0,'defaultlinelinewidth',2)
set(0,'defaultaxesfontsize',16)
set(0,'DefaultLegendFontSize',18)
set(0,'DefaultAxesTitleFontWeight','normal')

myblue = [0.2627 0.5765 0.7647];
mygreen = [0.1059 0.4706 0.2157];
mypurple = [0.4627 0.1647 0.5137];
myred = [0.6980 0.09412 0.1647];

% read data
reg = "FO";

if reg == "BE"
    % BOFFA EAST
    load('./Result/GIN/F1_BoffaEast/Elimination_MCMC_DetProj_M5_React0_GIN_F1_BoffaEast_DataGIN22[23]_Paras001_StratDef001.mat');                     
end
if reg == "BW"
    % BOFFA WEST
    load('./Result/GIN/F2_BoffaWest/Elimination_MCMC_DetProj_M5_React0_GIN_F2_BoffaWest_DataGIN22[23]_Paras001_StratDef001.mat');
end
if reg == "DU"
    % DUBREKA
    load('./Result/GIN/F3_Dubreka/Elimination_MCMC_DetProj_M5_React0_GIN_F3_Dubreka_DataGIN22[23]_Paras001_StratDef001.mat');
end
if reg == "FO"
    % FORECARIAH
    load('./Result/GIN/F4_Forecariah/Elimination_MCMC_DetProj_M5_React0_GIN_F4_Forecariah_DataGIN22[23]_Paras001_StratDef001.mat');
end

n = size(PostID,1);    % number of realisations
ns = size(YEOT,1)/n;   % number of samples per realisation
nstrat = size(PEOT,2); % number of strategies

figure(1)
for i = 1:nstrat
    plot(ElimByYears,PEOT(:,i))
    hold on
end
axis([min(ElimByYears) max(ElimByYears) 0 1.1])
xlabel('Year')
title('PEOT')

% get expected YEOT for each realisation
myYEOTs = zeros(n,1);
for samp = 1:n
    for strat = 1:nstrat
        myYEOTs(samp,strat) = mean(YEOT((samp-1)*ns+1:samp*ns,strat));
    end
end

% quantiles
YEOTQs = quantile(myYEOTs, [0.05 0.1 0.25 0.5 0.5 0.75 0.9 0.95]);

EOTyears = [min(ElimByYears):1:max(YEOTQs)];

f = figure(2);
f.Position = [205 450 1200 300];
%set(gcf,'Renderer','painters')

boxchart(YEOTQs,'BoxFaceColor',myblue,'BoxEdgeColor','k','BoxMedianLineColor','k','BoxFaceAlpha',0.5), view(90,90)
hold on
scatter([1, 2, 3, 4], mean(myYEOTs), 100, myred, "filled", "diamond")

for idx = 1:nstrat
    hold on
    scatter(idx.*ones(size(myYEOTs,1)),myYEOTs(:,idx),25,'k','filled')
end

ylim([min(EOTyears) max(EOTyears)+10])
xticklabels({'Mean uniform', 'Intensified uniform', 'Mean targeted', 'Intensified targeted'})
if reg == "FO"
    ylabel('Year')
end
if reg == "BE"
    title('Expected YEOT')
end
grid on

saveas(figure(2),strcat("YEOT_", reg, ".png"))