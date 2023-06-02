clear all

% Plotting preferences
set(0,'defaultlinelinewidth',2)
set(0,'defaultaxesfontsize',16)
set(0,'DefaultLegendFontSize',18)
set(0,'DefaultAxesTitleFontWeight','normal')

myblue = [0.2627 0.5765 0.7647];
mygreen = [0.1059 0.4706 0.2157];

% read data
% BOFFA EAST
BEloc = './Result/GIN/F1_BoffaEast/Posterior_MCMC_M5_GIN_F1_BoffaEast_DataGIN22_Paras001.mat';
% BOFFA WEST
BWloc = './Result/GIN/F2_BoffaWest/Posterior_MCMC_M5_GIN_F2_BoffaWest_DataGIN22_Paras001.mat';
% DUBREKA
DUloc = './Result/GIN/F3_Dubreka/Posterior_MCMC_M5_GIN_F3_Dubreka_DataGIN22_Paras001.mat';
% FORECARIAH
FOloc = './Result/GIN/F4_Forecariah/Posterior_MCMC_M5_GIN_F4_Forecariah_DataGIN22_Paras001.mat';

loc = BEloc;

load('FitPars.mat')
load(BEloc);

% use name of MCMC run as filename for plot images
[~,name,~] = fileparts(BEloc);
string = strcat(name,'.png');

% matrix of posteriors for ease of plotting
n = 1000;
post = [Posterior.u Posterior.rFW Posterior.rFP Posterior.rMY Posterior.rMW Posterior.rMP];
%post = post(n+1:2*n,:);
pars = ["u", "r^{FW}", "r^{FP}", "r^{MY}", "r^{MW}", "r^{MP}"];
nbins = 40;

f = figure(1);
f.Position = [105 430 600 800];
%sgtitle(strrep(name, '_', ' '))

for i = 1:6
    subplot(3,2,i)
    %yyaxis left
    histogram(post(:,i),nbins,'FaceColor',myblue,'FaceAlpha',0.9,Normalization='pdf')
    if mod(i,2) == 1
        ylabel('PDF')
    end
    xlim([0 max(post(:,i))])
    set(gca,'YColor','k')
    hold on
    %yyaxis right
    if i > 1
        pdf = gampdf([0:0.01:20]-FittedParameters.Lower(i+3),FittedParameters.Parameters{i+3}(1),FittedParameters.Parameters{i+3}(2));
        plot([0:0.01:20],pdf,'r','LineWidth',2)
    else
        pdf = betapdf([0:0.01:1],FittedParameters.Parameters{30}(1),FittedParameters.Parameters{30}(2));
        plot([0:0.01:1],pdf,'r','LineWidth',2)
    end
    %set(gca,'YColor','k','YTicklabel',[])
    xlabel(pars{i})
end


% matrix of posteriors for ease of plotting
% R0, gamma_H0, gamma_H, eta_H, specificity
%post2 = [Posterior.R0 Posterior.gamma_H0 Posterior.gamma_H Posterior.eta_H Posterior.specificity];
%pars2 = ["R0", "gamma H0", "gamma H", "eta H", "specificity"];


%f2 = figure(2);
%f2.Position = [1005 430 800 600];
%sgtitle(strrep(name, '_', ' '))
%histogram(Posterior.u,nbins,'Normalization','probability')
%xlim([min(Posterior.u) max(Posterior.u)])
%set(gca,'YColor','k')
%hold on
%yyaxis right
%if i == 14
%    plot([1:0.001:1.1],exppdf([1:0.001:1.1]-FittedParameters.Lower(1),FittedParameters.Parameters{1}(1)),'r','LineWidth',2)
%end
%title('u')

% save
saveas(figure(1),strcat("AGEFIT_",string))

% save
%saveas(figure(2),strcat("2",string))