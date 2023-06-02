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
% BOFFA EAST
load('./Result2/GIN/F1_BoffaEast/Output_MCMC_M5_GIN_F1_BoffaEast_DataGIN22_Paras001.mat');
% BOFFA WEST
%load('./Result3/GIN/F2_BoffaWest/Output_MCMC_M5_GIN_F2_BoffaWest_DataGIN22_Paras001.mat');
% DUBREKA
%load('./Result/GIN/F3_Dubreka/Output_MCMC_M5_GIN_F3_Dubreka_DataGIN22_Paras001.mat');
% FORECARIAH
%load('./Result3/GIN/F4_Forecariah/Output_MCMC_M5_GIN_F4_Forecariah_DataGIN22_Paras001.mat');
n = round(size(posterior,1)/2);

names = ["r^{FW}", "r^{FP}", "r^{MY}", "r^{MW}", "r^{MP}", "u"];

f = figure(1);
f.Position = [205 450 600 600];
semilogy([1:n],neg_log_likelihood(1:n),'color',myblue)
hold on
semilogy([n+1:2*n],neg_log_likelihood(n+1:2*n),'color',myred)
xlabel('Iterations')
title('Log Posterior')
legend('Chain 1','Chain 2')
grid on

f = figure(2);
f.Position = [205 450 1200 600];
for i = 2:6
    subplot(2,3,i-1)
    plot([1:n],posterior(1:n,i),'color',myblue)
    hold on
    plot([n+1:2*n],posterior(n+1:2*n,i),'color',myred)
    if i > 4
        xlabel('Iterations')
    end
    title(names{i-1})
    %legend('Chain 1','Chain 2')
    grid on
end

% u
subplot(2,3,6)
plot([1:n],posterior(1:n,11),'color',myblue)
hold on
plot([n+1:2*n],posterior(n+1:2*n,11),'color',myred)
xlabel('Iterations')
title(fitted_para_names{11})
%legend('Chain 1','Chain 2')
grid on

% save
%saveas(figure(1),strcat("BE_","AgeLogPost.png"))
%saveas(figure(2),strcat("BE_","AgeTraces.png"))