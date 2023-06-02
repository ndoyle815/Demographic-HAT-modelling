
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                     %
%   This code plots fitting results of the Warwick HAT model                                          %
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

Cloc = 'GIN'; % options: 'CIV', 'DRC', 'GIN', 'TCD', 'UGA'
Lv1loc = 1;
Lv2loc = 0;
Lv3loc = 0;
Model = 'M4';
DataStr = '18'; % full data
%ID = {'13-000', '13-101', '19-005'};
ParaStr = '006';
DispNames = {'Original model fit (2000-2013)', 'Updated model fit (2000-2013)', 'Updated model fit (2000-2019)'};
pos = [[0.1 0.8 0.854 0.165]; [0.107 0.571 0.843 0.182]; [0.107 0.347 0.843 0.183]; [0.1 0.125 0.85 0.18]];

Y1 = 19; % number of years with data
% Y1 = 0;
% Y = 14; % years without VC

MyBlue = [67 147 195]/255; % 80%
MyLightBlue = [103 169 207]/255; %[146 197 222]/255; % 60%
MyDarkBlue = [5 48 97]/255; % 90%
MyLightPurple = [175 141 195]/255;
MyGreen = [27 120 55]/255;
MyGrey = [0.6 0.6 0.6];
MyBlack = [0 0 0];
MyWhite = [1 1 1];

switch Cloc
    case {'CIV', 'GIN', 'TCD'}
        if Lv1loc ~= 0 % focus level simulation
            load(['../Data/', Cloc, DataStr, '/Data.mat']);
            f.name = strcat('F', num2str(Lv1loc), '_', CCLOC{Lv1loc});
            names = {f.name, Cloc};
            location = Lv1loc;
            Dir = ['../ResultRon/', Cloc, '/', f.name, '/'];
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
IDStr = ['_ID', Cloc, DataStr, '-', ParaStr];

% Data
Active = ACTIVE1(location, :) + ACTIVE2(location, :) + ACTIVENa(location, :);
Passive = PASSIVE1(location, :) + PASSIVE2(location, :) + PASSIVENa(location, :);
Screened = SCREENED(location, :);

% Fitted part
load([Dir, 'Fitted_', Model, '_', LocStr, IDStr, '.mat']);
ActiveQ1 = quantile(SampledActive1 + SampledActive2, [0.025 0.25 0.5 0.5 0.75 0.975]);
PassiveQ1 = quantile(SampledPassive1 + SampledPassive2, [0.025 0.25 0.5 0.5 0.75 0.975]);


        
% Figure setting
figure('Name', ['Fitting_', LocStr], 'NumberTitle', 'off')

% ===== Screening =====
subplot(3, 1, 1)
stairs(0 : Y1, [Screened Screened(end)], 'LineWidth', 1, 'Color', 'k')

plt = gca;
set(plt, 'children', flipud(get(gca, 'children')))
xticks(0.5 : 5 : Y1+0.5)
xticklabels({'2000', '2005', '2010', '2015'})
plt.XRuler.TickLabelGapOffset = -2;
plt.XAxis.TickDirection = 'out';
plt.XAxis.TickLength = [0.018 1];
plt.YAxis.TickLength = [0.012 1];
plt.XAxis.MinorTick = 'on';
plt.XAxis.MinorTickValues = (0 : Y1) + 0.5;
box off
    
xlim([0, Y1])
ylim([0, 1.15 * max(Screened)])
    
xlabel('Year')
ylabel('# screened')
title(CCLOC{location}, 'FontSize', 12)
ax2 = axes('Position', get(plt, 'Position'), 'Color', 'None', 'LineWidth', 0.5, 'XAxisLocation', 'top', 'XTick', [], 'YAxisLocation', 'right', 'YTick', []);

        
% ===== Active =====
subplot(3, 1, 2)
stairs(0 : Y1, [Active Active(end)], 'LineWidth', 1, 'Color', 'k')
hold on
ActiveBox1 = boxplot(ActiveQ1, 'symbol', '', 'Colors', MyBlue, 'positions', 0.5 : Y1-0.5);
set(ActiveBox1, {'linew'}, {1})

%get handles for boxplot
uw1 = findobj(ActiveBox1, 'tag', 'Upper Whisker');   % get handle to "Upper Whisker" line
uav1 = findobj(ActiveBox1, 'tag', 'Upper Adjacent Value');   %get handle to "Upper Adjacent Value" line
lw1 = findobj(ActiveBox1, 'tag', 'Lower Whisker');   % get handle to "Lower Whisker" line
lav1 = findobj(ActiveBox1, 'tag', 'Lower Adjacent Value');   %get handle to "Lower Adjacent Value" line
m1 = findobj(ActiveBox1, 'tag', 'Median');   %get handle to "Median" line
out1 = findobj(ActiveBox1, 'tag', 'Outliers');   %get handle to outliers
box1 = findobj(ActiveBox1, 'tag', 'Box');   %get handle to box

for i = 1 : Y1
    %Ensure whiskers are at 97.5% and 2.5% give solid whiskers
    uw1(i).YData(:) = [ActiveQ1(5,i) ActiveQ1(6,i)];
    uw1(i).LineStyle = '-';
    uav1(i).YData(:) = [ActiveQ1(6,i) ActiveQ1(6,i)];
    lw1(i).YData(:) = [ActiveQ1(1,i) ActiveQ1(2,i)];
    lw1(i).LineStyle = '-';
    lav1(i).YData(:) = [ActiveQ1(1,i) ActiveQ1(1,i)];

    %Fill box
    patch(get(box1(i), 'XData'), get(box1(i), 'YData'), MyBlue);
    m1(i).LineWidth = 1.5;
    m1(i).Color = 'w';
end

plt = gca;
set(plt, 'children', flipud(get(gca, 'children')))
xticks(0.5 : 5 : Y1+0.5)
xticklabels({'2000', '2005', '2010', '2015'})
plt.XRuler.TickLabelGapOffset = -2;
plt.XAxis.TickDirection = 'out';
plt.XAxis.TickLength = [0.018 1];
plt.YAxis.TickLength = [0.012 1];
plt.XAxis.MinorTick = 'on';
plt.XAxis.MinorTickValues = (0 : Y1) + 0.5;
box off
    
xlim([0, Y1])
ylim([0, 1.15 * max([Active ActiveQ1(6,:)])])
    
xlabel('Year')
ylabel('Active cases')
ax2 = axes('Position', get(plt, 'Position'), 'Color', 'None', 'LineWidth', 0.5, 'XAxisLocation', 'top', 'XTick', [], 'YAxisLocation', 'right', 'YTick', []);


% ===== Passive =====
subplot(3, 1, 3)
stairs(0 : Y1, [Passive Passive(end)], 'LineWidth', 1, 'Color', 'k')
hold on
PassiveBox1 = boxplot(PassiveQ1, 'symbol', '', 'Colors', MyBlue, 'positions', 0.5 : Y1-0.5);
set(PassiveBox1, {'linew'}, {1})

%get handles for boxplot
uw1 = findobj(PassiveBox1, 'tag', 'Upper Whisker');   % get handle to "Upper Whisker" line
uav1 = findobj(PassiveBox1, 'tag', 'Upper Adjacent Value');   %get handle to "Upper Adjacent Value" line
lw1 = findobj(PassiveBox1, 'tag', 'Lower Whisker');   % get handle to "Lower Whisker" line
lav1 = findobj(PassiveBox1, 'tag', 'Lower Adjacent Value');   %get handle to "Lower Adjacent Value" line
m1 = findobj(PassiveBox1, 'tag', 'Median');   %get handle to "Median" line
out1 = findobj(PassiveBox1, 'tag', 'Outliers');   %get handle to outliers
box1 = findobj(PassiveBox1, 'tag', 'Box');   %get handle to box

for i = 1 : Y1
    %Ensure whiskers are at 97.5% and 2.5% give solid whiskers
    uw1(i).YData(:) = [PassiveQ1(5,i) PassiveQ1(6,i)];
    uw1(i).LineStyle = '-';
    uav1(i).YData(:) = [PassiveQ1(6,i) PassiveQ1(6,i)];
    lw1(i).YData(:) = [PassiveQ1(1,i) PassiveQ1(2,i)];
    lw1(i).LineStyle = '-';
    lav1(i).YData(:) = [PassiveQ1(1,i) PassiveQ1(1,i)];

    %Fill box
    patch(get(box1(i), 'XData'), get(box1(i), 'YData'), MyBlue);
    m1(i).LineWidth = 1.5;
    m1(i).Color = 'w';
end

plt = gca;
set(plt, 'children', flipud(get(gca, 'children')))
xticks(0.5 : 5 : Y1+0.5)
xticklabels({'2000', '2005', '2010', '2015'})
plt.XRuler.TickLabelGapOffset = -2;
plt.XAxis.TickDirection = 'out';
plt.XAxis.TickLength = [0.018 1];
plt.YAxis.TickLength = [0.012 1];
plt.XAxis.MinorTick = 'on';
plt.XAxis.MinorTickValues = (0 : Y1) + 0.5;
box off
    
xlim([0, Y1])
ylim([0, 1.15 * max([Passive PassiveQ1(6,:)])])
    
xlabel('Year')
ylabel('Passive cases')
ax2 = axes('Position', get(plt, 'Position'), 'Color', 'None', 'LineWidth', 0.5, 'XAxisLocation', 'top', 'XTick', [], 'YAxisLocation', 'right', 'YTick', []);

        
        
        
     