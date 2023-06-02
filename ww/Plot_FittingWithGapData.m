
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

Cloc = 'CIV'; % options: 'CIV', 'DRC', 'GIN', 'TCD', 'UGA'
Lv1loc = 2;
Lv2loc = 0;
Lv3loc = 0;
Model = 'M4';
DataStr = '18'; % full data
ParaStr = '002';

Y1 = 19; % number of years with data
GapY = 3;
dGY = 3;
firstyear = 2000;

% pos1 = [0.2:0.4:1.8];
% year2 = repmat(cellstr(num2str([firstyear + Y0 + Y1 : firstyear + Y0 + Y1 + Y2 - 1]'))', 1, NumStrat);
% pos2 = sort([2.2:2.2+Y2-1 2.2+0.3:ceil(2.2+Y2-1) 2.2+0.6:ceil(2.2+Y2-1)]);
xmax = Y1 + GapY;
xall = 0.5 : xmax - 0.5;
yearall = firstyear : firstyear + Y1 + GapY - 1;
xpos = xall(mod(yearall, 5) == 0);%[0.4 0.8 3 4 8 9 13 14 18 19 23 24 28 29]; % need to modify manually
xlabels = num2cell(yearall(mod(yearall, 5) == 0));%{'2000','','       2020','','       2025','','       2030',''};
xminorpos = (0 : Y1 + GapY) + 0.5;%xall(mod(yearall, 5) > 1);%[0 1.2 1.6 2 5 6 7 10 11 12 15 16 17 20 21 22 25 26 27];

MyBlue = [67 147 195]/255;

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
Active = [ACTIVE1(location, :) + ACTIVE2(location, :) + ACTIVENa(location, :) GapActive(location, :)];
Passive = [PASSIVE1(location, :) + PASSIVE2(location, :) + PASSIVENa(location, :) GapPassive(location, :)];
Screened = [SCREENED(location, :) GapScreening{location, 1:GapY-dGY}];

% Fitted part
load([Dir, 'Fitted_', Model, '_', LocStr, '_ID', Cloc, DataStr, '-', ParaStr, '.mat']);
ActiveQ1 = quantile(SampledActive1 + SampledActive2, [0.025 0.25 0.5 0.5 0.75 0.975]);
PassiveQ1 = quantile(SampledPassive1 + SampledPassive2, [0.025 0.25 0.5 0.5 0.75 0.975]);


        
% Figure setting
figure('Name', ['Fitting_', LocStr], 'NumberTitle', 'off')

% ===== Screening =====
subplot(3, 1, 1)
stairs(0 : Y1 + GapY - dGY, [Screened Screened(end)], 'LineWidth', 1, 'Color', 'k')

plt = gca;
set(plt, 'children', flipud(get(gca, 'children')))
xticks(xpos)
xticklabels(xlabels)
plt.XRuler.TickLabelGapOffset = -2;
plt.XAxis.TickDirection = 'out';
plt.XAxis.TickLength = [0.018 1];
plt.YAxis.TickLength = [0.012 1];
plt.XAxis.MinorTick = 'on';
plt.XAxis.MinorTickValues = xminorpos;
box off
    
xlim([0, Y1])
ylim([0, 1.15 * max(Screened)])
    
xlabel('Year')
ylabel('# screened')
title(CCLOC{location}, 'FontSize', 12)
ax2 = axes('Position', get(plt, 'Position'), 'Color', 'None', 'LineWidth', 0.5, 'XAxisLocation', 'top', 'XTick', [], 'YAxisLocation', 'right', 'YTick', []);

        
% ===== Active =====
subplot(3, 1, 2)
stairs(0 : Y1 + GapY - dGY, [Active Active(end)], 'LineWidth', 1, 'Color', 'k')
hold on
ActiveBox1 = boxplot(ActiveQ1, 'symbol', '', 'Colors', MyBlue, 'positions', 0.5 : Y1 + GapY -0.5);
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
xticks(xpos)
xticklabels(xlabels)
plt.XRuler.TickLabelGapOffset = -2;
plt.XAxis.TickDirection = 'out';
plt.XAxis.TickLength = [0.018 1];
plt.YAxis.TickLength = [0.012 1];
plt.XAxis.MinorTick = 'on';
plt.XAxis.MinorTickValues = xminorpos;
box off
    
xlim([0, Y1])
ylim([0, 1.15 * max([Active ActiveQ1(6,:)])])
% xlim([10, Y1])
% ylim([0, 1.15 * max(ActiveQ1(6,11:end))])
    
xlabel('Year')
ylabel('Active cases')
ax2 = axes('Position', get(plt, 'Position'), 'Color', 'None', 'LineWidth', 0.5, 'XAxisLocation', 'top', 'XTick', [], 'YAxisLocation', 'right', 'YTick', []);


% ===== Passive =====
subplot(3, 1, 3)
stairs(0 : Y1 + GapY - dGY, [Passive Passive(end)], 'LineWidth', 1, 'Color', 'k')
hold on
PassiveBox1 = boxplot(PassiveQ1, 'symbol', '', 'Colors', MyBlue, 'positions', 0.5 : Y1 + GapY -0.5);
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
xticks(xpos)
xticklabels(xlabels)
plt.XRuler.TickLabelGapOffset = -2;
plt.XAxis.TickDirection = 'out';
plt.XAxis.TickLength = [0.018 1];
plt.YAxis.TickLength = [0.012 1];
plt.XAxis.MinorTick = 'on';
plt.XAxis.MinorTickValues = xminorpos;
box off
    
xlim([0, Y1])
ylim([0, 1.15 * max([Passive PassiveQ1(6,:)])])
% xlim([10, Y1])
% ylim([0, 1.15 * max(PassiveQ1(6,11:end))])
    
xlabel('Year')
ylabel('Passive cases')
ax2 = axes('Position', get(plt, 'Position'), 'Color', 'None', 'LineWidth', 0.5, 'XAxisLocation', 'top', 'XTick', [], 'YAxisLocation', 'right', 'YTick', []);

        
set(gcf,'Units','Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [8.5, 8.5/pos(3)*pos(4)])        
        
     