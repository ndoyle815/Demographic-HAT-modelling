%Kat box plot

%Take model output data for boxes and actual real-world data for lines
%MyColor is the colour of boxes
%pos, xpos, xlabels and xminor pos are about position of boxes, ticks and
%labels
function k = Katboxplot(BoxData, RealData, ShowYears, MyColor, pos, xpos, xlabels, xminorpos, width)

if nargin==8
    Box = boxplot(BoxData, 'symbol', '', 'Colors', MyColor, 'positions', pos);
elseif nargin==9
    Box = boxplot(BoxData, 'symbol', '', 'Colors', MyColor, 'positions', pos,'widths',width);
end
set(Box, {'linew'}, {1.5})

%get handles for boxplot
uw1 = findobj(Box, 'tag', 'Upper Whisker');   % get handle to "Upper Whisker" line
uav1 = findobj(Box, 'tag', 'Upper Adjacent Value');   %get handle to "Upper Adjacent Value" line
lw1 = findobj(Box, 'tag', 'Lower Whisker');   % get handle to "Lower Whisker" line
lav1 = findobj(Box, 'tag', 'Lower Adjacent Value');   %get handle to "Lower Adjacent Value" line
m1 = findobj(Box, 'tag', 'Median');   %get handle to "Median" line
out1 = findobj(Box, 'tag', 'Outliers');   %get handle to outliers
box1 = findobj(Box, 'tag', 'Box');   %get handle to box

for i = 1 : size(BoxData,2)
    %Ensure whiskers are at 97.5% and 2.5% give solid whiskers
    uw1(i).YData(:) = [BoxData(5,i) BoxData(6,i)];
    uw1(i).LineStyle = '-';
    uav1(i).YData(:) = [BoxData(6,i) BoxData(6,i)];
    lw1(i).YData(:) = [BoxData(1,i) BoxData(2,i)];
    lw1(i).LineStyle = '-';
    lav1(i).YData(:) = [BoxData(1,i) BoxData(1,i)];

    %Fill box
    k = patch(get(box1(i), 'XData'), get(box1(i), 'YData'), MyColor,'FaceAlpha',.4);
    m1(i).LineWidth = 1.5;
    m1(i).Color = 'k';
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
    
xlim([0, ShowYears])
ylim([0, 1.15 * max([RealData BoxData(6,:)])])