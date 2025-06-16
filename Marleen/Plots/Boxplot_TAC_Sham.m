% import TACstats
% TACstats.Mouse = cellstr(TACstats.Mouse);
% Make sure you give the table (Results) and y_variable. It takes "Group"
% in table automatically to group the boxplot so make sure the table
% contains a column labeled "Group"

function [fk] = Boxplot_TAC_Sham(Results, y_variable, x_variable, astile, tilenr)
% (Grouping, groups, overviewtable, ylimvalues, x, y)

if ~exist('astile', 'var') || astile == 0
    fk = figure('InvertHardcopy','off','Color',[1 1 1]);
elseif ~exist('tilenr', 'var')
    nexttile;
    fk = 1;
else
    nexttile(tilenr);
    fk = 1;
end

if ~exist('x_variable','var') || isempty(x_variable)
    Results.x = categorical(repmat({'x'}, size(Results,1),1));
else
    temp = find(matches(Results.Properties.VariableNames, x_variable));
    Results.Properties.VariableNames{temp} = 'x';
end

temp = find(matches(Results.Properties.VariableNames, y_variable));
Results.Properties.VariableNames{temp} = 'y';

Results.idx = grp2idx(Results.x); %this gives the groups based on alphabet, so sort the ROI labels as well:
labels = cellstr(unique(Results.x))';

groups = cellstr(unique(Results.Group));


% BOXPLOT
%dont display outliers because we will do scatter that will show them
b = boxchart(Results.x, Results.y, 'GroupByColor', Results.Group, 'LineWidth', 2, 'MarkerStyle', 'none');
hold on
    
% SCATTER
% give scatter index
xaxisstep = 1/length(groups);
xaxisplacement = 1 - 0.5*xaxisstep - (length(groups)/2-1)*xaxisstep - xaxisstep;

for ind_xgroups = 1:length(labels)
    currentxgroup = Results.x == labels{ind_xgroups};

    for indgroup = 1:length(groups)
        currentgroup = Results.Group == groups{indgroup};
        currentindex = currentgroup.*currentxgroup;
        Results.idx(currentindex == 1) = xaxisplacement + xaxisstep*indgroup;
    end
    xaxisplacement = xaxisplacement + 1;
end
clear ind_xgroups indgroup currentgroup currentxgroup currentindex ind indmouse tablemouse xaxisplacement group idx Mousegroup toplot

% go per mouse to get different markers
if any('Markertype' == string(Results.Properties.VariableNames))
    for indmouse = 1:size(Results,1)
        hs = scatter(Results.idx(indmouse), Results.y(indmouse), 70, "black", Results.Markertype{indmouse}, 'jitter','on','JitterAmount',0.05, 'LineWidth', 1);
    end
else
    hs = scatter(Results.idx, Results.y, 70, 'filled', 'jitter','on','JitterAmount',0.02);
end

% MAKE PRETTY
% hs.MarkerFaceColor = [0 0 0];
% hs.MarkerFaceAlpha = 0.3;
% % xticks(1:length(labels));
% % xticklabels(labels);
% % xlim([0.2 length(labels)+0.7])
% 
% axes1 = b.Parent;
% hold(axes1,'on');
% set(axes1,'FontSize',20,'FontWeight','bold','LineWidth',2);
legend(b, 'Location', 'northwest', 'NumColumns', 2);
%     ylim([8 12]);

% b(1).SeriesIndex = 7;
% b(2).SeriesIndex = 1;
% if size(b,1)>2
%     b(3).SeriesIndex = 2;
%     b(4).SeriesIndex = 6;
% end

t = b; %to give return

end
