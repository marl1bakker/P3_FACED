% options for variables are same as Results tabs, so Diameter,
% MeanVelocity, Pulsatility, Depth
% this is per acquisition. Things like nr of microbleeds and vessel wall
% thickness are per mouse. Make other function for that later. 

function [fk] = Scatter_TAC_Sham(Results, x_variable, y_variable, astile, tilenr)

if ~exist('astile', 'var') || astile == 0
    fk = figure;
elseif ~exist('tilenr', 'var')
    nexttile;
    fk = 1;
else
    nexttile(tilenr);
    fk = 1;
end


% % set up
Mice = unique(Results.Mouse);
Results.plotcolour(1) = {''};
Results(Results.Group=='TAC','plotcolour') = cellstr('red');
Results(Results.Group=='Sham','plotcolour') = cellstr('blue');

% go per mouse for different markers
for indmouse = 1:length(Mice) 

    MouseTable = Results(matches(Results.Mouse,Mice{indmouse}),:);
    MouseTable = MouseTable(MouseTable.UseAcq=="Yes",:);

    if MouseTable.Markerfilled{1} == 'Y'
        scatter(MouseTable.(x_variable), MouseTable.(y_variable), MouseTable.plotcolour{1}, 'Marker',MouseTable.Markertype{1}, 'filled', 'LineWidth',2);        
    else
        scatter(MouseTable.(x_variable), MouseTable.(y_variable), MouseTable.plotcolour{1}, 'Marker',MouseTable.Markertype{1}, 'LineWidth',2);
    end
    hold on
end

% get right labels
if contains(x_variable, 'velocity', 'IgnoreCase', true)
    xlabel('Velocity (mm/s)');
elseif contains(x_variable, 'diameter', 'IgnoreCase', true)
    xlabel('Diameter (um)');
else
    xlabel(x_variable);
end

if contains(y_variable, 'velocity', 'IgnoreCase', true)
    ylabel('Velocity (mm/s)');
elseif contains(y_variable, 'diameter', 'IgnoreCase', true)
    ylabel('Diameter (um)');
else
    ylabel(y_variable);
end


end