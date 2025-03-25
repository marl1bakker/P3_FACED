
function [fk] = Plot_Kymograph(kymoImg, frmRate, pixelsize, secs_to_plot, velocity, astile, tilenr)

% setup
if ~exist('pixelsize', 'var') || isempty(pixelsize)
    pixelsize = 1;
    label_y = 'micron (estimated)';
else
    label_y = 'micron';
end
if ~exist('secs_to_plot', 'var')
    secs_to_plot = 1;
end
if ~exist('astile', 'var') || astile == 0
    fk = figure;
elseif ~exist('tilenr', 'var')
    nexttile;
    fk = 1;
else
    nexttile(tilenr);
    fk = 1;
end
if exist('velocity', 'var') && isempty(velocity)
    clear velocity
end

if size(kymoImg, 2)<size(kymoImg,1) % make sure second axis is time
    kymoImg = kymoImg';
end

% if youre plotting velocity as well, make left and right y axis
if exist('velocity', 'var')
    yyaxis left
end

% Plot figure
imagesc(kymoImg(:,1:round(frmRate*secs_to_plot)));
colormap('gray');

% get right x axis
% marker every 1000/stepnr_ms milliseconds
if secs_to_plot <= 0.05 % marker every millisecond
    stepnr_ms = 1000;
elseif    secs_to_plot <= 0.1 % marker every 10 milliseconds
    stepnr_ms = 100;
elseif secs_to_plot <= 0.5 % marker every 20 milliseconds
    stepnr_ms = 50;
elseif secs_to_plot <=1 % marker every 100 milliseconds
    stepnr_ms = 10;
else % marker every 200 ms
    stepnr_ms = 5;
end
Acq_msec = size(kymoImg, 2)/frmRate*1000;
xticks(0:frmRate/stepnr_ms:size(kymoImg,2));
xticklabels(0:1000/stepnr_ms:Acq_msec);
xlabel('Milliseconds');

% get right y axis
ROIlength = pixelsize*size(kymoImg,1);

if ROIlength <= 10
    stepsize_um = 5;
elseif ROIlength <= 50
    stepsize_um = 10;
else
    stepsize_um = 25;
end

yticks(0:1/pixelsize*stepsize_um:size(kymoImg,1))
yticklabels(1:stepsize_um:pixelsize*size(kymoImg,1))
ylabel(label_y)

if exist('velocity', 'var')
    yyaxis right
    plot(velocity, 'Color', 'red', 'LineWidth', 2); ylabel('velocity mm/s')
end

end

