function Plot_Plunging(DataFolder)
% Only necessary for 2D scans. For a linescan, you do the same as the
% kymograph plotting.

if ~strcmp(DataFolder(end), filesep)
    DataFolder = [DataFolder filesep];
end

load([DataFolder 'AcqInfos.mat'], 'AcqInfoStream');

%% Linescan
if AcqInfoStream.height <= 1
    disp(['This acquisition is a linescan, so the corresponding kymograph-like', ...
        'figure made with Plot_Kymograph. Exited function.'])
    return
elseif ~contains(AcqInfoStream.Vessel, 'Plunging', 'IgnoreCase', true) ...
        && ~contains(AcqInfoStream.Comments, 'Plunging', 'IgnoreCase', true)
    disp('This is not a plunging vessel. Exited function.')
    return
end

%% 2D scan
% get data
fid = fopen([DataFolder 'faced.dat']);
dat = fread(fid, inf, '*single');
fclose(fid);

load([DataFolder 'morphology_img.mat'], 'av_image')
dat = reshape(dat, size(av_image,1), size(av_image,2), []);

num_of_im = 60;
im_per_row = 10;
num_rows = ceil(num_of_im/im_per_row);

f1 = figure;
imagesc(mean(dat, 3))
lims = clim;
newlims = [lims(1)+(0.5*(lims(2)-lims(1))), lims(2)]; % to see better
close(f1)

f1 = figure;
t = tiledlayout(num_rows,im_per_row, 'TileSpacing','tight');
colormap('gray')

for ind = 1:num_of_im
    nexttile;
    imagesc(dat(:,:,ind), newlims)
    axis('off')
end

title(t, [AcqInfoStream.Mouse ' ' AcqInfoStream.DatasetName ' Individual Frames'])


%% get fluorescent signal
Acq_msec = size(dat,3)/AcqInfoStream.FrameRateHz * 1000;
time = linspace(1, Acq_msec, size(dat, 3));

temp = reshape(dat, [], size(dat, 3));
fluo_signal = mean(temp, 1);
smoothed_sign = movmean(fluo_signal, 5);
% TF = islocalmin(smoothed_sign);

f2 = figure;
t = tiledlayout('Vertical');

nexttile
plot(time, fluo_signal)
hold on
plot(time, smoothed_sign)
% plot(time(TF), smoothed_sign(TF), '^')

xlabel('Time (ms)')
title([AcqInfoStream.Mouse ' ' AcqInfoStream.DatasetName ' Fluorescence signal - Whole Acq'])

nexttile
sec_to_plot = 0.3;
plotlength = round(AcqInfoStream.FrameRateHz*sec_to_plot);

% plot(time(1:plotlength), fluo_signal(1:plotlength))
hold on
plot(time(1:plotlength), smoothed_sign(1:plotlength))

TF = islocalmin(smoothed_sign(1:plotlength));
plot(time(TF), smoothed_sign(TF), '^', 'Color', 'black')

[~, inddips] = findpeaks(smoothed_sign(1:plotlength).*-1, 'MinPeakDistance', 3);
plot(time(inddips), smoothed_sign(inddips), '^', 'Color', 'red')

xlabel('Time (ms)')
title([AcqInfoStream.Mouse ' ' AcqInfoStream.DatasetName ' Fluorescence signal - ' num2str(sec_to_plot) ' seconds'])

nexttile
sec_to_plot = 0.1;
plotlength = round(AcqInfoStream.FrameRateHz*sec_to_plot);
% plot(time(1:plotlength), fluo_signal(1:plotlength))
hold on
plot(time(1:plotlength), smoothed_sign(1:plotlength))

[~, inddips] = findpeaks(smoothed_sign(1:plotlength).*-1, 'MinPeakDistance', 3);
plot(time(inddips), smoothed_sign(inddips), '^', 'Color', 'red')

xlabel('Time (ms)')
title([AcqInfoStream.Mouse ' ' AcqInfoStream.DatasetName ' Fluorescence signal - ' num2str(sec_to_plot) ' seconds'])


end





