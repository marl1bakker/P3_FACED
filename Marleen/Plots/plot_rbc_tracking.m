% Example_plots_single_acqs

function plot_rbc_tracking(DataFolder, give_input)
% Fig 1D Meng

% set up
if ~exist('DataFolder', 'var')
    DataFolder = 'D:\FACED\Data P3\M04\A20\';
elseif ~strcmp(DataFolder(end), filesep)
    DataFolder = [DataFolder filesep];
end

if ~exist('give_input', 'var')
    give_input = 0;
end

load([DataFolder 'AcqInfos.mat'], 'AcqInfoStream');
if ~isfield(AcqInfoStream, 'Coupled_acq')
    AddInfo(DataFolder, 'Coupled_acq');
    load([DataFolder 'AcqInfos.mat'], 'AcqInfoStream')
end

if AcqInfoStream.height < 2
    gave = 'line';
elseif AcqInfoStream.height >= 2
    gave = 'twoD';
end

coupled_acq_list = AcqInfoStream.Coupled_acq;
seps = strfind(DataFolder, filesep);
MouseFolder = DataFolder(1:seps(end-1));
clear AcqInfoStream

% start plotting
f1 = figure('Color', 'white');
tiledlayout(12,3, 'TileSpacing','tight', 'Padding', 'compact');

%% plot what you gave
if matches(gave, 'line')
    plot_rbc_tracking_linescan(DataFolder, give_input);
elseif matches(gave, 'twoD')
    plot_rbc_tracking_2dscan(DataFolder, give_input);
end


%% Get matching acq
if length(coupled_acq_list)>1 && matches(gave, 'line')
    best_c_acq = choose_best_kymo(coupled_acq_list, MouseFolder, 'twoD');
elseif length(coupled_acq_list)>1 && matches(gave, 'twoD')
    best_c_acq = choose_best_kymo(coupled_acq_list, MouseFolder, 'line');
else
    best_c_acq = coupled_acq_list{1};
end


%% plot coupled acq
if matches(gave, 'line')
    plot_rbc_tracking_2dscan([MouseFolder best_c_acq], give_input);
elseif matches(gave, 'twoD')
    plot_rbc_tracking_linescan([MouseFolder best_c_acq], give_input);
end


%% save
f1 = gcf;
% f1.Position = [1.8000    1.0000  591.2000  948.0000];
f1.Position = [-1.3342e+03 125 591.2 948]; % tweede scherm

saveas(f1, [DataFolder 'RBC_tracking_methods.svg'])

end


function plot_rbc_tracking_2dscan(DataFolder, give_input)
if ~exist('give_input', 'var')
    give_input = 0;
end
if ~strcmp(DataFolder(end), filesep)
    DataFolder = [DataFolder filesep];
end

% get data
gcf;
load([DataFolder 'AcqInfos.mat'], 'AcqInfoStream');
load([DataFolder 'morphology_img.mat'], 'av_image');
proportion_XY_ax = [1 AcqInfoStream.pxl_sz_x/AcqInfoStream.pxl_sz_y_avged 1];
frmRate = AcqInfoStream.FrameRateHz;
fid = fopen([DataFolder 'faced.dat']);
dat = fread(fid, inf, '*single');
fclose(fid);
dat = reshape(dat, size(av_image,1), size(av_image,2), []);


% Average image
nexttile(31, [2 2]); imagesc(av_image); title('average image')
ax = gca; ax.DataAspectRatio = proportion_XY_ax;
colormap('gray')
av_im_um_x = AcqInfoStream.width;
xticks(0:size(av_image,2)/av_im_um_x*10:size(av_image,2));
xticklabels(0:10:av_im_um_x);
xlabel('um');
av_im_um_y = AcqInfoStream.height;
yticks(0:size(av_image,1)/av_im_um_y*10:size(av_image,1));
yticklabels(fliplr(0:10:av_im_um_y));
ylabel('um');

if exist([DataFolder 'kymoROI_1.mat'],'file')
    ROI_exists = 1;
    load([DataFolder 'kymoROI_1.mat'], 'kymoImg', 'ROI_info');
else
    ROI_exists = 0;
end

% track rbc over time
if ROI_exists && give_input
    % load([DataFolder 'kymoROI_1.mat'], 'kymoImg');
    f3 = figure; imagesc(kymoImg(1:round(frmRate),:)'); colormap('gray');
    opts.WindowStyle = 'normal';
    frmstart = inputdlg({'What frame to start kymo?'},'input', [1 45], {''}, opts);
    close(f3);
    frmstart = str2double(frmstart{1});
else
    frmstart = 500;
end

frmend = frmstart+(0.1*frmRate); % 10 ms between images
frmskipped = (frmend-frmstart)/10; % to have 10 frames
tileind = 1;
for ind = frmstart : frmskipped : frmend-1
    nexttile(tileind);
    imagesc(mean(dat(:,:,round(ind):round(ind+frmRate*0.01)),3))
    hold on

    if ROI_exists % make arrows for rbc
        [~, rbcind] = min(mean(kymoImg(round(ind):round(ind+frmRate*0.01),:),1));
        plot(ROI_info.Coor(1,1)+rbcind, ROI_info.Coor(1,2)-5, 'v', 'Color', 'red', 'LineWidth', 1)
    end

    axis off
    title([num2str((ind/frmRate*1000)-(frmstart/frmRate*1000)) ' ms'])
    tileind = tileind+3;
end


if ROI_exists
    % load([DataFolder 'kymoROI_1.mat'], 'kymoImg', 'ROI_info');

    % show ROI on tile
    nexttile(31)
    title('ROI')
    line(ROI_info.start_end_points(:,1), ROI_info.start_end_points(:,2), 'Color', 'red', 'LineWidth', 2)

    % show kymo
    nexttile(2, [10 1])
    lengthkymo = round(2*frmRate); % 2 seconds
    imagesc(kymoImg(frmstart:frmstart+lengthkymo,:))
    Acq_msec = lengthkymo/frmRate*1000;
    yticks(0:frmRate/10:lengthkymo);
    yticklabels(0:100:Acq_msec);
    ylabel('Milliseconds');
    length_ROI = av_im_um_x/size(av_image,2)*size(kymoImg,2);
    xticks(0:size(kymoImg,2)/length_ROI*10:size(kymoImg,2));
    xticklabels(0:10:length_ROI);
    xlabel('um');
    % Acq_sec = lengthkymo/frmRate;
    % yticks(0:frmRate/2:lengthkymo);
    % yticklabels(0:0.5:Acq_sec);
    % ylabel('Seconds');
    hold on
    title('ROI kymograph')


    % show box on kymo for frames used in 'a'

end
end



function plot_rbc_tracking_linescan(DataFolder, give_input)
if ~strcmp(DataFolder(end), filesep)
    DataFolder = [DataFolder filesep];
end
if ~exist('give_input', 'var')
    give_input = 0;
end
if ~exist([DataFolder 'kymoROI_1.mat'], 'file')
    return
end

gcf;
load([DataFolder 'AcqInfos.mat'], 'AcqInfoStream');
frmRate = AcqInfoStream.FrameRateHzLinescan;
warning('off')
load([DataFolder 'kymoROI_1.mat'], 'kymoImg', 'PixelSize', 'Velocity_calc');
warning('on')


% plot kymograph
nexttile(3, [10 1])
lengthkymo = round(2*frmRate); % 2 seconds
imagesc(kymoImg(1:round(frmRate*2),:))
Acq_msec = lengthkymo/frmRate*1000;
yticks(0:frmRate/10:lengthkymo);
yticklabels(0:100:Acq_msec);
ylabel('Milliseconds');
title('Linescan kymograph')
if ~exist('PixelSize', 'var')
    return % if you dont have pixelsize, you didnt do velocity calc, you wont have skipamt
end
scan_um = size(kymoImg, 2)*PixelSize.pxlSize;
xticks(0:20/PixelSize.pxlSize:scan_um);
xticklabels(0:20:scan_um);
xlabel('um');


% Show two frames to xcorr later
if give_input
    % plot points on which frames you xcorr
    opts.WindowStyle = 'normal';
    frame_one = inputdlg({'Which frame first for xcorr (give ms)'},'input', [1 45], {''}, opts);
    frame_one = round(str2double(frame_one)*frmRate*0.001);
else
    frame_one = 1000;
end
frame_two = frame_one+Velocity_calc.skipamt;
hold on
plot(1,frame_one, 'Marker', '>', 'Color', 'green', 'MarkerSize', 10)
line([1 size(kymoImg,2)], [frame_one frame_one], 'Color', 'green', 'MarkerSize', 2)
plot(1,frame_two, 'Marker', '>', 'Color', 'blue', 'MarkerSize', 10)
line([1 size(kymoImg,2)], [frame_two frame_two], 'Color', 'blue', 'MarkerSize', 2)


% show the xcorr of those frames
nexttile(33, [2 1])
plot(kymoImg(frame_one,:), 'Color', 'green')
hold on
plot(kymoImg(frame_two,:), 'Color', 'blue')
ylabel('Intensity (a.u.)')
xticks(0:20/PixelSize.pxlSize:scan_um);
xticklabels(0:20:scan_um);
xlabel('um linescan kymograph');

end

function [best_kymo_acq] = choose_best_kymo(coupled_acq_list, MouseFolder, scantype)
if ~exist('scantype', 'var')
    scantype = 'line';
end

% For linescan/kymographs:
if matches(scantype, 'line')
    % plot all kymographs of coupled acqs.
    f2 = figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout('Vertical');
    secs_plotted =2;
    for c_acq_ind = 1:length(coupled_acq_list)
        if exist([MouseFolder coupled_acq_list{c_acq_ind} filesep 'kymoROI_1.mat'], 'file') && ...
                exist([MouseFolder coupled_acq_list{c_acq_ind} filesep 'AcqInfos.mat'], 'file')

            load([MouseFolder coupled_acq_list{c_acq_ind} filesep 'AcqInfos.mat'], 'AcqInfoStream');
            load([MouseFolder coupled_acq_list{c_acq_ind} filesep 'kymoROI_1.mat'], 'kymoImg', 'PixelSize');
            if isfield(AcqInfoStream, 'FrameRateHzLinescan')
                frmRate = AcqInfoStream.FrameRateHzLinescan;
            else
                frmRate = AcqInfoStream.FrameRateHz;
            end

            % if you didnt get the kymo for the full timespan, pad with nans
            if size(kymoImg,1)<round(frmRate*secs_plotted)
                kymoImg = [kymoImg; NaN(round(frmRate*secs_plotted)-size(kymoImg, 1), size(kymoImg,2))];
            end

            if exist('PixelSize', 'var')
                Plot_Kymograph(kymoImg, frmRate, PixelSize.pxlSize, secs_plotted, 1);
            else
                Plot_Kymograph(kymoImg, frmRate, [], secs_plotted, 1);
            end
            title([coupled_acq_list{c_acq_ind} ' ' AcqInfoStream.Acquisition_Type]);
            clear AcqInfoStream kymoImg PixelSize

        else
            nexttile
            %keep plot emtpy
        end
    end

    clear secs_plotted c_acq_ind
    [best_kymo_ind] = listdlg('PromptString', {'Which kymograph looks the best? (Cancel if none)'},...
        'ListString', coupled_acq_list);
    close(f2)

    best_kymo_acq = coupled_acq_list{best_kymo_ind};

elseif matches(scantype, 'twoD')
    % plot all av_ims of coupled acqs.
    f2 = figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout('Flow');
    for c_acq_ind = 1:length(coupled_acq_list)
        if exist([MouseFolder coupled_acq_list{c_acq_ind} filesep 'morphology_img.mat'], 'file') && ...
                exist([MouseFolder coupled_acq_list{c_acq_ind} filesep 'AcqInfos.mat'], 'file')

            load([MouseFolder coupled_acq_list{c_acq_ind} filesep 'AcqInfos.mat'], 'AcqInfoStream');
            load([MouseFolder coupled_acq_list{c_acq_ind} filesep 'morphology_img.mat'], 'av_image');

            nexttile
            imagesc(av_image)
            colormap('gray')

            title([coupled_acq_list{c_acq_ind}]);
            clear AcqInfoStream av_image

        else
            nexttile
            %keep plot emtpy
        end
    end

    [best_kymo_ind] = listdlg('PromptString', {'Which average image looks '...
        'the best? (Cancel if none)'},...
        'ListString', coupled_acq_list);
    close(f2)

    best_kymo_acq = coupled_acq_list{best_kymo_ind};

end

end
