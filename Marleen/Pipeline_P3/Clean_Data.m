% This function does the following:
% - Makes vessels bright and background dark
% - Divides by illumination profile

% For linescan:
% - Averages over 200 images
% - Saves kymoImg and new framerate
% - 

% For 2D scan:
% - Crops y-axis so Galvo return gets excluded
% - Crops x-axis so you only have the part with signal left
% - Averages y-axis so that one pixel is about one micron
% - Calculates pixel size in x and y directions
% - Saves the proportions of x and y to display the figure accurately, to
%       use like:
%       ax = gca;
%       ax.DataAspectRatio = proportion_XY_ax;

% Note: there's a lot of things called "range" in this function. Check per
% range in the comments what it is for.

function Clean_Data(DataFolder)

if ~strcmp(DataFolder(end), filesep)
    DataFolder = [DataFolder filesep];
end

%% load data
load([DataFolder 'AcqInfos.mat'], 'AcqInfoStream')
% check if you have a width
if isstring(AcqInfoStream.width)
    error('Width not known in infofile, or saved as string. Fix before running Clean_Data.')
end
datSize = [AcqInfoStream.nx, (AcqInfoStream.ny + AcqInfoStream.ny_extra)];


frameFormat = {'single', [double(datSize(1)), double(datSize(2))], 'imgj'}; % 

dat = memmapfile([DataFolder 'faced.dat'],...
    'Format', frameFormat, 'repeat', inf);

b = waitbar(0.3, 'Loading data...');
dat = dat.Data;
waitbar(0.6, b);
dat = reshape([dat.imgj], datSize(2), datSize(1), []);
waitbar(1, b);
close(b)
faced_cut = 0;


%% to make vessels bright and background dark:
dat = -1*dat;
lowest_value = min(dat, [], 'all');
dat = dat- lowest_value;

%% illumination profile
illumination_profile = mean(mean(dat,3),1)/max(mean(mean(dat,3),1));
dat = dat./illumination_profile;

%% in case of linescan:
if AcqInfoStream.height <= 1

    %% average over y axis
    % Linescan, so will get only 2D image that is immediately a kymograph.
    dat = permute(dat, [2,1,3]); 
    dat = reshape(dat, datSize(1), []);
    num_of_frames = size(dat,2);

    % figure;
    % tiledlayout('vertical');
    % nexttile
    % imagesc(dat(20:40,1:(AcqInfoStream.FrameRateHz*1064))); colormap('gray')
    
    % num_avs = 200;
    num_avs = 100; % new try: 27-2-25
    kymoImg = zeros(datSize(1), floor(num_of_frames/num_avs));

    line_avgs = 1;
    for ind = 1:num_avs:size(dat,2)-num_avs+1
        kymoImg(:,line_avgs) = mean(dat(:,ind:ind+num_avs-1),2);
        line_avgs = line_avgs +1;
    end

    % nexttile
    % imagesc(kymoImg(20:40,1:(AcqInfoStream.FrameRateHz*1064)/num_avs), [2000 5500]); colormap('gray')
    % title([num2str(num_avs) ' lines averaged, ' num2str(AcqInfoStream.FrameRateHz.*1064/num_avs) ' Hz'])

    maxframerate = (AcqInfoStream.ny+AcqInfoStream.ny_extra)*AcqInfoStream.FrameRateHz;
    AcqInfoStream.FrameRateHzLinescan = maxframerate/num_avs;
    % im_per_frame = (AcqInfoStream.ny+AcqInfoStream.ny_extra)/num_avs; % get x "images" per frame - increase framerate by x
    % AcqInfoStream.FrameRateHzLinescan = AcqInfoStream.FrameRateHz*im_per_frame;

    %% get borders
    showsec = 2;

    f1 = figure;
    tiledlayout('vertical');
    nexttile
    imagesc(kymoImg(:,1:showsec*round(AcqInfoStream.FrameRateHzLinescan)))
    colormap('gray');
    title(['original - ' num2str(showsec) ' seconds'])
    f1.Position = [f1.Position(1) f1.Position(2) f1.Position(3)*2 f1.Position(4)];
    y_range = [1, size(kymoImg,1)];

    answer_isgood = questdlg('Does it need to be cut?');
    nexttile
    while matches(answer_isgood, 'Yes')
        prompt = {'y begin', 'y end'};
        definput = {num2str(y_range(1)), num2str(y_range(2))};
        answer = inputdlg(prompt, 'Cropping', [1 45;1 45], definput);
        y_range = str2double(answer(1:2))';

        imagesc(kymoImg(y_range(1):y_range(2),1:round(showsec*AcqInfoStream.FrameRateHzLinescan)))
        colormap('gray')
        title(['cropped - ' num2str(showsec) ' seconds'])
        answer_isgood = questdlg('Does it need to be cut different?');
    end

    close(f1);
    kymoImg = kymoImg(y_range(1):y_range(2),:);
    kymoImg = kymoImg'; % so you save the same as other kymo's
    ROI_type = 'linescan';

    %% pixelsize
    % terminology of this is a bit confusing. pxl_sz_x is for the x-axis of
    % the original image, but the y-axis as the kymograph before was
    % depicted.
    AcqInfoStream.pxl_sz_y = 1;
    AcqInfoStream.pxl_sz_x = AcqInfoStream.width/(y_range(2)-y_range(1)+1);

    %% save
    datsize = size(dat);
    save([DataFolder 'AcqInfos.mat'], 'AcqInfoStream', '-append')
    CleanData.num_avs = num_avs;
    CleanData.y_range = y_range;
    CleanData.datsize = datsize;
    save([DataFolder 'AcqInfos.mat'], 'CleanData', '-append');

    ROI_info.calculate_velocity = 'Cancel';

    % if exist([DataFolder 'kymoROI_1.mat'], 'file')
    save([DataFolder 'kymoROI_1.mat'], 'kymoImg', 'ROI_type', 'ROI_info');
    % else
    %     save([DataFolder 'kymoROI_1.mat'], 'kymoImg', 'ROI_type','ROI_info');
    % end

    % because you normalized etc.
    fid = fopen([DataFolder 'faced.dat'],'w');
    fwrite(fid, single(dat), 'single');
    fclose(fid);

    return
end

%% in case of 2D scan:
%% find borders of field of view
%% Y-axis correction: get rid of galvo moving back
av_image_og = mean(dat, 3);
% since the return is a mirror image of the relevant data, check at which
% point they have the most correlation. This is the most likely point of
% return. The pixels needed to return is more or less at:
% height*0.5 + 100
% We will check in a range of this function +- 50 pixels, to see which
% return point matches the mirror image the best.

% return_range_start = AcqInfoStream.height*0.5 + 50;
% return_range_end = AcqInfoStream.height*0.5 + 150;
return_range_start = round((AcqInfoStream.ny+AcqInfoStream.ny_extra)*0.1);
return_range_end = round((AcqInfoStream.ny+AcqInfoStream.ny_extra)*0.5);

corr_values = zeros(1,(return_range_end-return_range_start+1));
for est_return = return_range_start:return_range_end
    new_im = av_image_og(est_return:end,:);
    cutoff_part = flipud(av_image_og(1:est_return,:));
    cutoff_part = imresize(cutoff_part, size(new_im));
    corr_values(est_return-return_range_start+1) = corr2(new_im, cutoff_part);
end
[~, indmax] = max(corr_values);
return_point = indmax+return_range_start-1;

% if show_plots
% verify:
new_im = av_image_og(return_point:end,:);
cutoff_part = flipud(av_image_og(1:return_point,:));
cutoff_part = imresize(cutoff_part, size(new_im));

f1 = figure();
t1 = tiledlayout(1,3);
colormap('gray')
nexttile(); imagesc(new_im); title('New image');
nexttile(); imagesc(cutoff_part); title('Part that is cut off');
nexttile(); imshowpair(new_im, cutoff_part); title('Comparison');
title(t1, 'Y-axis correction (galvo return)')
f1.Position = [1.8000   49.8000  766.4000  828.8000];
% end

y_range = [return_point, size(dat, 1)];


%% X-axis correction:
%% Note: this method is a bit too restrictive now. See if you can fix this later.
% x-axis is around 40-60 um, for 128 pixels. We measure the
% x-range by checking when a bead goes out of the image, so we only know
% the range in um for the part of the image where we have signal. Check
% which part of the image has signal in two ways:
signal_range = max(av_image_og, [], 1)-min(av_image_og, [], 1);

% 1. Get cut-off x-axis based on signal intensity range over x-axis, and
% take std of that as threshold
x_range_std = signal_range >= std(signal_range);
x_range_std = find(diff(x_range_std));
if length(x_range_std) > 2
    [~, maxind] = max(signal_range);
    dist_to_max = x_range_std - maxind;

    startind = dist_to_max;
    startind(startind > 0) = -100;
    [~, startind] = max(startind);
    new_x_range_std(1) = x_range_std(startind);

    endind = dist_to_max;
    endind(endind < 0) = 100;
    [~, endind] = min(endind);
    new_x_range_std(2) = x_range_std(endind);

    x_range_std = new_x_range_std;
    clear startind endind dist_to_max maxind new_x_range_std
end

% 2. Get it based on findchangepts matlab
x_range_change = findchangepts(signal_range, "MaxNumChanges",2, "Statistic","linear");

if length(x_range_change) == 2 && length(x_range_std) == 2
    % Take widest range of both, since they tend to be too restrictive
    x_range = [min(x_range_change(1), x_range_std(1)), max(x_range_change(2), x_range_std(2))];
elseif length(x_range_change) == 2
    x_range = x_range_change;
elseif length(x_range_std) == 2
    x_range = x_range_std;
else
    disp('cannot find borders on x-axis, do manually (code)')
end


%% Cropped image & pixel sizes:
answer_isgood = 'No';

while matches(answer_isgood, 'No')
    av_image_cropped = av_image_og(y_range(1):y_range(2),x_range(1):x_range(2));
    pxl_sz_x = AcqInfoStream.width/size(av_image_cropped,2); % 1 pixel = ... um in x direction
    pxl_sz_y = AcqInfoStream.height/size(av_image_cropped,1);
    proportion_XY_ax = [1 pxl_sz_x/pxl_sz_y 1];

    % if show_plots
    % check:
    f2 = figure;
    t2 = tiledlayout('flow');

    nexttile
    imagesc(av_image_og);
    title('Original image')

    nexttile
    yyaxis right; plot(illumination_profile)
    yyaxis left; plot(signal_range)
    hold on; xline(x_range(1)); xline(x_range(2))
    title('illumination profile/intensity range x-axis')

    nexttile
    imagesc(av_image_og)
    hold on;
    xline(x_range(1), 'LineWidth',1, 'Color', 'red');
    xline(x_range(2), 'LineWidth',1, 'Color', 'red');
    yline(y_range(1), 'LineWidth',1, 'Color', 'red');
    title('og image with ranges')

    nexttile
    imagesc(av_image_cropped);
    title('cropped image')

    nexttile([1 2])
    % um_size = [size(av_image,1), size(av_image,1)*(AcqInfoStream.width/AcqInfoStream.height)];
    % av_image_inproportion = imresize(av_image, um_size);
    % imagesc(av_image_inproportion); axis('image');
    % title('cropped image in proportion')
    % subtitle('imresize done, num of pixels differs')
    imagesc(av_image_cropped)
    ax = gca;
    ax.DataAspectRatio = proportion_XY_ax;
    title('cropped image in proportion')

    colormap('gray');
    title(t2, 'X and Y correction')
    f2.Position = [769.8000   49.8000  766.4000  828.8000];

    answer_isgood = questdlg('Does it make sense?');
    switch answer_isgood
        case {'No', 'Cancel'}

            prompt = {'x begin', 'x end', 'y begin', 'y end'};
            definput = {num2str(x_range(1)), num2str(x_range(2)), num2str(y_range(1)), num2str(y_range(2))};
            answer = inputdlg(prompt, 'Cropping', [1 45;1 45; 1 45;1 45], definput);
            x_range = str2double(answer(1:2))';
            y_range = str2double(answer(3:4))';

    end

    close(f2);

end
% end
close(f1)



%% average over y-axis to get 1 um per line
% this has to do with how precise our system can measure
% Apply prev cropping
dat = dat(y_range(1):y_range(2),x_range(1):x_range(2), :);

num_avs = floor(size(dat, 1)/AcqInfoStream.height); % round down
% warning: could be that some pixels are cut off because of this (at most
% num_avs-1). Could change it later so that it does take the average if
% there's more than half of the num_avs number of lines.
% dat = movmean(dat, num_avs, 1);

temp = zeros(AcqInfoStream.height, size(dat, 2), size(dat, 3));
line_avgs = 1;
for ind = 1:num_avs:size(dat, 1)-num_avs+1
    temp(line_avgs,:,:) = mean(dat(ind:ind+num_avs-1,:,:),1);
    line_avgs = line_avgs +1;
end
dat = temp;
av_image = mean(dat,3);
clear temp

% get new pixelsize that corresponds with that
AcqInfoStream.pxl_sz_y_avged = AcqInfoStream.height/size(av_image,1);


%% Save
% because you normalized etc.
fid = fopen([DataFolder 'faced.dat'],'w');
fwrite(fid, single(dat), 'single');
fclose(fid);

save([DataFolder 'morphology_img.mat'], 'av_image', 'av_image_cropped', 'proportion_XY_ax', 'av_image_og');

AcqInfoStream.pxl_sz_x = pxl_sz_x;
AcqInfoStream.pxl_sz_y = pxl_sz_y;
save([DataFolder 'AcqInfos.mat'], 'AcqInfoStream', '-append')

CleanData.faced_cut = faced_cut;
CleanData.illumination_profile = illumination_profile;
CleanData.faced_cut = faced_cut;
CleanData.y_range = y_range;
CleanData.x_range = x_range;
save([DataFolder 'AcqInfos.mat'], 'CleanData', '-append');

end
