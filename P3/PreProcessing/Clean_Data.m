% This function does the following:
% - Makes vessels bright and background dark
% - Divides by illumination profile
% - Crops y-axis so Galvo return gets excluded
% - Crops x-axis so you only have the part with signal left
% - Calculates pixel size in x and y directions
% - Saves the proportions of x and y to display the figure accurately, to
%       use like:
%       ax = gca;
%       ax.DataAspectRatio = proportion_XY_ax;

% Note: there's a lot of things called "range" in this function. Check per
% range in the comments what it is for.

function Clean_Data(DataFolder, show_plots)

if ~strcmp(DataFolder(end), filesep)
    DataFolder = [DataFolder filesep];
end

if ~exist('show_plots', 'var')
    show_plots = 1;
end

%% load data
load([DataFolder 'AcqInfos.mat'], 'AcqInfoStream')
datSize = [AcqInfoStream.nx, (AcqInfoStream.ny + AcqInfoStream.ny_extra)];

fid = fopen([DataFolder 'faced.dat']);
try
    dat = fread(fid, inf, '*single');
catch
    fid = fopen([DataFolder 'faced.dat']);
    nr_of_frames = 9000; 
    dat = fread(fid, datSize(1)*datSize(2)*nr_of_frames, '*single');
    disp('WARNING, FRAMES CUT TO 9000')
    disp('fix later')
end
fclose(fid);
dat = reshape(dat, datSize(2), datSize(1), []);

%% to make vessels bright and background dark:
dat = -1*dat;
lowest_value = min(dat, [], 'all');
dat = dat- lowest_value;

%% illumination profile
illumination_profile = mean(mean(dat,3),1)/max(mean(mean(dat,3),1));
dat = dat./illumination_profile;

%% Y-axis correction: get rid of galvo moving back
% since the return is a mirror image of the relevant data, check at which
% point they have the most correlation. This is the most likely point of
% return. The pixels needed to return is more or less at: 
% height*0.5 + 100
% We will check in a range of this function +- 50 pixels, to see which
% return point matches the mirror image the best. 

av_image_og = mean(dat, 3);
return_range_start = AcqInfoStream.height*0.5 + 50;
return_range_end = AcqInfoStream.height*0.5 + 150;

corr_values = zeros(1,(return_range_end-return_range_start+1));
for est_return = return_range_start:return_range_end 
    new_im = av_image_og(est_return:end,:);
    cutoff_part = flipud(av_image_og(1:est_return,:));
    cutoff_part = imresize(cutoff_part, size(new_im));
    corr_values(est_return-return_range_start+1) = corr2(new_im, cutoff_part);
end
[~, indmax] = max(corr_values);
return_point = indmax+return_range_start-1;

if show_plots
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
end

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
av_image_cropped = av_image_og(y_range(1):y_range(2),x_range(1):x_range(2));
pxl_sz_x = AcqInfoStream.width/size(av_image_cropped,2); % 1 pixel = ... um in x direction
pxl_sz_y = AcqInfoStream.height/size(av_image_cropped,1);
proportion_XY_ax = [1 pxl_sz_x/pxl_sz_y 1];

if show_plots
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

    answer = questdlg('Does it make sense?');
    switch answer
        case {'No', 'Cancel'}
            error('Cropping was not done properly, function exited')
    end

    close(f1, f2);
end


%% Save
save([DataFolder 'morphology_img.mat'], 'av_image_cropped', 'proportion_XY_ax', 'av_image_og');

% AcqInfoStream.new_ny = size(av_image,1);
% AcqInfoStream.new_nx = size(av_image,2);
AcqInfoStream.pxl_sz_x = pxl_sz_x;
AcqInfoStream.pxl_sz_y = pxl_sz_y; 
save([DataFolder 'AcqInfos.mat'], 'AcqInfoStream')

% dat = dat(y_range(1):y_range(2), x_range(1):x_range(2), :);
% fid = fopen([DataFolder 'faced.dat'],'w');
% fwrite(fid, single(dat), 'single');
% fclose(fid);

end


%% correlate per line -- does not change much :(
% datcor = zeros(size(dat));
% 
% for indframe = 1:size(dat,3)
%     datcor(1,:,indframe) = dat(1,:,indframe);
% 
%     for indy = 2:size(dat, 1)
%         corrvals = zeros(11,1);
% 
%         for shift = -5:1:5
%             corrvals(shift+6) = corr(datcor(indy,:,indframe)', circshift(dat(indy,:,indframe)', shift));
%         end
%         [~, indmaxcorr] = max(corrvals);
%         indmaxcorr = indmaxcorr - 6;
% 
%         datcor(indy,:,indframe) = circshift(dat(indy,:,indframe), indmaxcorr);
%     end
% end
