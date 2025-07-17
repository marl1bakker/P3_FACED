% Last alterations that make a lot of difference:
% 16-07-25

% - Needs a datafolder to work with. Will check all kymoROI_x.mat files in
%       there but you can also specify one (ROIname later)
% - show_plots can be 0 (don't show) or 1 (show)
% - ROIname for if you want a specific ROI (or ROIs), otherwise will
%       calculate all roi that were found
%       Give like: {'kymoROI_01'}; or {'ROI_1'; 'ROI_2'} (or, theoretically
%       {'roi_1', 'roi_2}), or you can say 'auto_list' so it does it
%       automatically
% - if overwrite is 1 it recalculates the velocity even though it was
%       already done and saved in the kymoroi.mat file
% - exclusionmethod is based on the fit itself ('rsquares') or based on
%       outlier detection ('outliers'). takes rsquares automatically.

% Script is based on previous publications (see below) but altered for the
% current settings and to give the best velocity.

% % 'Line-Scanning Particle Image Velocimetry: an Optical Approach for
% % Quantifying a Wide Range of Blood Flow Speeds in Live Animals'
% % by Tyson N. Kim, Patrick W. Goodwill, Yeni Chen, Steven M. Conolly, Chris
% % B. Schaffer, Dorian Liepmann, Rong A. Wang
% %
% % Adapted by Guanghan Meng. See details in:
% % 'Ultrafast two-photon fluorescence imaging of cerebral blood circulation
% % in the awake mouse brain in vivo' by G. Meng, et al.
% % 08/31/2021



function Linescan_Velocimetry(DataFolder, show_plots, ROIname, overwrite, secs, exclusionmethod)
%% set up
if ~strcmp(DataFolder(end), filesep)
    DataFolder = [DataFolder filesep];
end

if ~exist('overwrite', 'var')
    overwrite = 0;
end

if ~exist('exclusionmethod', 'var')
    exclusionmethod = 'rsquares';
end

if ~exist('show_plots', 'var')
    show_plots = 1;
end

load([DataFolder 'AcqInfos.mat'], 'AcqInfoStream')
if isfield(AcqInfoStream, 'FrameRateHzLinescan')
    frmRate = AcqInfoStream.FrameRateHzLinescan;
else
    frmRate = AcqInfoStream.FrameRateHz;
end

% Optimized parameters
maxGaussWidth = 100;
numstd = 3;
windowsize = frmRate/2;
shiftamt = 5;
numavgs = 100;

% Pre-calculate batch size for processing
batch_size = min(50, frmRate/10); % Process in batches to reduce memory usage

method = 'fft';

%% Get kymographs
% get list of ROI made
if exist('ROIname', 'var') && ~matches(ROIname, 'auto_list')
    kymograph_list = {ROIname};
else
    kymograph_list = dir([DataFolder 'kymoROI*.mat']);
    kymograph_list = struct2cell(kymograph_list);
    kymograph_list = kymograph_list(1,:);
end

if isempty(kymograph_list)
    error('No kymographs found, function exited.')
end

% go per kymograph
for ind_kymo = 1:length(kymograph_list)
    % clear stuff from last one
    clear kymoImg PixelSize velocity ROIname

    % load kymograph:
    ROIname = kymograph_list{ind_kymo};
    if ~strcmp(ROIname(end-3:end), '.mat')
        ROIname = [ROIname '.mat'];
    end
    warning('off');
    load([DataFolder ROIname], 'kymoImg', 'Velocity_calc', 'ROI_type', 'ROI_info', 'RBC');
    warning('on');



    % check if already done
    if  exist('Velocity_calc', 'var') && isfield(Velocity_calc, method) ...
            && overwrite == 0
        disp(['Velocity already calculated for ' ROIname])
        clear Velocity_calc
        continue
    elseif exist('Velocity_calc', 'var') && isfield(Velocity_calc, method) ...
            && overwrite == 1
        disp(['Velocity already calculated for ' ROIname])
        disp('OVERWRITING VELOCITY')
    end

    if matches(ROI_info.calculate_velocity, 'No')
        disp([ROIname ' not good enough for velocity calculation - skipped.']);
        disp('If this is wrong, change ROI_info.calculate_velocity.')
        continue
    end


    %% Both methods - Get pixelsize - once
    pxlSize = 1; % Calculated with beads calibration
    % note: if you have a differnt system, this will be different!!

% kymoImg = kymoImg(1:round(0.5*frmRate),:);
% kymoImg = fliplr(kymoImg);
    %% Clean up kymo
    % Pre-calculate moving average window
    mov_window = max(1, round(ROI_info.est_frames_to_cross_kymo/10));
    kymoImg = movmean(kymoImg, mov_window, 1);

    % Vectorized DC offset removal
    DCoffset = mean(kymoImg, 1);
    imageLinesDC = kymoImg - DCoffset;

    % Vectorized intensity correction
    profile = mean(abs(imageLinesDC), 1);
    profile_offset = mean(profile);
    correction_factor = (profile - profile_offset) * -1 + profile_offset;
    imageLinesDC = imageLinesDC .* correction_factor;


    %% start calculation
    % Pre-allocate arrays
    num_lines = size(imageLinesDC,1) - shiftamt;
    if exist('secs', 'var')
        num_lines_calc = round(frmRate*secs) - shiftamt;
    else
        num_lines_calc = num_lines;
    end
    index_vals = 1:num_lines_calc - numavgs;
    numpixels = size(imageLinesDC,2);
    % checked with artificial kymo, when width is 44, center should be 23!! in
    % this case, the velocity calculation is similar to that of an odd width
    % kymograph. With 22 or 22.5 as center, a forward kymograph will
    % underestimate and a backward kymo will overestimate the velocity.
    centerPixel = ceil((numpixels+1)/2);

    % Pre-allocate output arrays
    velocity = nan(length(index_vals), 1);
    amps = nan(length(index_vals), 1);
    sigmas = nan(length(index_vals), 1);
    goodness = nan(length(index_vals), 1);

    % Batch FFT processing
    b = waitbar(0, 'Computing velocity...');

    % Compute all FFTs at once
    scene_data = imageLinesDC(1:num_lines,:);
    test_data = imageLinesDC(shiftamt+1:num_lines+shiftamt,:);

    % % Vectorized FFT computation
    scene_fft = fft(scene_data, [], 2);
    test_fft = fft(test_data, [], 2);

    % Vectorized weight calculation
    W = 1 ./ sqrt(abs(scene_fft)) ./ sqrt(abs(test_fft));

    % Vectorized cross-correlation
    LSPIVresultFFT = scene_fft .* conj(test_fft) .* W;
    LSPIVresult = ifft(LSPIVresultFFT, [], 2);

    % Batch processing for peak fitting
    fprintf('Processing %d velocity points in batches...\n', length(index_vals));

    % Pre-calculate search bounds
    maxpxlshift = round(numpixels/2) - 1;
    searchStart = max(1, round(centerPixel - maxpxlshift));
    searchEnd = min(numpixels, round(centerPixel + maxpxlshift));

    % Pre-allocate fit options
    options = fitoptions('gauss1');
    options.Lower = [0, searchStart, 0, 0];
    options.Upper = [1e9, searchEnd, maxGaussWidth, 1];

    % Process in batches to reduce memory usage
    batch_indices = 1:batch_size:length(index_vals);

    for batch_start = batch_indices
        batch_end = min(batch_start + batch_size - 1, length(index_vals));
        batch_range = batch_start:batch_end;

        waitbar(find(batch_indices == batch_start)/length(batch_indices), b)

        % Vectorized averaging and peak finding
        for indbatch = 1:length(batch_range)
            index = batch_range(indbatch);
            actual_index = index_vals(index);

            % Vectorized averaging
            segment_indices = actual_index:(actual_index + numavgs);
            LSPIVresult_segment = LSPIVresult(segment_indices, :);
            LSPIVresult_AVG = fftshift(sum(LSPIVresult_segment, 1, 'omitnan'));
            LSPIVresult_AVG = LSPIVresult_AVG / max(LSPIVresult_AVG);

            % Efficient peak finding
            c = zeros(1, numpixels);
            c(searchStart:searchEnd) = LSPIVresult_AVG(searchStart:searchEnd);
            [~, maxindex] = max(c);

            % Optimized fitting with better initial guess
            options.StartPoint = [max(LSPIVresult_AVG), maxindex, 10, 0.1];

            warning('off', 'curvefit:fit:noStartPoint');
            try
                [q, good] = fit((1:length(LSPIVresult_AVG))', LSPIVresult_AVG', ...
                    'a1*exp(-((x-b1)/c1)^2) + d1', options);

                velocity(index) = (q.b1 - centerPixel) / shiftamt;
                amps(index) = q.a1;
                sigmas(index) = q.c1;
                goodness(index) = good.rsquare;

            catch
                velocity(index) = NaN;
                amps(index) = NaN;
                sigmas(index) = NaN;
                goodness(index) = NaN;
            end
            warning('on', 'curvefit:fit:noStartPoint');
        end
    end

    clear amps sigmas options batch_indices
    close(b)



    %% find bad fits
    % toc
    switch exclusionmethod
        case 'rsquares'
            % method 1: rsquares
            % take just 0.1 as cutoff of the good.rsquares that you obtained with
            % fitting:
            goodvals = find(goodness>=0.1);
            badvals = find(goodness<0.1);

        case 'outliers'
            % % method 2: outliers
            % % This is method used by Na Ji paper
            % % Find bad velocity points using a moving window
            pixel_windowsize = round(windowsize / shiftamt);
            badpixels = zeros(size(velocity));

            % Vectorized outlier detection
            for index = 1:(length(velocity) - pixel_windowsize)
                window_data = velocity(index:(index + pixel_windowsize - 1));
                pmean = mean(window_data, 'omitnan');
                pstd = std(window_data, 'omitnan');

                outliers = (window_data > pmean + pstd * numstd) | ...
                    (window_data < pmean - pstd * numstd);

                badpixels(index + find(outliers) - 1) = badpixels(index + find(outliers) - 1) + 1;
            end

            badvals = find(badpixels > 0);
            goodvals = find(badpixels == 0);
    end

    clear pixel_windowsize index goodness numstd pbadpts


    %% get actual velocity
    % velocity_mm = pxlSize/1000*frmRate*resample(velocity,1:length(velocity)); % gives the same as below
    % Change frames into seconds and velocity into mm/s, right now the units
    % of velocity are pixels per frame
    velocity = velocity * pxlSize; % um per frame
    velocity = velocity / 1000; % mm per frame
    velocity = velocity * frmRate; % mm per sec

    if mean(velocity, 'omitnan')<0
        velocity = -velocity;
    end

    % get mean and std but exclude bad fits
    meanvel  = mean(velocity(goodvals), 'omitnan'); %overall mean (exclude bad fits)
    stdvel   = std(velocity(goodvals), 'omitnan');  %overall std

    fprintf('Estimated velocity: %.2f mm/s, Calculated mean: %.2f mm/s\n', ...
        ROI_info.est_mm_per_sec, meanvel);


    %% show results
    if show_plots
        if exist('secs', 'var')
            [f1] = plot_vel_calc(kymoImg, imageLinesDC, velocity, ...
                numavgs, LSPIVresult, badvals, frmRate, secs);
        else
            [f1] = plot_vel_calc(kymoImg, imageLinesDC, velocity, ...
                numavgs, LSPIVresult, badvals, frmRate);
        end

        % plot one second
        [f2] = plot_vel_calc(kymoImg, imageLinesDC, velocity, ...
            numavgs, LSPIVresult, badvals, frmRate, 1);


        % title = ['Velocity calculation for ' num2str(seconds) ' seconds - ' method];
        % seps = strfind(DataFolder, filesep);
        % ROI_info_title = [DataFolder(seps(end-1)+1:seps(end)-1) ' ' ROIname(1:end-4) ' Vessel: ' AcqInfoStream.Vessel];
        % sgtitle({['{\bf\fontsize{14}' figuretitle '}'],ROI_info_title});

    end


    %% save stuff

    %% temp to check if it makes sense
    if exist('Velocity_calc', 'var') && isfield(Velocity_calc, 'fft') ...
            && length(Velocity_calc.fft.velocity) >= length(velocity)
        comp_to_old = sum(velocity - Velocity_calc.fft.velocity(1:length(velocity)));
        if comp_to_old == 0 
            disp('Old fft calc gave same values and is longer, new one not saved.')
            save([DataFolder ROIname], 'imageLinesDC', '-append');
            continue
        end

        if ~isfield(Velocity_calc, 'old_fft')
            % f10 = figure;
            % plot(Velocity_calc.fft.velocity);
            % hold on
            % plot(velocity)
            % legend({'old way', 'new way'})
            % % check if it makes sense

            % keep old calculation
            Velocity_calc.fft_old = Velocity_calc.fft;
            % close(f10)
        end
    else % this is the first time you calculate vel, but you do it correctly so store an empty "old" one
        Velocity_calc.fft_old = NaN;
    end
    %% end temp
    
    % save stuff
    Velocity_calc.(method).velocity = velocity;
    Velocity_calc.(method).goodvals = goodvals;
    Velocity_calc.(method).badvals = badvals;
    Velocity_calc.(method).meanvel = meanvel;
    Velocity_calc.(method).stdvel = stdvel;

    Velocity_calc.(method).shiftamt = shiftamt;
    Velocity_calc.(method).numavgs = numavgs;


    save([DataFolder ROIname], 'Velocity_calc', '-append');
    save([DataFolder ROIname], 'imageLinesDC', '-append');
    clear Velocity_calc  goodvals* badvals* meanvel stdvel time seconds perc_of_acq figuretitle ROIinfo seps


end
end




function [f1] = plot_vel_calc(kymoImg, imageLinesDC, velocity, ...
    numavgs, LSPIVresult, badvals, frmRate, seconds)


if ~exist('seconds', 'var')
    seconds = length(velocity)/frmRate;
    plotframes = length(velocity);
else
    plotframes = floor(seconds*frmRate - numavgs - 5); % 5 is for shiftamt, hardcoded
end

numpixels = size(kymoImg,2);

     % All
        f1 = figure;
        ax(1) = subplot(4,1,1);
        imagesc(kymoImg(1:plotframes,:)')
        title('Raw Data');
        ylabel('[pixels]');
        colormap(ax(1), 'gray');

        ax(2) = subplot(4,1,2);
        imagesc(imageLinesDC(1:plotframes,:)')
        title('Filtered Data');
        ylabel('[pixels]');

        ax(3) = subplot(4,1,3);
        imagesc(1:plotframes-numavgs,-numpixels/2:numpixels/2,fftshift(LSPIVresult(1:plotframes,:),2)');

        colormap(ax(3), 'jet');
        title('LSPIV xcorr');
        ylabel({'displacement'; '[pixels/scan]'});

        ax(4) = subplot(4,1,4);
        time = linspace(0, seconds, plotframes); % in sec
        plot(time, velocity(1:plotframes),'.');
        hold on
        plot(time(badvals), velocity(badvals), 'ro');
        hold off
        xlim([time(1) time(end)]);
        % ylim([meanvel-stdvel*4 meanvel+stdvel*4]);
        title('Fitted Pixel Displacement');
        ylabel({'velocity'; '[mm/sec]'});
        xlabel('seconds');

end