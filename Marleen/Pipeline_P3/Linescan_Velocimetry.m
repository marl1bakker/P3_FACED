% Last alterations that make a lot of difference:
% 06-2-25
% CHECK M03 A25 WHY WRONG

% - Needs a datafolder to work with. Will check all kymoROI_x.mat files in
%       there but you can also specify one (ROIname later)
% - method can be old/fft (which is largely taken from the original code)
%       or new/xcorr. Newer method works better, it adheres better to the
%       value that you would expect when checking by eye.
% - show_plots can be 0 (don't show) or 1 (show)
% - ROIname for if you want a specific ROI (or ROIs), otherwise will
%       calculate all roi that were found
%       Give like: {'kymoROI_01'}; or {'ROI_1'; 'ROI_2'} (or, theoretically
%       {'roi_1', 'roi_2}), or you can say 'auto_list' so it does it
%       automatically
% - if overwrite is 1 it recalculates the velocity even though it was
%       already done and saved in the kymoroi.mat file
% - exclusionmethod is based on the fit itself ('rsquares') or based on
%       outlier detection ('outliers')

% MB: Was LSPIV_parallel in original.
% Adapted by Marleen (MB) 15-8-2024. For more annotations see original.
% % For additional information, please see corresponding manuscript:
% %
% % 'Line-Scanning Particle Image Velocimetry: an Optical Approach for
% % Quantifying a Wide Range of Blood Flow Speeds in Live Animals'
% % by Tyson N. Kim, Patrick W. Goodwill, Yeni Chen, Steven M. Conolly, Chris
% % B. Schaffer, Dorian Liepmann, Rong A. Wang
% %
% % Adapted by Guanghan Meng. See details in:
% % 'Ultrafast two-photon fluorescence imaging of cerebral blood circulation
% % in the awake mouse brain in vivo' by G. Meng, et al.
% % 08/31/2021



function Linescan_Velocimetry(DataFolder, method, show_plots, ROIname, overwrite, exclusionmethod)
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


%% old method
if contains(method, {'old', 'fft', 'bla'})
    % 
    
    % Optimized parameters
    maxGaussWidth = 100;
    numstd = 3;
    windowsize = frmRate/2;
    shiftamt = 5; 
    numavgs = 100;
    
    % Pre-calculate batch size for processing
    batch_size = min(50, frmRate/10); % Process in batches to reduce memory usage
end


%% Both methods - Get kymographs
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

%% Both methods - go per kymograph
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

    %% temp 12-6-25 to speed up process
    if matches(method, 'fft') || matches(method, 'bla')
        % weird nr to make fit with pulsatility calc. Could be 1*shiftamt
        % kymoImg = kymoImg(1:round(2*frmRate)+numavgs+(2*shiftamt)+1,:);
        % kymoImg = kymoImg(1:round(frmRate*0.50),:);
    end

    % if matches(ROI_info.orientation, 'Forwards')
    %     kymoImg = fliplr(kymoImg);
    % else

    % end
    if matches(method, 'xcorr')
        skipamt = max(3, round(ROI_info.est_frames_to_cross_kymo/5));
    end

    % shiftamt = max(3, ROI_info.est_frames_to_cross_kymo/5);
    % numavgs = ROI_info.est_frames_to_cross_kymo/5;

    % if contains(ROI_type, 'perpendicular')
    if ~contains(ROI_type, 'line')
        disp([ROIname ' is a perpendicular ROI, velocity calculation skipped.'])
        continue
    elseif exist('Velocity_calc', 'var') && isfield(Velocity_calc, method) ...
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

    % %% Both methods - Clean up kymo - old
    % kymoImg = movmean(kymoImg, ROI_info.est_frames_to_cross_kymo/10, 1);
    % DCoffset = sum(kymoImg,1) / size(kymoImg,1);
    % imageLinesDC = kymoImg - repmat(DCoffset,size(kymoImg,1),1);
    % 
    % % get rid of differences in imaging strenght over the imageing line 21/4/25
    % profile = mean(abs(imageLinesDC),1);
    % profile_offset = mean(profile,2);
    % correction_factor = (profile-profile_offset)*-1 + profile_offset;
    % imageLinesDC = imageLinesDC .* correction_factor;
    % clear profile profile_offset correction_factor

    %% Clean up kymo - optimized
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


    %% Method old/1/fft
    if contains(method, {'old', 'fft'})
        method = 'fft';

        % Pre-allocate arrays
        num_lines = size(imageLinesDC,1) - shiftamt;
        index_vals = 1:num_lines - numavgs;
        numpixels = size(imageLinesDC,2);
% checked with artificial kymo, when width is 44, center should be 23!! in 
% this case, the velocity calculation is similar to that of an odd width
% kymograph. With 22 or 22.5 as center, a forward kymograph will
% underestimate and a backward kymo will overestimate the velocity. 
        centerPixel = ceil(numpixels+1)/2;

        % Pre-allocate output arrays
        velocity = nan(length(index_vals), 1);
        amps = nan(length(index_vals), 1);
        sigmas = nan(length(index_vals), 1);
        goodness = nan(length(index_vals), 1);

        % Batch FFT processing
        b = waitbar(0, 'Computing velocity...');
        % fprintf('Computing FFTs for %d lines...\n', num_lines);

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
        
        % bla = NaN(size(test_data,1), size(test_data,2)*2-1);
        % for index = 1:size(scene_data,1)
        %     bla(index, :) = xcorr(scene_data(index,:), test_data(index,:), 'normalized');
        % end

        % Batch processing for peak fitting
        fprintf('Processing %d velocity points in batches...\n', length(index_vals));
        
        % Pre-calculate search bounds
        maxpxlshift = round(numpixels/2) - 1;
        searchStart = max(1, round(centerPixel - maxpxlshift));
        % searchStart = 1;
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


                % % new 15-7                
                % % only take current peak of linecorr so that the fit makes sense
                % try
                %     [~, lowind] = findpeaks(LSPIVresult_AVG(1:maxindex)*-1); % find lows
                %     lowbefore = lowind(end);
                % catch
                %     lowbefore = 1;
                % end
                % try
                %     [~, lowind] = findpeaks(LSPIVresult_AVG(maxindex:end)*-1);
                %     lowafter = lowind(1)+maxindex-1;
                % catch
                %     lowafter = length(LSPIVresult_AVG);
                % end
                % 
                % linecorrtofit = [repmat(LSPIVresult_AVG(lowbefore), lowbefore-1,1)' ...
                %     LSPIVresult_AVG(lowbefore:lowafter) ...
                %     repmat(LSPIVresult_AVG(lowafter), length(LSPIVresult_AVG)-lowafter,1)'];

                
                warning('off', 'curvefit:fit:noStartPoint');
                try
                    [q, good] = fit((1:length(LSPIVresult_AVG))', LSPIVresult_AVG', ...
                        'a1*exp(-((x-b1)/c1)^2) + d1', options);
                    
                    velocity(index) = (q.b1 - centerPixel) / shiftamt;
                    amps(index) = q.a1;
                    sigmas(index) = q.c1;
                    goodness(index) = good.rsquare;

                    % [q, good] = fit((1:length(linecorrtofit))', linecorrtofit', ...
                    %     'a1*exp(-((x-b1)/c1)^2) + d1', options);
                    
                    % velocity_x(index) = (q.b1 - centerPixel) / shiftamt;
                    
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

        % velocity_x = velocity_x * pxlSize / 1000 * frmRate;

    elseif matches(method, {'new', 'xcorr'})
        method = 'xcorr';
        orientation = ROI_info.orientation;

        corr_values = nan(size(kymoImg,1)-skipamt,1);
        shift_amounts = nan(size(kymoImg,1)-skipamt,1);
        middle = ceil(size(imageLinesDC,2)/2);
        maxlag = round(size(imageLinesDC,2)/2)-1;
        LSPIV_xcorr = nan(size(kymoImg, 1)-skipamt, maxlag*2+1);
        velocity = nan(size(kymoImg, 1), 1);

              % Pre-allocate fit options
        options = fitoptions('gauss1');
        options.Lower = [0, 1, 0, 0];
        maxGaussWidth = 100;
        options.Upper = [1e9, size(imageLinesDC,2), maxGaussWidth, 1];

        if matches(orientation, 'Forwards')
            % peak should be at the first half of the xcorr
            for ind_corr = 1:size(kymoImg,1)-skipamt
                linecorr = xcorr(imageLinesDC(ind_corr,:), imageLinesDC(ind_corr+skipamt,:), 'normalized', maxlag);
                LSPIV_xcorr(ind_corr,:) = linecorr;

                % [max_corr_val, max_corr_ind] = max(linecorr);
                % if max_corr_ind >= middle
                %     continue
                % end
                [max_corr_val, max_corr_ind] = max(linecorr(:,1:middle));
                
                % Optimized fitting with better initial guess
                options.StartPoint = [max_corr_val, max_corr_ind, 10, 0.1];
                
                % only take current peak of linecorr so that the fit makes sense
                try
                [~, lowind] = findpeaks(linecorr(1:max_corr_ind)*-1); % find lows
                catch
                    disp('x')
                    continue
                end
                if ~isempty(lowind)
                    lowbefore = lowind(end);
                else 
                    lowbefore = 1;
                end
                try
                [~, lowind] = findpeaks(linecorr(max_corr_ind:end)*-1);
                catch
                    disp('xx')
                end
                if ~isempty(lowind)
                    lowafter = lowind(1)+max_corr_ind-1;
                else 
                    lowafter = length(linecorr);
                end                

                linecorrtofit = [repmat(linecorr(lowbefore), lowbefore-1,1)' ...
                    linecorr(lowbefore:lowafter) ...
                    repmat(linecorr(lowafter), length(linecorr)-lowafter,1)'];

                warning('off', 'curvefit:fit:noStartPoint');
                try
                    % [q, good] = fit((1:length(linecorr))', linecorr', ...
                    %     'a1*exp(-((x-b1)/c1)^2) + d1', options);
                    [q, good] = fit((1:length(linecorrtofit))', linecorrtofit','a1*exp(-((x-b1)/c1)^2) + d1', options);

                    velocity(ind_corr) = (middle-q.b1) / skipamt;
                    
                catch
                    velocity(index) = NaN;
                    % amps(index) = NaN;
                    % sigmas(index) = NaN;
                    % goodness(index) = NaN;
                end
                warning('on', 'curvefit:fit:noStartPoint');

                shift_amounts(ind_corr) = middle-max_corr_ind;
                corr_values(ind_corr) = max_corr_val;
            end
        elseif matches(orientation, 'Backwards')
            % peak should be after middle of xcorr
            for ind_corr = 1:size(kymoImg,1)-skipamt

                linecorr = xcorr(imageLinesDC(ind_corr,:), imageLinesDC(ind_corr+skipamt,:), 'normalized', maxlag);
                LSPIV_xcorr(ind_corr,:) = linecorr;

                [max_corr_val, max_corr_ind] = max(linecorr(:,middle:end));

                shift_amounts(ind_corr) = -(max_corr_ind-1); %because if corr is highest at middle, you need to get 0
                corr_values(ind_corr) = max_corr_val;
            end
        end


        % figure;tiledlayout('vertical')
        % nexttile; plot(linecorr);
        % nexttile; imagesc(imageLinesDC(ind_corr,:)); title(['frame ' num2str(ind_corr)])
        % nexttile; imagesc(imageLinesDC(ind_corr+skipamt,:)); title(['frame ' num2str(ind_corr+skipamt)])
        % 
        % nexttile; imagesc(circshift(imageLinesDC(ind_corr,:), shift_amounts(ind_corr))); title(['frame ' num2str(ind_corr) ' shifted'])


        % numavgs = 0; % tried up to 5, doesnt make a big difference.
        % for ind_corr = 1:size(kymoImg,1)-skipamt-numavgs
        % linecorr = xcorr(mean(imageLinesDC(ind_corr:ind_corr+numavgs,:),1), mean(imageLinesDC(ind_corr+skipamt:ind_corr+numavgs+skipamt,:),1), 'normalized', maxlag);
        % LSPIV_xcorr(ind_corr,:) = linecorr;
        %
        % [max_corr_val, max_corr_ind] = max(linecorr(:,middle:end));
        %
        %     shift_amounts(ind_corr) = max_corr_ind;
        %     corr_values(ind_corr) = max_corr_val;
        % end

        % velocity = shift_amounts*pxlSize*(frmRate/skipamt)/1000; % in mm/sec
        velocity = velocity * pxlSize / 1000 * frmRate;

        % index_vals = 1:size(LSPIV_xcorr,1)-numavgs;
        threshold_corr = 0.5;
        index_vals = 1:size(LSPIV_xcorr,1);
        badvals = find(corr_values<threshold_corr);
        goodvals = find(corr_values>=threshold_corr);

        raw_velocity = velocity;
        velocity(badvals) = NaN;
        velocity = movmean(velocity, 50, 'omitnan');


        %% X correlation with range
    elseif matches(method, {'xcorr_range'})
        orientation = ROI_info.orientation;

        % Set-up of settings
        if matches(orientation, 'Forwards')
            % flip the kymograph.
            kymoImg = fliplr(kymoImg);
            imageLinesDC = fliplr(imageLinesDC);
        end

        % If speed is low, the max corr peak will lie in the middle, and to
        % recognize that as peak, you need to take the middle
        % into account. Could make speed slightly negative
        % if ROI_info.est_mm_per_sec<5 % changed 17-4-25
            if matches(orientation, 'Backwards') || matches(orientation, 'Forwards')
                playroom_middle = round(size(kymoImg,2)/10);
            elseif matches(orientation, 'Varies')
                playroom_middle = ceil(size(kymoImg,2)/5);
            else
                playroom_middle = ceil(size(kymoImg,2)/5);
            end
        % else
        %     playroom_middle = 0;
        % end

        % skiprange: compare multiple upcoming frames to frame of
        % interest.
        skipamt = round(0.5*ROI_info.est_frames_to_cross_kymo);% couldve put this at calc-skipamt but already calculated with 0.25, this is easier
        % skipamt_range_start = skipamt - 50;
        if skipamt>100 % 21-4-25
            skipamt_range_start = skipamt - 50;
        else
            skipamt_range_start = 1;
        end

        velocity = NaN(size(kymoImg-skipamt, 1), 1);
        weights = NaN(size(kymoImg-skipamt, 1), 1);
        middle = ceil(size(imageLinesDC,2)/2);
        maxlag = round(size(imageLinesDC,2)/2)-1;


        b = waitbar(0, 'Finding best fit per line...');
        progress_bar_points = 1:frmRate:size(kymoImg,1)-skipamt-1;
        progress_bar_points = round(progress_bar_points);

        % Start going per frame
        for ind_corr = 1:size(kymoImg,1)-skipamt
            % maxskipamt_reached = [];
            max_corr_val_prev = 1;
            velocityline = NaN(skipamt, 1);
            % figure;hold on; % to plot skipamts

            for ind_skip_amt = skipamt_range_start:skipamt
                % if ~isempty(maxskipamt_reached)
                %     break
                % end

                linecorr = xcorr(imageLinesDC(ind_corr,:), imageLinesDC(ind_corr+ind_skip_amt,:), 'normalized', maxlag);
                [max_corr_val, max_corr_ind] = findpeaks(linecorr(:,middle-playroom_middle:end));

                if isempty(max_corr_val) % skip if there is no peak
                    continue
                else % if there's multiple peaks, take the first
                    max_corr_val = max_corr_val(1);
                    max_corr_ind = max_corr_ind(1)-playroom_middle;
                end
                if max_corr_ind == 1 % if peak is by not shifting the frame, it doesnt give any info, so skip
                    continue
                end

                if max_corr_val < 0.5 % skip if correlation is low (this is end of this line xcorr also, correlation should not get better further away)
                    maxskipamt_reached = ind_skip_amt;
                    break
                elseif max_corr_val > max_corr_val_prev+0.1 % skip if correlation is higher than previous line, since you're then likely onto the next RBC
                    maxskipamt_reached = ind_skip_amt;
                    break
                end

                % % to visualize: do only one line (ind_corr) and make
                % % figure with hold on before doing this.
                % plot(linecorr);
                % plot(max_corr_ind+middle-1, max_corr_val, '*')

                max_corr_val_prev = max_corr_val;
                velocityline(ind_skip_amt) = (max_corr_ind-1)*pxlSize * (frmRate/ind_skip_amt) / 1000; % in mm/sec
            end
            velocity(ind_corr) = mean(velocityline, 1, 'omitnan');
            weights(ind_corr) = sum(~isnan(velocityline));

            if ~isempty(find(progress_bar_points == ind_corr, 1))
                waitbar(ind_corr/size(kymoImg,1), b)
            end
        end

        close(b)
        disp('x')

        badvals = find(isnan(velocity));
        goodvals = ones(size(velocity));
        goodvals(badvals) = 0;
        goodvals = find(goodvals);

        raw_velocity = velocity;
        velocity = movmean(velocity, skipamt, 'omitnan');

        % f1 = figure;imagesc(kymoImg(1:round(frmRate),:)')
        if raw_velocity>round(frmRate)
            f2 = figure;plot(raw_velocity(1:round(frmRate))); hold on; plot(velocity(1:round(frmRate)))
        % close(f1, f2)
        end

        f5 = figure;
        yyaxis left; imagesc(kymoImg'); colormap('gray'); ylabel('kymograph')
        %      Acq_msec = size(kymoImg,1)/frmRate*1000;
        %     xticks(0:frmRate/10:size(kymoImg,1)); % this times next line needs to be 1000
        %     xticklabels(0:100:Acq_msec);
        %     xlabel('Milliseconds');
        Acq_sec = size(kymoImg,1)/frmRate;
        xticks(0:frmRate/2:size(kymoImg,1)); % this times next line needs to be 1
        xticklabels(0:0.5:Acq_sec);
        xlabel('Seconds');

        yyaxis right; plot(velocity, 'Color', 'red', 'LineWidth', 2); ylabel('velocity mm/s')

        close(f5)
                if raw_velocity>round(frmRate)
close(f2)
                end

    % elseif matches(orientation, 'varies')
        % to do: code for if the rbc's also go backwards sometimes.
        % for example: check where peak lies for three frames further,
        % or take peak around middle instead of on one side.

        % end % orientation

    elseif matches(method, 'bla')
        % maxlag = round(size(imageLinesDC,2)/2)-1;
        % skipamt = max(3, round(ROI_info.est_frames_to_cross_kymo/5));
        % 
        % velocity = nan(size(imageLinesDC,2),1);
        % for ind_corr = 1:size(kymoImg,1)-skipamt
        % 
        %     [vel, fit_quality, peakinfo] = improved_peak_fitting(imageLinesDC, ind_corr, skipamt, maxlag, 'PreviousVelocity', ROI_info.est_mm_per_sec);
        %     velocity(ind_corr) = vel;
        % end
    end %methods



    %% Both methods - get mean and std but exclude bad fits
    meanvel  = mean(velocity(goodvals), 'omitnan'); %overall mean (exclude bad fits)
    stdvel   = std(velocity(goodvals), 'omitnan');  %overall std

    fprintf('Estimated velocity: %.2f mm/s, Calculated mean: %.2f mm/s\n', ...
        ROI_info.est_mm_per_sec, meanvel);

    %% show results
    if show_plots

        % All
        figure
        ax(1) = subplot(3,1,1);
        % imagesc(kymoImg(:, startColumn:endColumn)')
        imagesc(kymoImg')
        title('Raw Data');
        ylabel('[pixels]');
        colormap(ax(1), 'gray');

        ax(2) = subplot(3,1,2);
        if matches(method, {'old', 'fft'})
            imagesc(index_vals,-numpixels/2:numpixels/2,fftshift(LSPIVresult(:,:),2)');
        elseif matches(method, {'new', 'xcorr'})
            imagesc(LSPIV_xcorr', [0 1])
        end

        colormap(ax(2), 'jet');
        title('LSPIV xcorr');
        ylabel({'displacement'; '[pixels/scan]'});

        ax(3) = subplot(3,1,3);
        timeacq = size(imageLinesDC,1)/frmRate;
        time = linspace(0, timeacq, length(velocity)); % in sec
        plot(time, velocity,'.');
        hold on
        plot(time(badvals), velocity(badvals), 'ro');
        hold off
        xlim([time(1) time(end)]);
        ylim([meanvel-stdvel*4 meanvel+stdvel*4]);
        title('Fitted Pixel Displacement');
        ylabel({'velocity'; '[mm/sec]'});
        xlabel('seconds');

        figuretitle = ['Velocity calculation for entire acquisition - ' method];
        seps = strfind(DataFolder, filesep);
        ROI_info_title = [DataFolder(seps(end-1)+1:seps(end)-1) ' ' ROIname(1:end-4) ' Vessel: ' AcqInfoStream.Vessel];
        sgtitle({['{\bf\fontsize{14}' figuretitle '}'],ROI_info_title});


        % One second
        seconds = 1; % so you can change it

        figure
        ax(1) = subplot(3,1,1);
        imagesc(kymoImg(1:round(seconds*frmRate), :)');
        % imagesc(kymoImg(1:round(seconds*frmRate), startColumn:endColumn)');
        title('Raw Data');
        ylabel('pixels along ROI');
        colormap(ax(1), 'gray');
        xlabel('frames')

        ax(2) = subplot(3,1,2);
        if matches(method, {'old', 'fft'})
            imagesc(index_vals,-numpixels/2:numpixels/2,fftshift(LSPIVresult(1:round(seconds*frmRate),:),2)');
        elseif matches(method, {'new', 'xcorr'})
            imagesc(LSPIV_xcorr(1:round(seconds*frmRate),:)', [0 1])
        end
        colormap(ax(2), 'jet')
        title('LSPIV xcorr');
        ylabel({'displacement'; '(pixels/scan)'});
        xlabel('frames')

        ax(3) = subplot(3,1,3);
        perc_of_acq = (seconds*frmRate)/size(kymoImg,1); % let's say this gives 0.1044, so 10% of acq. Take 10% of velocity vector. Have to do it this way because they have different sizes and it will change if you cahnge parameters
        velocity_sec = velocity(1:round(length(velocity)*perc_of_acq));
        time = linspace(0, seconds, length(velocity_sec)); % in sec
        plot(time, velocity_sec, '.');
        badvals_sec = badvals(badvals<(length(velocity_sec)));
        hold on
        plot(time(badvals_sec), velocity_sec(badvals_sec), 'ro');
        hold off
        ylim([meanvel-stdvel*4 meanvel+stdvel*4]);
        title('Velocity');
        ylabel('mm/seconds');
        xlabel('time (seconds)');

        if seconds == 1
            figuretitle = ['Velocity calculation for ' num2str(seconds) ' second'];
        else
            figuretitle = ['Velocity calculation for ' num2str(seconds) ' seconds'];
        end
        sgtitle({['{\bf\fontsize{14}' figuretitle '}'],ROI_info_title});

    end
 

    %% save stuff
    % savefig(f5, [DataFolder 'Kymograph_overview.fig']);

    %% temp to check if it makes sense
    if isfield(Velocity_calc, 'fft') && ~isfield(Velocity_calc, 'old_fft')
        f10 = figure;
        plot(Velocity_calc.fft.velocity);
        hold on
        plot(velocity)
        legend({'old way', 'new way'})
        % check if it makes sense

        % keep old calculation
        Velocity_calc.fft_old = Velocity_calc.fft;
        close(f10)
    end

    Velocity_calc.(method).velocity = velocity;
    Velocity_calc.(method).goodvals = goodvals;
    Velocity_calc.(method).badvals = badvals;
    Velocity_calc.(method).meanvel = meanvel;
    Velocity_calc.(method).stdvel = stdvel;

    if matches(method, 'xcorr')
        Velocity_calc.(method).corr_values = corr_values;
        Velocity_calc.(method).raw_velocity = raw_velocity;
    elseif matches(method, 'xcorr_range')
        Velocity_calc.(method).raw_velocity = raw_velocity;
        Velocity_calc.(method).weights = weights;
    elseif matches(method, 'fft')
        Velocity_calc.(method).shiftamt = shiftamt;
        Velocity_calc.(method).numavgs = numavgs;        
    end


    save([DataFolder ROIname], 'Velocity_calc', '-append');

    clear Velocity_calc  goodvals* badvals* meanvel stdvel time seconds perc_of_acq figuretitle ROIinfo seps


end
end


