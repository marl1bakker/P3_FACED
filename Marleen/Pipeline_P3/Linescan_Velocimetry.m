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
if contains(method, {'old', 'fft'})
    method = 'fft';
    % Parameters to improve fits
    maxGaussWidth = 100;  % maximum width of peak during peak fitting

    % Judge correctness of fit
    numstd        = 3;  %num of stdard deviation from the mean before flagging
    windowsize    = frmRate/2; % windowsize is only used for detecting bad fits
    % In the original code, this was the description:
    %"in # scans, this will be converted to velocity points
    %if one scan is 1/2600 s, then windowsize=2600 means
    %a 1 second moving window.  Choose the window size
    %according to experiment."
    % But as far as I (MB) can see it's only determining the windowsize of the
    % bad fit detection and does not change the nr of velocity points.


    % added 12/6/2025
    shiftamt = 5;
    numavgs = 100;

    % if matches(AcqInfoStream.Vessel, 'Capillary')
    %     % est_um_per_sec = 1000; % um/sec estimated speed, can also be around 8 mm/s (fig 4 meng) or 0.2 (fig 3 meng)
    % 
    %     shiftamt      = 5;
    %     % numavgs       = 25;
    %     numavgs = 50;
    %     % skipamt       = 25;   %if it is 2, it skips every other point.  3 = skips 2/3rds of points, etc.
    % elseif matches(AcqInfoStream.Vessel, 'Artery') || matches(AcqInfoStream.Vessel, 'Arteriole')
    %     % 30 mm/s (fig 2 meng), fig 5 shows 6 mm/s or 8 mm/s though
    %     % ??10-35 cm/s carotid -- In Vivo MRI Assessment of Blood Flow in Arteries and Veins from Head-to-Toe Across Age and Sex in C57BL/6 Mice
    %     % est_um_per_sec = 30000;
    %     % est_um_per_sec = 10000;
    %     numavgs       = 50;  %up to 100 (or more) for noisy or slow data
    %     % skipamt       = 25;
    %     shiftamt      = 2; 
    %     skipamt = 2;
    % elseif matches(AcqInfoStream.Vessel, 'Vein') || matches(AcqInfoStream.Vessel, 'Venule')
    %     % est_um_per_sec = 5000; % 5 mm/s (fig 2, 3, 5 meng)
    %     numavgs       = 100;  %up to 100 (or more) for noisy or slow data
    %     % skipamt       = 25;
    %     shiftamt      = 1;
    % end
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
    load([DataFolder ROIname], 'kymoImg', 'PixelSize', 'Velocity_calc', 'ROI_type', 'ROI_info', 'RBC');
    warning('on');

    %% temp 12-6-25 to speed up process
    if matches(method, 'fft')
        % weird nr to make fit with pulsatility calc. Could be 1*shiftamt
        kymoImg = kymoImg(1:round(2*frmRate)+numavgs+(2*shiftamt)+1,:);
    end

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
    % moved to make_kymograph
    % if ~exist('PixelSize', 'var')
    %     [PixelSize] = Pixelsize_ROI(DataFolder, ROIname);
    % end
    % pxlSize = PixelSize.pxlSize;
    % error('think of pixelsize')
    pxlSize = 1; % TRY 
    clear PixelSize


    %% Both methods - Spot where there are no RBC - once
    % if ~exist('RBC', 'var')
    %     [RBC] = RBC_presence(DataFolder, ROIname, 1);
    % end
    % no_rbc = RBC.no_rbc;
    % segments = RBC.segments;
    % clear RBC


    %% Both methods - Clean up kymo
    % kymoImg = movmean(kymoImg, 20, 1); % old: changed 17-4-25. With very fast vessels this blurs instead of makes more clear
    kymoImg = movmean(kymoImg, ROI_info.est_frames_to_cross_kymo/10, 1);
    DCoffset = sum(kymoImg,1) / size(kymoImg,1);
    imageLinesDC = kymoImg - repmat(DCoffset,size(kymoImg,1),1);

    % get rid of differences in imaging strenght over the imageing line 21/4/25
    profile = mean(abs(imageLinesDC),1);
    profile_offset = mean(profile,2);
    correction_factor = (profile-profile_offset)*-1 + profile_offset;
    imageLinesDC = imageLinesDC .* correction_factor;
    clear profile profile_offset correction_factor

    % X = [ones(size(imageLinesDC,2),1), mean(imageLinesDC,1)'];
    % B = X\imageLinesDC';
    % bla = imageLinesDC' - X'*B;

    % MB: get rid of slow freq variations [Dxx,Dxy,Dyy] = Hessian2D(I,Sigma)
    % if linescan is not long enough it will give weird results
    % if size(kymoImg,2)>20 && ROI_info.est_mm_per_sec < 5 % changed 17-4-25
    %     Sigma = 6/pxlSize; % RBC is about 6 micron
    %     [Dxx, Dxy, Dyy] = Hessian2D(imageLinesDC', Sigma);
    %     imageLinesDC = (-Dyy-Dxy)';
    %     clear Dxx Dxy Dyy DCoffset
    % else
    %     % disp('CHECK WHAT TO DO')
    % end

    % %% Both methods - get number of frames to skip for correlation and
    % %% identify orientation -- moved to make_kymograph
    % % Do this after cleaning up kymo so it's easier to see
    % % options of calculations are "whole" for estimating how long one rbc
    % % takes to cross the entire kymo, or "part" where you can take two
    % % points.
    % 
    % if exist('Velocity_calc', 'var') && isfield(Velocity_calc, 'skipamt')
    %     skipamt = Velocity_calc.skipamt;
    %     est_mm_per_sec = Velocity_calc.est_mm_per_sec;
    %     orientation = Velocity_calc.orientation;
    %     try
    %         est_frames_to_cross_kymo = Velocity_calc.est_frames_to_cross_kymo;
    %     catch % temp patch
    %         est_frames_to_cross_kymo = 4*skipamt;
    %         Velocity_calc.est_frames_to_cross_kymo = est_frames_to_cross_kymo;
    %         save([DataFolder ROIname], 'Velocity_calc',  '-append')
    %     end
    % else
    %     [skipamt, est_frames_to_cross_kymo, ~, est_mm_per_sec] = calc_skip_amount(kymoImg, frmRate, pxlSize, 'part');
    %     % [skipamt, ~, ~, est_mm_per_sec] = calc_skip_amount(imageLinesDC, frmRate, pxlSize, 'part');
    %     Velocity_calc.skipamt = skipamt;
    %     Velocity_calc.est_mm_per_sec = est_mm_per_sec;
    %     [orientation] = identify_orientation(kymoImg, frmRate);
    %     Velocity_calc.orientation = orientation;
    %     Velocity_calc.est_frames_to_cross_kymo = est_frames_to_cross_kymo;
    %     save([DataFolder ROIname], 'Velocity_calc',  '-append')
    % end


    %% Method old/1/fft
    if contains(method, {'old', 'fft'})
        startColumn = 1;
        endColumn = size(kymoImg,2);

        scene_fft  = fft(imageLinesDC(1:end-shiftamt,:),[],2);
        test_img   = zeros(size(scene_fft));
        test_img(:,startColumn:endColumn)   = imageLinesDC(shiftamt+1:end, startColumn:endColumn);
        test_fft   = fft(test_img,[],2);
        W      = 1./sqrt(abs(scene_fft)) ./ sqrt(abs(test_fft)); % phase only

        LSPIVresultFFT      = scene_fft .* conj(test_fft) .* W;
        LSPIVresult         = ifft(LSPIVresultFFT,[],2);


        % find shift amounts
        index_vals = 1:size(LSPIVresult,1)-numavgs;
        numpixels = size(LSPIVresult,2);

        velocity  = nan(size(index_vals));
        amps      = nan(size(index_vals));
        sigmas    = nan(size(index_vals));
        goodness  = nan(size(index_vals));

        % iterate through
        b = waitbar(0, 'Finding best fit per line...');

        minGaussWidth = 0;
        maxpxlshift = round(size(kymoImg,2)/2)-1;

        for index = 1:length(index_vals)

            if mod(index_vals(index),100) == 0
                waitbar(index/length(index_vals), b)
            end


            % LSPIVresult_AVG   = fftshift(sum(LSPIVresult(index_vals(index):index_vals(index)+numavgs,:),1)) ...
            %     / max(sum(LSPIVresult(index_vals(index):index_vals(index)+numavgs,:),1)); % problem: sometimes a NaN line in LSPIVresult
            LSPIVresult_AVG   = fftshift(sum(LSPIVresult(index_vals(index):index_vals(index)+numavgs,:),1, 'omitnan')) ...
                / max(sum(LSPIVresult(index_vals(index):index_vals(index)+numavgs,:),1, 'omitnan'));


            % find a good guess for the center
            c = zeros(1, numpixels);
            c(round(numpixels/2-maxpxlshift):round(numpixels/2+maxpxlshift)) = ...
                LSPIVresult_AVG(round(numpixels/2-maxpxlshift):round(numpixels/2+maxpxlshift));
            [~, maxindex] = max(c);

            % fit a guassian to the xcorrelation to get a subpixel shift
            options = fitoptions('gauss1');
            % options.Lower      = [0    numpixels/2-maxpxlshift   0            0];
            options.Lower      = [0    numpixels/2-maxpxlshift   minGaussWidth            0];
            options.Upper      = [1e9  numpixels/2+maxpxlshift  maxGaussWidth 1];
            options.StartPoint = [1 maxindex 10 .1];
            warning('off')
            try
                [q,good] = fit((1:length(LSPIVresult_AVG))',LSPIVresult_AVG','a1*exp(-((x-b1)/c1)^2) + d1',options);
            catch
                disp('x')
            end
            warning('on')

            %save the data
            velocity(index)  = (q.b1 - size(LSPIVresult,2)/2 - 1)/shiftamt;

            amps(index)      = q.a1;
            sigmas(index)    = q.c1;
            goodness(index)  = good.rsquare;

        end
        close(b)
        clear amps sigmas q index minGaussWidth maxpxlshift options b

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
                % method 2: outliers
                % This is method used by Na Ji paper
                % Find bad velocity points using a moving window
                pixel_windowsize = round(windowsize / skipamt);

                badpixels = zeros(size(velocity));
                for index = 1:1:length(velocity)-pixel_windowsize
                    pmean = mean(velocity(index:index+pixel_windowsize-1)); %partial window mean
                    pstd  = std(velocity(index:index+pixel_windowsize-1));  %partial std

                    pbadpts = find((velocity(index:index+pixel_windowsize-1) > pmean + pstd*numstd) | ...
                        (velocity(index:index+pixel_windowsize-1) < pmean - pstd*numstd));

                    badpixels(index+pbadpts-1) = badpixels(index+pbadpts-1) + 1; %running sum of bad pts
                end
                badvals  = find(badpixels > 0); % turn pixels into indicies
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

    elseif matches(method, {'new', 'xcorr'})
        method = 'xcorr';
        orientation = ROI_info.orientation;

        corr_values = nan(size(kymoImg,1)-skipamt,1);
        shift_amounts = nan(size(kymoImg,1)-skipamt,1);
        middle = ceil(size(imageLinesDC,2)/2);
        maxlag = round(size(imageLinesDC,2)/2)-1;
        LSPIV_xcorr = nan(size(kymoImg, 1)-skipamt, maxlag*2+1);

        if matches(orientation, 'Forwards')
            % peak should be at the first half of the xcorr
            for ind_corr = 1:size(kymoImg,1)-skipamt
                linecorr = xcorr(imageLinesDC(ind_corr,:), imageLinesDC(ind_corr+skipamt,:), 'normalized', maxlag);
                LSPIV_xcorr(ind_corr,:) = linecorr;

                [max_corr_val, max_corr_ind] = max(linecorr(:,1:middle));

                shift_amounts(ind_corr) = middle-max_corr_ind;
                corr_values(ind_corr) = max_corr_val;
            end
        elseif matches(orientation, 'Backwards')
            % peak should be after middle of xcorr
            for ind_corr = 1:size(kymoImg,1)-skipamt

                linecorr = xcorr(imageLinesDC(ind_corr,:), imageLinesDC(ind_corr+skipamt,:), 'normalized', maxlag);
                LSPIV_xcorr(ind_corr,:) = linecorr;

                [max_corr_val, max_corr_ind] = max(linecorr(:,middle:end));

                shift_amounts(ind_corr) = max_corr_ind-1; %because if corr is highest at middle, you need to get 0
                corr_values(ind_corr) = max_corr_val;
            end
        end

        % verify:
        if matches(orientation, 'Forwards')
            shift_amt = middle-max_corr_ind; %fwd - shift forward to "middle"
        else
            shift_amt = -(max_corr_ind-1); % bkwd - minus bc shift to the left
        end
        figure;tiledlayout('vertical')
        nexttile; plot(linecorr);
        nexttile; imagesc(imageLinesDC(ind_corr,:), [-1 1]); title(['frame ' num2str(ind_corr)])
        nexttile; imagesc(imageLinesDC(ind_corr+skipamt,:), [-1 1]); title(['frame ' num2str(ind_corr+skipamt)])
        nexttile; imagesc(circshift(imageLinesDC(ind_corr,:), shift_amt), [-1 1]); title(['frame ' num2str(ind_corr) ' shifted'])

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

        velocity = shift_amounts*pxlSize*(frmRate/skipamt)/1000; % in mm/sec

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

    elseif matches(orientation, 'varies')
        % to do: code for if the rbc's also go backwards sometimes.
        % for example: check where peak lies for three frames further,
        % or take peak around middle instead of on one side.

        % end % orientation

    end %methods



    %% Both methods - get mean and std but exclude bad fits
    meanvel  = mean(velocity(goodvals), 'omitnan'); %overall mean (exclude bad fits)
    stdvel   = std(velocity(goodvals), 'omitnan');  %overall std

    disp(['Estimated velocity was ' num2str(ROI_info.est_mm_per_sec), ...
        ', mean calculated velocity was ' num2str(meanvel)])


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

%% show all velocity calcs
if show_plots
    kymograph_list = dir([DataFolder 'kymoROI*.mat']);
    kymograph_list = struct2cell(kymograph_list);
    kymograph_list = kymograph_list(1,:);

    figure
    hold on
    legendnames = {};

    % go per kymograph
    for ind_kymo = 1:length(kymograph_list)
        ROIname = kymograph_list{ind_kymo};
        warning('off');
        load([DataFolder ROIname], 'Velocity_calc', 'ROI_type');
        warning('on');

        if contains(ROI_type, 'perpendicular')
            continue
        elseif ~exist('Velocity_calc', 'var')
            continue
        end

        plot(Velocity_calc.xcorr_range.velocity)
        legendnames = [legendnames; ROIname];
        clear Velocity_calc ROI_type
    end

    legend(legendnames)

end

end


% 
% 
% function [skipframes, est_frames_to_cross_kymo, est_sec_to_cross_kymo, est_mm_per_sec] = calc_skip_amount(kymoImg, framerate, pixelsize, option)
% 
% if exist('option', 'var') && matches(option, 'whole')
% 
%     %% skipamt calculation:
%     % if it is 2, it skips every other point.  3 = skips 2/3rds of points, etc
%     % estimate speed
%     f1 = figure; imagesc(kymoImg(1:round(framerate),:)'); colormap('gray');
%     opts.WindowStyle = 'normal';
%     est_frames_to_cross_kymo = inputdlg({'How many frames for rbc to cross the kymograph?'},'input', [1 45], {''}, opts);
%     close(f1);
%     est_frames_to_cross_kymo = str2double(est_frames_to_cross_kymo{1});
%     est_sec_to_cross_kymo = est_frames_to_cross_kymo/framerate;
%     length_um_ROI = size(kymoImg,2)*pixelsize;
%     est_um_per_sec = (1/est_sec_to_cross_kymo)*length_um_ROI;
%     est_mm_per_sec = est_um_per_sec/1000;
%     disp(['Estimated speed = ' num2str(est_mm_per_sec) ' mm per sec.'])
% 
%     if est_frames_to_cross_kymo < 1
%         warning(['It is estimated that RBC''s leave the kymo within one frame. ' ...
%             'Either the kymograph ROI is too short, or the acquisition is too ' ...
%             'slow to make an accurate speed calculation.']);
%         skipframes = 1;
%         return
%     end
% 
% else
%     %% on a part of the RBC trajectory
%     if size(kymoImg, 1)<round(framerate)
%         f1 = figure; imagesc(kymoImg(:,:)); colormap('gray');
%     else
%         f1 = figure; imagesc(kymoImg(1:round(framerate),:)); colormap('gray');
%     end
%     opts.WindowStyle = 'normal';
%     dlgtitle = 'Estimate speed';
%     prompt = {'Frames over x-axis (space)', 'Frames over y-axis (time)'};
%     fieldsize = [1 45; 1 45];
%     answers = inputdlg(prompt,dlgtitle, fieldsize, {'', ''}, opts);
% 
%     close(f1);
%     shift_amt = str2double(answers{1});
%     skipamt = str2double(answers{2});
%     est_um_per_sec = shift_amt*pixelsize*(framerate/skipamt);
%     est_mm_per_sec = est_um_per_sec/1000;
% 
%     length_um_ROI = size(kymoImg,2)*pixelsize;
%     est_sec_to_cross_kymo = length_um_ROI/est_um_per_sec;
%     est_frames_to_cross_kymo = (size(kymoImg,2)/shift_amt)*skipamt;
% 
%     disp(['Estimated speed = ' num2str(est_mm_per_sec) ' mm per sec.'])
% 
%     if est_frames_to_cross_kymo < 1
%         warning(['It is estimated that RBC''s leave the kymo within one frame. ' ...
%             'Either the kymograph ROI is too short, or the acquisition is too ' ...
%             'slow to make an accurate speed calculation.']);
%         skipframes = 1;
%         return
%     end
% 
% end
% 
% %% calculate number of frames you should skip
% skipframes = round(0.25*est_frames_to_cross_kymo);
% 
% if skipframes < 1
%     warning(['Calculated frames to skip was smaller than 1. Adjusted to'...
%         ' 4 but RBCs may be too fast to accurately calculate velocity.']);
%     skipframes = 4;
% elseif skipframes < 5
%     warning('Skipframes was smaller than 5, so kept at 5.')
%     skipframes = 5;
% end
% 
% end
% 
% 
% 
% function [orientation] = identify_orientation(kymoImg, framerate)
% 
% f1 = figure('Position', [40 60 1000 800]);
% tiledlayout(6,2)
% 
% nexttile(1, [1 2]);
% imagesc(kymoImg(:,:)'); colormap('gray');
% title('all')
% 
% if size(kymoImg,1)>round(framerate*2)
%     nexttile(3, [1 2]);
%     imagesc(kymoImg(1:round(framerate*2),:)'); colormap('gray');
%     title('2 sec')
% end
% 
% if size(kymoImg,1)>round(framerate/2)
%     nexttile(5, [1 2]);
%     imagesc(kymoImg(1:round(framerate/2),:)'); colormap('gray');
%     title('0.5 sec')
% end
% 
% nexttile(7, [1 2]);
% imagesc(kymoImg(1:100,:)'); colormap('gray');
% title('100 frames')
% 
% nexttile(9, [2 1]);
% plot([1 2], [2 1]);
% title('Forwards')
% 
% nexttile(10, [2 1]);
% plot([1 2], [1 2]);
% title('Backwards')
% 
% opts.Interpreter = 'none';
% opts.Default = 'Cancel';
% opts.WindowStyle = 'normal';
% orientation = questdlg('Are the RBC moving forwards or backwards?', 'Orientation', 'Forwards', 'Backwards', 'Cancel', opts);
% 
% if matches(orientation, 'Cancel')
%     return
% end
% 
% close(f1);
% end
% 


