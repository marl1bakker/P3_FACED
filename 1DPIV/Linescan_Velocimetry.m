%MB: Was LSPIV_parallel in original, calls upon vesselAngleD in original
%(Vessel_Diameter_Pixelsize here).
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


% DataFolder = 'D:\FACED\T2_matlab\2024-07-18-10';
% DataFolder = 'D:\FACED\T2_matlab\M99\A07';


% - ROIname for if you want a specific ROI (or ROIs), otherwise will
%       calculate all roi that were found
%       Give like: {'kymoROI_01'}; or {'ROI_1'; 'ROI_2'} (or, theoretically
%       {'roi_1', 'roi_2}), or you can say 'auto_list' so it does it
%       automatically
% - show_plots can be 0 (don't show) or 1 (show)
% - manualCrop can be 0 (take whole kymograph) or 1 (take only a part)
% - if overwrite is 1 it recalculates the velocity even though it was
%       already done and saved in the kymoroi.mat file

function Linescan_Velocimetry(DataFolder, method, show_plots, ROIname, overwrite)

% method = 1;
% have two methods: original one and one where i tried to simplify the
% original one. should give the same end results.

%% set up
if ~strcmp(DataFolder(end), filesep)
    DataFolder = [DataFolder filesep];
end

if ~exist('show_plots', 'var')
    show_plots = 1;
end

load([DataFolder 'AcqInfos.mat'], 'AcqInfoStream')
frmRate = AcqInfoStream.FrameRateHz;

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


% skipamt calculation: 
  %if it is 2, it skips every other point.  3 = skips 2/3rds of points, etc.
% Skipampt is which frame it correlates to which other frame. if you
% are very fast in your scanning, or the RBC's go too slowly, the frame
% next to your frame of interest might be the same (the RBC did not
% have time to move 1 pixel yet). Adjust accordingly.
% skipamt is what determines the number of velocity points at the end. 

% Capillary blood speed is 0.5-2 mm/s (goes up to 10 in Meng et al though).
% Size of kymograph: let's say 200 points as example
% rbc should be halfway through kymo at 

if matches(AcqInfoStream.Vessel, 'Capillary')
    % speedSetting == 1   % CAPILLARY SETTING
    numavgs       = 25;  %up to 100 (or more) for noisy or slow data
    % numavgs = round(frmRate/100);
    % skipamt       = 25;   %if it is 2, it skips every other point.  3 = skips 2/3rds of points, etc.
    skipamt = round(frmRate/100);
    shiftamt      = 5;
elseif matches(AcqInfoStream.Vessel, 'Artery')
    % speedSetting == 2   % ARTERY SETTING
    % numavgs       = 50;  %up to 100 (or more) for noisy or slow data
    skipamt       = 25;   %if it is 2, it skips every other point.  3 = skips 2/3rds of points, etc.
    shiftamt      = 1;
    % elseif speedSetting == 3   % USER SETTING
    %     disp('settings are hard coded in the script, see script!');
    %     numavgs       = 100;  %up to 200 (or more) for troublesome data. However
    %                           %you will lose some of the info in the peaks and
    %                           %troughs
    %     skipamt       = 10;   %if it is 2, it skips every other point.  3 = skips 2/3rds of points, etc.
    %     shiftamt      = 1;
end

if ~exist('overwrite', 'var')
    overwrite = 0;
end

%% Get kymographs
% get list of ROI made
if exist('ROIname', 'var') && ~matches(ROIname, 'auto_list')
    kymograph_list = ROIname;
else
    kymograph_list = dir([DataFolder 'kymoROI*.mat']);
    kymograph_list = struct2cell(kymograph_list);
    kymograph_list = kymograph_list(1,:);
end

if isempty(kymograph_list)
    error('No kymographs found, function exited.')
end


%% go per kymograph
for ind_kymo = 1:length(kymograph_list)
    % clear stuff from last one
    clear kymoImg PixelSize velocity ROIname

    % load kymograph:
    ROIname = kymograph_list{ind_kymo};
    warning('off');
    load([DataFolder ROIname], 'kymoImg', 'PixelSize', 'Velocity_calc', 'ROI_type');
    warning('on');

    if contains(ROI_type, 'perpendicular')
        disp([ROIname ' is a perpendicular ROI, velocity calculation skipped.'])
        continue
    elseif exist('Velocity_calc', 'var') && overwrite == 0
        disp(['Velocity already calculated for ' ROIname])
        clear Velocity_calc
        continue
    elseif exist('Velocity_calc', 'var') && overwrite == 1
        disp(['Velocity already calculated for ' ROIname])
        disp('OVERWRITING VELOCITY')
    end

    if size(kymoImg,2) < 20
        answer = questdlg('WARNING: The size of the ROI in the x-axis is less than 20. Try to check if you can make the ROI longer.', 'WARNING - short ROI', 'Okay I''ll change it', 'Don''t tell me what to do!', 'Okay I''ll change it');
        if matches(answer, 'Okay I''ll change it')
            continue
            % note to self: make code so you can change it directly?
        end
        % disp('Warning: ROI not long enough to make calculations on velocity!')
        % continue
    end

    if ~exist('PixelSize', 'var')
        [PixelSize] = Pixelsize_ROI(DataFolder, ROIname);
        save([DataFolder ROIname], 'PixelSize', '-append');
    end
    pxlSize = PixelSize.pxlSize;
    clear PixelSize


    %% do LSPIV correlation
    % make sure you are centered around 0 (otherwise correlation will do
    % weird things)
    DCoffset = sum(kymoImg,1) / size(kymoImg,1);
    imageLinesDC = kymoImg - repmat(DCoffset,size(kymoImg,1),1);

    %% do LSPIV correlation
    startColumn = 1;
    endColumn = size(kymoImg,2);
    % windowsize = frmRate/2;
    % numavgs = round(3/pxlSize);
    % maxGaussWidth = 100;
    % numstd = 3;

    if matches(method, 'old')
        scene_fft  = fft(imageLinesDC(1:end-shiftamt,:),[],2);
        test_img   = zeros(size(scene_fft));
        test_img(:,startColumn:endColumn)   = imageLinesDC(shiftamt+1:end, startColumn:endColumn);
        test_fft   = fft(test_img,[],2);
        W      = 1./sqrt(abs(scene_fft)) ./ sqrt(abs(test_fft)); % phase only

        LSPIVresultFFT      = scene_fft .* conj(test_fft) .* W;
        LSPIVresult         = ifft(LSPIVresultFFT,[],2);
    elseif matches(method, 'new')
        % LSPIV_xcorr = nan(size(kymoImg,1)-skipamt, size(kymoImg,2)*2-1);
        LSPIV_xcorr = nan(size(kymoImg,1)-skipamt, size(kymoImg,2)-1);
        for ind_corr = 1:size(kymoImg,1)-skipamt
            linecorr = xcorr(imageLinesDC(ind_corr,:), imageLinesDC(ind_corr+skipamt,:), 'normalized', size(imageLinesDC,2)/2-1);
            LSPIV_xcorr(ind_corr,:) = linecorr;
        end
    end

    % disp('LSPIV complete');
    % toc

    %% find shift amounts
    % disp('Find the peaks');
    velocity = [];
    maxpxlshift = round(size(kymoImg,2)/2)-1;

    if matches(method, 'old')
        index_vals = skipamt:skipamt:(size(LSPIVresult,1) - numavgs);
        numpixels = size(LSPIVresult,2);
    elseif matches(method, 'new')
        index_vals = skipamt:skipamt:(size(LSPIV_xcorr,1) - numavgs);
        numpixels = size(LSPIV_xcorr,2);
    end
    velocity  = nan(size(index_vals));
    amps      = nan(size(index_vals));
    sigmas    = nan(size(index_vals));
    goodness  = nan(size(index_vals));

    % note to self: might be better to do movmean instead of just a
    % mean over x nr of frames and move on to the next ones. Maybe
    % change later.

    %% iterate through
    b = waitbar(0, 'Finding best fit per line...');

    for index = 1:length(index_vals)
        if mod(index_vals(index),100) == 0
            waitbar(index/length(index_vals), b)
        end

        if matches(method, 'old')
            LSPIVresult_AVG   = fftshift(sum(LSPIVresult(index_vals(index):index_vals(index)+numavgs,:),1)) ...
                / max(sum(LSPIVresult(index_vals(index):index_vals(index)+numavgs,:),1));
        elseif matches(method, 'new')
            LSPIVresult_AVG = mean(LSPIV_xcorr(index_vals(index):(index_vals(index)+numavgs),:), 1);
            % LSPIVresult_AVG = sum(LSPIV_xcorr(index_vals(index):index_vals(index)+numavgs,:),1) ...
            %     / max(sum(LSPIV_xcorr(index_vals(index):index_vals(index)+numavgs,:),1));
        end

        % find a good guess for the center
        c = zeros(1, numpixels);
        % c(numpixels/2-maxpxlshift:numpixels/2+maxpxlshift) = ...
        %     LSPIVresult_AVG(numpixels/2-maxpxlshift:numpixels/2+maxpxlshift);
        c(round(numpixels/2-maxpxlshift):round(numpixels/2+maxpxlshift)) = ...
            LSPIVresult_AVG(round(numpixels/2-maxpxlshift):round(numpixels/2+maxpxlshift));
        [~, maxindex] = max(c);

        % fit a guassian to the xcorrelation to get a subpixel shift
        options = fitoptions('gauss1');
        options.Lower      = [0    numpixels/2-maxpxlshift   0            0];
        options.Upper      = [1e9  numpixels/2+maxpxlshift  maxGaussWidth 1];
        options.StartPoint = [1 maxindex 10 .1];
        warning('off')
        [q,good] = fit((1:length(LSPIVresult_AVG))',LSPIVresult_AVG','a1*exp(-((x-b1)/c1)^2) + d1',options);
        warning('on')

        %save the data
        if matches(method, 'old')
            velocity(index)  = (q.b1 - size(LSPIVresult,2)/2 - 1)/shiftamt;
        elseif matches(method, 'new')
            velocity(index)  = (q.b1 - size(LSPIV_xcorr,2)/2 - 1)/shiftamt;
        end
        amps(index)      = q.a1;
        sigmas(index)    = q.c1;
        goodness(index)  = good.rsquare;
    end

    close(b)

    %% find possible bad fits
    % toc

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

    % meanvel  = mean(velocity(goodvals)); %overall mean
    % stdvel   = std(velocity(goodvals));  %overall std

    %     case 2
    %
    %
    %         %% new - cross correlation
    %         % smooth data: do movmean over 3 micron (about half a red blood cell)
    %         % this takes the noise of the galvo out a bit
    %         imageLinesDC = movmean(imageLinesDC, round(3/pxlSize), 2);
    %
    %         % Go per frame, correlate to a frame skipamt frames from your frame of
    %         % interest. You will get results that are twice as long as your line
    %         % ROI, since you "slide" them over each other to see which shift amount
    %         % has the highest correlation.
    %         LSPIV_xcorr = nan(size(kymoImg,1)-skipamt, size(kymoImg,2)*2-1);
    %         for ind_corr = 1:size(kymoImg,1)-skipamt
    %             linecorr = xcorr(imageLinesDC(ind_corr,:), imageLinesDC(ind_corr+skipamt,:), 'normalized');
    %             LSPIV_xcorr(ind_corr,:) = linecorr;
    %         end
    %
    %         %% find shift amounts
    %         % numavgs = round(3/pxlSize);
    %         % disp('Find the peaks');
    %         maxpxlshift = round(size(kymoImg,2)/2)-1;
    %         index_vals = skipamt:skipamt:(size(LSPIV_xcorr,1) - numavgs);
    %         numpixels = size(LSPIV_xcorr,2);
    %         velocity  = nan(size(index_vals));
    %
    %         %% iterate through
    %         maxGaussWidth = 100;
    %         b = waitbar(0, 'Finding best fit per line...');
    %
    %         for index = 1:length(index_vals)
    %
    %             if mod(index_vals(index),100) == 0
    %                 waitbar(index/length(index_vals), b)
    %             end
    %
    %             LSPIVresult_AVG = mean(LSPIV_xcorr(index_vals(index):(index_vals(index)+numavgs),:), 1);
    %
    %             % find a good guess for the center
    %             c = zeros(1, numpixels);
    %             c(round(numpixels/2-maxpxlshift) : round(numpixels/2+maxpxlshift)) = ...
    %                 LSPIVresult_AVG(round(numpixels/2-maxpxlshift) : round(numpixels/2+maxpxlshift));
    %             [~, maxindex] = max(c);
    %             maxindex = maxindex - size(kymoImg,2);
    %
    %             % fit a guassian to the xcorrelation to get a subpixel shift
    %             options = fitoptions('gauss1');
    %             options.Lower      = [0    numpixels/2-maxpxlshift   0            0];
    %             options.Upper      = [1e9  numpixels/2+maxpxlshift  maxGaussWidth 1];
    %             options.StartPoint = [1 maxindex 10 .1];
    %             warning('off')
    %             [q,~] = fit((1:length(LSPIVresult_AVG))',LSPIVresult_AVG','a1*exp(-((x-b1)/c1)^2) + d1',options);
    %             warning ('on')
    %
    %             %save the data
    %             velocity(index)  = (q.b1 - size(LSPIV_xcorr,2)/2 - 1)/shiftamt;
    %
    %         end
    %         clear options c maxindex maxGaussWidth LSPIVresult_AVG q numpixels maxpxlshift numavgs
    %         close(b)
    %
    %         %% find possible bad fits
    %         % Find bad velocity points using a moving window
    %         windowsize = frmRate/2; % windowsize is half a second
    %         pixel_windowsize = round(windowsize / skipamt);
    %         numstd = 3;
    %
    %         badpixels = zeros(size(velocity));
    %         for index = 1:length(velocity)-pixel_windowsize
    %             pmean = mean(velocity(index:index+pixel_windowsize-1)); %partial window mean
    %             pstd  = std(velocity(index:index+pixel_windowsize-1));  %partial std
    %
    %             pbadpts = find((velocity(index:index+pixel_windowsize-1) > pmean + pstd*numstd) | ...
    %                 (velocity(index:index+pixel_windowsize-1) < pmean - pstd*numstd));
    %
    %             badpixels(index+pbadpts-1) = badpixels(index+pbadpts-1) + 1; %running sum of bad pts
    %         end
    %         badvals  = find(badpixels > 0); % turn pixels into indicies
    %         goodvals = find(badpixels == 0);
    %
    %         % Above detects outliers. From what I understand, you should get the same
    %         % results if you do isoutlier with a movmean of the pixel windowsize, but
    %         % it gives different results. No idea why... Will keep the way of the
    %         % original paper, since that seems to detect more.
    %         % outliers = isoutlier(velocity, 'movmean', pixel_windowsize);
    %         % plot(outliers)
    %         % hold on
    %         % plot(badpixels+2)
    %
    %         clear pmean pstd pbadpts numstd windowsize pixel_windowsize numstd badpixels
    %
    %
    % end


    %% get actual velocity
    % velocity_mm = pxlSize/1000*frmRate*resample(velocity,1:length(velocity)); % gives the same as below
    % Change frames into seconds and velocity into mm/s, right now the units
    % of velocity are pixels per frame
    velocity = velocity * pxlSize; % um per frame
    velocity = velocity / 1000; % mm per frame
    velocity = velocity * frmRate; % mm per sec

    meanvel  = mean(velocity(goodvals)); %overall mean (exclude bad fits)
    stdvel   = std(velocity(goodvals));  %overall std
    if meanvel < 0 % if the mean is negative, the flow is from right to left insted of left to right. Of course, the actual velocity is not negative, so flip data.
        velocity = -1*velocity;
        meanvel = -1*meanvel;
    end

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
        if matches(method, 'old')
            imagesc(index_vals,-numpixels/2:numpixels/2,fftshift(LSPIVresult(:,:),2)');
        elseif matches(method, 'new')
            imagesc(LSPIV_xcorr', [0 1])
        end

        colormap(ax(2), 'jet');
        title('LSPIV xcorr');
        ylabel({'displacement'; '[pixels/scan]'});

        ax(3) = subplot(3,1,3);
        plot(index_vals, velocity,'.');
        hold on
        plot(index_vals(badvals), velocity(badvals), 'ro');
        hold off
        xlim([index_vals(1) index_vals(end)]);
        ylim([meanvel-stdvel*4 meanvel+stdvel*4]);
        title('Fitted Pixel Displacement');
        ylabel({'displacement'; '[pixels/scan]'});
        xlabel('index [pixel]');

        figuretitle = 'Velocity calculation for entire acquisition';
        seps = strfind(DataFolder, filesep);
        ROI_info = [DataFolder(seps(end-1)+1:seps(end)-1) ' ' ROIname(1:end-4) ' Vessel: ' AcqInfoStream.Vessel];
        sgtitle({['{\bf\fontsize{14}' figuretitle '}'],ROI_info});


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
        if matches(method, 'old')
            imagesc(index_vals,-numpixels/2:numpixels/2,fftshift(LSPIVresult(1:round(seconds*frmRate),:),2)');
        elseif matches(method, 'new')
            imagesc(LSPIV_xcorr(1:round(seconds*frmRate),:)', [0 1])
        end
        colormap(ax(2), 'jet')
        title('LSPIV xcorr');
        ylabel({'displacement'; '(pixels/scan)'});
        xlabel('frames')

        ax(3) = subplot(3,1,3);
        perc_of_acq = (seconds*frmRate)/size(kymoImg,1); % let's say this gives 0.1044, so 10% of acq. Take 10% of velocity vector. Have to do it this way because they have different sizes and it will change if you cahnge parameters
        velocity_sec = velocity(1:round(size(velocity,2)*perc_of_acq));
        time = linspace(0, seconds, size(velocity_sec,2)); % in sec
        plot(time, velocity_sec, '.');
        badvals_sec = badvals(badvals<(size(velocity_sec,2)));
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
        sgtitle({['{\bf\fontsize{14}' figuretitle '}'],ROI_info});

    end

    %% save stuff
    Velocity_calc.velocity = velocity;
    Velocity_calc.goodvals = goodvals;
    Velocity_calc.meanvel = meanvel;
    Velocity_calc.stdvel = stdvel;

    save([DataFolder ROIname], 'Velocity_calc', '-append');

    clear Velocity_calc velocity* goodvals* badvals* meanvel stdvel time seconds perc_of_acq figuretitle ROIinfo seps

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
        end

        plot(Velocity_calc.velocity)
        legendnames = [legendnames; ROIname];

    end

    legend(legendnames)

end




































%% add?

%     medV = median(velocity);
%     velocity = medV/abs(medV)*velocity; % flip it upside down
%
%
%     %% save results and input parameters;
%     time = (index_vals + numavgs/2 )/frmRate*1000;
%     time_bad = (index_vals(badvals)+numavgs/2)/frmRate*1000;
%     if b_medianfilt
%         velocity_m = velocity;
%         if setNaN
%             velocity_m(badvals) = nan;
%         end
%         velocity_m = medfilt1(velocity_m,3);
%         %fill in the 'nan' values;
%         velocity_mm = pxlSize/1000*frmRate*resample(velocity_m,1:length(velocity_m));
%     else
%         velocity_mm = pxlSize/1000*frmRate*resample(velocity,1:length(velocity));
%     end
%
%     N = length(velocity_mm);
%     T = size(kymograph,1)/frmRate;
%     Fs = N/T;
%     f = Fs*(0:N/2)/N;
%     Y1 = abs(fft(velocity_mm,N))/N;
%     P1 = Y1(1:floor(N/2)+1);
%
%     save(fullfile(outputdir,['velocity_' ROIname]), 'velocity','velocity_mm', ...
%         'time','time_bad','index_vals','startColumn','endColumn','windowsize','numstd','f','P1');
%     % save(fullfile(outputdir,['velocity_',saveStr,'.mat']), 'velocity','velocity_mm', ...
%     %     'time','time_bad','index_vals','startColumn','endColumn','windowsize','numstd','f','P1');
%
%     if show_plots
%         %%  show the velocity;
%         hg1 = figure(111);
%         clf;
%         set(hg1,'Units','Normalized','Position',[0.2641    0.2063    0.4871    0.2590]);
%         plot(time, velocity_mm,'-','Marker','.','MarkerSize',5);
%         hold all;
%         plot(time_bad, velocity_mm(badvals), 'ro');
%         for i = 1:1:length(badvals)
%             plot( [time_bad(i),time_bad(i)],[velocity_mm(badvals(i)) velocity(badvals(i))*pxlSize/1000*frmRate],'color',[1 0.3 0.3] );
%         end
%         xlim([0 size(LSPIVresult,1)/frmRate*1000]);
%         ylim( [min([0,min(velocity_mm)]), max([nanmax(velocity_mm), nanmean(velocity_mm) + 4*nanstd(velocity_mm)])] );
%         ylim ( [-inf inf]);
%         h = line([time(1) time(end)], [nanmean(velocity_mm), nanmean(velocity_mm)] );
%         set(h, 'LineStyle','--','Color','k');
%         hold off;
%         [~,name,~] = fileparts(fname);
%         strpos = regexp(name,'ROI');
%         tpos1 = regexp(name,'_\d\d_');
%         tpos2 = regexp(name,'D');
%         tpos3 = regexp(name(tpos2:end),'_')+tpos2 - 1;
%         stackIndxStr = [name(tpos1+1 : tpos1+2),' ',name(tpos2: tpos3-1), ' '];
%         roiIndxStr = [stackIndxStr, name(strpos:end)];
%         if ~isempty(regexp(roiIndxStr,'_','once'))
%             roiIndxStr(regexp(roiIndxStr,'_')) = ' ';
%         end
%         ylabel({'Velocity'; '(mm/s)'}, 'FontSize',14);
%         xlabel('Time(ms)','FontSize',14);
%         titleStr = sprintf('%s avg%02d skip%d shift%d, mean: %.2f mm/s', roiIndxStr, numavgs,skipamt,shiftamt, mean(velocity_mm));
%         title(titleStr, 'FontSize',16);
%         print(hg1,'-r600','-dpng',fullfile(outputdir,[saveStr,'v.png']));
%
%         %% show the freqency spectrum
%         hg2 = figure(110);
%         Nyquist = frmRate/numavgs/2;
%         fLim = min(Nyquist,40);
%         subplot(1,2,1);
%         P1_norm =  mat2gray( smoothdata( P1,'gaussian',3) );
%         plot(f,P1_norm,'LineWidth',2);
%         axis([0 fLim 0 inf]);
%         ylabel('Amplitude(a.u.)', 'FontSize',12)
%         xlabel('f(Hz)', 'FontSize',12);
%         title('Frequency spectrum','FontSize',14);
%         subplot(1,2,2);
%         semilogy(f,P1_norm,'LineWidth',2);
%         axis([0 fLim 10^-3 inf]);
%         ylabel('Amplitude(a.u.)', 'FontSize',12)
%         xlabel('f(Hz)', 'FontSize',12);
%         title('Frequency spectrum - semilog','FontSize',14);
%     end
% end
%
% %% note to self
% % make sure you save per roi the velocity and diameter in the general
% % overview or the checklist or something like that
%
%
% end