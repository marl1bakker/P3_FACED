% TODO SAVE PULSATILITY IN OVERVIEW.MAT RESULTS.PULSATILITY

% Pulsatility index defined as difference between peak systolic flow and
% minimum diastolic flow / mean velocity
%
% Could not find original code of Meng et al to calculate this, so made
% from scratch.
% Marleen Bakker 18-10-24 // 23-01-25

function [pulsatility_output] = Pulsatility_Index(DataFolder, show_plots, ROIname, overwrite)


%% set up
if ~strcmp(DataFolder(end), filesep)
    DataFolder = [DataFolder filesep];
end

if ~exist('show_plots', 'var')
    show_plots = 1;
end

if ~exist('overwrite', 'var')
    overwrite = 0;
end

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
else % make table for output
    pulsatility_output = table('Size', [length(kymograph_list), 5], ...
        'VariableNames', {'Name', 'Type', 'Pulsatility_mean', 'Pulsatility_median', 'Heartbeat'},...
        'VariableTypes', {'string', 'string', 'single', 'single', 'single'});
end

% naming stuff
load([DataFolder 'AcqInfos.mat'], 'AcqInfoStream')
seps = strfind(DataFolder, filesep);
Mouse = DataFolder(seps(end-2)+1:seps(end-1)-1);
Acq = AcqInfoStream.DatasetName; % or: DataFolder(seps(end-1)+1:end-1);
clear seps



%% go per kymograph
for ind_kymo = 1:length(kymograph_list)

    %% load kymograph:
    ROIname = kymograph_list{ind_kymo};
    warning('off');
    load([DataFolder ROIname], 'Velocity_calc', 'kymoImg', 'ROI_type', 'Pulsatility_calc');
    warning('on');

    pulsatility_output.Name(ind_kymo) = ROIname;
    pulsatility_output.Type(ind_kymo) = ROI_type;

    %% check if already done
    if ~exist('Velocity_calc', 'var')
        error(['Velocity not found for ' ROIname ' of ' Mouse ', ' Acq '. Run Linescan_Velocimetry first.']);

    elseif exist('Pulsatility_calc', 'var') && overwrite == 0
        disp(['Pulsatility already calculated for ' ROIname ' ' Mouse ' ' Acq '. ROI skipped.'])
        continue

    elseif exist('Pulsatility_calc', 'var') && overwrite == 1
        disp(['Pulsatility already calculated for ' ROIname ' ' Mouse ' ' Acq '. OVERWRITING PULSATILITY CALCULATIONS.'])

    end

    %% load parameters
    vel = Velocity_calc.velocity;
    if isfield(AcqInfoStream, 'FrameRateHzLinescan')
        frmRate = AcqInfoStream.FrameRateHzLinescan;
    else
        frmRate = AcqInfoStream.FrameRateHz;
    end
    seconds = size(kymoImg, 1)/frmRate; % sec of kymoROI in question
    % vel_freq = size(vel,2)/seconds; %freq of velocity calc

    %% Mask outliers --- CHECK
    % for all the velocity values that were "bad", take an average of
    % values before and after so that there are no nans
    vel_mask = NaN(size(vel));
    vel_mask(Velocity_calc.goodvals) = 1;
    % make mask broader at the front
    mask_ind = find(isnan(vel_mask));
    vel_mask(mask_ind-1) = NaN;
    vel = vel.*vel_mask;
    clear mask_ind

    %% get heartbeats
    % make vector so you can work with it. Fit curve.
    time = linspace(0,seconds,length(vel));
    opts = fitoptions('Method', 'SmoothingSpline');
    opts.SmoothingParam = 0.99999; % 5 times 9 seems best
    [curve, ~] = fit(time',vel', "smoothingspline", opts);
    yFitted = feval(curve, time');

    % get peaks
    [pks, indx] = findpeaks(yFitted);
    av_heart_beat = length(pks)/seconds*60;

    % If heartbeat is off:
    if av_heart_beat < 300 || av_heart_beat > 550
        disp(['Heartbeat of ' num2str(av_heart_beat) ' detected, which is outside the expected range.'])
        return
    else
        pulsatility_output.Heartbeat(ind_kymo) = av_heart_beat;
    end

    % get pulsatility index
    PI = NaN(size(pks));
    for ind = 1:size(pks)-1
        % take peak, find low, calculate average
        peak = pks(ind);
        heartcycle = yFitted(indx(ind):indx(ind+1)-1);
        low = min(heartcycle);
        av_heartcycle = mean(heartcycle);

        PI(ind) = (peak-low)/av_heartcycle;
    end

    if show_plots
        startframe = 1;
        secsplot = 2;
        endframe = round(startframe+secsplot*frmRate);

        f2 = figure;
        tiledlayout('vertical')

        nexttile
        title([Mouse ' '])
        imagesc(kymoImg')
        xlim([startframe;endframe])

        colormap('gray')

        nexttile
        scatter(time', vel', '.')
        xlim([startframe/frmRate;endframe/frmRate])

        nexttile;
        plot(time,yFitted)
        hold on
        scatter(time(indx), pks)
        xlim([startframe/frmRate;endframe/frmRate])

        title(['Average heart beat = ' num2str(av_heart_beat) ...
            ' - Average Pulsatility Index = ' num2str(mean(PI, 'omitnan'))])
    end



    %% Save
    % as table to fit to results.mat
    pulsatility_output.Pulsatility_mean(ind_kymo) = mean(PI, 'omitnan');
    pulsatility_output.Pulsatility_median(ind_kymo) = median(PI, 'omitnan');

    % in kymoROI file:
    Pulsatility_calc.yFitted = yFitted;
    Pulsatility_calc.pulsatility = PI;
    Pulsatility_calc.meanpuls = mean(PI, 'omitnan');
    Pulsatility_calc.peakinds = indx;
    Pulsatility_calc.peaks = pks;

    save([DataFolder ROIname], 'Pulsatility_calc', '-append');


    %% FFT - freq spectrum
    y = fft(vel);
    n = length(vel);
    vel_freq = size(vel,2)/seconds; %freq of velocity calc
    f = (0:n-1)*(vel_freq/n); % freq range
    power = abs(y).^2/n; % power of DFT

    y0 = fftshift(y);
    f0 = (-n/2:n/2-1)*(vel_freq/n);
    power0 = abs(y0).^2/n;

    if show_plots
        f1 = figure;
        t = tiledlayout('flow');
        title(t, [Mouse ' ' Acq ' ' ROIname], 'Interpreter', 'none')

        nexttile
        plot(f(2:end), power(2:end)) % first one is insanely high so plot from second one
        xlabel('Frequency')
        ylabel('Power')

        nexttile
        plot(f0, power0)
        xlabel('Frequency')
        ylabel('Power')
        title('shifted')

        % close(f1)
    end


    clear peak* kymoImg Velocity_calc ROI_type Pulsatility_calc vel vel_freq vel_mask seconds pulsatility ROIname
end
end


% function plot_pulsatility(velocity, pulsatility)
% [peakvals, peakinds] = findpeaks(velocity);
%
% yyaxis left
% plot(velocity)
% hold on
% plot(peakinds, peakvals, 'LineStyle', 'none', 'Marker','o')
% ylabel('Velocity (mm/s)')
% yyaxis right
% plot(peakinds, pulsatility(1:length(peakinds)), 'LineStyle', 'none', 'Marker', '*')
% ylabel('Pulsatility')
%
% subtitle(['Pulsatility = ' num2str(mean(pulsatility(1:length(peakinds)), 'omitnan'))])
% end




% %% old
% % get heartbeats (hb)
%
% % get peak, and compare to the dip that comes after (could also take
% % the one that comes before). Don't take the unreliable values of the
% % velocity calculation (vel_mask).
% [peakvals, peakinds] = findpeaks(vel);
% % sys_dias_diff = NaN(length(peakinds),1);
% pulsatility = NaN(length(peakinds),1);
% for indhb = 1:length(peakinds)-1
%     hb_cycle = vel(peakinds(indhb):peakinds(indhb+1)-1); %dont take the second peak
%
%     if anynan(hb_cycle)
%         disp([num2str(indhb) ' contains nan, skipped'])
%         continue
%     end
%
%     % check if heartbeat cycle is an expected length of time
%     hb_time = (length(hb_cycle))/vel_freq;
%     hb_freq = 1/hb_time;
%
%     % let's take mouse heartbeat as between 300 and 700 beats per minute
%     % (anesthetized slower). Thats (300/60=) 5 to (700/60=) 12 beats per sec
%     if hb_freq < 5
%         % warning('Heartbeat of the mouse seems too slow!')
%         disp([num2str(indhb) ' Heartbeat of the mouse seems too slow! Skipped'])
%         continue
%     elseif  hb_freq > 12
%         % warning('Heartbeat of the mouse seems too fast!')
%         % disp(indhb)
%         disp([num2str(indhb) ' Heartbeat of the mouse seems too fast! Skipped'])
%         continue
%     end
%
%     hb_peak = peakvals(indhb);
%     hb_dip = min(hb_cycle);
%     % sys_dias_diff(indhb) = hb_peak-hb_dip;
%     hb_pulsatility = (hb_peak-hb_dip)/mean(hb_cycle);
%
%     % if pulsatility is negative or very high (can be for example when
%     % velocity is close to 0
%     if hb_pulsatility<0 || hb_pulsatility>5
%         disp([num2str(indhb) ' Pulsatility odd, skipped'])
%     else
%         pulsatility(indhb) = hb_pulsatility;
%     end
% end
% clear hb_cycle hb_dip hb_freq hb_peak hb_pulsatility hb_time indhb
%
% if show_plots
%
%     f1 = figure;
%     t = tiledlayout(2,1);
%
%     nexttile
%     % plot all
%     plot_pulsatility(vel, pulsatility)
%     title('Whole acquisition')
%
%     nexttile
%     % plot 2 sec of velocity like in paper
%     plot_pulsatility(vel(1:round(size(vel,2)/seconds*2)), pulsatility)
%     title('2 seconds')
%
%     title(t, [Mouse ' ' Acq ' ' ROIname], 'Interpreter', 'none')
%     subtitle(t, ROI_type, 'Interpreter', 'none')
% end


