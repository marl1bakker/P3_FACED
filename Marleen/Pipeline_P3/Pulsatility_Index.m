% TODO SAVE PULSATILITY IN OVERVIEW.MAT RESULTS.PULSATILITY

% Pulsatility index defined as difference between peak systolic flow and
% minimum diastolic flow / mean velocity
%
% Could not find original code of Meng et al to calculate this, so made
% from scratch.
% Marleen Bakker 18-10-24 // 23-01-25 // 20-3-25 (PI close to/below 0)

% "method" is the method from linescan_velocimetry that you want to use.

function [pulsatility_output] = Pulsatility_Index(DataFolder, method, show_plots, ROIname, overwrite)


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

if ~exist('method', 'var')
    method = 'xcorr';
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
        'VariableNames', {'Name', 'Type', 'PI_mean', 'PI_median', 'Heartbeat'},...
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
    load([DataFolder ROIname], 'Velocity_calc', 'kymoImg', 'ROI_type', 'Pulsatility_calc', 'no_rbc', 'ROI_info', 'PixelSize');
    warning('on');

    pulsatility_output.Name(ind_kymo) = ROIname;
    pulsatility_output.Type(ind_kymo) = ROI_type;

    %% check if already done
    if ~exist('Velocity_calc', 'var') || ~isfield(Velocity_calc, method)
        error(['Velocity not found for ' ROIname ' of ' Mouse ', ' Acq '. Run Linescan_Velocimetry first.']);
        % continue

    elseif ~contains(ROI_type, 'line') || matches(ROI_info.calculate_velocity, 'No')
        % skip if it's a block ROI (didnt make code for that yet) or if you
        % indicated the kymo is not good enough for velocity calculation
        continue

    elseif exist('Pulsatility_calc', 'var') && isfield(Pulsatility_calc, method) ...
            && overwrite == 0
        disp(['Pulsatility already calculated for ' ROIname ' ' Mouse ' ' Acq '. ROI skipped.'])

        pulsatility_output.PI_mean(ind_kymo) = mean(Pulsatility_calc.(method).PI, 'omitnan');
        pulsatility_output.PI_median(ind_kymo) = median(Pulsatility_calc.(method).PI, 'omitnan');
        pulsatility_output.Heartbeat(ind_kymo) = Pulsatility_calc.(method).heartbeat;
        continue

    elseif exist('Pulsatility_calc', 'var') && isfield(Pulsatility_calc, method) ...
            && overwrite == 1
        disp(['Pulsatility already calculated for ' ROIname ' ' Mouse ' ' Acq '. OVERWRITING PULSATILITY CALCULATIONS.'])

    end

    %% load parameters
    if matches(method, 'xcorr')
        vel = Velocity_calc.(method).raw_velocity;
    else
        vel = Velocity_calc.(method).velocity;
    end

    if size(vel,1) < size(vel,2)
        vel = vel';
    end

    if mean(vel, 'omitnan') < 0
        vel = vel * -1;
    end

    if isfield(AcqInfoStream, 'FrameRateHzLinescan')
        frmRate = AcqInfoStream.FrameRateHzLinescan;
    else
        frmRate = AcqInfoStream.FrameRateHz;
    end
    % vel_freq = size(vel,2)/seconds; %freq of velocity calc

    %% Exclude bad fits & parts without RBC
    % Already removed the bad vals for xcorr before doing movmean, so take
    % the nans as bad values. Badvals are correlations with low value.
    if matches(method, 'old') || matches(method, 'fft')
        badvals = Velocity_calc.(method).badvals;
    elseif matches(method, 'new') || matches(method, 'xcorr') || matches(method, 'xcorr_range')
        badvals = find(isnan(vel));
    end



    %% Fit curve
    % disp('Fitting curve over velocity...');
    % seconds = size(kymoImg, 1)/frmRate;
    % time = linspace(0,seconds,length(vel))';
    % opts = fitoptions('Method', 'SmoothingSpline');
    % % opts.Exclude = find(badvals);
    % opts.Exclude = badvals;
    % opts.SmoothingParam = 0.99999; % 5 times 9 seems best
    % if matches(method, 'xcorr')
    %     w = Velocity_calc.(method).corr_values;
    %     w(w<0) = 0;
    %     opts.Weights = w;
    % % elseif matches(method, 'xcorr_range') % less noise when ignoring weights.
    % %     w = Velocity_calc.(method).weights;
    % %     w(isnan(w)) = 0;
    % %     opts.Weights = w;
    % end
    % [curve, ~] = fit(time,vel, "SmoothingSpline", opts);
    % yFitted = feval(curve, time);

    % Cannot fit all 20 seconds in one go. Go in periods
    disp('Fitting curve over velocity...');
    opts = fitoptions('Method', 'SmoothingSpline');
    opts.SmoothingParam = 0.99999; % 5 times 9 seems best
    if matches(method, 'xcorr')
        w = Velocity_calc.(method).corr_values;
        w(w<0) = 0;
        opts.Weights = w;
    end

    b = waitbar(0, 'Fitting curve...');
    periodsecs = 2; % take 2 sec to fit curve on
    periodfrms = frmRate*periodsecs;
    start_frame = 1;
    yFitted = NaN(size(vel));
    for indperiod = 1:floor(length(vel)/periodfrms)
        end_frame = start_frame+round(periodfrms)-1;

        periodvel = vel(start_frame:end_frame);
        time = linspace(start_frame/frmRate, (start_frame/frmRate)+periodsecs, ...
            length(periodvel))';
        periodbadvals = badvals(badvals>=start_frame);
        periodbadvals = periodbadvals(periodbadvals<=end_frame)-start_frame+1;
        opts.Exclude = periodbadvals;

        [curve, ~] = fit(time,periodvel, "SmoothingSpline", opts);
        yFitted(start_frame:end_frame) = feval(curve, time);

        start_frame = round(start_frame+periodfrms);

        waitbar(indperiod/(floor(length(vel)/periodfrms)), b);
    end
    close(b);

    secs = length(vel)/frmRate;
    time = linspace(0,secs,length(vel))';

    % get peaks
    [pks, indx] = findpeaks(yFitted, 'MinPeakDistance',1/(700/60)*frmRate); % assume heartbeat < 700 (700 would already be way too high for anesthetized mouse)
    av_heart_beat = length(pks)/secs*60;

    % % check
    % figure;plot(time, vel);
    % hold on; plot(time,yFitted);
    % scatter(time(indx), pks, 'Color', 'black', 'Marker', 'v', 'MarkerEdgeColor', 'black')
    % xlabel('seconds')
    % ylabel('velocity mm/sec')

    % If heartbeat is off:
    if av_heart_beat < 300 || av_heart_beat > 550
        disp(['Heartbeat of ' num2str(av_heart_beat) ' detected, which is outside the expected range.'])
        % return
        saveresult = 0;
    else
        pulsatility_output.Heartbeat(ind_kymo) = av_heart_beat;
        saveresult = 1;
    end

    % get pulsatility index
    PI = NaN(size(pks));
    if size(pks)<1
        error('Peaks not detected')
    end
    for ind = 1:length(pks)-1
        % take peak, find low, calculate average
        peak = pks(ind);
        heartcycle = yFitted(indx(ind):indx(ind+1)-1);
        low = min(heartcycle);
        av_heartcycle = mean(heartcycle);

        if peak<0 || av_heartcycle<0.1 % if the peak is below 0, the entire heartbeat should be reversed so this doesnt make sense
            continue
        end

        % if av_heartcycle<0.1 % skip if too close to zero, will give weird values
        %     continue
        % else
        PI(ind) = (peak-low)/abs(av_heartcycle);
        % end
    end



    %% FFT - freq spectrum
    [P1, f] = Frequency_spectrum_velocity(vel, frmRate);

    %% Save
    % as table to fit to results.mat
    pulsatility_output.PI_mean(ind_kymo) = mean(PI, 'omitnan');
    pulsatility_output.PI_median(ind_kymo) = median(PI, 'omitnan');

    if saveresult == 0
        continue
    end
    % in kymoROI file:
    Pulsatility_calc.(method).yFitted = yFitted;
    Pulsatility_calc.(method).PI = PI;
    Pulsatility_calc.(method).peakinds = indx;
    Pulsatility_calc.(method).peaks = pks;
    Pulsatility_calc.(method).heartbeat = av_heart_beat;

    save([DataFolder ROIname], 'Pulsatility_calc', '-append');




    %% Visualize
    if show_plots
        startframe_plot = 1;
        secsplot = 2;
        endframe_plot = round(startframe_plot+secsplot*frmRate);

        f2 = figure('Color', 'white');
        tiledlayout('vertical');

        Plot_Kymograph(kymoImg, frmRate, PixelSize.pxlSize, secsplot, [], 1, 1);
        title([Mouse ' ' Acq])

        nexttile
        if isfield(Velocity_calc, 'xcorr_range')
            vel = abs(Velocity_calc.xcorr_range.velocity);
            if size(vel,1) < size(vel,2)
                vel = vel';
            end
            time = linspace(0,secs,length(vel))';

            plot(time, vel, '.')
            xlim([startframe_plot/frmRate;endframe_plot/frmRate])
            xlabel('Seconds')
            ylabel('Velocity (mm/s)')
            title('Velocity xcorr_range')
        end

        nexttile;
        if isfield(Pulsatility_calc, 'xcorr_range')
            plot(time,Pulsatility_calc.xcorr_range.yFitted)
            hold on
            scatter(time(Pulsatility_calc.xcorr_range.peakinds), Pulsatility_calc.xcorr_range.peaks)
            xlim([startframe_plot/frmRate;endframe_plot/frmRate])
            title(['Average heart beat = ' num2str(Pulsatility_calc.xcorr_range.heartbeat) ...
                ' - Average Pulsatility Index = ' num2str(mean(Pulsatility_calc.xcorr_range.PI, 'omitnan'))])
            xlabel('Seconds'); ylabel('Fitted velocity xcorr_range');
        end

        nexttile
        if isfield(Velocity_calc, 'fft')
            vel = abs(Velocity_calc.fft.velocity);
            if size(vel,1) < size(vel,2)
                vel = vel';
            end
            time = linspace(0,secs,length(vel))';

            plot(time, vel, '.')
            xlim([startframe_plot/frmRate;endframe_plot/frmRate])
            xlabel('Seconds')
            ylabel('Velocity (mm/s)')
            title('fft')
        end

        nexttile;
        if isfield(Pulsatility_calc, 'fft')
            plot(time,Pulsatility_calc.fft.yFitted)
            hold on
            scatter(time(Pulsatility_calc.fft.peakinds), Pulsatility_calc.fft.peaks)
            xlim([startframe_plot/frmRate;endframe_plot/frmRate])
            title(['Average heart beat = ' num2str(Pulsatility_calc.fft.heartbeat) ...
                ' - Average Pulsatility Index = ' num2str(mean(Pulsatility_calc.fft.PI, 'omitnan'))])
            xlabel('Seconds'); ylabel('Fitted velocity fft/old');
        end

        nexttile;
        % semilogy(f, P1);
        plot(f, P1);
        ylim([0.001 1]); ylabel('Amplitude');
        xlim([1 100]); xlabel('f (Hz)');
        axis([0 100 10^-3 inf]);
        title('FFT on velocity');

        f2.Position = [500 60 700 780];
        pause(2);
        savefig(gcf, [DataFolder 'Velocity.fig']);
        close(f2)
    end

    clear peak* kymoImg Velocity_calc ROI_type Pulsatility_calc vel vel_freq vel_mask secs pulsatility ROIname
end
end

%%
function [P1, f] = Frequency_spectrum_velocity(velocity, framerate)

if size(velocity,1)>size(velocity,2)
    velocity = velocity';
end

velocity = velocity(1:find(~isnan(velocity),1,'last'));
if sum(isnan(velocity))>0
    disp('No FFT - still nan in velocity, fix!')
end

velocity = velocity - mean(velocity); %??
L = 2^(nextpow2(length(velocity))-2); % N in meng et al
f = framerate/L*(0:(L/2)); %Fs = L/T; f = Fs*(0:L/2)/L;
% Nyquist = frmRate/numavgs/2;
% fLim = min(Nyquist,40);
% fLim = 100;

% regular fft:
Y = fft(velocity, L);
P2 = abs(Y/L); % two sided spectrum
P1 = P2(1:floor(L/2)+1); % one sided spectrum
P1(2:end-1) = 2*P1(2:end-1);  %Multiply the spectrum in the positive frequencies by 2. You do not need to multiply P1(1) and P1(end) by 2 because these amplitudes correspond to the zero and Nyquist frequencies, respectively, and they do not have the complex conjugate pairs in the negative frequencies.

% f1 = figure;
% semilogy(f, P1);
% ylim([0.001 1]); ylabel('Amplitude');
% xlim([1 100]); xlabel('f (Hz)');
% axis([0 fLim 10^-3 inf]);
% title('fft on velocity')

% hold on
% P1_norm =  mat2gray( smoothdata( P1,'gaussian',3) );
% P1_smoothed = movmedian(P1, 2);
% semilogy(f, P1_norm);
% semilogy(f, P1_smoothed);
% legend({'P1', 'P1 norm', 'P1 smoothed'})
%



% %%
%
% % hanning window fft:
% f3 = figure;
% % hanvelocity = hanning(length(velocity))' .* velocity;
% hanvelocity = hann(length(velocity), 'periodic')' .* velocity;
% Y = fft(hanvelocity, L);
% P2 = abs(Y/L); % two sided spectrum
% P1 = P2(1:floor(L/2)+1); % one sided spectrum
% P1(2:end-1) = 2*P1(2:end-1);  %Multiply the spectrum in the positive frequencies by 2. You do not need to multiply P1(1) and P1(end) by 2 because these amplitudes correspond to the zero and Nyquist frequencies, respectively, and they do not have the complex conjugate pairs in the negative frequencies.
%
% P1_norm =  mat2gray( smoothdata( P1,'gaussian',3) );
% P1_smoothed = movmedian(P1, 2);
%
% semilogy(f, P1);
% ylim([0.001 1]); ylabel('Aplitude (a.u.)');
% xlim([1 100]); xlabel('f (Hz)');
% % hold on
% % semilogy(f, P1_norm);
% % semilogy(f, P1_smoothed);
%
% legend({'P1', 'P1 norm', 'P1 smoothed'})
% axis([0 fLim 10^-3 inf]);
% title('Hanning window fft')
%
% velocity_fft = P1;

end

