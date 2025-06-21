%% P3_Pipeline
% to do: M07 A73 opnieuw, get raw files, same for A16, A68, A72

%M09 A52 & M13 A25 gave problems, noted as 01-01-2050

%% One time calibration of x-axis
% CalibrationFolder = 'D:\FACED\Raw Data\M03\calibration\';
% ReadRawFACEDFiles(CalibrationFolder, CalibrationFolder, 1, 1, 0, 0, 10); % take av over 10 frames


%% Add recording to overview (have to do "manually" per mouse)
% SaveDirectory = 'D:\FACED\Data P3\';
% BackupDirectory = 'E:\PhD\P3 - Ultrafast two photon';
SaveDirectory = '/media/mbakker/data1/P3/Data P3/';
BackupDirectory = '/home/mbakker/Documents/';
% for ind = 1:size(Overview)
%     Overview.RawDataFolder(ind) = {['/media/mbakker/data1/P3/Raw Data/' Overview.Mouse{ind} filesep Overview.Acq{ind}]};
%     Overview.DataFolder(ind) = {['/media/mbakker/data1/P3/Data P3/' Overview.Mouse{ind} filesep Overview.Acq{ind} filesep]};
% end

% AddRecording('D:\FACED\Raw Data\M16\', SaveDirectory);
load([SaveDirectory 'Overview.mat'], 'Overview', 'Pipeline_check')
% save([BackupDirectory filesep 'Overview_backup_' num2str(day(datetime)) '-' num2str(month(datetime)) '-' num2str(year(datetime)) '.mat'], 'Overview', 'Pipeline_check')

%% Go per Recording

% for ind = 87:size(Overview,1)
for ind = 1:size(Overview,1)
% for ind = 1:390
     Mouse = Overview.Mouse{ind};
    Acq = Overview.Acq{ind};
    check_ind = find(matches(Pipeline_check.Mouse, Mouse).*matches(Pipeline_check.Acq, Acq));
    DataFolder = Overview.DataFolder{ind};
    RawDataFolder = Overview.RawDataFolder{ind};

    fprintf('\n');    fprintf('\n');
    disp('-------------------------------');
    disp(['Pipeline for ' Mouse ' ' Acq]);
      
    % Make faced.dat file
    % current_step = 'ReadRawFACEDFiles';
    exist_step = any(strcmp('ReadRawFACEDFiles',Pipeline_check.Properties.VariableNames));
    if ~exist_step || isnat(Pipeline_check.ReadRawFACEDFiles(check_ind)) ...
            || Pipeline_check.ReadRawFACEDFiles(check_ind) < datetime('9-Mar-2025')
        try
            % temp
             if ~exist([DataFolder filesep 'faced.dat'], 'file') 
                ReadRawFACEDFiles(RawDataFolder, DataFolder, 1, 1, 0);
             else
                s = dir([DataFolder filesep 'faced.dat']);
                if s.bytes == 0
                    ReadRawFACEDFiles(RawDataFolder, DataFolder, 1, 1, 0,1);
                end
                clear s
            end

            Pipeline_check.ReadRawFACEDFiles(check_ind) = datetime;
            save([SaveDirectory, 'Overview.mat'], 'Pipeline_check', '-append')
            disp('ReadRawFACEDFiles done!')
        catch
            disp('ReadRawFACEDFiles Error!')
            continue
        end
    else
        disp('ReadRawFACEDFiles already done.')
    end

    % continue

    % Add information to the info file
    % current_step = 'Add_to_AcqInfos';
    exist_step = any(strcmp('Add_to_AcqInfos',Pipeline_check.Properties.VariableNames));
    if ~exist_step || isnat(Pipeline_check.Add_to_AcqInfos(check_ind))
        try
            AddInfo(DataFolder, 'Mouse', 'Change', Mouse);
            AddInfo(DataFolder, 'Comments', 'Add');
            AddInfo(DataFolder, 'Coupled_acq');

            Pipeline_check.Add_to_AcqInfos(check_ind) = datetime;
            save([SaveDirectory, 'Overview.mat'], 'Pipeline_check', '-append')
            disp('Add_to_AcqInfos done!')
        catch
            disp('Add_to_AcqInfos Error!')
        end
    else
        disp('Add_to_AcqInfos already done.')
    end

% continue
    % Take off the part where the galvo has to reset and get average image
    % current_step = 'Clean_Data';
    exist_step = any(strcmp('Clean_Data',Pipeline_check.Properties.VariableNames));
    if ~exist_step || isnat(Pipeline_check.Clean_Data(check_ind)) || ...
            Pipeline_check.Clean_Data(check_ind) < Pipeline_check.ReadRawFACEDFiles(check_ind) 
        try
            % temp
            x = who('-file', [DataFolder filesep 'AcqInfos.mat']);
            if ~ismember('CleanData', x)
                Clean_Data(DataFolder)
            end
            clear x 

            Pipeline_check.Clean_Data(check_ind) = datetime;
            save([SaveDirectory, 'Overview.mat'], 'Pipeline_check', '-append')
            disp('Clean_Data done!')
        catch
            disp('Clean_Data Error!')
            continue
        end
    else
        disp( 'Clean_Data already done.')
        % patch_code(DataFolder)
    end
% continue

    % % vessel diameter
    % current_step = 'Vessel_Diameter';
    exist_step = any(strcmp('Vessel_Diameter',Pipeline_check.Properties.VariableNames));
    if ~exist_step || isnat(Pipeline_check.Vessel_Diameter(check_ind))
        try
            Diameter_Vessel(DataFolder, 0)

            Pipeline_check.Vessel_Diameter(check_ind) = datetime;
            save([SaveDirectory, 'Overview.mat'], 'Pipeline_check', '-append')
            disp('Vessel_Diameter done!')
        catch
            disp('Vessel_Diameter Error!')
        end
    else
        disp('Vessel_Diameter already done.')
    end


    % continue 

    % Make ROI for Kymograph
    % current_step = 'ROI_Kymograph';
    exist_step = any(strcmp('ROI_Kymograph',Pipeline_check.Properties.VariableNames));
    if ~exist_step || isnat(Pipeline_check.ROI_Kymograph(check_ind)) || ...
            Pipeline_check.ROI_Kymograph(check_ind) < Pipeline_check.Clean_Data(check_ind) % hardcoded
        try
            %temp
            if ~exist([DataFolder filesep 'kymoROI_1.mat'], 'file')
                ROI_Kymograph(DataFolder);
            end

            Pipeline_check.ROI_Kymograph(check_ind) = datetime;
            save([SaveDirectory, 'Overview.mat'], 'Pipeline_check', '-append')
            disp('ROI_Kymograph done!')
        catch
            disp('ROI_Kymograph Error!')
            % continue
        end
    else
        disp('ROI_Kymograph already done.')
    end


    % Plot Kymograph
    % current_step = 'Plot_Kymograph';
    exist_step = any(strcmp('Plot_Kymograph',Pipeline_check.Properties.VariableNames));
    if ~exist_step || isnat(Pipeline_check.Plot_Kymograph(check_ind)) || ...
            Pipeline_check.Plot_Kymograph(check_ind) < Pipeline_check.ROI_Kymograph(check_ind) ...% hardcoded
            || Pipeline_check.Plot_Kymograph(check_ind) < datetime('25-Mar-2025') %temp, because changed est velocity to be in ROI_info
        try
            Make_Kymograph(DataFolder, 1);

            Pipeline_check.Plot_Kymograph(check_ind) = datetime;
            save([SaveDirectory, 'Overview.mat'], 'Pipeline_check', '-append')
            disp('Plot_Kymograph done!')
        catch
            disp('Plot_Kymograph Error!')
            % continue
        end
    else
        disp('Plot_Kymograph already done.')
    end

% continue

    % % Check where there are red blood cells to measure
    % current_step = 'Detect_RBC';
    % exist_step = any(strcmp(current_step,Pipeline_check.Properties.VariableNames));
    % if ~exist_step || isnat(Pipeline_check.(current_step)(check_ind))
    %     try
    %         RBC_presence(DataFolder, 'auto_list', 1); 
    % 
    %         Pipeline_check.(current_step)(check_ind) = datetime;
    %         save([SaveDirectory, 'Overview.mat'], 'Pipeline_check', '-append')
    %         disp([current_step ' done!'])
    %     catch
    %         disp([current_step ' Error!'])
    %     end
    % else
    %     disp([current_step ' already done.'])
    % end    

    % continue


    % only do velocimetry/pulsatility for acquisitions that make sense...
    

    % % % Calculate speed of blood cells
    % % % current_step = 'Linescan_Velocimetry';
    % exist_step = any(strcmp('Linescan_Velocimetry',Pipeline_check.Properties.VariableNames));
    % if ~exist_step || isnat(Pipeline_check.Linescan_Velocimetry(check_ind)) ||...
    %      Pipeline_check.Linescan_Velocimetry(check_ind) < Pipeline_check.Plot_Kymograph(check_ind) % hardcoded
    %     try
    %         % Linescan_Velocimetry(DataFolder, 'old', 0, 'auto_list', 1); 
    %         Linescan_Velocimetry(DataFolder, 'xcorr_range', 0, 'auto_list', 0); 
    %         Linescan_Velocimetry(DataFolder, 'fft', 0, 'auto_list', 0); 
    % 
    %         Pipeline_check.Linescan_Velocimetry(check_ind) = datetime;
    %         save([SaveDirectory, 'Overview.mat'], 'Pipeline_check', '-append')
    %         disp( 'Linescan_Velocimetry done!')
    % 
    %     catch
    %         disp('Linescan_Velocimetry Error!')
    %     end
    % else
    %     disp('Linescan_Velocimetry already done.')
    % end
    % 
    % % Calculate pulsatility
    % % current_step = 'Pulsatility';
    % exist_step = any(strcmp('Pulsatility',Pipeline_check.Properties.VariableNames));
    % if ~exist_step || isnat(Pipeline_check.Pulsatility(check_ind)) ||...
    %         Pipeline_check.Pulsatility(check_ind) < Pipeline_check.Linescan_Velocimetry(check_ind) % hardcoded
    %     try
    %         [pulsatility_output] = Pulsatility_Index(DataFolder, 'xcorr_range', 1, 'auto_list', 0);
    %         [pulsatility_output] = Pulsatility_Index(DataFolder, 'fft', 1, 'auto_list', 1);
    % 
    %         Pipeline_check.Pulsatility(check_ind) = datetime;
    %         save([SaveDirectory, 'Overview.mat'], 'Pipeline_check', '-append')
    %         disp('Pulsatility done!')
    %     catch
    %         disp('Pulsatility Error!')
    %     end
    % else
    %     disp('Pulsatility already done.')
    % end
    % 

    % % continue
    % % Plot of transitioning cells through a 2D plunging vessel scan
    % current_step = '2D_plunging';
    % exist_step = any(strcmp(current_step,Pipeline_check.Properties.VariableNames));
    % if ~exist_step || isnat(Pipeline_check.(current_step)(check_ind))
    %     try
    %         Plot_Plunging(DataFolder); 
    % 
    %         Pipeline_check.(current_step)(check_ind) = datetime;
    %         save([SaveDirectory, 'Overview.mat'], 'Pipeline_check', '-append')
    %         disp([current_step ' done!'])
    %     catch
    %         disp([current_step ' Error!'])
    %     end
    % else
    %     disp([current_step ' already done.'])
    % end


    % Determine vessel type
    

end


%% For all Acquisitions:
try 
    UpdateOverview(Overview, SaveDirectory)
    load([SaveDirectory 'Overview.mat'], 'Overview')
catch
    disp('CHECK IF YOU HAVE ENOUGH SPACE!!!!! PROBABLY NOT!!')
    save([BackupDirectory filesep 'Overview_backup_' num2str(day(datetime)) '-' num2str(month(datetime)) '-' num2str(year(datetime)) '_XXX.mat'], 'Overview', 'Pipeline_check')
end

% update results per mouse
% UpdateResults('D:\FACED\Raw Data\M01\', SaveDirectory);
% UpdateResults('D:\FACED\Raw Data\M02\', SaveDirectory);
% [Results] = UpdateResults(SaveDirectory, {'M03'});
% [Results] = UpdateResults(SaveDirectory, {'M14'});
[Results] = UpdateResults(SaveDirectory, unique(Overview.Mouse));


%% Plot things
Overview_Kymographs(Overview)

% single acq rbc tracking:
% plot_rbc_tracking('D:\FACED\Data P3\M04\A21', 1) % for methods fig

%% Results.mat
% SaveDirectory = 'D:\FACED\Data P3\';
load([SaveDirectory 'Results.mat'], 'Results')
Results = Results(Results.UseAcq == 'Yes',:);

% make groups based on size
Results.Vesselgroup(Results.Diameter<10) = repmat({'Small'}, 1, size(Results(Results.Diameter<10,:),1));
Results.Vesselgroup(Results.Diameter>10) = repmat({'Large'}, 1, size(Results(Results.Diameter>10,:),1));
Results.Vesselgroup = categorical(Results.Vesselgroup);

% scatters
figure('Color',  'white'); tiledlayout('flow')
Scatter_TAC_Sham(Results, 'Diameter', 'MeanVelocity_fft', 1, 1); % as tile, tilenr
Scatter_TAC_Sham(Results, 'Diameter', 'MeanVelocity_xcorr', 1, 2); % as tile, tilenr
Scatter_TAC_Sham(Results, 'MeanVelocity_fft', 'Pulsatility_fft', 1, 3);
Scatter_TAC_Sham(Results, 'MeanVelocity_xcorr', 'Pulsatility_xcorr', 1, 4);
Scatter_TAC_Sham(Results, 'Diameter', 'Pulsatility_fft', 1, 5);
Scatter_TAC_Sham(Results, 'Diameter', 'Pulsatility_xcorr', 1, 6);

% boxplots
figure('Color',  'white'); tiledlayout('flow')
Boxplot_TAC_Sham(Results, 'Pulsatility_fft','Vesselgroup', 1, 1);
Boxplot_TAC_Sham(Results, 'Pulsatility_xcorr','Vesselgroup', 1, 2);
Boxplot_TAC_Sham(Results, 'MeanVelocity_fft','Vesselgroup', 1, 3);
Boxplot_TAC_Sham(Results, 'MeanVelocity_xcorr','Vesselgroup', 1, 4);



%% proof surgery worked
% % import TACstats
% TACstats.Mouse = cellstr(TACstats.Mouse);

% Get right markers
% note: M18 and M19 have no imaging because of bleeding. Give '*' marker and make gray
% SaveDirectory = 'D:\FACED\Data P3\';
load([SaveDirectory 'Results.mat'], 'Results')
load([SaveDirectory 'TACstats.mat'], 'TACstats')

for indMouse = 1:size(TACstats,1)
    try
        tempind = find(ismember(Results.Mouse, TACstats.Mouse{indMouse}));
        tempind = tempind(1);
        TACstats.Markertype(indMouse) = Results.Markertype(tempind);
    catch
        TACstats.Markertype(indMouse) = {'*'};
    end
end

% save
save('D:\FACED\Data P3\TACstats.mat', 'TACstats');

load('D:\FACED\Data P3\TACstats.mat');


figure('Color',  'white'); tiledlayout('flow')
% Boxplot_TAC_Sham(TACstats, 'HW_BW', 'Sex', 1, 1);  % as tile, tilenr
% ylabel('heartweight (mg)/bodyweight (g)')

f = figure('Color',  'white'); tiledlayout('flow')
Boxplot_TAC_Sham(TACstats, 'HW_BW', [], 1, 2);  % as tile, tilenr
ylabel('hw/bw in mg/g')
title('TAC surgery verification')
f.Position = [800  450  500  350];

% Stats
[h, p] = ttest2(TACstats(TACstats.Group == 'TAC',:).HW_BW, TACstats(TACstats.Group == 'Sham',:).HW_BW);
subtitle(['p = ' num2str(p)])

saveas(gcf, 'E:\PhD\P3 - Ultrafast two photon\Article\Figures and images\TAC_hw_bw.svg');