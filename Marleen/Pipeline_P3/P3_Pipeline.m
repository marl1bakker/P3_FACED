%% P3_Pipeline
%% One time calibration of x-axis
% CalibrationFolder = 'D:\FACED\Raw Data\M03\calibration\';
% ReadRawFACEDFiles(CalibrationFolder, CalibrationFolder, 1, 1, 0, 0, 10); % take av over 10 frames


%% Add recording to overview (have to do "manually" per mouse)
SaveDirectory = 'D:\FACED\Data P3\';
% AddRecording('D:\FACED\Raw Data\M01\', SaveDirectory);
% AddRecording('D:\FACED\Raw Data\M02\', SaveDirectory);
% AddRecording('D:\FACED\Raw Data\M04\', SaveDirectory);

%% Go per Recording
load([SaveDirectory 'Overview.mat'], 'Overview', 'Pipeline_check')

% for ind = 87:size(Overview,1)
for ind = 129:size(Overview,1)
    Mouse = Overview.Mouse{ind};
    Acq = Overview.Acq{ind};
    check_ind = find(matches(Pipeline_check.Mouse, Mouse).*matches(Pipeline_check.Acq, Acq));
    DataFolder = Overview.DataFolder{ind};
    RawDataFolder = Overview.RawDataFolder{ind};

    disp(['Pipeline for ' Mouse ' ' Acq]);
      
    % Make faced.dat file
    current_step = 'ReadRawFACEDFiles';
    exist_step = any(strcmp(current_step,Pipeline_check.Properties.VariableNames));
    if ~exist_step || isnat(Pipeline_check.(current_step)(check_ind)) ...
            || Pipeline_check.(current_step)(check_ind) < datetime('9-Mar-2025')
        try
            ReadRawFACEDFiles(RawDataFolder, DataFolder, 1, 1, 0);

            Pipeline_check.(current_step)(check_ind) = datetime;
            save([SaveDirectory, 'Overview.mat'], 'Pipeline_check', '-append')
            disp([current_step ' done!'])
        catch
            disp([current_step ' Error!'])
            continue
        end
    else
        disp([current_step ' already done.'])
    end

    % Add information to the info file
    current_step = 'Add_to_AcqInfos';
    exist_step = any(strcmp(current_step,Pipeline_check.Properties.VariableNames));
    if ~exist_step || isnat(Pipeline_check.(current_step)(check_ind))
        try
            AddInfo(DataFolder, 'Mouse', 'Change', Mouse);
            AddInfo(DataFolder, 'Comments', 'Add');
            AddInfo(DataFolder, 'Coupled_acq');
            
            Pipeline_check.(current_step)(check_ind) = datetime;
            save([SaveDirectory, 'Overview.mat'], 'Pipeline_check', '-append')
            disp([current_step ' done!'])
        catch
            disp([current_step ' Error!'])
        end
    else
        disp([current_step ' already done.'])
    end


    % Take off the part where the galvo has to reset and get average image
    current_step = 'Clean_Data';
    exist_step = any(strcmp(current_step,Pipeline_check.Properties.VariableNames));
    if ~exist_step || isnat(Pipeline_check.(current_step)(check_ind)) || ...
            Pipeline_check.Clean_Data(check_ind) < Pipeline_check.ReadRawFACEDFiles(check_ind) % hardcoded
        try
            Clean_Data(DataFolder)

            Pipeline_check.(current_step)(check_ind) = datetime;
            save([SaveDirectory, 'Overview.mat'], 'Pipeline_check', '-append')
            disp([current_step ' done!'])
        catch
            disp([current_step ' Error!'])
            continue
        end
    else
        disp([current_step ' already done.'])
        patch_code(DataFolder)
    end


    % vessel diameter
    current_step = 'Vessel_Diameter';
    exist_step = any(strcmp(current_step,Pipeline_check.Properties.VariableNames));
    if ~exist_step || isnat(Pipeline_check.(current_step)(check_ind))
        try
            Diameter_Vessel(DataFolder, 0)

            Pipeline_check.(current_step)(check_ind) = datetime;
            save([SaveDirectory, 'Overview.mat'], 'Pipeline_check', '-append')
            disp([current_step ' done!'])
        catch
            disp([current_step ' Error!'])
        end
    else
        disp([current_step ' already done.'])
    end
    

    % Make ROI for Kymograph
    current_step = 'ROI_Kymograph';
    exist_step = any(strcmp(current_step,Pipeline_check.Properties.VariableNames));
    if ~exist_step || isnat(Pipeline_check.(current_step)(check_ind)) || ...
            Pipeline_check.ROI_Kymograph(check_ind) < Pipeline_check.Clean_Data(check_ind) % hardcoded
        try
            ROI_Kymograph(DataFolder);

            Pipeline_check.(current_step)(check_ind) = datetime;
            save([SaveDirectory, 'Overview.mat'], 'Pipeline_check', '-append')
            disp([current_step ' done!'])
        catch
            disp([current_step ' Error!'])
            % continue
        end
    else
        disp([current_step ' already done.'])
    end


    % Plot Kymograph
    current_step = 'Plot_Kymograph';
    exist_step = any(strcmp(current_step,Pipeline_check.Properties.VariableNames));
    if ~exist_step || isnat(Pipeline_check.(current_step)(check_ind)) || ...
            Pipeline_check.Plot_Kymograph(check_ind) < Pipeline_check.ROI_Kymograph(check_ind) % hardcoded
        try
            Make_Kymograph(DataFolder, 1);

            Pipeline_check.(current_step)(check_ind) = datetime;
            save([SaveDirectory, 'Overview.mat'], 'Pipeline_check', '-append')
            disp([current_step ' done!'])
        catch
            disp([current_step ' Error!'])
            % continue
        end
    else
        disp([current_step ' already done.'])
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
    

    % Calculate speed of blood cells
    current_step = 'Linescan_Velocimetry';
    exist_step = any(strcmp(current_step,Pipeline_check.Properties.VariableNames));
    if ~exist_step || isnat(Pipeline_check.(current_step)(check_ind)) ||...
         Pipeline_check.Linescan_Velocimetry(check_ind) < Pipeline_check.Plot_Kymograph(check_ind) % hardcoded
        try
            % Linescan_Velocimetry(DataFolder, 'old', 0, 'auto_list', 1); 
            Linescan_Velocimetry(DataFolder, 'xcorr_range', 0, 'auto_list', 0); 

            Pipeline_check.(current_step)(check_ind) = datetime;
            save([SaveDirectory, 'Overview.mat'], 'Pipeline_check', '-append')
            disp([current_step ' done!'])

        catch
            disp([current_step ' Error!'])
        end
    else
        disp([current_step ' already done.'])
    end

    % Calculate pulsatility
    current_step = 'Pulsatility';
    exist_step = any(strcmp(current_step,Pipeline_check.Properties.VariableNames));
    if ~exist_step || isnat(Pipeline_check.(current_step)(check_ind)) ||...
            Pipeline_check.Pulsatility(check_ind) < Pipeline_check.Linescan_Velocimetry(check_ind) % hardcoded
        try
            [pulsatility_output] = Pulsatility_Index(DataFolder, 'xcorr_range', 1, 'auto_list', 1);

            Pipeline_check.(current_step)(check_ind) = datetime;
            save([SaveDirectory, 'Overview.mat'], 'Pipeline_check', '-append')
            disp([current_step ' done!'])
        catch
            disp([current_step ' Error!'])
        end
    else
        disp([current_step ' already done.'])
    end

 
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



end


%% For all Acquisitions:
UpdateOverview(Overview, SaveDirectory)
load([SaveDirectory 'Overview.mat'], 'Overview')

% update results per mouse
% UpdateResults('D:\FACED\Raw Data\M01\', SaveDirectory);
% UpdateResults('D:\FACED\Raw Data\M02\', SaveDirectory);
[Results] = UpdateResults(SaveDirectory, {'M03'});
[Results] = UpdateResults(SaveDirectory, {'M04'});


%% Plot things
Overview_Kymographs(Overview)

% single acq rbc tracking:
plot_rbc_tracking('D:\FACED\Data P3\M04\A21', 1) % for methods fig

SaveDirectory = 'D:\FACED\Data P3\';
load([SaveDirectory 'Results.mat'], 'Results')


figure('Color',  'white'); tiledlayout('flow')
Scatter_TAC_Sham(Results, 'Diameter', 'MeanVelocity', 1, 1); % as tile, tilenr
Scatter_TAC_Sham(Results, 'MeanVelocity', 'Pulsatility', 1, 2);
Scatter_TAC_Sham(Results, 'Diameter', 'Pulsatility', 1, 3);

Boxplot_TAC_Sham(Results)