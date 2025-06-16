% Give Mice like {'M04', 'M05}
% This function always checks for velocity/pulsatility/depth etc. and
% overwrites the latest value, regardless of the value of the variable
% overwrite. The variable overwrite is for choosing the best acquisition in
% the coupled_acqs. 
% This automatically takes kymoROI_1! that's hardcoded.

function [Results] = UpdateResults(SaveDirectory, Mice, overwrite)

% set-up
if ~strcmp(SaveDirectory(end), filesep)
    SaveDirectory = [SaveDirectory filesep];
end

if exist([SaveDirectory 'Results.mat'], 'file')
    load([SaveDirectory 'Results.mat'], 'Results');
else
    varNames = {'Mouse', 'Group', 'Acq', 'ScanType', 'Vessel', 'Diameter', 'MeanVelocity_xcorr', 'Pulsatility_xcorr', 'MeanVelocity_fft', 'Pulsatility_fft', 'Depth', 'UseAcq', 'Markertype'};
    sz = [0 size(varNames,2)];
    varTypes = {'cell', 'categorical' 'cell', 'cell', 'cell', 'double', 'double', 'double', 'double', 'double', 'double', 'categorical', 'cell'};
    Results = table('size', sz, 'VariableTypes', varTypes, 'VariableNames',varNames);
    clear sz varTypes varNames
end

if ~exist('Mice', 'var') % if you didnt specify mice, take all of them
    Mice = struct2table(dir([SaveDirectory, 'M*']));
    Mice = Mice(Mice.isdir,:);
    Mice = Mice.name;
end

if ~exist('overwrite', 'var')
    overwrite = 0;
%     answeroverwrite = 'No';
% elseif overwrite == 1
%     answeroverwrite = questdlg('Do you want to recheck which acquisition was the best?');
end

breakpointanswer = 'No';


%% Go per mouse
for indMouse = 1:length(Mice)
    MouseFolder = [SaveDirectory Mice{indMouse} filesep];
    Acquisitions = struct2table(dir([MouseFolder, 'A*']));
    Acquisitions = Acquisitions(Acquisitions.isdir,:);
    Acquisitions = Acquisitions.name;

    % get mouse results if they existed
    if any(contains(Results.Mouse, Mice{indMouse})) && overwrite == 0
        MouseResults = Results(matches(Results.Mouse, Mice{indMouse}),:);
    else
        varNames = {'Mouse', 'Group', 'Acq', 'ScanType', 'Vessel', 'Diameter',  'MeanVelocity_xcorr', 'Pulsatility_xcorr', 'MeanVelocity_fft', 'Pulsatility_fft','Depth', 'UseAcq', 'Markertype'};
        sz = [0 size(varNames,2)];
        varTypes = {'cell', 'categorical' 'cell', 'cell', 'cell', 'double', 'double', 'double', 'double', 'double', 'double','categorical', 'cell'};
        MouseResults = table('size', sz, 'VariableTypes', varTypes, 'VariableNames',varNames);
        clear sz varTypes varNames
    end


    %% Go per Acq.
    for indAcq = 1:length(Acquisitions)

        % If the acq was already in results, update de info 
        if any(strcmp(MouseResults.Acq, Acquisitions(indAcq))) 
            tablewriteind = find(strcmp(MouseResults.Acq, Acquisitions(indAcq)));
        % If you didnt start any of the analysis yet, skip it
        elseif ~exist([MouseFolder Acquisitions{indAcq} filesep 'AcqInfos.mat'], 'file')
            continue
        % otherwise, add acquisition to the end of the table
        else
            tablewriteind = size(MouseResults,1)+1;
        end

        % updates results
        try
            [MouseResults] = AddAcqToTable(MouseFolder, Acquisitions{indAcq}, MouseResults, tablewriteind);
        catch
            disp([Mice{indMouse} ' ' Acquisitions{indAcq} ' not done yet, skipped.'])
            continue
        end

        % if useacq is already there, skip
        if ~isundefined(MouseResults.UseAcq(tablewriteind))
            continue
        end

        %% See if current acquisition is useable
        load([MouseFolder Acquisitions{indAcq} filesep 'AcqInfos.mat'], 'AcqInfoStream')

        % if there is no ROI, skip
        if ~exist([MouseFolder Acquisitions{indAcq} filesep 'kymoROI_1.mat'], 'file')
            continue
        else 
            load([MouseFolder Acquisitions{indAcq} filesep 'kymoROI_1.mat'], 'ROI_info');
        end

        % if it is not indicated whether velocity can be calculated, skip
        if ~isfield(ROI_info, 'calculate_velocity')
            continue
         % if acq not useable for pulsatility/velocity calc, skip         
        elseif matches(ROI_info.calculate_velocity, 'No')
            MouseResults.UseAcq(tablewriteind) = 'No';
            continue
        % if there are no coupled acqs, use this one
        elseif ~isfield(AcqInfoStream, 'Coupled_acq') || isempty(AcqInfoStream.Coupled_acq)
            MouseResults.UseAcq(tablewriteind) = 'Yes';
            continue
        end
        

        %% start with coupled acquisitions
        coupled_acq_list = [AcqInfoStream.Coupled_acq; AcqInfoStream.DatasetName];
        Coupled_table = MouseResults(contains(MouseResults.Acq, coupled_acq_list),:);
                
        % Check if all coupled acq already have a UseAcq label (ie you already chose). If so, skip
        if all(contains(coupled_acq_list, MouseResults.Acq)) && ...
            ~any(ismissing(Coupled_table.UseAcq))
            continue
        end

        % If ROI_info indicates "no" for velocity calc, take away from
        % coupled list and add to results with useacq 'No'
        new_coupled_acq_list = coupled_acq_list;
        for c_acq_ind = 1:length(coupled_acq_list)
            if exist([MouseFolder coupled_acq_list{c_acq_ind} filesep 'kymoROI_1.mat'], 'file')
                load([MouseFolder coupled_acq_list{c_acq_ind} filesep 'kymoROI_1.mat'], 'ROI_info');

                if isfield(ROI_info, 'calculate_velocity') && matches(ROI_info.calculate_velocity, 'No')
                    if any(strcmp(MouseResults.Acq, coupled_acq_list(c_acq_ind)))
                        tablewriteind_coupled = find(strcmp(MouseResults.Acq, coupled_acq_list(c_acq_ind)));
                    else
                        tablewriteind_coupled = size(MouseResults,1)+1;
                        [MouseResults] = AddAcqToTable(MouseFolder, coupled_acq_list{c_acq_ind}, MouseResults, tablewriteind_coupled);
                    end
                    MouseResults.UseAcq(tablewriteind_coupled) = 'No';

                    new_coupled_acq_list(c_acq_ind) = {'del'};
                end
            end
        end
        new_coupled_acq_list(matches(new_coupled_acq_list, 'del')) = [];
        coupled_acq_list = new_coupled_acq_list;
        clear Coupled_table AcqInfoStream c_acq_ind ROI_info tablewriteind_coupled new_coupled_acq_list

        % if youre only left with one good acquisition, use that one
        if length(coupled_acq_list) == 1
            if any(strcmp(MouseResults.Acq, coupled_acq_list(1)))
                tablewriteind = find(strcmp(MouseResults.Acq, coupled_acq_list(1)));
            else
                tablewriteind = size(MouseResults,1)+1;
                [MouseResults] = AddAcqToTable(MouseFolder, coupled_acq_list{1}, MouseResults, tablewriteind);
            end
            MouseResults.UseAcq(tablewriteind) = 'Yes';
            continue
        end


        %% plot - compare acquisitions of same vessel
        % plot all kymographs of coupled acqs.
        f1 = figure('units','normalized','outerposition',[0 0 1 1]);
        tiledlayout('Vertical');
        secs_plotted =0.5;
        for c_acq_ind = 1:length(coupled_acq_list)
            if exist([MouseFolder coupled_acq_list{c_acq_ind} filesep 'kymoROI_1.mat'], 'file') && ...
                    exist([MouseFolder coupled_acq_list{c_acq_ind} filesep 'AcqInfos.mat'], 'file')
                
                load([MouseFolder coupled_acq_list{c_acq_ind} filesep 'kymoROI_1.mat'], 'kymoImg', 'PixelSize', 'ROI_info');
                load([MouseFolder coupled_acq_list{c_acq_ind} filesep 'AcqInfos.mat'], 'AcqInfoStream');
                if isfield(AcqInfoStream, 'FrameRateHzLinescan')
                    frmRate = AcqInfoStream.FrameRateHzLinescan;
                else
                    frmRate = AcqInfoStream.FrameRateHz;
                end
                
                % if you didnt get the kymo for the full timespan, pad with nans
                if size(kymoImg,1)<round(frmRate*secs_plotted)
                    kymoImg = [kymoImg; NaN(round(frmRate*secs_plotted)-size(kymoImg, 1), size(kymoImg,2))];
                end

                if exist('PixelSize', 'var')
                    Plot_Kymograph(kymoImg, frmRate, PixelSize.pxlSize, secs_plotted, [], 1);
                else
                    Plot_Kymograph(kymoImg, frmRate, [], secs_plotted, [], 1);
                end
                
                if isfield(ROI_info, 'calculate_velocity')
                    title([coupled_acq_list{c_acq_ind} ' - Calculate velocity: ' ROI_info.calculate_velocity]);
                else
                    title([coupled_acq_list{c_acq_ind} ' - Calculate velocity: unknown']);
                end
                clear AcqInfoStream kymoImg PixelSize ROI_info

            else
                nexttile
                %keep plot emtpy
            end
        end

        clear secs_plotted c_acq_ind
        [best_kymo_ind] = listdlg('PromptString', {'Which kymograph looks the best? (Cancel if none)'},...
            'ListString', coupled_acq_list);
        close(f1)


        % Add codes to table
        UseAcqCodes = zeros(1,length(coupled_acq_list));
        UseAcqCodes(best_kymo_ind) = 1;
        CodeOptions = {'No', 'Yes'}; 

        for c_acq_ind = 1:length(coupled_acq_list)
            if any(strcmp(MouseResults.Acq, coupled_acq_list(c_acq_ind)))
                tablewriteind = find(strcmp(MouseResults.Acq, coupled_acq_list(c_acq_ind)));
            else
                tablewriteind = size(MouseResults,1)+1;
                [MouseResults] = AddAcqToTable(MouseFolder, coupled_acq_list{c_acq_ind}, MouseResults, tablewriteind);
            end
            MouseResults.UseAcq(tablewriteind) = CodeOptions(UseAcqCodes(c_acq_ind)+1);
        end
        clear c_acq_ind best_kymo_ind UseAcqCodes CodeOptions frmRate tablewriteind coupled_acq_list f1


        % give breakpoint for yourself
        breakpointanswer = questdlg('Do you want to break?');
        if matches(breakpointanswer, 'Yes')
            break
        end

    end

    % Delete old Mouse results from table, replace with new
    Results(contains(Results.Mouse, Mice{indMouse}),:) = [];
    Results = [Results; MouseResults];
save([SaveDirectory 'Results_checkpoint.mat'], 'Results');
    if matches(breakpointanswer, 'Yes')
        break
    end

end

%% Get marker type per mouse
% note: have to do another mice loop because you need all of them here, not
% just the one you want to update
markerstyles = {"o","+","*",".","x","_","|","square","diamond","^","v",">","<","pentagram","hexagram"};
Mice = sort(unique(Results.Mouse));
for indmouse = 1:length(Mice) 
    Results(matches(Results.Mouse, Mice{indmouse}),'Markertype') = cellstr(markerstyles(indmouse));
end

%% save
save([SaveDirectory 'Results.mat'], 'Results')

% save backup
save(['E:\PhD\P3 - Ultrafast two photon' filesep 'Results_backup_' num2str(day(datetime)) '-' num2str(month(datetime)) '-' num2str(year(datetime)) '.mat'], 'Results')

end




function [CurrentTable] = AddAcqToTable(MouseFolder, Acq, CurrentTable, tablewriteind)
% Don't give Acq as cell

if ~exist("tablewriteind", 'var')
    tablewriteind = size(CurrentTable,1)+1;
end

load([MouseFolder Acq filesep 'AcqInfos.mat'], 'AcqInfoStream')

warning('off')
CurrentTable.Mouse(tablewriteind) = {AcqInfoStream.Mouse};
CurrentTable.Group(tablewriteind) = AcqInfoStream.Group;
CurrentTable.Acq(tablewriteind) = {Acq};
CurrentTable.ScanType(tablewriteind) = {AcqInfoStream.Acquisition_Type};
CurrentTable.Depth(tablewriteind) = AcqInfoStream.Depth;
CurrentTable.Vessel(tablewriteind) = {AcqInfoStream.Vessel};

% Add Diameter if it's calculated
if isfield(AcqInfoStream, 'Diameter') && length(AcqInfoStream.Diameter) == 1
    CurrentTable.Diameter(tablewriteind) = AcqInfoStream.Diameter;
elseif isfield(AcqInfoStream, 'Diameter')
    CurrentTable.Diameter(tablewriteind) = 1000;
else
    CurrentTable.Diameter(tablewriteind) = NaN;
end

if exist([MouseFolder Acq filesep 'kymoROI_1.mat'], 'file')
    load([MouseFolder Acq filesep 'kymoROI_1.mat'], 'Velocity_calc', 'Pulsatility_calc')

    % Velocity xcorr
    if any(matches(who('-file', [MouseFolder Acq filesep 'kymoROI_1.mat']), 'Velocity_calc')) &&...
            isfield(Velocity_calc, 'xcorr_range')
        CurrentTable.MeanVelocity_xcorr(tablewriteind) = Velocity_calc.xcorr_range.meanvel;
    else
        CurrentTable.MeanVelocity_xcorr(tablewriteind) = NaN;
    end

    % Velocity fft
    if any(matches(who('-file', [MouseFolder Acq filesep 'kymoROI_1.mat']), 'Velocity_calc')) &&...
            isfield(Velocity_calc, 'fft')
        CurrentTable.MeanVelocity_fft(tablewriteind) = Velocity_calc.fft.meanvel;
    else
        CurrentTable.MeanVelocity_fft(tablewriteind) = NaN;
    end

    % Pulsatility xcorr
    if any(matches(who('-file', [MouseFolder Acq filesep 'kymoROI_1.mat']), 'Pulsatility_calc')) &&...
            isfield(Pulsatility_calc, 'xcorr_range')
        CurrentTable.Pulsatility_xcorr(tablewriteind) = mean(Pulsatility_calc.xcorr_range.PI, 'omitnan'); 
    else
        CurrentTable.Pulsatility_xcorr(tablewriteind) = NaN;
    end

    % Pulsatility fft
    if any(matches(who('-file', [MouseFolder Acq filesep 'kymoROI_1.mat']), 'Pulsatility_calc')) &&...
            isfield(Pulsatility_calc, 'fft')
        CurrentTable.Pulsatility_fft(tablewriteind) = mean(Pulsatility_calc.fft.PI, 'omitnan'); 
    else
        CurrentTable.Pulsatility_fft(tablewriteind) = NaN;
    end

else
    CurrentTable.MeanVelocity_fft(tablewriteind) = NaN;
    CurrentTable.Pulsatility_fft(tablewriteind) = NaN;
    CurrentTable.MeanVelocity_xcorr(tablewriteind) = NaN;
    CurrentTable.Pulsatility_xcorr(tablewriteind) = NaN;
end

warning('on')

end