% Give rawdata folder of Mouse, not acquisition! 
% This function is made to go per mouse. If you add more acquisitions for a
% mouse later, but you already have this mouse in the overviewtable, the
% function will not add the new ones! You'll have to do that manually. 

% RawDataFolder = 'D:\FACED\Rawdata\T2\M99'; 
% SaveDirectory = 'D:\FACED\T2_matlab\';

function AddRecording(RawDataFolder, SaveDirectory)

%% set up
if ~strcmp(RawDataFolder(end), filesep)
    RawDataFolder = [RawDataFolder filesep];
end

if ~strcmp(SaveDirectory(end), filesep)
    SaveDirectory = [SaveDirectory filesep];
end

%% Get information
Acquisitions = dir([RawDataFolder, 'A*']);
varNames = {'Mouse', 'Acq', 'Rec', 'Vessel', 'Type', 'Group', 'Depth', 'Comments', 'RawDataFolder', 'DataFolder'};
sz = [size(Acquisitions,1) size(varNames,2)];
varTypes = repmat({'cell'}, 1, size(varNames,2));
Mousetable = table('size', sz, 'VariableTypes', varTypes, 'VariableNames',varNames);

% Mouse
seps = strfind(RawDataFolder, filesep);
Mouse = RawDataFolder(seps(end-1)+1:seps(end)-1);
Mousetable.Mouse = repmat({Mouse}, sz(1), 1);

%% check if already done
if exist([SaveDirectory 'Overview.mat'], 'file')
    load ([SaveDirectory 'Overview.mat'], 'Overview')
    if contains(Overview.Mouse, Mouse)
        disp('Mouse already added to overview, function exited.')
        return
    end
end

%% Add per acquisition
for ind = 1:size(Acquisitions,1)

    % Acq
    Acq = Acquisitions(ind).name;
    Mousetable.Acq(ind) = {Acq};

    %% temp
    infofile_patch([RawDataFolder Acq]);
    
    % RawDataFolder
    Mousetable.RawDataFolder(ind) = {[RawDataFolder Acq]};

    % DataFolder
    SaveFolder = [SaveDirectory Mouse filesep Acq filesep];
    if ~exist(SaveFolder, 'dir')
        mkdir(SaveFolder)
    end
    Mousetable.DataFolder(ind) = {SaveFolder};

    % get a column per relevant parameter in the info file
    infofile = readlines([RawDataFolder Acq filesep 'info.txt']);

    % Note: only takes first line of Comments. If you hit enter in the
    % comment section, it will not show in overview. 
    info_parameters = {'Vessel', 'Type', 'Group', 'Depth', 'Comments'};
    for inddropboxes = 1:size(info_parameters,2)
        line= infofile(contains(infofile, [info_parameters{inddropboxes} ': ']));
        line = line(1,:); % get first occurence (in case you put "type" in comments)
        line = convertStringsToChars(line);
        line = line(regexp(line, ':', 'once')+2:end); % get only after :
        Mousetable.(info_parameters{inddropboxes})(ind) = {line};
    end

end

%% Attach to overview table if relevant
if exist([SaveDirectory 'Overview.mat'], 'file')
    load ([SaveDirectory 'Overview.mat'], 'Overview', 'Pipeline_check')
    Overview = [Overview; Mousetable];

    warning('off')
    startindex = size(Pipeline_check, 1)+1;
    endindex = startindex+size(Acquisitions,1)-1;
    Pipeline_check.Mouse(startindex:endindex) = Mousetable.Mouse;
    Pipeline_check.Acq(startindex:endindex) = Mousetable.Acq;
    Pipeline_check.AddRecording(startindex:endindex) = repmat(datetime, size(Mousetable,1),1);
    % Pipeline_check = [Pipeline_check; Pipeline_check_mouse];
    warning('on')
else
    % make pipeline check table
    Pipeline_check_mouse = Mousetable(:,1:2);
    Pipeline_check_mouse.AddRecording = repmat(datetime, size(Mousetable,1),1);
    Overview = Mousetable;
    Pipeline_check = Pipeline_check_mouse;
end

%% Save
save([SaveDirectory 'Overview.mat'], 'Overview', 'Pipeline_check');

end