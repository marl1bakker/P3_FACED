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
% varNames = {'Mouse', 'Acq', 'Rec', 'Vessel', 'Type', 'Group', 'Depth', 'Comments', 'RawDataFolder', 'DataFolder'};
varNames = {'Mouse', 'Acq', 'Comments', 'Vessel', 'Group', 'Depth', 'nx', 'ny', 'ny_extra', 'height', 'width', 'Acq Freq', 'FrameRateHz', 'RawDataFolder', 'DataFolder'};
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
    if any(contains(Overview.Mouse, Mouse))
        compare_table = Overview(contains(Overview.Mouse, Mouse),:);
        % disp('Mouse already added to overview, function exited.')
        % return
    end
end

%% Add per acquisition
for ind = 1:size(Acquisitions,1)

    Acq = Acquisitions(ind).name;
    
    % check if done
    if exist("compare_table",'var') && any(matches(compare_table.Acq, Acq))
       Mousetable.Acq(ind) = {'dummy'};
       % Mousetable(ind,:) = compare_table(matches(compare_table.Acq, Acq),:);
       continue
    end

    % Acq
    Mousetable.Acq(ind) = {Acq};

    %% temp
    % infofile_patch([RawDataFolder Acq]);
    
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
    % info_parameters = {'Vessel', 'Type', 'Group', 'Depth', 'Comments'};
    info_parameters = {'Comments', 'Vessel', 'Group', 'Depth', 'nx', 'ny', 'ny_extra', 'height', 'width', 'Acq Freq', 'FrameRateHz'};
    for inddropboxes = 1:size(info_parameters,2)
        line= infofile(contains(infofile, [info_parameters{inddropboxes} ': ']));
        line = line(1,:); % get first occurence (in case you put "type" in comments)
        line = convertStringsToChars(line);
        line = line(regexp(line, ':', 'once')+2:end); % get only after :
        Mousetable.(info_parameters{inddropboxes})(ind) = {line};
    end

end

Mousetable(matches(Mousetable.Acq, 'dummy'),:) = [];

%% Attach to overview table if relevant
if size(Mousetable, 1) == 0
    disp('Mouse already added to overview, function exited.')
    return
elseif exist([SaveDirectory 'Overview.mat'], 'file')
    load ([SaveDirectory 'Overview.mat'], 'Overview', 'Pipeline_check')
    Overview = [Overview; Mousetable];

    warning('off')
    startindex = size(Pipeline_check, 1)+1;
    % endindex = startindex+size(Acquisitions,1)-1;
    endindex = startindex+size(Mousetable,1)-1;
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



function infofile_patch(RawDataFolder)

if ~strcmp(RawDataFolder(end), filesep)
    RawDataFolder = [RawDataFolder filesep];
end

infofile = readlines([RawDataFolder 'info.txt']);
writelines(infofile, [RawDataFolder 'info_copy.txt']);

% If comments are on other line than Comments: (put first line next to
% Comments: , put other lines with an enter in between below that
indline = find(matches(infofile, 'Comments: '));
if indline < size(infofile, 1)
    infofile(indline,:) = [];
    infofile(indline,:) = append("Comments: ", infofile(indline,:));
end

% Pos_to_brain to Depth
if sum(contains(infofile, 'Pos_to_brain: '))
    indline = find(contains(infofile, 'Pos_to_brain: '));
    line= infofile(indline);
    line = convertStringsToChars(line);
    line = ['Depth: ' line(15:end)];
    infofile(indline) = line;
end

% Get rid of () in ny extra
if sum(contains(infofile, 'n_extra (in ny): '))
    indline = find(contains(infofile, 'n_extra (in ny): '));
    line= infofile(indline);
    line = convertStringsToChars(line);
    line = ['ny_extra: ' line(18:end)];
    infofile(indline) = line;
end

% Get rid of s after height
indline = find(contains(infofile, 'height: '));
line = infofile(indline);
if contains(line, 's')
    line = convertStringsToChars(line);
    line = line(1:end-1);
    infofile(indline) = line;
end

% If no Vessel/Group/Type:
info_parameters = {'Vessel', 'Acquisition_Type', 'Group', 'Depth', 'width'};
for inddropboxes = 1:size(info_parameters,2)
    if ~sum(find(contains(infofile, [info_parameters{inddropboxes} ': '])))
        infofile(end+1) = [info_parameters{inddropboxes} ': unknown'];
    end
end

% Just for mouse M01:
seps = strfind(RawDataFolder, filesep);
Mouse = RawDataFolder(seps(end-2)+1:seps(end-1)-1);
if matches(Mouse, 'M01') && ~sum(find(contains(infofile, 'FrameRateHz_old: ')))
    infofile(contains(infofile, 'width:')) = 'width: 60';
    
    % the FrameRateHz is half the value!
    indline = find(contains(infofile, 'FrameRateHz: '));
    line = infofile(indline);
    line = convertStringsToChars(line);
    hz_old = str2double(line(14:end));
    hz_new = hz_old/2;
    infofile(indline) = ['FrameRateHz: ' num2str(hz_new)];
    infofile(end+1) = ['FrameRateHz_old: ' num2str(hz_old)];

end



%% save
writelines(infofile, [RawDataFolder 'info.txt']);

end
