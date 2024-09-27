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



