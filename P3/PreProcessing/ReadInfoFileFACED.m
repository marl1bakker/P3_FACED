function out = ReadInfoFileFACED(FolderPath, varargin)
% This function parses the "info.txt" file and saves the data to a
% structure.
% Inputs:
% FolderPath (char): Path to folder containing the "info.txt" file.
% infoFile (char): Optional Parameter. Name of the .TXT file containing the
% acquisition information. Use this parameter to read a file with a
% different name as "info".

% Read the info.txt file:
if nargin == 1
    txt = readcell(fullfile(FolderPath, 'info.txt'), 'Delimiter', ':', 'NumHeaderLines',1);
else
    [~,infoFile,ext] = fileparts(varargin{:});
    if isempty(ext)
        ext = '.txt';
    end        
    txt = readcell(fullfile(FolderPath, [infoFile, ext]), 'Delimiter', ':', 'NumHeaderLines',1);
end
    
% Rebuild strings that were split by the delimiter:
b_hasMissingVals = cellfun(@(x) isa(x,'missing'), txt);
if size(txt,2)>2
    for ind = 1:size(txt,1)
        if all(~b_hasMissingVals(ind,2:end))            
            txt{ind,2} = strjoin(cellfun(@num2str,txt(ind,2:end), 'UniformOutput',false),': ');
        end
    end
end

% Remove Parameters with missing values:
txt(b_hasMissingVals(:,2),:) = [];
% Replace white spaces in Parameters column by underscores:
txt(:,1) = cellfun(@(x) strrep(x, ' ', '_'),txt(:,1), 'UniformOutput', false);
% Create structure from text:
out = struct();
for ind = 1:size(txt,1)
    Param = txt{ind,1};
    Value = txt{ind,2};

    if ischar(Value)
        if ( all(isstrprop(erase(Value, ' '), 'digit')) )
            Value = str2num(Value);%#ok. "str2double" would return NaN in this case.
        end
    end

    if contains(Param, '(') || contains(Param, ')')
        Param = strrep(Param, '(', '_');
        Param = strrep(Param, ')', '_');
    end
    % Save parameter value to structure:
    out.(Param) = Value;
end
% Save MultiCam index:
out.MultiCam = 0;
end