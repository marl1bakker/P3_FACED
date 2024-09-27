function ReadRawFACEDFiles(DataFolder, SaveFolder, BinningSpatial, BinningTemp, b_SubROI)
% function ImagesClassification(DataFolder, SaveFolder, BinningSpatial, BinningTemp, b_SubROI)
%%% 
% Adapted from ImagesClassification to work for FACED. Adapted by Marleen
% 16-2-2024. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channels classification for Labeotech IOS systems.
% (see manual page 26 for more details)
%
% At the end, 1 .dat file and 1 .mat file will be generated.
% The first file (.dat) contains the raw images in chronological order.
% The second file (.mat) contains all the informations about the
% acquisition (Freq., Stimulation vector, ROI, etc.)
%
%%% Input Parameters:
% 1- DataFolder Path:
%  Path contaning dataset from Labeo's system
% 2- SaveFolder Path:
%  Path where to save
% 3- Spatial Binning:
%  Set to 1 for no binning; 2 for a 2x2 binning; 4 for a 4x4 binning and so on...
% 4- Temporal Binning:
%  Set to 1 for no binning; Otherwise, enter de number of frames to be combined
%  for exemple: 4 will merge the images by group of 4. So, images #1,2,3,4
%  will become image #1 after avering them together. Images #5, 6, 7, 8 will
%  become frame #2, etc.
% 5- Region of Interest (ROI)
%  this parameter is a boolean (0 or 1) to tell the software if
%  we want to keep the whole image or if we want to select a smaller ROI


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set-up
if(nargin < 5)
    b_SubROI = 0;
end

AcqInfoStream = ReadInfoFileFACED(DataFolder);

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = strcat(DataFolder, filesep);
end
if( ~strcmp(SaveFolder(end), filesep) )
    SaveFolder = strcat(SaveFolder, filesep);
end

% Create save folder if it does not exist.
if ~exist(SaveFolder,'dir')
    mkdir(SaveFolder)
end
save([SaveFolder 'AcqInfos.mat'],'AcqInfoStream'); %To keep this information and avoid multiple reading of the txt file.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Format and Header Information:

hWima = 5; %wat is dit? MB
nx = AcqInfoStream.nx;
ny = AcqInfoStream.ny;
ny = ny + AcqInfoStream.ny_extra;
% extra_ny = AcqInfoStream.n_extra__in_ny_;
% ny = ny + extra_ny;

%Header format for each individual image:
%frameFormat = {'uint64', 3, 'framej';'uint16', [double(nx), double(ny)], 'imgj'}; % wat is dit MB
frameFormat = {'uint16', [double(nx), double(ny)], 'imgj'}; % 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select ROI if applicable
LimX = [1 nx];
LimY = [1 ny];

Rx = round((LimX(2) - LimX(1) + 1)/BinningSpatial); % nieuwe pixel grootte? MB
Ry = round((LimY(2) - LimY(1) + 1)/BinningSpatial);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% How many colors and in which order?
% only one, is faced. deleted "if( AcqInfoStream.MultiCam )"
imgFilesList = dir([DataFolder 'faced_*.bin']);
% matH = ChannelsSort(imgFilesList); %matH is to make files of same lenght.
% only have one file so dont need matH.
ChannelsSort(imgFilesList);
fprintf('Done. \n');

    function fColor = ChannelsSort(fList)
        % CHANNELSORT performs the classification of the raw data into the
        % existing colors and saves the data into separate .dat files.
        hTag = 'faced.mat';
        dTag = 'faced.dat';

        if( exist([SaveFolder hTag], 'file') )
            delete([SaveFolder hTag]);
        end

        %matfile:
        %MB%
        % fColor = matfile([SaveFolder hTag], 'Writable', true);
        % fColor.datFile = dTag;
        % fColor.datSize = [Ry, Rx];
        % fColor.datLength = 0;
        % fColor.FirstDim = 'y';
        % fColor.Datatype = 'single';
        % fColor.datName = 'data';
        % fColor.dim_names = {'Y', 'X', 'T'};
        % fColor.Freq = (AcqInfoStream.FrameRateHz)/BinningTemp; %MB%

        % fColor.tExposure = colors(1).Exposure;
        fid = fopen([SaveFolder dTag],'w');

        % Opening Images Files:
        for indF = 1:size(fList,1)
            fprintf('Sorting %s.', fList(indF).name);
            data = memmapfile([DataFolder fList(indF).name],...
                'Offset', hWima*4, 'Format', frameFormat,...
                'repeat', inf);
            data = data.Data;
            % hData = reshape([data.framej], 3, []);
            iData = reshape([data.imgj], nx, ny, []);
            Images = permute(iData,[2 1 3]);
            clear data;

            if isempty(Images)
                % In cases where the temporal binning causes the overflow
                % of all frames. This should happen only in the last ".bin"
                % file.
                break
            end

            if( any(sum(sum(Images,1),2) == 0) )
                idx = find(sum(sum(Images,1),2) > 1);
                Images = interp1(idx, single(reshape(Images(:,:,idx),[], length(idx)))', 1:size(Images,3),'linear','extrap');
                Images = reshape(Images', ny, nx, []);
            end

            %SubROI
            if( b_SubROI )
                Images = Images(round(LimY(1)):round(LimY(2)),round(LimX(1)):round(LimX(2)),:);
            end
            %Temporal Binning
            if( BinningTemp > 1 )
                Images = imresize3(Images, [size(Images,1), size(Images,2),...
                    size(Images,3)/BinningTemp], 'linear');
            end
            %Spatial Binning
            if( BinningSpatial > 1 )
                Images = imresize(Images,1/BinningSpatial);
            end

            % Remove returning pmt (MB)
            % Images = 

            fwrite(fid, single(Images), 'single');
            % fColor.datLength = fColor.datLength + size(Images,3); %MB%
            fprintf('\n');
        end

        fclose(fid);
    end
end


function out = ReadInfoFileFACED(FolderPath, varargin)
% This function parses the "info.txt" file and saves the data to a
% structure.
% Inputs:
% FolderPath (char): Path to folder containing the "info.txt" file.
% infoFile (char): Optional Parameter. Name of the .TXT file containing the
% acquisition information. Use this parameter to read a file with a
% different name as "info".

infofile_patch(FolderPath);

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

    if contains(Param, '(') || contains(Param, ')') % for ny extra, redundant now
        Param = strrep(Param, '(', '_');
        Param = strrep(Param, ')', '_');
    end
    % Save parameter value to structure:
    out.(Param) = Value;
end
% Save MultiCam index:
out.MultiCam = 0;
end