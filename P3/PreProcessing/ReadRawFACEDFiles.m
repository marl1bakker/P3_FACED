% DataFolder = 'C:\Users\marle\OneDrive\Documenten\PhD\P3 - Ultrafast two photon\FACED dummy data\dummy5';
% SaveFolder = 'C:\Users\marle\OneDrive\Documenten\PhD\P3 - Ultrafast two photon\FACED dummy data\Dummy5Save';

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
imgFilesList = dir([DataFolder 'faced*.bin']);
nx = AcqInfoStream.nx;
ny = AcqInfoStream.ny;

%Header format for each individual image:
frameFormat = {'uint64', 3, 'framej';'uint16', [double(nx), double(ny)], 'imgj'}; % wat is dit MB
% ImRes_XY = [nx, ny]; % replaced by changing ImRes_XY(1) by nx etc.
SizeImage = nx*ny*2 + 3*8; % waarom +3*8 MB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select ROI if applicable
% SubROI...
if( b_SubROI )
    fprintf('Redefining Region Of Interest post-process: \n');
    %Dialog (there are different options to determine the new ROI):
    ButtonName = questdlg('Would you like to use a pre-defined ROI?', ...
        'ROI', ...
        'Pre-defined', 'Draw', 'Cancel', 'Draw');
    switch ButtonName  %Depending on user choice:
        case 'Pre-defined' %Used a ROI from an other acquisition:
            [filename, pathname] = uigetfile('*.mat', 'Select ROI file');
            if isequal(filename,0) || isequal(pathname,0)
                disp('User pressed cancel')
                Pos = [1 1 nx ny];
            else
                load([pathname filesep filename]);%#ok
            end
        case 'Draw' %Select ROI directly on a frame:
            dat = memmapfile([DataFolder...
                imgFilesList(1).name],...
                'Offset', hWima*4 + 5*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
            dat = dat.Data.imgj;
            fig = figure; imagesc(dat);
            h = drawrectangle();
            wait(h);
            Pos = h.Position;
            close(fig);
        case 'Cancel' %User Changed is mind and want to use the original ROI
            disp('User pressed cancel')
            Pos = [1 1 nx ny];
    end

    LimX = [round(Pos(1)) round(Pos(1)+Pos(3))];
    LimY = [round(Pos(2)) round(Pos(2)+Pos(4))];
    save([SaveFolder 'ROI.mat'],'Pos'); %Save region of interest in a .mat file
else
    LimX = [1 nx];
    LimY = [1 ny];
end
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
        fColor = matfile([SaveFolder hTag], 'Writable', true);
        fColor.datFile = dTag;
        fColor.datSize = [Ry, Rx];
        fColor.datLength = 0;
        fColor.FirstDim = 'y';
        fColor.Datatype = 'single';
        fColor.datName = 'data';
        fColor.dim_names = {'Y', 'X', 'T'};
        fColor.Freq = (AcqInfoStream.FrameRateHz)/BinningTemp;
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

            fwrite(fid, single(Images), 'single');
            fColor.datLength = fColor.datLength + size(Images,3);
            fprintf('\n');
        end

        fclose(fid);
    end
end
