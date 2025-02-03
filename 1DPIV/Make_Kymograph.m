% Code adjusted from getKymo
% MB: for original code, see bottom (commented)

% ROIname for if you want a specific ROI (or ROIs), otherwise will
% calculate all roi that were found
% Give like: {'kymoROI_01'} or {'ROI_1'; 'ROI_2'} (or, theoretically
% {'roi_1', 'roi_2}) or 'auto_list' (for all roi)
% plotoverview can be 1 (give overview) or 0 (dont)
% overwrite can be 1 (ignore existing things and overwrite them) or 0 (dont)

% example:
% DataFolder = 'C:\Users\marle\OneDrive\Documenten\MATLAB\P3\FACED2PFM-vessel-main\FACED2PFM-vessel-main\1DPIV\sample data\';
% own data:
% DataFolder = 'D:\FACED\T2_matlab\2024-07-18-10';

function Make_Kymograph(DataFolder, PlotOverview, ROIname, overwrite)

%% set up
if ~strcmp(DataFolder(end), filesep)
    DataFolder = [DataFolder filesep];
end

% get list of ROI made
if exist('ROIname', 'var') && ~matches(ROIname, 'auto_list')
    ROI_list = ROIname;
else
    ROI_list = dir([DataFolder 'kymoROI*.mat']);
    ROI_list = struct2cell(ROI_list);
    ROI_list = ROI_list(1,:);
end

if isempty(ROI_list)
    error('No ROI found, function exited.')
end

if ~exist("PlotOverview", 'var')
    PlotOverview = 0;
end

if ~exist("overwrite", 'var')
    overwrite = 0;
end

%% load data
load([DataFolder 'AcqInfos.mat'], 'AcqInfoStream');

if AcqInfoStream.height == 1
    disp('Plot_kymograph already done within Clean_Data because it is a linescan')
    
    load([DataFolder 'kymoROI_1.mat'], 'kymoImg');
    f2 = figure;
    lengthkymo = round(AcqInfoStream.FrameRateHz/2);
    imagesc(kymoImg(1:lengthkymo,:)')

    % get right x axis
    % Acq_sec = lengthkymo/AcqInfoStream.FrameRateHz;
    % xticks(0:AcqInfoStream.FrameRateHz:lengthkymo);
    % xticklabels(0:Acq_sec);
    % xlabel('Seconds');
    Acq_msec = lengthkymo/AcqInfoStream.FrameRateHz*1000;
    xticks(0:AcqInfoStream.FrameRateHz/10:lengthkymo); % this and next line needs to be 1000
    xticklabels(0:100:Acq_msec);
    xlabel('Milliseconds');

    colormap('gray');
    f2.Position = [20 40 f2.Position(3)*2 f2.Position(4)];
    savefig([DataFolder 'Kymograph_overview.fig']);
    pause(5)
    close(f2)

    return
end

datSize = [AcqInfoStream.ny+AcqInfoStream.ny_extra AcqInfoStream.nx];

% % with orignal data (incl mirror)
% load([DataFolder 'morphology_img.mat'], 'av_image_og', 'proportion_XY_ax');
% av_image = av_image_og;

% with cropped data averaged over y axis
load([DataFolder 'morphology_img.mat'], 'av_image');
load([DataFolder 'AcqInfos.mat'], 'AcqInfoStream');
proportion_XY_ax = [1 AcqInfoStream.pxl_sz_x/AcqInfoStream.pxl_sz_y_avged 1];

fid = fopen([DataFolder 'faced.dat']);
dat = fread(fid, inf, '*single');
fclose(fid);
dat = reshape(dat, size(av_image,1), size(av_image,2), []);

% %% Filter -- temp CHECK
% for ind = 1:size(dat,3)
%     dat(:,:,ind) = imbilatfilt(mat2gray(dat(:,:,ind)),0.1,2);
% end

%% Go per roi
for ind_ROI = 1:length(ROI_list)
    currentROI = ROI_list{ind_ROI};
    warning('off')
    load([DataFolder currentROI], 'kymoImg', 'ROI_type');
    warning('on')

    if exist('kymoImg', 'var') && overwrite == 0
        disp('ROI kymoimg already made, skipped.')
        clear Coor kymoImg ROI_type 
        continue
    end

    % Get diameter, pixelsize, theta. Also saves this in currentROI
    % [theta, pxlSize] = Vessel_Diameter_Pixelsize(currentROI, av_image, pxlX,pxlY,AcqInfoStream.Vessel,1,0,[1 2]);
    % [theta, pxlSize] = Vessel_Diameter_Pixelsize(currentROI, av_image, AcqInfoStream.pxl_sz_x,AcqInfoStream.pxl_sz_y,AcqInfoStream.Vessel,1,0,[1 2]);    

    if ~exist('ROI_type', 'var') %means it's done by fiji
        ROI_type = 'Fiji';
    end

    switch ROI_type
        % case 'linescan'
        %     shapedat = size(dat);
        %     dat = permute(dat, [1, 3, 2]);
        %     dat = reshape(dat, [], shapedat(2));
        %     kymoImg = dat;
        %     save([DataFolder currentROI], 'shapedat', '-append')


        case {'line', 'perpendicular_line', 'linescan'}
            load([DataFolder currentROI], 'Coor', 'roi_pixels');

            % kymoImg = zeros(size(dat,3), size(Coor,1));
            kymoImg = zeros(size(dat,3), length(roi_pixels));
            for frmIndx = 1:size(dat,3)
                current_frame = dat(:,:,frmIndx);
                kymoImg(frmIndx,:) = current_frame(roi_pixels);
                % for pxlIndx = 1:size(Coor,1)
                %     kymoImg(frmIndx, pxlIndx) = current_frame(Coor(pxlIndx,2), Coor(pxlIndx,1));
                % end
            end
            
            clear current_frame frmIndx pxlIndx 

            f1 = figure;
            lengthkymo = round(AcqInfoStream.FrameRateHz/2);
            imagesc(kymoImg(1:lengthkymo,:)')
            colormap('gray');

            % get right x axis
            Acq_sec = lengthkymo/AcqInfoStream.FrameRateHz;
            xticks(0:AcqInfoStream.FrameRateHz:lengthkymo);
            xticklabels(0:Acq_sec);
            xlabel('Seconds');

            close(f1)

        case 'automatic'
            load([DataFolder currentROI], 'mask');

            mask(mask == 0) = NaN;
            kymoImg = dat.*mask;
            kymoImg = mean(kymoImg, 1, 'omitnan'); % mean over y axis
            kymoImg = reshape(kymoImg, datSize(2), []);
            kymoImg = kymoImg(~isnan(kymoImg(:,1)),:);
            kymoImg = kymoImg'; % to make sure all roi types are the same
            f1 = figure;
            imagesc(kymoImg)
            colormap('gray')
            close(f1)

        case {'block', 'block_fixed_height'}
            load([DataFolder currentROI], 'roiY', 'roiX', 'block_height_um');

            kymoImg = mean(dat(roiY(1):roiY(2), roiX(1):roiX(2),:), 1);
            kymoImg = reshape(kymoImg, [], size(dat, 3))';

            f1 = figure;
            imagesc(kymoImg(1:300,:)');
            colormap('gray');
            title('block roi')
            subtitle([num2str(roiY(2)-roiY(1)) ' points in y, over ' num2str(block_height_um) ' um'])
            close(f1)
            
        case 'perpendicular_block'
            load([DataFolder currentROI], 'roiY', 'roiX');

            kymoImg = mean(dat(roiY(1):roiY(2), roiX(1):roiX(2),:), 2);
            kymoImg = reshape(kymoImg, [], size(dat, 3))';

            f1 = figure;
            imagesc(kymoImg(1:300,:)');
            colormap('gray');
            title('perpendicular block roi')
            % subtitle([num2str(roiY(2)-roiY(1)) ' points in y, over ' num2str(block_height_um) ' um'])
            close(f1)

        case 'Fiji'
            % Open ImageJ/Fiji
            % Load the sum.tiff or sum.png file
            % Hold the line button down to make a line ROI by hand (freestyle)
            % Draw the ROI along the vessel
            % File - Save As - Selection - name.roi

            % ROIlist = dir(fullfile(ROIpath, [datName,'*.zip']));
            % if isempty (ROIlist)
            %     ROIlist = dir(fullfile(ROIpath, '*.roi'));
            % end
            % 
            % if isempty(ROIlist) %MB: if your roi list is empty
            %     disp('no roi found')
            %     return
            %     % write code to get ROI
            % end
            % 
            % ROIname = ROIlist(1).name;
            % sROI = ReadImageJROI( fullfile(ROIpath,ROIname)); % load imageJ ROIs
            % if strcmp(ROIname(end-3:end),'.roi')
            %     tmpROI = sROI;
            %     sROI = cell(1,1);
            %     sROI{1} = tmpROI;
            % end
            % 
            % % MB: sROI are structures with information about the ROI, but
            % % only need Coor.
            % for ROIindx = 1:length(sROI)
            %     Coor = sROI{ROIindx}.mnCoordinates;
            %     Coor = Coor( Coor(:,1)>0 & Coor(:,2)>0, :);
            % 
            %     kymoImg = zeros(size(dat,3), size(Coor,1));
            %     for frmIndx = 1:1:size(dat,3)
            %         current_frame = dat(:,:,frmIndx);
            %         for pxlIndx = 1:1:size(Coor,1)
            %             kymoImg(frmIndx, pxlIndx) = current_frame(Coor(pxlIndx,2), Coor(pxlIndx,1));
            %         end
            %     end
            % 
            %     % MB: save, take name of mouse/acq, add _kymoROI and then
            %     % the roi index, so that each roi has a seperate file
            %     % example name: '20200302_47_D430_kymoROI01.tif'

            % end
    end

    % Save kymograph
    % kymoName = sprintf('%s_kymoROI%02d.tif',datName,ROIindx);
    % imwrite(kymoImg,fullfile(ROIpath, kymoName));
    save([DataFolder currentROI], 'kymoImg', '-append')

    disp(['Kymograph made for ' currentROI])
    clear Coor kymoImg ROI_type roiY roiX

end

% get angle/speed
% see Linescan_Velocimetry()

if PlotOverview
    %% All together
    % have to redo the ROI_list, because if you specified one roi your list
    % will be only one
    seps = strfind(DataFolder, filesep);
    ROI_list = dir([DataFolder 'kymoROI*.mat']); 
    ROI_list = struct2cell(ROI_list);
    ROI_list = ROI_list(1,:);

    f2 = figure;
    t = tiledlayout(3,size(ROI_list,2));
    nexttile([1,size(ROI_list, 2)])
    imagesc(av_image)
    ax = gca;
    ax.DataAspectRatio = proportion_XY_ax;
    colormap('gray')
    title([DataFolder(seps(end-2)+1:seps(end-1)-1) ' ' DataFolder(seps(end-1)+1:seps(end)-1)])

    for ind_ROI = 1:length(ROI_list)
        warning('off')
        % load([DataFolder 'kymoROI_' num2str(ind_ROI) '.mat'], 'kymoImg', 'ROI_type', 'Coor', 'roi_pixels')
        load([DataFolder ROI_list{ind_ROI}], 'kymoImg', 'ROI_type', 'Coor', 'roi_pixels', 'roiY', 'block_height_um')
        warning('on')

        % show ROI
        nexttile(size(ROI_list,2)+ind_ROI)
        ROI_image = av_image;
        if contains(ROI_type, 'block')
            ROI_image(Coor(:,2), Coor(:,1)) = 0;
        elseif contains(ROI_type, 'line') || matches(ROI_type, 'automatic')
            ROI_image(roi_pixels) = 0;
        end
        imagesc(ROI_image)

        % show kymograph of roi
        nexttile(size(ROI_list,2)*2+ind_ROI)
        % imagesc(kymoImg(1:round(AcqInfoStream.FrameRateHz/2),:)')
        lengthkymo = round(AcqInfoStream.FrameRateHz/2);
        imagesc(kymoImg(1:lengthkymo,:)')

        % get right x axis
        % Acq_sec = lengthkymo/AcqInfoStream.FrameRateHz;
        % xticks(0:AcqInfoStream.FrameRateHz:lengthkymo);
        % xticklabels(0:Acq_sec);
        % xlabel('Seconds');
        Acq_msec = lengthkymo/AcqInfoStream.FrameRateHz*1000;
        xticks(0:AcqInfoStream.FrameRateHz/10:size(kymoImg,1)); % this plus next line needs to be 1000
        xticklabels(0:100:Acq_msec);
        xlabel('Milliseconds');

        colormap('gray');
        title(ROI_type, 'Interpreter','none')


        if contains(ROI_type, 'block')
            % load([DataFolder 'kymoROI_' num2str(ind_ROI) '.mat'], 'roiY')
            % load([DataFolder ROI_list{ind_ROI}], 'roiY')
            % subtitle([num2str(roiY(2)-roiY(1)) ' points in y'])
            subtitle([num2str(roiY(2)-roiY(1)) ' points in y, over ' num2str(block_height_um) ' um'], 'Interpreter','none')

        end
    end

    % saveas(f2, [DataFolder 'Kymograph_overview.png'], 'png');
    savefig([DataFolder 'Kymograph_overview.fig']);
    pause(5)
    close(f2)
    % waitfor(f2)
end

end



% function [ind, label] = drawline(p1,p2,image_size)
% %DRAWLINE Returns the geometric space (matrix indices) occupied by a line segment
% %in a MxN matrix.  Each line segment is defined by two endpoints.
% %
% %   IND = DRAWLINE(P1, P2, IMAGE_SIZE) returns the matrix indices
% %   of the line segment with endpoints p1 and p2.
% %   If both points are out of the image boundary no line is drawn and an error will appear.
% %   If only one of the endpoints is out of the image boundary a line is still drawn.
% %
% %       ARGUMENT DESCRIPTION:
% %                       P1 - set of endpoints (Nx2). ([row column; ...])
% %                       P2 - set of endpoints that connect to p1 (Nx2). ([row column; ...])
% %               IMAGE_SIZE - vector containing image matrix dimensions,
% %                            where the first element is the number of rows
% %                            and the second element is the number of
% %                            columns.
% %
% %       OUTPUT DESCRIPTION:
% %                      IND - matrix indices occupied by the line segments.
% %                    LABEL - label tag of each line drawn (from 1 to N).
% %
% %
% %       1
% %    1 _|_ _ _ _> COLUMNS
% %       |_|_|_|_
% %       |_|_|_|_
% %       |_|_|_|_
% %       V
% %      ROWS
% %
% %   Example
% %   -------------
% %   BW = zeros(250,250);
% %   p1 = [10 10; 23 100; -14 -40];
% %   p2 = [50 50; 90 100;  50  50];
% %   [ind label] = drawline(p1,p2,[250 250]); % OR ...drawline(p2,p1,...
% %   BW(ind) = label;
% %   figure, imshow(BW,[])
% %
% % See also line, ind2sub.
% % Credits:
% % Daniel Simoes Lopes
% % ICIST
% % Instituto Superior Tecnico - Universidade Tecnica de Lisboa
% % danlopes (at) civil ist utl pt
% % http://www.civil.ist.utl.pt/~danlopes
% %
% % June 2007 original version.
% % Input verification.
% if max(size(p1) ~= size(p2))
%     error('The number of points in p1 and p2 must be the same.')
% end
% if length(size(image_size)) ~= 2
%     error('Image size must be bi-dimensional.')
% end
% % Cicle for each pair of endpoints.
% ind = [];
% label = [];
% for line_number = 1:size(p1,1)
% 
%     % Point coordinates.
%     p1r = p1(line_number,1);   p1c = p1(line_number,2);
%     p2r = p2(line_number,1);   p2c = p2(line_number,2);
% 
%     % Image dimension.
%     M = image_size(1); % Number of rows.
%     N = image_size(2); % Number of columns.
% 
%     % Boundary verification.
%     % A- Both points are out of range.
%     if  ((p1r < 1 || M < p1r) || (p1c < 1 || N < p1c)) && ...
%             ((p2r < 1 || M < p2r) || (p2c < 1 || N < p2c)),
%         error(['Both points in line segment n� ', num2str(line_number),...
%             ' are out of range. New coordinates are requested to fit',...
%             ' the points in image boundaries.'])
%     end
% 
%     % Reference versors.
%     % .....r..c.....
%     eN = [-1  0]';
%     eE = [ 0  1]';
%     eS = [ 1  0]';
%     eW = [ 0 -1]';
% 
%     % B- One of the points is out of range.
%     if (p1r < 1 || M < p1r) || (p1c < 1 || N < p1c) || ...
%             (p2r < 1 || M < p2r) || (p2c < 1 || N < p2c),
%         % ....Classify the inner and outer point.
%         if     (p1r < 1 || M < p1r) || (p1c < 1 || N < p1c)
%             out = [p1r; p1c];  in = [p2r; p2c];
%         elseif (p2r < 1 || M < p2r) || (p2c < 1 || N < p2c)
%             out = [p2r; p2c];  in = [p1r; p1c];
%         end
%         % Vector defining line segment.
%         v = out - in;
%         aux = sort(abs(v)); aspect_ratio = aux(1)/aux(2);
%         % Vector orientation.
%         north = v'*eN;
%         west  = v'*eW;              east  = v'*eE;
%         south = v'*eS;
%         % Increments.
%         deltaNS = [];
%         if north > 0, deltaNS = -1; end
%         if south > 0, deltaNS =  1; end
%         if isempty(deltaNS), deltaNS = 0; end
%         deltaWE = [];
%         if east > 0, deltaWE =  1; end
%         if west > 0, deltaWE = -1; end
%         if isempty(deltaWE), deltaWE = 0; end
%         % Matrix subscripts occupied by the line segment.
%         if abs(v(1)) >= abs(v(2))
%             alpha(1) = in(1); beta(1) = in(2);
%             iter = 1;
%             while (1 <= alpha(iter)) && (alpha(iter) <= M) && ...
%                     (1 <=  beta(iter)) && (beta(iter)  <= N),
%                 alpha(iter+1) = alpha(iter) + deltaNS;              % alpha grows throughout the column direction.
%                 beta(iter+1)  = beta(iter)  + aspect_ratio*deltaWE; % beta grows throughout the row direction.
%                 iter = iter + 1;
%             end
%             alpha = round(alpha(1:end-1)); beta = round(beta(1:end-1));
%             ind = cat(2,ind,sub2ind(image_size,alpha,beta));
%             label = cat(2,label,line_number*ones(1,max(size(alpha))));
%         end
%         % ...
%         if abs(v(1)) < abs(v(2))
%             alpha(1) = in(2); beta(1) = in(1);
%             iter = 1;
%             while (1 <= alpha(iter)) && (alpha(iter) <= N) &&...
%                     (1 <=  beta(iter)) && (beta(iter)  <= M),
%                 alpha(iter+1) = alpha(iter) + deltaWE;              % alpha grows throughout the row direction.
%                 beta(iter+1)  = beta(iter)  + aspect_ratio*deltaNS; % beta grows throughout the column direction.
%                 iter = iter + 1;
%             end
%             alpha = round(alpha(1:end-1)); beta = round(beta(1:end-1));
%             ind = cat(2,ind,sub2ind(image_size,beta,alpha));
%             label = cat(2,label,line_number*ones(1,max(size(alpha))));
%         end
%         clear alpha beta
%         continue
%     end
%     % C- Both points are in range.
%     in = [p1r; p1c];  out = [p2r; p2c]; % OR in = p2; out = p1;
%     % Vector defining line segment.
%     v = out - in;
%     aux = sort(abs(v)); aspect_ratio = aux(1)/aux(2);
%     % Vector orientation.
%     north = v'*eN;
%     west  = v'*eW;              east  = v'*eE;
%     south = v'*eS;
%     % Increments.
%     deltaNS = [];
%     if north > 0, deltaNS = -1; end
%     if south > 0, deltaNS =  1; end
%     if isempty(deltaNS), deltaNS = 0; end
%     deltaWE = [];
%     if east > 0, deltaWE =  1; end
%     if west > 0, deltaWE = -1; end
%     if isempty(deltaWE), deltaWE = 0; end
%     % Matrix subscripts occupied by the line segment.
%     row_range = sort([p1r p2r]);
%     col_range = sort([p1c p2c]);
%     if abs(v(1)) >= abs(v(2))
%         alpha(1) = in(1); beta(1) = in(2);
%         iter = 1;
%         while (row_range(1) <= alpha(iter)) && (alpha(iter) <= row_range(2)) && ...
%                 (col_range(1) <=  beta(iter)) && (beta(iter)  <= col_range(2)),
%             alpha(iter+1) = alpha(iter) + deltaNS;              % alpha grows throughout the column direction.
%             beta(iter+1)  = beta(iter)  + aspect_ratio*deltaWE; % beta grows throughout the row direction.
%             iter = iter + 1;
%         end
%         alpha = round(alpha(1:end-1)); beta = round(beta(1:end-1));
%         ind = cat(2,ind,sub2ind(image_size,alpha,beta));
%         label = cat(2,label,line_number*ones(1,max(size(alpha))));
%     end
%     % ...
%     if abs(v(1)) < abs(v(2))
%         alpha(1) = in(2); beta(1) = in(1);
%         iter = 1;
%         while (col_range(1) <= alpha(iter)) && (alpha(iter) <= col_range(2)) &&...
%                 (row_range(1) <=  beta(iter)) && (beta(iter)  <= row_range(2)),
%             alpha(iter+1) = alpha(iter) + deltaWE;              % alpha grows throughout the row direction.
%             beta(iter+1)  = beta(iter)  + aspect_ratio*deltaNS; % beta grows throughout the column direction.
%             iter = iter + 1;
%         end
%         alpha = round(alpha(1:end-1)); beta = round(beta(1:end-1));
%         ind = cat(2,ind,sub2ind(image_size,beta,alpha));
%         label = cat(2,label,line_number*ones(1,max(size(alpha))));
%     end
%     clear alpha beta
%     continue
% end
% 
% end
%% commented "raw" version:

% % % MB: get all tiff files in folder
% tifStr = '*.tif'; % tif file name
% tifList = dir(fullfile(DataFolder, tifStr));

% %% MB: go per tif file
% for tifIndx = 1:1:numel(tifList)
%     tifName = tifList(tifIndx).name;
%     stackFullName = fullfile(DataFolder,tifName); % MB: same as [tifPath, filesep, tifName]
%     imgInfo = imfinfo(stackFullName);
%
%     sliceNum = length( imgInfo ); % number of frames in the stacks -- MB: very big structure, with "information" about each image in the file. Just need this to get the number of frames
%     %MB: check if it's more than one frame, if not, skip to next .tif file
%
%     if sliceNum > 1
%
%         %% MB: get ROI
%         % MB: regexp looks for certain expressions. In this case, it's the
%         % same as using strfind(tifName, 'D'). Takes name of acq.
%         ROIstr = tifName(1:regexp(tifName,'D')+3); % the char string used to
%         % identify the matched ROI files. Change this line if the ROI name is
%         % different
%         if isempty(ROIstr)
%             ROIstr = tifName(1:end-4);
%         end
%         ROIlist = dir(fullfile(ROIpath, [ROIstr,'*.zip'])); %MB: roi is zipfile...?
%         if isempty (ROIlist)
%             ROIlist = dir(fullfile(ROIpath, [ROIstr,'*.roi'])); %MB: roi is .roifile????
%         end
%
%         if ~isempty(ROIlist) %MB: if your roi list is not empty
%             ROIname = ROIlist(1).name;
%             sROI = ReadImageJROI( fullfile(ROIpath,ROIname)); % load imageJ ROIs
%             %MB: sROI is a cell-array with structures with a lot of info
%             if strcmp(ROIname(end-3:end),'.roi') % MB: if you have .roi file
%                 tmpROI = sROI;
%                 sROI = cell(1,1);
%                 sROI{1} = tmpROI;
%             end
%
%             %% MB: read images
%             img1 = imread(stackFullName,1); %MB: get first image
%             imgHeight = size(img1,1); %MB: get height & width
%             imgWidth = size(img1,2);
%
%             dat = zeros( imgHeight, imgWidth, sliceNum); %create zero file with x, y, time
%             % read in the images;
%             disp(['reading frames from ',tifName,'...']);
%             tic;
%             TiffLink = Tiff(stackFullName,'r'); %MB: creates a Tiff object for read access to the TIFF file filename.
%             h_waitBar = waitbar(0,['Loading stack ', tifName,'...']);
%             for sliceIndx=1:1:sliceNum %MB: extra 1: is unecessary, just ind = 1:sliceNum
%                 TiffLink.setDirectory(sliceIndx);
%                 dat(:,:,sliceIndx) = TiffLink.read(); %MB: read the directory you just set i guess?
%                 waitbar(sliceIndx/sliceNum,h_waitBar);
%             end
%             close(h_waitBar);
%             TiffLink.close();
%             toc;
%             % MB: now you have a 3D matrix of x y time in double.
%             % Difference with our .dat files is the values in the file
%             % (here 0-120ish, ours around 30 000) and datatype (here
%             % double, ours single)
%
%             %MB: saving stuff in best format?
%             if max(dat(:)) < 1.1
%                 dtypeSave = 'double';
%             elseif max(dat(:)) > 256
%                 dtypeSave = 'uint16';
%             else
%                 dtypeSave = 'uint8';
%             end
%
%             % save morphology image - temporally averaged/summed
%             av_image = mean(dat, 3);
%             if strcmp(dtypeSave,'uint8')
%                 if max(dat(:))<15 % MB: if theres not a lot of "activation"/signal...
%                     av_image = sum(dat,3); %MB: make it the sum instead of the mean
%                 end
%             end
%             %MB: normalisation! of just the morphImg/average image
%             av_image = (av_image - min(av_image(:)))./(max(av_image(:)) - min(av_image(:))); % normalize the morphology image;
%             imwrite(av_image,fullfile(DataFolder,[tifName(1:end-4),'_sum.tif']))
%
%             % MB: sROI are structures with information about the ROI, but
%             % only need Coor.
%             for ROIindx = 1:length(sROI)
%                 % crrtROI = sROI{ROIindx}; %old
%                 % Coor = crrtROI.mnCoordinates; %old
%                 Coor = sROI{ROIindx}.mnCoordinates; % MB: new
%                 Coor = Coor( Coor(:,1)>0 & Coor(:,2)>0, :); %MB: make sure nothing falls out of frame (?)
%
%                 % % MB: Coor are X, Y in two columns. It's a line
%                 % % that follows a vessel. example:
%                 % bla = av_image;
%                 % for ind = 1:size(Coor,1)
%                 %     bla(Coor(ind,2), Coor(ind,1)) = 2;
%                 % end
%                 % figure
%                 % imagesc(bla)
%
%                 kymoImg = zeros(size(dat,3), size(Coor,1));
%                 for frmIndx = 1:1:size(dat,3)
%                     current_frame = dat(:,:,frmIndx);
%                     for pxlIndx = 1:1:size(Coor,1)
%                         kymoImg(frmIndx, pxlIndx) = current_frame(Coor(pxlIndx,2), Coor(pxlIndx,1));
%                     end
%                 end
%
%                 % MB: save, take name of mouse/acq, add _kymoROI and then
%                 % the roi index, so that each roi has a seperate file
%                 % example name: '20200302_47_D430_kymoROI01.tif'
%                 kymoName = sprintf('%s_kymoROI%02d.tif',tifName(1:end-4),ROIindx);
%                 switch dtypeSave
%                     case 'double'
%                         imwrite(kymoImg,fullfile(ROIpath, kymoName));
%                     case 'uint16'
%                         imwrite(uint16(kymoImg), fullfile(ROIpath, kymoName));
%                     case 'uint8'
%                         imwrite(uint8(kymoImg), fullfile(ROIpath, kymoName));
%                 end
%             end
%
%         else
%             disp(['no ROI found for ',tifName])
%         end
%     else
%         disp(['not a stack: ', tifName]);
%     end
%
% end
