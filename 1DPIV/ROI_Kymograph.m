% Add can be 1 or 0, if it's 1 it adds a ROI if there already is a ROI of
% that type. If it's 0 it skips.
%

function ROI_Kymograph(DataFolder, Add)
% Makes ROI for kymograph

%% set up
if ~strcmp(DataFolder(end), filesep)
    DataFolder = [DataFolder filesep];
end

if ~exist('Add', 'var')
    Add = 1;
end

if ~isempty(dir([DataFolder 'kymoROI*.mat']))
    if Add == 0
        disp('ROI already made for this Acq., add was set to 0. Function exited.')
        return
    end
    nr_of_roi = size(dir('kymoROI*.mat'), 1);
else
    nr_of_roi = 0;
end

%% load data
load([DataFolder 'morphology_img.mat'], 'av_image_og', 'proportion_XY_ax');
av_image = av_image_og;
clear av_image_og;

%% Start doing ROI, keep looping until you have all the roi you want.
more_to_add = 1;

while more_to_add
    nr_of_roi = nr_of_roi+1;

    f1 = figure;
    imagesc(av_image)
    colormap('gray')
    f1.Position= [50 50 500 500];
    ax = gca;
    ax.DataAspectRatio = proportion_XY_ax;
    ROI_type_list = {'line', 'block', 'block_fixed_height', 'automatic', ...
        'perpendicular_line', 'perpendicular_block'};
    [indx, ~] = listdlg('PromptString', 'Choose a ROI type', 'ListString', ROI_type_list);
    ROI_type = ROI_type_list{indx};

    f10 = figure;
    imagesc(av_image)
    colormap('gray')

    switch ROI_type
        case {'line', 'perpendicular_line'}
            roi_pixels = drawpolyline;
            roi_pixels.Position = round(roi_pixels.Position);
            start_end_points = roi_pixels.Position;

            p1 = [];
            p2 = [];
            for indlinepoints = 1:size(roi_pixels.Position,1)-1
                p1(indlinepoints,:) = roi_pixels.Position(indlinepoints,:);
                p2(indlinepoints,:) = roi_pixels.Position(indlinepoints+1,:);
            end

            [roi_pixels, ~] = draw_line(fliplr(p1), fliplr(p2), size(av_image)); % not matlab function (drawline)
            [row, col] = ind2sub(size(av_image), roi_pixels);
            % Note: this is not perfect, sometimes the line is more than one
            % pixel wide.
            Coor = [col', row']; % make sure it's the same as Coor from roi of fiji

            % mask_roi = zeros(size(av_image));
            % for pxlIndx = 1:size(Coor,1)
            %     mask_roi(Coor(pxlIndx,2), Coor(pxlIndx,1)) = 1;
            % end
            % mask_roi = bwmorph(mask_roi, 'thin', inf);

            close(f1, f10);
            savename = sprintf('kymoROI_%d.mat', nr_of_roi);
            save([DataFolder savename], 'Coor', 'ROI_type', 'start_end_points', 'roi_pixels');

        case 'automatic'
            load([DataFolder 'morphology_img.mat'], 'av_image_og', 'av_image_cropped');
            av_image = av_image_og;
            clear av_image_og

            % Get skeleton out of average image to automatically do region of interest
            % for kymograph

            % note: you can only use this for small vessels where there is
            % only a single rbc at a time, and there cannot be any
            % branching. 
            answer = questdlg({'WARNING: In order to be able to do an automatic ROI detection, the following conditions must be met: ', '', ...
                '- The vessel is small, only one RBC can pass at a time', ...
                '- There is no branching in the vessel', ...
                '- The vessel is horizontal', ...
                '', 'Do you wish to continue with the automatic segmentation of this vessel?'});
            if matches(answer, 'No')
                continue
            end

            % based on the vessel being brighter than the rest of the
            % image:
            min_val = min(av_image, [], 'all');
            max_val = max(av_image, [], 'all');
            cutoff = ((max_val-min_val)*0.75) + min_val; % take top 75 % of your values as "vessel"

            
            roi_verification = av_image;
            roi_verification(roi_verification>cutoff) = NaN;
            roi_verification(1:size(av_image,1) - size(av_image_cropped,1),:) = 0;

            f2 = figure;
            imagesc(roi_verification)
            f2.Position = [20 60 f2.Position(3) f2.Position(4)];
            title('Detected Vessel')
            
            % Get Mask
            roi_pixels = find(isnan(roi_verification));
            [row, col] = find(isnan(roi_verification));
            Coor = [col, row];

            mask = zeros(size(av_image));
            mask(roi_pixels) = 1;

            % make sure that there is no pixels skipped vertically, and
            % that single pixels are removed:
            for ind_col = 1:size(mask, 2)
                column = mask(:, ind_col);
                first_pix = find(column,1,'first');
                last_pix = find(column,1,'last');
                column(first_pix:last_pix) = 1;
                if first_pix == last_pix
                    column = zeros(size(mask,1),1);
                end
                mask(:,ind_col) = column;
            end

            f3 = figure;
            imagesc(mask)
            f3.Position = [700 60 f2.Position(3) f2.Position(4)];

            answer = questdlg('Does it make sense?');
            if matches(answer, 'Yes')
                roi_pixels = find(mask);
            else
                disp('pls fix')
                pause
            end

            close(f1, f2, f3, f10)
            savename = sprintf('kymoROI_%d.mat', nr_of_roi);
            save([DataFolder savename], 'Coor', 'ROI_type', 'roi_pixels', 'mask');

        case 'block'
            roi_pixels = drawrectangle;
            roiX = round([roi_pixels.Position(1), roi_pixels.Position(1)+roi_pixels.Position(3)]);
            roiY = round([roi_pixels.Position(2), roi_pixels.Position(2)+roi_pixels.Position(4)]);

            Coorx = (roiX(1):roiX(2));
            Coory = (roiY(1):roiY(2));
            Coor = [];
            for indy = 1:size(Coory, 2)
                Cooradd = [Coorx' repmat(Coory(indy), size(Coorx,2), 1)];
                Coor = [Coor; Cooradd];
            end

            % Coor is saved as x, y. Safety check:
            % bla(Coor(:,2), Coor(:,1)) = 0;
            % figure
            % imagesc(bla)

            load([DataFolder 'AcqInfos.mat'], 'AcqInfoStream')
            block_height_um = ( roiY(2)-roiY(1) ) * AcqInfoStream.pxl_sz_y;

            close(f1, f10)
            savename = sprintf('kymoROI_%d.mat', nr_of_roi);
            save([DataFolder savename], 'Coor', 'roiX', 'roiY', 'ROI_type', 'block_height_um');


        case 'block_fixed_height'
            block_height_um = 3; %um            
            % one eurytrocyte mouse is about 6 um, so make box 3 um as default

            roi_pixels = drawpolyline; % note: can also do with drawline but then you have to hold the mouse
            roi_pixels.Position = round(roi_pixels.Position);
            start_end_points = roi_pixels.Position;

            roiX = [start_end_points(1), start_end_points(2)];
            Coorx = roiX(1):roiX(2);
                
            y_height = round(mean(start_end_points([3, 4])));
            load([DataFolder 'AcqInfos.mat'], 'AcqInfoStream')
            roiY = [round(y_height - (block_height_um/AcqInfoStream.pxl_sz_y) ), ...
                round(y_height + (block_height_um/AcqInfoStream.pxl_sz_y) )];
            Coory = roiY(1):roiY(2);
  
            Coor = [];
            for indy = 1:size(Coory, 2)
                Cooradd = [Coorx' repmat(Coory(indy), size(Coorx,2), 1)];
                Coor = [Coor; Cooradd];
            end

            f2 = figure;
            roi_verification = av_image;
            roi_verification(Coor(:,2), Coor(:,1)) = 0;
            imagesc(roi_verification)
            colormap('gray')
            f2.Position = [20 60 f2.Position(3) f2.Position(4)];
            
            answer = questdlg('Does it make sense?');
            if matches(answer, 'Yes')
                close(f2)
            else
                disp('pls fix')
                pause
            end

            close(f1, f10)
            savename = sprintf('kymoROI_%d.mat', nr_of_roi);
            save([DataFolder savename], 'Coor', 'roiX', 'roiY', 'ROI_type', 'block_height_um');

        case 'perpendicular_block'
            roi_pixels = drawrectangle;
            roiX = round([roi_pixels.Position(1), roi_pixels.Position(1)+roi_pixels.Position(3)]);
            roiY = round([roi_pixels.Position(2), roi_pixels.Position(2)+roi_pixels.Position(4)]);

            Coorx = (roiX(1):roiX(2));
            Coory = (roiY(1):roiY(2));
            Coor = [];
            for indy = 1:size(Coory, 2)
                Cooradd = [Coorx' repmat(Coory(indy), size(Coorx,2), 1)];
                Coor = [Coor; Cooradd];
            end

            % Coor is saved as x, y. Safety check:
            % bla(Coor(:,2), Coor(:,1)) = 0;
            % figure
            % imagesc(bla)

            load([DataFolder 'AcqInfos.mat'], 'AcqInfoStream')
            block_height_um = ( roiY(2)-roiY(1) ) * AcqInfoStream.pxl_sz_y;

            close(f1, f10)
            savename = sprintf('kymoROI_%d.mat', nr_of_roi);
            save([DataFolder savename], 'Coor', 'roiX', 'roiY', 'ROI_type', 'block_height_um');

        case 'polygon'
            close(f1)
            %     f1 = figure;
            %     imagesc(av_image)
            %     colormap('gray')
            %     roi = drawpolygon;
            %     close(f1)
            %     roi.Position = round(roi.Position);
            %     roiX = round([roi.Position(1), roi.Position(1)+roi.Position(3)]);
            %     roiY = round([roi.Position(2), roi.Position(2)+roi.Position(4)]);
            %     kymo = mean(dat(roiY(1):roiY(2), roiX(1):roiX(2),:), 1);
            %     kymo = reshape(kymo, [], size(dat, 3));
            %     f1 = figure;
            %     imagesc(kymo')
            %     close(f1)


        case 'Fiji'
            close(f1)
            % Open ImageJ/Fiji
            % Load the sum.tiff or sum.png file
            % Hold the line button down to make a line ROI by hand (freestyle)
            % Draw the ROI along the vessel
            % File - Save As - Selection - name.roi

            ROIlist = dir(fullfile(ROIpath, [datName,'*.zip']));
            if isempty (ROIlist)
                ROIlist = dir(fullfile(ROIpath, '*.roi'));
            end

            if isempty(ROIlist) %MB: if your roi list is empty
                disp('no roi found')
                return
                % write code to get ROI
            end

            ROIname = ROIlist(1).name;
            sROI = ReadImageJROI( fullfile(ROIpath,ROIname)); % load imageJ ROIs
            if strcmp(ROIname(end-3:end),'.roi')
                tmpROI = sROI;
                sROI = cell(1,1);
                sROI{1} = tmpROI;
            end

            % MB: sROI are structures with information about the ROI, but
            % only need Coor.
            for ROIindx = 1:length(sROI)
                Coor = sROI{ROIindx}.mnCoordinates;
                Coor = Coor( Coor(:,1)>0 & Coor(:,2)>0, :);

                % kymoImg = zeros(size(dat,3), size(Coor,1));
                % for frmIndx = 1:1:size(dat,3)
                %     current_frame = dat(:,:,frmIndx);
                %     for pxlIndx = 1:1:size(Coor,1)
                %         kymoImg(frmIndx, pxlIndx) = current_frame(Coor(pxlIndx,2), Coor(pxlIndx,1));
                %     end
                % end
            end

    end

    %% save, check if you want to do more
    % savename = sprintf('kymoROI_%d.mat', nr_of_roi);
    % save([DataFolder savename], 'Coor', 'ROI_type');

    more_to_add_answer = questdlg('More ROI to add?');
    if matches(more_to_add_answer, 'No')
        more_to_add = 0;
    end

end
end

%%
function [ind, label] = draw_line(p1,p2,image_size)
%DRAWLINE Returns the geometric space (matrix indices) occupied by a line segment
%in a MxN matrix.  Each line segment is defined by two endpoints.
%
%   IND = DRAWLINE(P1, P2, IMAGE_SIZE) returns the matrix indices
%   of the line segment with endpoints p1 and p2.
%   If both points are out of the image boundary no line is drawn and an error will appear.
%   If only one of the endpoints is out of the image boundary a line is still drawn.
%
%       ARGUMENT DESCRIPTION:
%                       P1 - set of endpoints (Nx2). ([row column; ...])
%                       P2 - set of endpoints that connect to p1 (Nx2). ([row column; ...])
%               IMAGE_SIZE - vector containing image matrix dimensions,
%                            where the first element is the number of rows
%                            and the second element is the number of
%                            columns.
%
%       OUTPUT DESCRIPTION:
%                      IND - matrix indices occupied by the line segments.
%                    LABEL - label tag of each line drawn (from 1 to N).
%
%
%       1
%    1 _|_ _ _ _> COLUMNS
%       |_|_|_|_
%       |_|_|_|_
%       |_|_|_|_
%       V
%      ROWS
%
%   Example
%   -------------
%   BW = zeros(250,250);
%   p1 = [10 10; 23 100; -14 -40];
%   p2 = [50 50; 90 100;  50  50];
%   [ind label] = drawline(p1,p2,[250 250]); % OR ...drawline(p2,p1,...
%   BW(ind) = label;
%   figure, imshow(BW,[])
%
% See also line, ind2sub.
% Credits:
% Daniel Simoes Lopes
% ICIST
% Instituto Superior Tecnico - Universidade Tecnica de Lisboa
% danlopes (at) civil ist utl pt
% http://www.civil.ist.utl.pt/~danlopes
%
% June 2007 original version.
% Input verification.
if max(size(p1) ~= size(p2))
    error('The number of points in p1 and p2 must be the same.')
end
if length(size(image_size)) ~= 2
    error('Image size must be bi-dimensional.')
end
% Cicle for each pair of endpoints.
ind = [];
label = [];
for line_number = 1:size(p1,1)

    % Point coordinates.
    p1r = p1(line_number,1);   p1c = p1(line_number,2);
    p2r = p2(line_number,1);   p2c = p2(line_number,2);

    % Image dimension.
    M = image_size(1); % Number of rows.
    N = image_size(2); % Number of columns.

    % Boundary verification.
    % A- Both points are out of range.
    if  ((p1r < 1 || M < p1r) || (p1c < 1 || N < p1c)) && ...
            ((p2r < 1 || M < p2r) || (p2c < 1 || N < p2c)),
        error(['Both points in line segment n� ', num2str(line_number),...
            ' are out of range. New coordinates are requested to fit',...
            ' the points in image boundaries.'])
    end

    % Reference versors.
    % .....r..c.....
    eN = [-1  0]';
    eE = [ 0  1]';
    eS = [ 1  0]';
    eW = [ 0 -1]';

    % B- One of the points is out of range.
    if (p1r < 1 || M < p1r) || (p1c < 1 || N < p1c) || ...
            (p2r < 1 || M < p2r) || (p2c < 1 || N < p2c),
        % ....Classify the inner and outer point.
        if     (p1r < 1 || M < p1r) || (p1c < 1 || N < p1c)
            out = [p1r; p1c];  in = [p2r; p2c];
        elseif (p2r < 1 || M < p2r) || (p2c < 1 || N < p2c)
            out = [p2r; p2c];  in = [p1r; p1c];
        end
        % Vector defining line segment.
        v = out - in;
        aux = sort(abs(v)); aspect_ratio = aux(1)/aux(2);
        % Vector orientation.
        north = v'*eN;
        west  = v'*eW;              east  = v'*eE;
        south = v'*eS;
        % Increments.
        deltaNS = [];
        if north > 0, deltaNS = -1; end
        if south > 0, deltaNS =  1; end
        if isempty(deltaNS), deltaNS = 0; end
        deltaWE = [];
        if east > 0, deltaWE =  1; end
        if west > 0, deltaWE = -1; end
        if isempty(deltaWE), deltaWE = 0; end
        % Matrix subscripts occupied by the line segment.
        if abs(v(1)) >= abs(v(2))
            alpha(1) = in(1); beta(1) = in(2);
            iter = 1;
            while (1 <= alpha(iter)) && (alpha(iter) <= M) && ...
                    (1 <=  beta(iter)) && (beta(iter)  <= N),
                alpha(iter+1) = alpha(iter) + deltaNS;              % alpha grows throughout the column direction.
                beta(iter+1)  = beta(iter)  + aspect_ratio*deltaWE; % beta grows throughout the row direction.
                iter = iter + 1;
            end
            alpha = round(alpha(1:end-1)); beta = round(beta(1:end-1));
            ind = cat(2,ind,sub2ind(image_size,alpha,beta));
            label = cat(2,label,line_number*ones(1,max(size(alpha))));
        end
        % ...
        if abs(v(1)) < abs(v(2))
            alpha(1) = in(2); beta(1) = in(1);
            iter = 1;
            while (1 <= alpha(iter)) && (alpha(iter) <= N) &&...
                    (1 <=  beta(iter)) && (beta(iter)  <= M),
                alpha(iter+1) = alpha(iter) + deltaWE;              % alpha grows throughout the row direction.
                beta(iter+1)  = beta(iter)  + aspect_ratio*deltaNS; % beta grows throughout the column direction.
                iter = iter + 1;
            end
            alpha = round(alpha(1:end-1)); beta = round(beta(1:end-1));
            ind = cat(2,ind,sub2ind(image_size,beta,alpha));
            label = cat(2,label,line_number*ones(1,max(size(alpha))));
        end
        clear alpha beta
        continue
    end
    % C- Both points are in range.
    in = [p1r; p1c];  out = [p2r; p2c]; % OR in = p2; out = p1;
    % Vector defining line segment.
    v = out - in;
    aux = sort(abs(v)); aspect_ratio = aux(1)/aux(2);
    % Vector orientation.
    north = v'*eN;
    west  = v'*eW;              east  = v'*eE;
    south = v'*eS;
    % Increments.
    deltaNS = [];
    if north > 0, deltaNS = -1; end
    if south > 0, deltaNS =  1; end
    if isempty(deltaNS), deltaNS = 0; end
    deltaWE = [];
    if east > 0, deltaWE =  1; end
    if west > 0, deltaWE = -1; end
    if isempty(deltaWE), deltaWE = 0; end
    % Matrix subscripts occupied by the line segment.
    row_range = sort([p1r p2r]);
    col_range = sort([p1c p2c]);
    if abs(v(1)) >= abs(v(2))
        alpha(1) = in(1); beta(1) = in(2);
        iter = 1;
        while (row_range(1) <= alpha(iter)) && (alpha(iter) <= row_range(2)) && ...
                (col_range(1) <=  beta(iter)) && (beta(iter)  <= col_range(2)),
            alpha(iter+1) = alpha(iter) + deltaNS;              % alpha grows throughout the column direction.
            beta(iter+1)  = beta(iter)  + aspect_ratio*deltaWE; % beta grows throughout the row direction.
            iter = iter + 1;
        end
        alpha = round(alpha(1:end-1)); beta = round(beta(1:end-1));
        ind = cat(2,ind,sub2ind(image_size,alpha,beta));
        label = cat(2,label,line_number*ones(1,max(size(alpha))));
    end
    % ...
    if abs(v(1)) < abs(v(2))
        alpha(1) = in(2); beta(1) = in(1);
        iter = 1;
        while (col_range(1) <= alpha(iter)) && (alpha(iter) <= col_range(2)) &&...
                (row_range(1) <=  beta(iter)) && (beta(iter)  <= row_range(2)),
            alpha(iter+1) = alpha(iter) + deltaWE;              % alpha grows throughout the row direction.
            beta(iter+1)  = beta(iter)  + aspect_ratio*deltaNS; % beta grows throughout the column direction.
            iter = iter + 1;
        end
        alpha = round(alpha(1:end-1)); beta = round(beta(1:end-1));
        ind = cat(2,ind,sub2ind(image_size,beta,alpha));
        label = cat(2,label,line_number*ones(1,max(size(alpha))));
    end
    clear alpha beta
    continue
end

end
