% TODO SAVE DIAMETER IN OVERVIEW.MAT RESULTS.DIAMETER

% altered from vesselAngleD
% note: the old way of doing the vessel diameter calculation does not work
% in our case because of the big difference between the x and y axis.

% FIX, the diameter as calculated before seems off...

function Diameter_Vessel(DataFolder, ROIname, overwrite)

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

% if isempty(ROI_list)
%     error('No ROI found, function exited.')
% end

shiftFrmMax = 5;% the location to calculate the diameter is determined
% as 'shiftFrmMax' pixels away from the pixel with maximum intensity
% along the line ROI
tailflatten = [5 5];% the number of points on both ends to be set as 0
% this is to deal with close-by blood vessels causing a large gradient
% on the ends

%% load data
load([DataFolder 'morphology_img.mat'], 'av_image_cropped', 'av_image_og', 'proportion_XY_ax')
load([DataFolder 'AcqInfos.mat'], 'AcqInfoStream')

size_y = size(av_image_og, 1)*proportion_XY_ax(1); %usually stays same, *1
size_x = size(av_image_og, 2)*proportion_XY_ax(2);
av_im_proportion = imresize(av_image_og, [size_y, size_x]);
rescale_factor_y = size(av_im_proportion, 1)/size(av_image_og,1);
rescale_factor_x = size(av_im_proportion, 2)/size(av_image_og,2);
pxlY = AcqInfoStream.pxl_sz_y/rescale_factor_y;
pxlX = AcqInfoStream.pxl_sz_x/rescale_factor_x;

% if pxlY ~= pxlX % does not work
if abs(pxlY-pxlX) > 0.0001
    error('The Rescaling went wrong somewhere! Fix.')
end

switch AcqInfoStream.Vessel
    case 'Capillary'
        % estVesselDia = round(3/pxlX);
        % estVesselDia = 12;% unit: pixel
        est_diameter_um = 4;
        % est_diameter_pix =
    case {'Arteriole', 'Venule'}
        % estVesselDia = 10/pxlX;
        % estVesselDia = 70;
        est_diameter_um = 10;
end
estVesselDia = round(est_diameter_um/pxlX);


%% Go per roi - get pxlSize of orthogonal line (pxlSizeD1)
for ind_ROI = 1:length(ROI_list)
    currentROI = ROI_list{ind_ROI};
    load([DataFolder currentROI],  'ROI_type', 'Coor');

    switch ROI_type
        case 'line'
            load([DataFolder currentROI], 'start_end_points');
            segments = size(start_end_points,1) -1;
            if segments>2
                temp_start = start_end_points(round(segments/2),:);
                temp_end = start_end_points(round(segments/2)+1,:);
                start_end_points = [temp_start; temp_end];
            end

            start_end_points([1,2]) = start_end_points([1,2])*rescale_factor_x;
            start_end_points([3,4]) = start_end_points([3,4])*rescale_factor_y;

            slope = (start_end_points(4)-start_end_points(3)) / (start_end_points(2) - start_end_points(1) );
            slopeD1 = -1/slope;
            theta = atan(abs(slope));

            if theta == 0 || theta == pi/2 % line is completely horizontal or if angle is 90 degrees
                pxlSizeD1 = pxlY;
            else
                if tan(theta)<1
                    % pxlSize = sqrt((pxlY*tan(theta))^2 + pxlX^2); % unit: um.
                    pxlSizeD1 = sqrt(pxlY^2 + (pxlX*tan(theta))^2);
                else
                    % pxlSize = sqrt(pxlY^2 + (pxlX/tan(theta))^2);
                    pxlSizeD1 = sqrt((pxlY/tan(theta))^2 + pxlX^2);
                end
            end

            % Coor(:,1) = round(Coor(:,1)*rescale_factor_x);
            % Coor(:,2) = round(Coor(:,2)*rescale_factor_y);
            p1 = start_end_points(1,:);
            p2 = start_end_points(2,:);
            [roi_pixels, ~] = draw_line(fliplr(p1), fliplr(p2), size(av_im_proportion)); % not matlab function (drawline)
            [row, col] = ind2sub(size(av_im_proportion), roi_pixels);
            Coor = [col', row'];

        case 'automatic'
            load([DataFolder currentROI], 'mask');
            mask = imresize(mask, [size_y, size_x]);
            skeleton = bwskel(logical(mask));


        case 'perpendicular_line'
            pxlSizeD1 = pxlY;

        case 'perpendicular_block'
            pxlSizeD1 = pxlY;

    end

    %% get brightness along ROI line
    av_image2 = av_im_proportion;
    pxlValue = zeros(size(Coor,1),1); % the pixel brightness along the ROI line
    for pxlIndx = 1:1:size(Coor,1)
        av_image2(Coor(pxlIndx,2), Coor(pxlIndx,1)) = min(av_im_proportion(:));
        pxlValue(pxlIndx) = av_im_proportion(Coor(pxlIndx,2), Coor(pxlIndx,1));
    end

    %% calculate diamter
    % MB: av_image2 is the same as av_image, but the parts of the ROI are now the
    % same as the minimum value of the entire image (in this case 0...). The
    % pxlValue on the other hand takes the values of the pixels of the roi.
    pxlValue2 = pxlValue;
    % itr = 0;
    % diameter = [];
    % while isempty(diameter) %% MB: why is this a while function? doesnt ever seem to not be a diameter after the first try...
    %     itr = itr+1;
    [~,maxI] = max(pxlValue2); % MB: index of max value of line roi, 55 in example
    maxI = maxI + shiftFrmMax; % MB:add 5...?
    %MB:why max -2? in example, below pxlValue2(58:62) = 0
    pxlValue2(max(maxI-2,1):min(length(pxlValue2),maxI+2)) = 0; % put 5 maximum pixels to 0
    % if maxI < size(Coor, 1)/4 || maxI > ( 3/4* size(Coor, 1)) % MB: if the max value is too close to either end of the ROI, take the middle of the ROI
    %     maxI = floor( size(Coor, 1)/2 );
    % end
    while maxI < size(Coor, 1)/5 || maxI > ( 4/5* size(Coor, 1)) % MB: if the max value is too close to either end of the ROI, choose new max index
        pxlValue2(max(maxI-2,1):min(length(pxlValue2),maxI+2)) = 0;
        [~,maxI] = max(pxlValue2);
    end

    slopeVD = slopeD1;
    pxlSizeVD = pxlSizeD1;
    target_x = Coor(maxI,1);
    target_y = Coor(maxI,2);
    if abs(slopeVD) > 1
        row = max(target_y-estVesselDia,1):min(target_y+estVesselDia,size(av_im_proportion,1));
        column = round(1/slopeVD*(row - target_y) + target_x);
        vesselD =  pxlSizeVD*( max(-estVesselDia,1-target_y):1:min(size(av_im_proportion,1)-target_y,estVesselDia) );
    else
        column = max(target_x-estVesselDia,1):min(target_x+estVesselDia,size(av_im_proportion,2)); %21:45
        if slopeVD == 0
            row = target_y*ones(size(column));
        else
            row = round( slopeVD*(column - target_x) + target_y );
        end
        vesselD =  pxlSizeVD*( max(-estVesselDia,1-target_x):1:min(size(av_im_proportion,2)-target_x,estVesselDia) );
    end
    validIndx = (row<size(av_im_proportion,1))& (row>0) & (column>1)&(column<size(av_im_proportion,2));
    row = row(validIndx); % values for line crossection vessel
    column = column(validIndx);
    vesselD = vesselD(validIndx);
    lineIndx = size(av_im_proportion,1)*(column - 1) + row;
    lineImg = av_im_proportion(:);
    vesselB = lineImg(lineIndx); % vesselB is brightness along the ortogonal line

    % interpolate the brightness profile along the ROI line (MB = diff1)
    interpF = 10;
    vesselD_interp = interp(vesselD,interpF);
    vesselB_interp = interp(vesselB,interpF);
    diff1 = double( vesselB_interp(interpF+1:end ) - vesselB_interp(1:end-interpF) );
    % the brighness should have the maximum drop/increase at the
    % vessel wall.
    idx = ( 1:1:length(diff1) ).';
    diff1( (idx <tailflatten(1)*interpF | idx>max(idx)-tailflatten(2)*interpF) ) = 0;
    [~,indx1] = max(abs(diff1)); %% checked until here TODO
    [~,maxBIndx] = max(vesselB_interp);
    while ( abs(vesselB_interp(indx1) - min(vesselB))>0.45*(max(vesselB)-min(vesselB)) ) ||...
            ( (maxBIndx - indx1)* diff1(indx1) < 0 ) %% MB if certain condition (not sure what), put the found max at 0 and try again
        diff1(indx1) = 0;
        [~,indx1] = max(abs(diff1));
    end
    pk1 = diff1(indx1);
    indx2 = indx1;
    % in case vessel wall and the vessel tube has more brightness
    % drop than the difference between other wall and background
    if pk1>0
        diff1(1:indx1) = 0;
        diff1(diff1>0) = 0;
    elseif pk1<0
        diff1(indx1+1:end) = 0;
        diff1(diff1<0) = 0;
    end
    diff1bk = diff1;
    while (abs(indx2 - indx1)/interpF < estVesselDia/3) || (abs(vesselB_interp(indx2) - min(vesselB))>0.45*(max(vesselB)-min(vesselB)))
        diff1(indx2) = 0;
        [~,indx2] = max(abs(diff1));
        if indx2 < 10 || (indx2>length(diff1) - 10)
            diff1 = diff1bk;
            indx2 = indx1;
            while (abs(indx2 - indx1)/interpF < estVesselDia/2)
                [~,indx2] = max(abs(diff1));
                diff1(indx2) = 0;
                if indx2 < 10 || (indx2>length(diff1) - 10)
                    break
                end
            end
            break;
        end
    end

    indx1 = indx1 + 1;
    indx2 = indx2 + 1;
    crIndx = round( (indx1 + indx2)/2/10 );

    %% plot
    % image with roi and orthoganal line
    figure;
    t = tiledlayout('flow');
    title(t, currentROI, 'Interpreter', 'none')

    nexttile
    imagesc(av_image2)
    % axis('image')
    hold on;
    plot(column,row,'color','r','LineWidth',2);
    plot(column (crIndx),row (crIndx),'color','r','Marker','o','MarkerSize', 10);
    hold off

    % plot brightness of orthogonal line
    nexttile
    % h2 = figure(2);
    % set(h2,'Units','Normalized','Position',[0.6430 0.3271 0.2734 0.3646]); clf;
    plot(vesselD,vesselB,'color','k','LineWidth',2,...
        'Marker','.','MarkerEdgeColor','m','MarkerSize',20); hold on;
    yLim = get(gca,'YLim');
    plot( [vesselD_interp(indx1) vesselD_interp(indx1)], yLim,'color',[0.5 0.5 0.5],'LineWidth',1.5); % plot cut off lines
    plot( [vesselD_interp(indx2) vesselD_interp(indx2)], yLim,'color',[0.5 0.5 0.5],'LineWidth',1.5);
    hold off;
    diameter = abs( vesselD_interp(indx2) - vesselD_interp(indx1) ); % diameter is the distance between cutoff lines
    text( min(vesselD_interp(indx1), vesselD_interp(indx2)), mean(yLim),...
        sprintf('vessel diameter: %.2f %s',diameter, '\mum'),'FontSize',16 );

    % end
end

% close(h1, h2)
end



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





%% old way

% function Diameter_Vessel(DataFolder, ROIname, overwrite)
% % To calculate blood vessel diameter:
% % The function will first draw a line ideally perpendicular to the blood
% % vessel axis, and plot the intensity along the line. The line was determined
% % as the orthogonal intersect of either the line ROI used to extract the
% % kymograph ROI, or the vessel edges detected by a Radon transform. The
% % intensity profile along the line was then upsampled by a factor of 10,
% % after which the edges of the blood vessels were determined as the two
% % positions with slopes of maximal amplitude and opposite signs. The diameter
% % was finally calculated as the distance between these two edges.
%
% % pxlX: pixel size along x dimension in the morphology image.
% % pxlY: pixel size along y dimension in the morphology image;
% % ROIindx: since a single imageJ ROI file could contain multiple ROIs for
% %          different blood vessel segments in the morphology image, ROI
% %          index is required to determine the vessel segment to analyze.
% % vesselType: 'capillary' or 'large'. The line to determine the blood
% %          vessel diameter is longer for 'large' vessels than 'capillary'
% % cal_vesselD: if calculate blood vessel diameter; if not, the function
% %          will only calibrate the pixel size.
% % outputdir: the folder to save the results in.
% % method_input: the method to determine the line to plot the intensity profile.
% %             if method_input == 1, the line was determined
% %             as the orthogonal intersect of the line ROI
% %             if method_input == 2, the line was determined
% %             as the orthogonal intersect of the blood vessel edges
% %             detected by Radon transform.
%
% %% set up
% if ~strcmp(DataFolder(end), filesep)
%     DataFolder = [DataFolder filesep];
% end
%
% % get list of ROI made
% if exist('ROIname', 'var') && ~matches(ROIname, 'auto_list')
%     ROI_list = ROIname;
% else
%     ROI_list = dir([DataFolder 'kymoROI*.mat']);
%     ROI_list = struct2cell(ROI_list);
%     ROI_list = ROI_list(1,:);
% end
%
% if isempty(ROI_list)
%     error('No ROI found, function exited.')
% end
%
% if ~exist("overwrite", 'var')
%     overwrite = 0;
% end
%
% shiftFrmMax = 5;% the location to calculate the diameter is determined
% % as 'shiftFrmMax' pixels away from the pixel with maximum intensity
% % along the line ROI
% tailflatten = [5 5];% the number of points on both ends to be set as 0
% % this is to deal with close-by blood vessels causing a large gradient
% % on the ends
%
% if ~exist('method','var')
%     method_input = [1, 2]; % can only be 1 (determine orghogonal intersect with ROI line),
%     % 2 (determine orthogonal intersect with edges detected by Radon Transform
%     % or [1, 2](do calculations using both methods).
% end
%
% %% load data
% load([DataFolder 'AcqInfos.mat'], 'AcqInfoStream');
% % datSize = [AcqInfoStream.ny+AcqInfoStream.ny_extra AcqInfoStream.nx];
% load([DataFolder 'morphology_img.mat'], 'av_image_og', 'proportion_XY_ax')
% av_image = av_image_og;
% % av_image_proportional = imresize(av_image, [proportion_XY_ax(1)*size(av_image,1) proportion_XY_ax(2)*size(av_image,2)]);
% % pxlSize = AcqInfoStream.pxl_sz_x;
%
% switch AcqInfoStream.Vessel
%     case 'Capillary'
%         estVesselDia = 12;% unit: pixel
%         % est_diameter_um = 4;
%         % est_diameter_pix =
%     case {'Arteriole', 'Venule'}
%         estVesselDia = 70;
%         % est_diameter_um = 10;
% end
%
% pxlX = AcqInfoStream.pxl_sz_x;
% pxlY = AcqInfoStream.pxl_sz_y;
%
% %% Go per roi
% for ind_ROI = 1:length(ROI_list)
%     currentROI = ROI_list{ind_ROI};
%     load([DataFolder currentROI],  'Coor', 'PixelSize', 'Slope');
%     pxlSize = PixelSize.pxlSize;
%
%     % patch, can delete later
%     if ~exist('Slope', 'var')
%         Pixelsize_ROI(DataFolder, currentROI);
%         load([DataFolder currentROI], 'Slope');
%     end
%
%     av_image2 = av_image;
%     pxlValue = zeros(size(Coor,1),1); % the pixel brightness along the ROI line
%     for pxlIndx = 1:1:size(Coor,1)
%         av_image2(Coor(pxlIndx,2), Coor(pxlIndx,1)) = min(av_image(:));
%         pxlValue(pxlIndx) = av_image(Coor(pxlIndx,2), Coor(pxlIndx,1));
%     end
%
%     checkpoint = Coor(round(size(Coor,1)/2),:);
%
%
%     % MB: av_image2 is the same as av_image, but the parts of the ROI are now the
%     % same as the minimum value of the entire image (in this case 0...). The
%     % pxlValue on the other hand takes the values of the pixels of the roi.
%     pxlValue2 = pxlValue;
%     itr = 0;
%     diameter = [];
%     while isempty(diameter) %% MB: why is this a while function? doesnt ever seem to not be a diameter after the first try...
%         itr = itr+1;
%         [~,maxI] = max(pxlValue2); % MB: index of max value of line roi, 55 in example
%         maxI = maxI + shiftFrmMax; % MB:add 5...?
%         %MB:why max -2? in example, below pxlValue2(58:62) = 0
%         pxlValue2(max(maxI-2,1):min(length(pxlValue2),maxI+2)) = 0; % MB:max(bla,1) takes max of either bla or 1 if 1 is max
%         % extent = 10; % the number of pixels along the line ROI to determine the line slope
%         if maxI < size(Coor, 1)/4 || maxI > ( 3/4* size(Coor, 1)) % MB: ?? what's this?
%             maxI = floor( size(Coor, 1)/2 );
%         end
%
%
%     %% Determine the orthogonal intersect of the vessel length using the ROI line itself.
%     % CoorX = Coor(max(1, maxI-extent) : maxI+extent,1); %MB: example: 50:70,1
%     % CoorY = Coor(max(1, maxI-extent) : maxI+extent,2);
%     % !!!!!!!!!! FIX PER SEGMENT
%     % CoorX = Coor(start_end_points(1,1): start_end_points(2,1), 1); % FIXXX
%     % CoorY = Coor(start_end_points(1,2): start_end_points(2,2), 2);
%
%     %% get slope of orthogonal line
%     CoorX = Coor(:,1);
%     CoorY = Coor(:,2);
%
%     slope = Slope.slope_segments;
%     % slope = (CoorY(end)-CoorY(1)) / (CoorX(end)-CoorX(1));
%     slopeD1 = -1/slope;
%     theta = atan(abs(slope));
%
%     % sanity check:
%     f1 = figure;
%     plot(CoorX, -1*CoorY) %-1 is because image plots from top to bottom
%     hold on
%     x = CoorX(1):CoorX(end);
%     x1 = CoorX(1);
%     y = slope*(x-x1) + CoorY(1);
%     plot(x,-y) % line over ROI. -1 is because image plots from top to bottom
%     x = CoorX(1):CoorX(end);
%     x1 = CoorX(1);
%     y = slopeD1*(x-x1) + CoorY(1);
%     plot(x,-y) % orthogonal line. -1 is because image plots from top to bottom
%     axis('equal');
%     close(f1)
%     clear x x1 y
%
%     if theta == 0 % line is completely horizontal
%         % pxlSize = pxlX;
%         pxlSizeD1 = pxlY;
%     elseif theta == pi/2 % MB: if angle is 90 degrees
%         % pxlSize = pxlY;
%         pxlSizeD1 = pxlX;
%     else
%         if tan(theta)<1
%             % pxlSize = sqrt((pxlY*tan(theta))^2 + pxlX^2); % unit: um.
%             pxlSizeD1 = sqrt(pxlY^2 + (pxlX*tan(theta))^2);
%         else
%             % pxlSize = sqrt(pxlY^2 + (pxlX/tan(theta))^2);
%             pxlSizeD1 = sqrt((pxlY/tan(theta))^2 + pxlX^2);
%         end
%     end
%     % MB:pxlSize is for ROI_line, pxlSizeD1 is for orthogonal line
%
%
%
%     % %% Determine the orthogonal intersect of the vessel length using Radon transform.
%     % % take line ROI and enlarge by estimated vessel diameter so you get a
%     % % cropped square with the vessel in it
%     % cropImg = av_image(max(Coor(maxI,2)-estVesselDia,1):min(Coor(maxI,2)+estVesselDia, size(av_image,1)), ...
%     %     max(Coor(maxI,1)-estVesselDia,1): min(Coor(maxI,1) +estVesselDia,size(av_image,2)));
%     % edgeThr = 0.7; % the threshold for Radon transform
%     % edges = edge(cropImg,'Canny',edgeThr,2); % use Radon transform to detect
%     % % the blood vessel edges. MB: edges is logical mask of size cropImg
%     % % with vessel walls
%     % thetavector = 0:1:180;
%     % R = radon(edges,thetavector);
%     %
%     % h2 = figure(2); set(h2,'Units','Normalized','Position',[0.0730 0.4375 0.2734 0.3646]);
%     % clf; ha = subplot(2,2,1);
%     % set(ha, 'Units','Normalized','Position',[0.25 0.53 0.4 0.4]);
%     % imshow(cropImg,[]);
%     % title([['method 2 ',ROIname(tpos1+1 : tpos1+2),' ',ROIname(tpos2+1: tpos2+4),' kymo ROI'], num2str(ROIindx)],'FontSize',15)
%     % hb = subplot(2,2,3);
%     % set(hb,'Units','Normalized','Position',[0.05 0.05 0.4 0.4]);
%     % imshow(edges);
%     % title('Detected Edges','FontSize',15)
%     % hc = subplot(2,2,4);
%     % set(hc,'Units','Normalized','Position',[0.5 0.05 0.4 0.4]);
%     % imshow(R,[]);colormap('jet');
%     % title(['Radon Transform(0-180' char(176) ')'],'FontSize',15)
%     % xlabel(['\Theta (', char(176),')']);
%     %
%     % [~,thetaM] = max(max(R)); % MB get index of max row
%     % slopeD2 = -tan(thetavector(thetaM)/180*pi); % MB: slope of orthogonal line calculated in second method
%     % thetaMr = thetaM/180*pi; % MB: theta max radian
%     % if thetaMr ==0
%     %     pxlSizeD2 = pxlX;
%     % elseif thetaMr ==pi/2
%     %     pxlSizeD2 = pxlY;
%     % else
%     %     if abs( tan(thetaMr) )<1
%     %         pxlSizeD2 = sqrt((pxlY*tan(thetaMr))^2 + pxlX^2); % unit: um.
%     %     else
%     %         pxlSizeD2 = sqrt(pxlY^2 + (pxlX/tan(thetaMr))^2);
%     %     end
%     % end
%     % slopeD = [slopeD1 slopeD2];
%     % pxlSizeD = [pxlSizeD1 pxlSizeD2]; % combine the calculations from two methods together
%
%    slopeD = [slopeD1];
%     pxlSizeD = [pxlSizeD1];
%
%     %% calculate the blood vessel diameter using the intensity profile along the orthogonal intersect decided in the above sections
%     for method = method_input(1):1:max(1,length(method_input))
%         % titleStr = ['method ', num2str(method), ' ',ROIname(tpos1+1 : tpos1+2),' ',ROIname(tpos2+1: tpos2+4),' kymo ROI'];
%         % titleStr(regexp(titleStr,'_')) = ' ';
%         slopeVD = slopeD(method);
%         pxlSizeVD = pxlSizeD(method);
%         if abs(slopeVD) > 1
%             row= max(Coor(maxI,2)-estVesselDia,1):min(Coor(maxI,2)+estVesselDia,size(av_image,1));
%             column =round(1/slopeVD*(row - Coor(maxI,2)) + Coor(maxI,1));
%             vesselD =  pxlSizeVD*( max(-estVesselDia,1-Coor(maxI,2)):1:min(size(av_image,1)-Coor(maxI,2),estVesselDia) );
%         else
%             column = max(Coor(maxI,1)-estVesselDia,1):min(Coor(maxI,1)+estVesselDia,size(av_image,2)); %21:45
%             if slopeVD == 0
%                 row = Coor(maxI,2)*ones(size(column));
%             else
%                 row = round( slopeVD*(column - Coor(maxI,1)) + Coor(maxI,2) );
%             end
%             vesselD =  pxlSizeVD*( max(-estVesselDia,1-Coor(maxI,1)):1:min(size(av_image,2)-Coor(maxI,1),estVesselDia) );
%         end
%         validIndx = (row<size(av_image,1))& (row>0) & (column>1)&(column<size(av_image,2));
%         row = row(validIndx); % values for line crossection vessel
%         column = column(validIndx);
%         vesselD = vesselD(validIndx);
%         lineIndx = size(av_image,1)*(column - 1) + row;
%         lineImg = av_image(:);
%         vesselB = lineImg(lineIndx);
%
%         h1 = figure(11);
%         set(gcf,'Units','Normalized','Position',[0.4    0.2    0.2  0.56]);
%         cla;
%         set(gca,'Units','Normalized','Position',[0.05 0.05 0.9 0.9]);
%         imshow(av_image2,[]);
%         hold on;
%         plot(column,row,'color','r','LineWidth',2);
%         hold off;
%         % title([titleStr, num2str(ROIindx)],'FontSize',16);
%
%         % fileNameStr = sprintf('%s_kymoROI%02d_vesselD_method%d',ROIname(1:pos-1),ROIindx,method);
%
%         % interpolate the brightness profile along the ROI line (MB = diff1)
%         interpF = 10;
%         vesselD_interp = interp(vesselD,interpF);
%         vesselB_interp = interp(vesselB,interpF);
%         diff1 = double( vesselB_interp(interpF+1:end ) - vesselB_interp(1:end-interpF) );
%         % the brighness should have the maximum drop/increase at the
%         % vessel wall.
%         idx = ( 1:1:length(diff1) ).';
%         diff1( (idx <tailflatten(1)*interpF | idx>max(idx)-tailflatten(2)*interpF) ) = 0;
%         [~,indx1] = max(abs(diff1)); %% checked until here TODO
%         [~,maxBIndx] = max(vesselB_interp);
%         while ( abs(vesselB_interp(indx1) - min(vesselB))>0.45*(max(vesselB)-min(vesselB)) ) ||...
%                 ( (maxBIndx - indx1)* diff1(indx1) < 0 ) %% MB if certain condition (not sure what), put the found max at 0 and try again
%             diff1(indx1) = 0;
%             [~,indx1] = max(abs(diff1));
%         end
%         pk1 = diff1(indx1);
%         indx2 = indx1;
%         % in case vessel wall and the vessel tube has more brightness
%         % drop than the difference between other wall and background
%         if pk1>0
%             diff1(1:indx1) = 0;
%             diff1(diff1>0) = 0;
%         elseif pk1<0
%             diff1(indx1+1:end) = 0;
%             diff1(diff1<0) = 0;
%         end
%         diff1bk = diff1;
%         while (abs(indx2 - indx1)/interpF < estVesselDia/3) || (abs(vesselB_interp(indx2) - min(vesselB))>0.45*(max(vesselB)-min(vesselB)))
%             diff1(indx2) = 0;
%             [~,indx2] = max(abs(diff1));
%             if indx2 < 10 || (indx2>length(diff1) - 10)
%                 diff1 = diff1bk;
%                 indx2 = indx1;
%                 while (abs(indx2 - indx1)/interpF < estVesselDia/2)
%                     [~,indx2] = max(abs(diff1));
%                     diff1(indx2) = 0;
%                     if indx2 < 10 || (indx2>length(diff1) - 10)
%                         break
%                     end
%                 end
%                 break;
%             end
%         end
%
%         indx1 = indx1 + 1;
%         indx2 = indx2 + 1;
%         crIndx = round( (indx1 + indx2)/2/10 );
%         h3 = figure(3);
%         set(h3,'Units','Normalized','Position',[0.6430 0.3271 0.2734 0.3646]); clf;
%         plot(vesselD,vesselB,'color','k','LineWidth',2,...
%             'Marker','.','MarkerEdgeColor','m','MarkerSize',20); hold on;
%         yLim = get(gca,'YLim');
%         plot( [vesselD_interp(indx1) vesselD_interp(indx1)], yLim,'color',[0.5 0.5 0.5],'LineWidth',1.5)
%         plot( [vesselD_interp(indx2) vesselD_interp(indx2)], yLim,'color',[0.5 0.5 0.5],'LineWidth',1.5);
%         % title([titleStr, num2str(ROIindx)],'FontSize',16);
%         hold off;
%         diameter = abs( vesselD_interp(indx2) - vesselD_interp(indx1) );
%         text( min(vesselD_interp(indx1), vesselD_interp(indx2)), mean(yLim),...
%             sprintf('vessel diameter: %.2f %s',diameter, '\mum'),'FontSize',16 );
%         % save(fullfile(outputdir,[fileNameStr,'.mat']), 'diameter','vesselD', 'vesselB','pxlSizeVD',...
%         %     'column','row','crIndx','shiftFrmMax');
%         %
%         %         end
%         h1 = figure(11);
%         hold on;
%         plot(column (crIndx),row (crIndx),'color','r','Marker','o','MarkerSize', 10);
%         % print(h1,'-dpng','-r600',fullfile(outputdir,[fileNameStr,'lines.png']));
%         % print(h3,'-dpng','-r600',fullfile(outputdir,[fileNameStr,'.png']));
%     end
%     print(h2,'-dpng','-r600',fullfile(outputdir,[fileNameStr,'edges.png']));
% end
% end
% end
%
% function [slope,b] = LineFit(x, y)
% lineF =@(f) f(1)*x + f(2) - y;
% f0 = [ (y(end) -y(1))/(x(end)-x(1)) 0];
% [f,~,~] = lsqnonlin(lineF,f0);
% slope = f(1);
% b = f(2);
% end




