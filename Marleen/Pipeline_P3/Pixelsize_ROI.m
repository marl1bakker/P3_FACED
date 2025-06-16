% always calculates with 2 seperate methods

function [PixelSize] = Pixelsize_ROI(DataFolder, ROIname, saveresult)

%% set up
if ~strcmp(DataFolder(end), filesep)
    DataFolder = [DataFolder filesep];
end

if ~exist('saveresult', 'var')
    saveresult = 1;
end

% get list of ROI made
if exist('ROIname', 'var') && ~iscell(ROIname)
    ROI_list = {ROIname};
elseif ~exist('ROIname', 'var')
    ROI_list = dir([DataFolder 'kymoROI*.mat']); % fix
    ROI_list = struct2cell(ROI_list);
    ROI_list = ROI_list(1,:);
end

if isempty(ROI_list)
    error('No ROI found, function exited.')
end

load([DataFolder 'AcqInfos.mat'], "AcqInfoStream");
pxlX = AcqInfoStream.pxl_sz_x;
% pxlY = AcqInfoStream.pxl_sz_y;
if AcqInfoStream.height <= 1 || matches(AcqInfoStream.Acquisition_Type, 'Linescan') % second one because you did linescan of 10 um once (can delete later)
    pxlY = 1;
else
    pxlY = AcqInfoStream.pxl_sz_y_avged;
end

% Go per roi
for ind_ROI = 1:size(ROI_list,2)
    ROIname = ROI_list{ind_ROI};
    load([DataFolder ROIname], 'ROI_type');

    switch ROI_type
        case 'linescan'
            PixelSize.pxlSize = AcqInfoStream.pxl_sz_x;
            Slope.slope = NaN;

        case 'line'
            % My way (MB):
            load([DataFolder ROIname], 'ROI_info');
            % TEMP
            if ~exist('ROI_info', 'var')
                patch_code(DataFolder)
                load([DataFolder ROIname], 'ROI_info');
            end
            start_end_points = ROI_info.start_end_points;
            Coor = ROI_info.Coor;

            nr_of_segments = size(start_end_points, 1)-1;
            length_segments_um = nan(nr_of_segments, 1);
            length_segments_pix = nan(nr_of_segments, 1);
            pxlSize_segments = nan(nr_of_segments, 1);
            slope_segments = nan(nr_of_segments, 1);

            for ind_segm = 1:nr_of_segments
                start_point = start_end_points(ind_segm,:);
                end_point = start_end_points(ind_segm+1,:);

                % new:
                % slope
                slope_segments(ind_segm) = (end_point(2)-start_point(2)) / (end_point(1)-start_point(1));

                % lengths
                a_pix = abs(end_point(2)-start_point(2))+1;
                b_pix = abs(end_point(1)-start_point(1))+1;
                % c_pix = sqrt( a_pix^2 + b_pix^2); %length of roiline in pix

                if nr_of_segments == 1
                    c_pix = length(Coor);
                else
                    % find c pixels number
                    c_start = find(Coor(:,1) == start_point(1) & Coor(:,2) == start_point(2), 1, 'first');
                    c_end = find(Coor(:,1) == end_point(1) & Coor(:,2) == end_point(2), 1, 'first');
                    c_pix = c_end-c_start+1;
                end

                % sometimes it does not include endpoint in coor, so add:
                % TEMP patch, fixed in ROI_Kymograph. Can delete if you ran
                % all ROI once
                if isempty(c_pix) && ind_segm == nr_of_segments
                    Coor(end+1,:) = end_point;
                    % load([DataFolder ROIname], 'roi_pixels');
                    load([DataFolder 'morphology_img.mat'], "av_image");
                    ROI_info.roi_pixels = [ROI_info.roi_pixels, sub2ind(size(av_image), end_point(2), end_point(1))];
                    ROI_info.Coor = Coor;
                    save([DataFolder ROIname], 'ROI_info', '-append');
                    clear roi_pixels av_image 

                    c_start = find(Coor(:,1) == start_point(1) & Coor(:,2) == start_point(2), 1, 'first');
                    c_end = find(Coor(:,1) == end_point(1) & Coor(:,2) == end_point(2), 1, 'first');
                    c_pix = c_end-c_start+1;
                    % warning(['redo the kymograph ' DataFolder ' ' ROIname ])
                    % Plot_Kymograph(DataFolder, 0, {ROIname}, 1)
                    Make_Kymograph(DataFolder, 0, {ROIname}, 1)
                end

                a_um = a_pix * pxlY;
                b_um = b_pix * pxlX;
                c_um = sqrt( a_um^2 + b_um^2 ); % length of roiline in um

                length_segments_um(ind_segm) = c_um;
                length_segments_pix(ind_segm) = c_pix;
                pxlSize_segments(ind_segm) = c_um/c_pix;
            end

            % PixelSize.length_segments_pix = length_segments_pix;
            PixelSize.length_segments_um = length_segments_um;
            PixelSize.pxlSize_segments = pxlSize_segments;
            PixelSize.pxlSize = sum(length_segments_um)/sum(length_segments_pix);

            Slope.slope_segments = slope_segments;

            %% Method 2 Na Ji
            pxlSize_segments_method_2 = nan(nr_of_segments, 1);
            slope_segments_method_2 = nan(nr_of_segments, 1);

            for ind_segm = 1: nr_of_segments
                start_point = start_end_points(ind_segm,:);
                end_point = start_end_points(ind_segm+1,:);

                try % cant do this if x doesnt change (divide by 0)
                    % slope = (end_point(2)-start_point(2)) / (end_point(1)-start_point(1)); % in degrees
                    indsegstart = find(Coor(:,1) == start_point(1) & Coor(:,2) == start_point(2), 1, 'first');
                    indsegend = find(Coor(:,1) == end_point(1) & Coor(:,2) == end_point(2), 1, 'first');
                    slope = LineFit(Coor(indsegstart:indsegend,1), Coor(indsegstart:indsegend,2));
                    clear insegstart indsegend
                    % slope_D1 = -1/slope; % the slope of the orthogonal intersect of the ROI line
                    theta = atan(slope); % theta is in radians
                catch % line is completely vertical:
                    % slope is infinite
                    theta = pi/2;
                    slope = NaN;
                    % slopeD1 = 0;
                end

                % calibrate the pixel size using the slope (i.e. theta) calculated above
                if theta == 0 % line is completely horizontal
                    pxlSize_2 = pxlX;
                    % pxlSizeD1 = pxlY;
                elseif theta == pi/2 % line is completely vertical
                    pxlSize_2 = pxlY;
                    % pxlSizeD1 = pxlX;
                else % line is angled
                    if tan(theta)<1 %MB: if slope is smaller than 45 degrees??
                        pxlSize_2 = sqrt((pxlY*tan(theta))^2 + pxlX^2); % unit: um.
                        % pxlSizeD1 = sqrt(pxlY^2 + (pxlX*tan(theta))^2);
                    else
                        pxlSize_2 = sqrt(pxlY^2 + (pxlX/tan(theta))^2);
                        % pxlSizeD1 = sqrt((pxlY/tan(theta))^2 + pxlX^2);
                    end
                end

                slope_segments_method_2(ind_segm) = slope;
                pxlSize_segments_method_2(ind_segm) = pxlSize_2;

            end
            PixelSize.pxlSize_segments_method_2 = pxlSize_segments_method_2;
            Slope.slope_segments_method_2 = slope_segments_method_2;


            if any(PixelSize.pxlSize_segments./PixelSize.pxlSize_segments_method_2 < 0.95) || ...
                    any(PixelSize.pxlSize_segments./PixelSize.pxlSize_segments_method_2 > 1.05)
                warning('something is up, pixel sizes of different methods vary a lot. Check!')
                disp(pxlSize_segments)
                disp(pxlSize_segments_method_2)
            end

            
        case {'block', 'block_fixed_height'}
            % a block is always straight, and since we average over the
            % y-axis, the "line" covers the distance of the pixel size over
            % the x axis
            load([DataFolder ROIname], 'ROI_info');
            PixelSize.length_segments_um = (ROI_info.roiX(2)-ROI_info.roiX(1)) * pxlX;
            PixelSize.pxlSize = pxlX;
            PixelSize.pxlSize_segments = pxlX;

            Slope.slope = 0;

        case 'automatic'
            load([DataFolder ROIname], 'mask')
            skeleton = bwskel(logical(mask));
            [row, col] = find(skeleton, 1, 'first');
            start_point = [row, col];
            [row, col] = find(skeleton, 1, 'last');
            end_point = [row, col];

            % slope
            slope_segments = (end_point(2)-start_point(2)) / (end_point(1)-start_point(1));

            a_pix = end_point(1)-start_point(1); % pay attention, this is different than case: line, becuase the first one is y and the second one x
            b_pix = end_point(2)-start_point(2);
            c_pix = sqrt( a_pix^2 + b_pix^2); %length of roiline in pix

            a_um = a_pix * pxlY;
            b_um = b_pix * pxlX;
            c_um = sqrt( a_um^2 + b_um^2 ); % length of roiline in um

            PixelSize.length_segments_um = c_um;
            PixelSize.pxlSize = c_um/c_pix;

            Slope.slope = slope_segments;
    end

    if saveresult
        save([DataFolder ROIname], 'PixelSize', 'Slope', '-append');
    end

end

end


function [slope,b] = LineFit(x, y)
lineF =@(f) f(1)*x + f(2) - y;
f0 = [ (y(end) -y(1))/(x(end)-x(1)) 0];
[f,~,~] = lsqnonlin(lineF,f0);
slope = f(1);
b = f(2);
end