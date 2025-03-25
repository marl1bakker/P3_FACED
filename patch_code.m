function patch_code(DataFolder)

if ~strcmp(DataFolder(end), filesep)
    DataFolder = [DataFolder filesep];
end

if exist('ROIname', 'var') && ~matches(ROIname, 'auto_list')
    kymograph_list = {ROIname};
else
    kymograph_list = dir([DataFolder 'kymoROI*.mat']);
    kymograph_list = struct2cell(kymograph_list);
    kymograph_list = kymograph_list(1,:);
end

if isempty(kymograph_list)
    disp('No kymographs found, function exited.')
return
end

% load([DataFolder 'AcqInfos.mat'])
% if exist('CleanData', 'var')
%     disp('doesn''t need to be patched')
%     return
% end

%% Both methods - go per kymograph
for ind_kymo = 1:length(kymograph_list)
 
    % load kymograph:
    ROIname = kymograph_list{ind_kymo};
    if ~strcmp(ROIname(end-3:end), '.mat')
        ROIname = [ROIname '.mat'];
    end
    warning('off');
    load([DataFolder ROIname]);
    warning('on');

    if exist('ROI_info', 'var')
        disp('Already patched')
        continue
    end

    load([DataFolder 'AcqInfos.mat'])


    if matches(ROI_type, 'linescan') && ~exist('CleanData', 'var')
          CleanData.num_avs = num_avs;
          CleanData.y_range = y_range;
          CleanData.datsize = datsize;
          save([DataFolder 'AcqInfos.mat'], 'CleanData', '-append');
          clear num_avs y_range datsize CleanData

    elseif ~exist('CleanData', 'var')
        load([DataFolder 'morphology_img.mat'])
        save([DataFolder 'morphology_img.mat'], 'av_image', 'av_image_cropped', 'proportion_XY_ax', 'av_image_og');
        clear av_image av_image_cropped proportion_XY_ax av_image_og

        CleanData.faced_cut = faced_cut;
        CleanData.illumination_profile = illumination_profile;
        CleanData.y_range = y_range;
        CleanData.x_range = x_range;
        save([DataFolder 'AcqInfos.mat'], 'CleanData', '-append');
        clear faced_cut illumination_profile y_range x_range CleanData
    end

        switch ROI_type
            case 'linescan'
                ROI_info = [];
            case {'line', 'perpendicular_line'}
                ROI_info.Coor = Coor;
                ROI_info.start_end_points = start_end_points;
                ROI_info.roi_pixels = roi_pixels;
                clear Coor start_end_points roi_pixels
            case 'automatic'
                ROI_info.Coor = Coor;
                ROI_info.start_end_points = start_end_points;
                ROI_info.roi_pixels = roi_pixels;
                ROI_info.mask = mask;
                clear Coor start_end_points roi_pixels mask
            case {'block', 'block_fixed_height'}
                ROI_info.Coor = Coor;
                ROI_info.roiX = roiX;
                ROI_info.roiY = roiY;
                ROI_info.block_height_um = block_height_um;
                clear Coor roiX roiY block_height_um
        end

    if exist('Velocity_calc', 'var') && exist('skipamt', 'var')
        Velocity_calc.skipamt = skipamt;
        Velocity_calc.est_mm_per_sec = est_mm_per_sec;
        if isfield(Velocity_calc, 'orientation')
            Velocity_calc.orientation = orientation;
        end
        clear skipamt est_mm_per_sec est_sec_to_cross_kymo est_frames_to_cross_kymo orientation
        save([DataFolder ROIname], 'kymoImg', 'ROI_type', 'ROI_info',  'PixelSize', 'Slope', 'Velocity_calc')
        clear kymoImg ROI_type ROI_info PixelSize Slope Velocity_calc

    else
        save([DataFolder ROIname], 'kymoImg', 'ROI_type','ROI_info' )
        clear kymoImg ROI_type ROI_info
    end


end
end
