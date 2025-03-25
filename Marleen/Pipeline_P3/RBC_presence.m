% Check if there are RBC present over the whole time of the kymograph. If
% not, it will detect periods that are longer than 0.5 sec and save them as
% seperate kymoROI.

function [RBC] = RBC_presence(DataFolder, ROIname, saveresult)

%% set up
if ~strcmp(DataFolder(end), filesep)
    DataFolder = [DataFolder filesep];
end

if ~exist('saveresult', 'var')
    saveresult = 1;
end

% get list of ROI made
if exist('ROIname', 'var') && ~iscell(ROIname) && ~matches(ROIname, 'auto_list')
    ROI_list = {ROIname};
elseif ~exist('ROIname', 'var') || matches(ROIname, 'auto_list')
    ROI_list = dir([DataFolder 'kymoROI*.mat']); % fix
    ROI_list = struct2cell(ROI_list);
    ROI_list = ROI_list(1,:);
end
if isempty(ROI_list)
    error('No ROI found, function exited.')
end

load([DataFolder 'AcqInfos.mat'], "AcqInfoStream");
if isfield(AcqInfoStream, 'FrameRateHzLinescan')
    frmRate = AcqInfoStream.FrameRateHzLinescan;
else
    frmRate = AcqInfoStream.FrameRateHz;
end

% Go per roi
for ind_ROI = 1:size(ROI_list,2)
    ROIname = ROI_list{ind_ROI};
    load([DataFolder ROIname], 'kymoImg');

    % make between 1 and 0
    kymoImg = rescale(kymoImg);


    %% Mask for RBC parts
    % find parts without RBC
    no_rbc = min(kymoImg, [], 2); % if there is a rbc, the minimum is lower than in case of no rbc
    no_rbc(no_rbc>0.4) = 1; %more or less arbitrary cutoff, based on trial & error
    no_rbc(no_rbc<=0.4) = 0;

    % make less noisy
    window = round(0.01*frmRate);
    no_rbc = [zeros(window, 1); no_rbc; zeros(window,1)];  % make buffer
    tmp = zeros(size(no_rbc));
    nrbc = find(no_rbc);
    for ind_nrbc = 1:length(nrbc)
        % Before target
        if sum(no_rbc(nrbc(ind_nrbc)-window:nrbc(ind_nrbc),:))>1
            tmp(nrbc(ind_nrbc)-window:nrbc(ind_nrbc)) = 1;
        end
        % After target
        if sum(no_rbc(nrbc(ind_nrbc):nrbc(ind_nrbc)+window,:))>1
            tmp(nrbc(ind_nrbc):nrbc(ind_nrbc)+window) = 1;
        end
    end
    no_rbc = tmp(window+1:end-window);


    %% Find useful periods (>0.5 sec)
    if sum(no_rbc)>1

        % find longest stretch with RBC flow
        no_rbc = [1; no_rbc; 1]; % start and end with "bad" value
        min_stretch_length = round(0.5*frmRate);
        start_end_length = zeros(1,3);
        indstretch = 0;

        for ind_nrbc = 2:length(no_rbc)-1
            if no_rbc(ind_nrbc) == 0 && ...
                    no_rbc(ind_nrbc-1) == 1 && ...
                    no_rbc(ind_nrbc+1) == 0 % find beginning

                nbc_stretch = find(no_rbc(ind_nrbc+1:end), 1, 'first')-1; % find end

                if nbc_stretch > min_stretch_length
                    indstretch = indstretch+1;
                    start_end_length(indstretch, 1) = ind_nrbc-1;% because you padded no_rbc
                    start_end_length(indstretch, 2) = ind_nrbc-1+nbc_stretch;
                    start_end_length(indstretch, 3) = nbc_stretch;
                end
            end
        end
        no_rbc = no_rbc(2:end-1);

    else
        start_end_length = [1 size(kymoImg,1) size(kymoImg,1)-1];
    end


    %% Save
    no_rbc = find(no_rbc); % to have same as badvals
    RBC.no_rbc = no_rbc;
    RBC.segments = start_end_length;

    if saveresult
        save([DataFolder ROIname], 'RBC', '-append')
    end

    clear ind_nrbc indstretch min_stretch_length nbc_stretch nrbc tmp window RBC no_rbc


    %% Save the kymographs seperately
    alfabet = 'a':'z';

    if size(start_end_length,1) > 1
        for seg_ind = 1:size(start_end_length, 1)
            ind_start = start_end_length(seg_ind,1);
            ind_end = start_end_length(seg_ind,2);
            RBC.ind_start = ind_start;
            RBC.ind_end = ind_end;
            RBC.no_rbc = [];

            % fix: only load cleandata with linescan
            load([DataFolder ROIname], 'ROI_type', 'kymoImg', 'ROI_info')
            kymoImg = kymoImg(ind_start:ind_end,:);

            ROIsubname = [ROIname(1:end-4) alfabet(seg_ind) '.mat'];

            save([DataFolder ROIsubname], 'ROI_type', 'kymoImg', 'ROI_info', 'RBC');

            
        end
    end

end
end
