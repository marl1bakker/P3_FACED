% When doing acquisitions, it is not always clear straight away what the
% vessel type is. Often, I write things down like "plunging vessel on the
% left" or "connects to big vessel on the right". The direction of the
% blood flow is not discernible at high velocities when scanning though, so
% the type of vessel needs to be determined ad hoc for certain cases. This
% function will show the kymograph, and thus the direction of the vessel,
% and the comments associated with the acquisition. This way, you can tell
% what the type of vessel is. 
% 'Vessel' is already in AcqInfos, but likely doesnt make sense (did not
% change it during acquisition). Will save new variable in AcqInfos named
% 'Vessel_type'

function Vessel_Type(DataFolder, overwrite)

% set up
if ~strcmp(DataFolder(end), filesep)
    DataFolder = [DataFolder filesep];
end

load([DataFolder 'AcqInfos.mat'], 'AcqInfoStream')

if isfield(AcqInfoStream, 'Vessel_type') && overwrite == 0
    disp('Vessel type already determined. Function exited.')
    return
elseif isfield(AcqInfoStream, 'Vessel_type') && overwrite == 1
    disp('OVERWRITING VESSEL TYPE')
end


% Check if coupled acquisitions are already categorized
if overwrite == 0
    cpld_acqs = AcqInfoStream.Coupled_acq;
    clear AcqInfoStream
    seps = strfind(DataFolder, filesep);
    MouseFolder = DataFolder(1:seps(end-1));
    for ind_coupled = 1:size(cpld_acqs,1)
        load([MouseFolder cpld_acqs{ind_coupled} filesep 'AcqInfos.mat'], 'AcqInfoStream')

        if isfield(AcqInfoStream, 'Vessel_type') % if you find a vessel type for a coupled acquisition, youre done
            vesseltype = AcqInfoStream.Vessel_type;
            disp(['Vessel type found in coupled acquisition. Type was ' vesseltype])

            load([DataFolder 'AcqInfos.mat'], 'AcqInfoStream')
            AcqInfoStream.Vessel_type = vesseltype;
            save([DataFolder 'AcqInfos.mat'], 'AcqInfoStream', '-append')
            return
        end

        clear AcqInfoStream seps
    end
end



% Show kymograph of this acq, of coupled acq, and comments. 

    ROI_list = dir([DataFolder 'kymoROI*.mat']);
    ROI_list = struct2cell(ROI_list);


promptstring = {'Comments' AcqInfoStream.comments};

vessel_type = listdlg('PromptString', promptstring,'ListString', {'Artery', 'Vein', 'Capillary', 'Unknown'});




end