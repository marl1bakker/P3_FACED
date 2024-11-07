% Alteration can be 'Add' or 'Change'. Always give new_info between ' '.
% Field can be an existing field or a new one
% Changes the AcqInfos.mat. To chagne the overview after that, do
% UpdateOverview. 

function AddInfo(DataFolder, Field, Alteration, new_info)

if ~strcmp(DataFolder(end), filesep)
    DataFolder = [DataFolder filesep];
end

load([DataFolder 'AcqInfos.mat'], 'AcqInfoStream')

if ~exist('Field', 'var')
    Field = 'Comments';
end

if ~isfield(AcqInfoStream, Field)
    disp([Field ' did not exist, made field.'])
    Alteration = 'NewField';
    AcqInfoStream.(Field) = 'empty';
end

if ~exist('Alteration', 'var')
    Alteration = 'Add';
end

if ~exist('new_info', 'var')
    new_info = inputdlg([Alteration ' ' Field ' ' DataFolder newline ...
        'Old line was: ' AcqInfoStream.(Field)]);
elseif ~matches(class(new_info), 'cell')
    new_info = {new_info};
end



switch Alteration
    case 'Add'
        if ~ischar(AcqInfoStream.(Field))
            error('Can only add to a char field.')
        end
        AcqInfoStream.(Field) = [AcqInfoStream.(Field) ' -- ' new_info{1}];

    case 'Change'
        datatype = class(AcqInfoStream.(Field));
        if matches(datatype, 'double')
            AcqInfoStream.(Field) = str2double(new_info{1});
        elseif matches(datatype, 'char')
            AcqInfoStream.(Field) = [new_info{1}];
        else
            error('Datatype not recognized.')
        end

    case 'NewField'
        AcqInfoStream.(Field) = [new_info{1}];
end

save([DataFolder 'AcqInfos.mat'], 'AcqInfoStream');


end


