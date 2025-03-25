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
    if ~matches(Field, 'Coupled_acq')
        AcqInfoStream.(Field) = 'empty';
    end
elseif ~exist('Alteration', 'var')
    Alteration = 'Add';
end

%% coupled acq
% acq of same vessel
if matches(Field, 'Coupled_acq')
    % current_acq = AcqInfoStream.DatasetName;
    seps = strfind(DataFolder, filesep);

    % if you have one already
    if isfield(AcqInfoStream, 'Coupled_acq') && ~matches(Alteration, 'Change')
        add_acqs = questdlg(['Mouse: ' AcqInfoStream.Mouse '. '...
            'Found following coupled acquisitions for ' ...
            AcqInfoStream.DatasetName ': ' AcqInfoStream.Coupled_acq{:} ...
            ', do you want to add?']);
    elseif isfield(AcqInfoStream, 'Coupled_acq') && matches(Alteration, 'Change')
        change_acqs = questdlg(['Mouse: ' AcqInfoStream.Mouse '. '...
            'Found following coupled acquisitions for ' ...
            AcqInfoStream.DatasetName ': ' AcqInfoStream.Coupled_acq{:} ...
            ', do you want to CHANGE? Old couplings will be deleted.']);
        if matches(change_acqs, 'Yes')
            AcqInfoStream = rmfield(AcqInfoStream,'Coupled_acq');
            add_acqs = 'Yes';
        end
    else
        add_acqs = 'Yes';
    end

    if matches(add_acqs, 'Yes')
        acq_list = struct2table(dir(DataFolder(1:seps(end-1))));
        acq_list = acq_list(startsWith(acq_list.name, 'A'),:);
        acq_list = acq_list(acq_list.isdir, :);

        if isfield(AcqInfoStream, 'Coupled_acq')
            coupled_acq = listdlg('PromptString', [{'Select coupled acquisition(s) ',...
                'of mouse: '} AcqInfoStream.Mouse {'Scan: '} AcqInfoStream.DatasetName ...
                {' Already coupled are: '} AcqInfoStream.Coupled_acq{:}],...
                'ListString', acq_list.name);
            AcqInfoStream.Coupled_acq = [AcqInfoStream.Coupled_acq; acq_list.name(coupled_acq)];
        else
            coupled_acq = listdlg('PromptString', [{'Select coupled acquisition(s) ',...
                'of mouse: '} AcqInfoStream.Mouse {'Scan: '} AcqInfoStream.DatasetName ...
                {' Already coupled are: None'}],...
                'ListString', acq_list.name);
            AcqInfoStream.Coupled_acq = [acq_list.name(coupled_acq)];
        end

        save([DataFolder 'AcqInfos.mat'], 'AcqInfoStream', '-append');
    end

    % check if the other acq(s) has coupled_acq, and includes this one
    coupled_acq_list = [AcqInfoStream.Coupled_acq; AcqInfoStream.DatasetName];
    clear AcqInfoStream add_acqs
    for indcoupled = 1:length(coupled_acq_list)-1 % last one you just did
        if ~exist([DataFolder(1:seps(end-1)) coupled_acq_list{indcoupled} filesep 'AcqInfos.mat'], 'file')
            continue
        end

        load([DataFolder(1:seps(end-1)) coupled_acq_list{indcoupled} filesep 'AcqInfos.mat'], 'AcqInfoStream');
        if ~isfield(AcqInfoStream, 'Coupled_acq')
            AcqInfoStream.Coupled_acq = coupled_acq_list(~matches(coupled_acq_list, AcqInfoStream.DatasetName),:);
            save([DataFolder(1:seps(end-1)) coupled_acq_list{indcoupled} filesep 'AcqInfos.mat'], 'AcqInfoStream', '-append');
        else
            list_for_acq = coupled_acq_list(~matches(coupled_acq_list, AcqInfoStream.DatasetName));
            add_to_coupled_list = list_for_acq(~contains(list_for_acq, AcqInfoStream.Coupled_acq),:);
            if isempty(add_to_coupled_list)
                continue
            end
            AcqInfoStream.Coupled_acq = [AcqInfoStream.Coupled_acq; add_to_coupled_list];
            save([DataFolder(1:seps(end-1)) coupled_acq_list{indcoupled} filesep 'AcqInfos.mat'], 'AcqInfoStream', '-append');
            clear list_for_acq add_to_coupled_list AcqInfoStream

        end
    end
    return
end


%% other info
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
        AcqInfoStream.(Field) = [AcqInfoStream.(Field) new_info{1}];

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

save([DataFolder 'AcqInfos.mat'], 'AcqInfoStream', '-append');


end


