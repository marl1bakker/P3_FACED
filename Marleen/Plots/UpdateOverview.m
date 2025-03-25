% load('D:\FACED\Data P3\Overview.mat', 'Overview')
% if you want to add a column, simply make sure you add it to overview and
% it will fill it out automatically if it's in acqinfostream. For example:
% Overview.ny = cellstr(repmat('x', size(Overview,1), 1);

function UpdateOverview(Overview, SaveDirectory)

% set up
if ~exist('SaveDirectory', 'var')
    SaveDirectory = 'D:\FACED\Data P3\';
elseif ~strcmp(SaveDirectory(end), filesep)
    SaveDirectory = [SaveDirectory filesep];
end

if ~exist('Overview', 'var')
    load([SaveDirectory 'Overview.mat'], 'Overview');
end

updated_fields = 0;
OverviewColumns = Overview.Properties.VariableNames;

% go per acquisition
for indacq = 1:size(Overview,1)

    % check if you already made acqinfos for this acquisition
    if exist([Overview.DataFolder{indacq} 'AcqInfos.mat'], 'file')
        load([Overview.DataFolder{indacq} 'AcqInfos.mat'], 'AcqInfoStream');
    else
        continue
    end

    % go per variable
    for indcolumn = 1:size(OverviewColumns, 2) 

        % if AcqInfoStream doesnt contains that information 
        if ~isfield(AcqInfoStream, OverviewColumns{indcolumn}) 
            continue
        end

        % All types are cell in overview. To compare if fields are already
        % the same and to update field, acqinfostream fields need to be
        % transformed to cell.
        
        current_field = AcqInfoStream.(OverviewColumns{indcolumn});
        type_of_var = class(current_field);
            
        switch type_of_var
            case {'double', 'single'}
                % check if it differs less than 1 (to make sure you dont
                % overwrite depth when unecessary)
                if str2double(Overview.(OverviewColumns{indcolumn})(indacq)) - current_field == 0
                    continue
                end
                current_field = cellstr(num2str(current_field));

            case 'char'
                current_field = cellstr(current_field);
        end
        
        if ~matches(Overview.(OverviewColumns{indcolumn})(indacq), current_field) % if not already the same
       
            % % TEMP for m01
            % if matches(OverviewColumns(indcolumn), {'FrameRateHz'})
            %     answer = questdlg([Overview.Acq{indacq} ' change FrameRateHz ' Overview.(OverviewColumns{indcolumn}){indacq} ' to ' current_field{1} '?']);
            % 
            %     if matches(answer, 'No')
            %         AddInfo(Overview.DataFolder{indacq}, 'FrameRateHz', 'Change', Overview.(OverviewColumns{indcolumn}){indacq});
            %         continue
            %     end
            % end

            Overview.(OverviewColumns{indcolumn})(indacq) = current_field;
            updated_fields = updated_fields+1;
        
        end
    end
end

disp(['Overview update succesfull. Number of updated fields is ' num2str(updated_fields)])

% save
save([SaveDirectory 'Overview.mat'], 'Overview', '-append');


end





