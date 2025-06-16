% to delete faced files, to save room. Could be done after Clean Data
% theoretically

function Delete_faced_files(Overview)

for ind = 1:size(Overview, 1)
    datafolder = Overview.DataFolder{ind};
    seps = strfind(datafolder, filesep);

    if ~exist([datafolder 'Kymograph_overview.fig'], 'file')
        disp(['No kymograph overview for ' datafolder(seps(end-2)+1:end-1)]);
        continue
    elseif ~exist([datafolder 'faced.dat'], 'file')
        disp(['Faced.dat not made or already deleted for ' datafolder(seps(end-2)+1:end-1)]);
        continue
    end
    
    f1 = openfig([datafolder 'Kymograph_overview.fig']);
    f1.Position = [1.8000   49.8000  766.4000  828.8000];
    answer = questdlg('Can faced.dat be removed?');
    close(f1)
    
    if matches(answer, 'Yes')
        delete([datafolder 'faced.dat'])
        disp(['Faced.dat deleted for ' datafolder(seps(end-2)+1:end-1)]);
    elseif matches(answer, 'No')
        disp(['Not deleted ' datafolder(seps(end-2)+1:end-1)])
        continue
    elseif matches(answer, 'Cancel')
        error('Pipeline canceled.')
    end
    
end

end