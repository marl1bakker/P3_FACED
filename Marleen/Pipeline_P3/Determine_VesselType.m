function Determine_VesselType(Results, Overview)

if ~any('VesselType' ==  string(Results.Properties.VariableNames))
    Results.VesselType = categorical(repmat({''}, size(Results,1),1));
end

% Results = Results(Results.UseAcq == 'Yes', :);

% go per acquisition in results
for ind = 1:size(Results,1)

    % skip if it's not useable
    if Results.UseAcq(ind) == 'No'
        continue
     
    % skip if already done
    elseif ~isundefined(Results.VesselType(ind))
        continue

    % if diameter is <7 it's a capillary
    elseif Results.Diameter < 7
        Results.Vesseltype(ind) = 'Capillary';
        continue
    end

    %% determine vessel type
    % get right datafolder
    temptable = Overview(matches(Overview.Mouse, Results.Mouse(ind)),:);
    temptable = mousetable(matches(mousetable.Acq, Results.Acq(ind)),:);
    DataFolder = temptable.DataFolder{1};
    if ~strcmp(DataFolder(end), filesep)
        DataFolder = [DataFolder filesep];
    end
    seps = strfind(DataFolder, filesep);
    StemFolder = DataFolder(1:seps(end-1));
    clear temptable seps

    load([DataFolder 'AcqInfos.mat'], 'AcqInfoStream')

    % show kymograph overviews
    f1 = openfig([DataFolder 'Kymograph_overview.fig']);
    for ind_coupledacq = 1:length(AcqInfoStream.Coupled_acq)
        openfig([StemFolder AcqInfoStream.Coupled_acq{ind_coupledacq} filesep 'Kymograph_overview.fig']);
    end

    % show comments of current and coupled acqs
    b = figure;plot(1:10)

    text(1,9, 'Backwards');
    


end