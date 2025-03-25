% DataFolder = pwd;
% DataFile = 'm04_heart.nii';

function [slice, lv_thickness] = Lightsheet_Heart_Analysis(DataFolder, filename)
% note: this is very rough. just takes the middle of the scan and draws a
% line on one of the axes. Think of more sophisticated way to do it...

if ~strcmp(DataFolder(end), filesep)
    DataFolder = [DataFolder filesep];
end

% load data
b = waitbar(0.1, 'Opening nifti file (might take a while)...');
V = niftiread([DataFolder filename]);
info = niftiinfo([DataFolder filename]);

% get slices in middle
b = waitbar(0.8, 'Getting slices of file');

axlist = {'x', 'y', 'z'};
szlist = [2,3; 1,3; 1,2];
pxlsz = info.PixelDimensions;
datsz = info.ImageSize;


slice.x = single(V(datsz(1)/2,:,:));
slice.y = single(V(:,datsz(2)/2,:));
slice.z = single(V(:,:,datsz(3)/2));

f1 = figure; tiledlayout(1, 3)
for axind = 1:length(axlist)

    slice.(axlist{axind}) = reshape(slice.(axlist{axind}), datsz(szlist(axind, 1)), datsz(szlist(axind,2)));

    nexttile; imagesc(slice.(axlist{axind})*-1); axis('image');
    yticks(0:1/pxlsz(szlist(axind,1))*1000:datsz(szlist(axind,1)))
    yticklabels(0:1:pxlsz(szlist(axind,1))*datsz(szlist(axind,1))/1000)
    ylabel('mm')
    xticks(0:1/pxlsz(szlist(axind,2))*1000:datsz(szlist(axind,2)))
    xticklabels(0:1:pxlsz(szlist(axind,2))*datsz(szlist(axind,2))/1000)
    xlabel('mm')
    title([axlist{axind} ' axis'])
end

colormap('gray')


% Make plot 
% choose the axis to draw roi on
% roiaxis = listdlg('ListString', axlist);
close(f1)

f1 = figure;
% plot figure, make roi
figure;imagesc(flipud(slice.z)*-1); colormap('gray'); axis('image');

roi_pixels = drawpolyline;
thickness_slice = slice.z(size(slice.z,1)-round(mean(roi_pixels.Position(:,2))),:);
figure;
plot(thickness_slice)
hold on

% get thickness
change_points = findchangepts(thickness_slice, 'MaxNumChanges', 4, 'Statistic', 'linear');
scatter(change_points, ones(1, 4));

% answer = questdlg('Does it make sense?');

lv_thickness.left_wall = (change_points(2)-change_points(1))*pxlsz(3)/1000; % in mm
lv_thickness.right_wall = (change_points(4)-change_points(3))*pxlsz(3)/1000;

lv_thickness.average_wall = mean(left_wall_lv, right_wall_lv);

% save
save([DataFolder 'Middle_slices.mat'], 'slice', 'lv_thickness')

end






% slice.x = single(V(info.ImageSize(1)/2,:,:));
% slice.x = reshape(slice.x, info.ImageSize(2), info.ImageSize(3));
% 
% 
% slice.y = single(V(:,info.ImageSize(2)/2,:));
% slice.y = reshape(slice.y, info.ImageSize(1), info.ImageSize(3));
% 
% 
% slice.z = single(V(:,:,info.ImageSize(3)/2));
% slice.z = reshape(slice.z, info.ImageSize(1), info.ImageSize(2));
% 
% close(b)
% 
% % plot slices
% f1 = figure; tiledlayout(1, 3)
% nexttile; imagesc(x_slice*-1); axis('image');
% yticks(0:1/info.PixelDimensions(2)*1000:size(x_slice,1))
% yticklabels(0:1:info.PixelDimensions(2)*size(x_slice,1)/1000)
% ylabel('mm')
% xticks(0:1/info.PixelDimensions(3)*1000:size(x_slice,1))
% xticklabels(0:1:info.PixelDimensions(3)*size(x_slice,1)/1000)
% xlabel('mm')
% title('x axis')
% 
% nexttile; imagesc(y_slice*-1); axis('image');
% yticks(0:1/info.PixelDimensions(1)*1000:size(y_slice,1))
% yticklabels(0:1:info.PixelDimensions(1)*size(y_slice,1)/1000)
% ylabel('mm')
% xticks(0:1/info.PixelDimensions(3)*1000:size(y_slice,1))
% xticklabels(0:1:info.PixelDimensions(3)*size(y_slice,1)/1000)
% xlabel('mm')
% title('y axis')
% 
% nexttile; imagesc(z_slice*-1); axis('image');
% yticks(0:1/info.PixelDimensions(1)*1000:size(z_slice,1))
% yticklabels(0:1:info.PixelDimensions(1)*size(z_slice,1)/1000)
% ylabel('mm')
% xticks(0:1/info.PixelDimensions(2)*1000:size(z_slice,1))
% xticklabels(0:1:info.PixelDimensions(2)*size(z_slice,1)/1000)
% xlabel('mm')
% title('z axis')
% 
% colormap('gray')