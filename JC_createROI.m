
% NeuroElf library
n = neuroelf;

% find ROI files
rois = n.findfiles(FOLDER_WITH_ROI_NIFTIS, 'r???_*.nii');

% create new VOI
voi = xff('new:voi');

% import each of the files into the VOI file
for c=1:numel(rois)
    nvoi.ImportClusters(rois{c}, struct('sepclus', false, 'sort', 'value'));
end

% rename?
