%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%		CREATED BY:	JASON CRAGGS
%		CREATED ON:	2017-12-13
%    	MODIFIED ON:	 2018_02_16
%
%		USAGE:		CONVERT AND COMBINE .NII BRAIN MASKS INTO A
%                        SINGLE .VOI FILE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   February 16, 2018
%                   THE RESLICED AAL IMAGES MAY BE PROBLEMATIC
%                   I AM BUILDING A VOI WITH THE UN-RESLICED NII FILES TO SEE IF THIS
%                   CORRECTS THE PROBLEM
%

% NeuroElf library
n = neuroelf;

maskPath = '/Volumes/Data/Imaging/R01/preprocessed/Test_files/RegionMasks/';

% find ROI files
%    THESE ARE REGIONAL MASKS IN THE .nii FORMAT
%rois = n.findfiles(FOLDER_WITH_ROI_NIFTIS, 'r???_*.nii');
%rois = n.findfiles(maskPath, 'r???_*.nii');
%rois = n.findfiles(maskPath, 'r*.nii');
rois = n.findfiles(maskPath, 'MNI*.nii');

% create new VOI
voi = xff('new:voi');

% import each of the .nii files into the VOI file
for c=1:numel(rois)
    voi.ImportClusters(rois{c}, struct('sepclus', false, 'sort', 'value'));
end

% rename?
voiName  = 'AALtest2.voi';

voi.SaveAs(voiName);
