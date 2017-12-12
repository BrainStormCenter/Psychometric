%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%		CREATED BY:		JASON CRAGGS
%		CREATED ON:		2017-12-07
%
%
%         USAGE: SETUP FOR THE JC_FCcalc SCRIPT
%              IDENTIFIES THE FUNCTIONAL AND ROI FILES AND PASSES THEM OFF
%              FOR THE FUNCTIONAL CONNECTIVITY ANALYSES
%
%              THE OUTPUT FILE IS A CORRELATION MAP REPRESENTING THE CORRELATION
%              OF THAT ROI WITH ALL OTHER VOXELS IN THE BRAIN
%
%         MODIFIED ON:	  2017_12_08
%    	MODIFIED ON:	 2017_12_12
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%              SETUP FOR BATCH PROCESSING
%              THIS MUST BE RUN FROM THE ROOT FOLDER
%
%              fcswa = preprocessed/cleaned functional datasets
%              roifiles = masks of brain regions being studied
%                        THE ROI FILES MUST BE IN EACH SUBJECT'S FOLDER 
%

n = neuroelf;
%
%rootpath = '/Volumes/Data/Imaging/R01/preprocessed/_Jason_0/';
rootpath = '/Volumes/Data/Imaging/R01/preprocessed/';
subpattern = 'Sub*_v*';
%roipath = '/Volumes/Data/Imaging/R01/preprocessed/_Jason_0/resliced_masks/'

% find subjects in root folder
dirinfo = dir([rootpath subpattern]);
subjlist = {dirinfo.name};
%
% pick subject according to job number
for sc = 1:numel(subjlist)

    % set primary path
    primary_path = [rootpath subjlist{sc} filesep];
    cd(primary_path);
    fcswa = n.findfiles([pwd], 'fcsw*.nii', '-d1')
    fcswa = n.findfiles([pwd], 'fcsw*.nii', '-d1');
    roifiles = n.findfiles([pwd], 'r???_*.nii', '-d1');
    %roifiles = n.findfiles(roipath, 'r???_*.nii', '-d1');

    % %              PASS THE VARIABLES TO THE FUNCTION BELOW
     JC_FCcalc(fcswa, roifiles)
end

% %fcswa = n.findfiles([pwd '/Sub*'], 'fcsw*.nii', '-d1')
% fcswa = n.findfiles([pwd], 'fcsw*.nii', '-d1');
% roifiles = n.findfiles(pwd, 'r???_*.nii', '-d1');
%
% %              PASS THE VARIABLES TO THE FUNCTION BELOW
% JC_FCcalc(fcswa, roifiles)
% %
%
%
%              END OF SCRIPT
%
