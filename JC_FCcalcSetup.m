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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%              SETUP FOR A SINGLE SUBJECT
%              THIS MUST BE RUN FROM THEIR FOLDER
%
%              fcswa = preprocessed/cleaned functional datasets
%              roifiles = masks of brain regions being studied
%
%
n = neuroelf;
%fcswa = n.findfiles([pwd '/Sub*'], 'fcsw*.nii', '-d1')
fcswa = n.findfiles([pwd], 'fcsw*.nii', '-d1');
roifiles = n.findfiles(pwd, 'r???_*.nii', '-d1');

%              PASS THE VARIABLES TO THE FUNCTION BELOW
JC_FCcalc(fcswa, roifiles)
%
%
%
%              END OF SCRIPT
%
