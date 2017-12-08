%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%		CREATED BY:		JASON CRAGGS
%
%		CREATED ON:		2017-12-07
%
%
%		USAGE:			SETUP FOR THE JC_FCcalc SCRIPT
%               IDENTIFIES THE FUNCTIONAL AND ROI FILES AND PASSES THEM OFF
%               FOR THE FUNCTIONAL CONNECTIVITY ANALYSES
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%           SETUP FOR A SINGLE SUBJECT - RUN FROM THEIR FOLDER

n = neuroelf;
%fcswa = n.findfiles([pwd '/Sub*'], 'fcsw*.nii', '-d1')
fcswa = n.findfiles([pwd], 'fcsw*.nii', '-d1');
hs = n.findfiles(pwd, 'r0*.nii', '-d1');

%         JOCHEN EXAMPLE FROM MATLAB HISTORY
% fcswa = n.findfiles([pwd '/Sub*'], 'fcsw*.nii', '-d1')
% hs = n.findfiles(pwd, 'r*.nii', '-d1');





JC_FCcalc(fcswa, hs)
%
%
%
%
%
