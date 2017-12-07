%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%		CREATED BY:		JASON CRAGGS
%
%		CREATED ON:		2017-12-07
%
%
%		USAGE:			CALCULATING THE FUNCTIONAL DATASETS FOR THE 
%                       FUNCTIONAL CONNECTIVITY ANALYSES
%                       
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%           SETUP FOR A SINGLE SUBJECT - RUN FROM THEIR FOLDER

n = neuroelf;
rpfiles = n.findfiles([pwd], 'rp*.txt', '-d1');
wc1 = n.findfiles([pwd], 'wc1*.nii', '-d1');
wc2 = n.findfiles([pwd], 'wc2*.nii', '-d1');
wc3 = n.findfiles([pwd], 'wc3*.nii', '-d1');
swa = n.findfiles([pwd], 'sw*.nii', '-d1');



JC_FCcalc(fcswa, hs)