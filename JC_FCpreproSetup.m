%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%		CREATED BY:		JASON CRAGGS
%		CREATED ON:		2017-12-07
%
%		USAGE:			PREPARING THE FUNCTIONAL DATASETS FOR THE
%                       FUNCTIONAL CONNECTIVITY ANALYSES
%
% 	MODIFIED ON:	  2017_12_08
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           SETUP FOR A SINGLE SUBJECT
%           THIS MUST BE RUN FROM THEIR FOLDER
%
%           rpfiles = rp-*.txt files = motion parameters
%                     Columns 1-3 are X/Y/Z movement in mm.
%                     Columns 4-6 are rotation around X/Y/Z in radians.
%           wc* = TISSUE CLASSES FROM THE STRUCTURAL DATA
%                   c1/c2/c3 images are the first three "tissue classes"
%                   (GM/WM/CSF) from the segmentation
%                   r = realigned
%                   w = normalized
%           swa = THE FUNCTIONAL DATA FROM EACH RUN
%                   a = slice time corrected
%                   s = smoothed
%                   w = normalized
%
%              OUTPUT FILES
%                   The function will make a copy of (each of) the input file(s)
%                   with a prefix of "fc" (e.g. "fcswraRUN.nii").
%
n = neuroelf;
rpfiles = n.findfiles([pwd], 'rp*.txt', '-d1');
wc1 = n.findfiles([pwd], 'wc1*.nii', '-d1');
wc2 = n.findfiles([pwd], 'wc2*.nii', '-d1');
wc3 = n.findfiles([pwd], 'wc3*.nii', '-d1');
swa = n.findfiles([pwd], 'sw*.nii', '-d1');
%
%           THE WC FILES NEED TO BE REFORMATTED TO MAINTAIN THE SAME
%           DIMENSIONS AS THE ARRAY IN THE RP FILE
wc1 = repmat(wc1', 4, 1);
wc2 = repmat(wc2', 4, 1);
wc3 = repmat(wc3', 4, 1);
%
%           RUNNING THE JC_FCprepro SCRIPT USING THE ABOVE VARIABLES
 for c=1:4,  %   THE NUMBER OF FUNCTIONAL RUNS
     JC_FCprepro(swa{c}, rpfiles{c}, 120/2.8, {wc1{c}; wc2{c}; wc3{c}});
end
%           SETUP FOR MULTIPLE SUBJECTS
% n = neuroelf;
% rpfiles = n.findfiles([pwd '/Sub004*'], 'rp*.txt', '-d1');
% wc1 = n.findfiles([pwd '/Sub*'], 'wc1*.nii', '-d1');
% wc2 = n.findfiles([pwd '/Sub*'], 'wc2*.nii', '-d1');
% wc3 = n.findfiles([pwd '/Sub*'], 'wc3*.nii', '-d1');
%
%
%
%
%
%         END OF SCRIPT
