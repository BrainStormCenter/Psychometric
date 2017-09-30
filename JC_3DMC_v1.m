% check SPM version

%%
%       CREATED BY:     JASON CRAGGS
%       CREATED ON:     2017_09_29
%
%       USAGE:          GET INFORMATION ABOUT MOTION CORRECTION ESTIMATES
%                       AND COMPILE THEM INTO A SINGLE FILE
%



% configure root path and subject pattern, as well as file patterns
%rootpath = '/cluster/folder/craggs/study/preprocessed/';
%rootpath = '/Volumes/Data/Imaging/R01/preprocessed/';
rootpath = '/Volumes/Data/Imaging/R01/preprocessed/_Jason_0/';
subpattern = 'Sub*_v*';
anatpattern = 'T1_*.nii';
funcpattern = 'RSrun*.nii';
numruns = 4;

% set variables, number of volumes and functional slices, TR
%nvols = 125;
%nslices = 42;
%TR = 2.8;

swa = n.findfiles([pwd '/Sub*'], 'sw*.nii', '-d1');
%
rpfiles = n.findfiles([pwd '/Sub*'], 'rp*.txt', '-d1');
%
wc1 = n.findfiles([pwd '/Sub*'], 'wc1*.nii', '-d1');
wc2 = n.findfiles([pwd '/Sub*'], 'wc2*.nii', '-d1');
wc3 = n.findfiles([pwd '/Sub*'], 'wc3*.nii', '-d1');
JC_FCprepro(swa{c}, rpfiles{c}, 120/2.8, {wc1{c}; wc2{c}; wc3{c}})
