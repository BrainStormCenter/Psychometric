% check SPM version

%%
%       CREATED BY:     JASON CRAGGS
%       CREATED ON:     2017_09_29
%
%       USAGE:          GET INFORMATION ABOUT MOTION CORRECTION ESTIMATES
%                       AND COMPILE THEM INTO A SINGLE FILE

%		MODIFIED BY:	JASON CRAGGS
%		MODIFIED ON:	2017_10_04
%		REASON:			NEED TO WRITE OUT RESULTS TO FILE


%

% libraries
n = neuroelf;
x = xff;

% configure root path and subject pattern, as well as file patterns
%rootpath = '/cluster/folder/craggs/study/preprocessed/';
%rootpath = '/Volumes/Data/Imaging/R01/preprocessed/';
rootpath = '/Volumes/Data/Imaging/R01/preprocessed/_Jason_0/';
subpattern = 'Sub*_v*';
%anatpattern = 'T1_*.nii';
funcpattern = 'RSrun*.nii';
%numruns = 4;


% find subjects in root folder
dirinfo = dir([rootpath subpattern]);
subjlist = {dirinfo.name};

% pick subject according to job number
for sc = 1:numel(subjlist)

    % set primary path
    primary_path = [rootpath subjlist{sc} filesep];
    cd(primary_path);
end

%rpfiles = n.findfiles([pwd '/Sub*'], 'rp*.txt', '-d1');
%rpfiles = n.findfiles(['rp*.txt'], '-d1');
rpfiles = n.findfiles(['rp*.txt']);

%   MOTION PARAMETERS ARE STORED IN THE rp*.txt files
%   THE FIRST STEP IS TO LOAD THE TXT FILES
%   THEN PLOT THEM, BEING SURE TO CONVERTING RADIANS INTO
%   DEGREES TO EASE VISUAL INSPECTION
%rp=load('rp_aRSrun1_v1.txt');
rp=load(rpfiles{1,1});
mot = (rp * diag([1, 1, 1, 180 ./ [pi, pi, pi]]));
figure(1), plot(mot, 'DisplayName', 'mot');

%   FRAME BY FRAME MOTION ESTIMATES
fbf_motion = abs(diff(rp * diag([1, 1, 1, 180 ./ [pi, pi, pi]])));
%figure(2), plot(fbf_motion * diag([1, 1, 1, 180 ./ [pi, pi, pi]]));
figure(2), plot(fbf_motion);

cd ../

% plot(fbf_motion,'DisplayName','fbf_motion')
% total_motion = sum(fbf_motion);
% bar(total_motion)
% plot(rp,'DisplayName','rp')
% plot(rp * diag([1, 1, 1, 180 ./ [pi, pi, pi]]));
% figure, plot(rp * diag([1, 1, 1, 180 ./ [pi, pi, pi]]));
% bar(total_motion)
% total_motion = sum(fbf_motion);
% max_motion = max(fbf_motion);
% bar(max_motion)
% figure2, bar(max_motion)


% set variables, number of volumes and functional slices, TR
%nvols = 125;
%nslices = 42;
%TR = 2.8;
%
% swa = n.findfiles([pwd '/Sub*'], 'sw*.nii', '-d1');
%
% rpfiles = n.findfiles([pwd '/Sub*'], 'rp*.txt', '-d1');
% %
% % wc1 = n.findfiles([pwd '/Sub*'], 'wc1*.nii', '-d1');
% % wc2 = n.findfiles([pwd '/Sub*'], 'wc2*.nii', '-d1');
% % wc3 = n.findfiles([pwd '/Sub*'], 'wc3*.nii', '-d1');
% %JC_FCprepro(swa{c}, rpfiles{c}, 120/2.8, {wc1{c}; wc2{c}; wc3{c}})
% %
% %
% % THIS COMMAND WORKS FOR A SINGLE SUBJECT
% %JC_FCprepro(swa,rpfiles,120/2.8, [wc1; wc2; wc3])
