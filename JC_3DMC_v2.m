% check SPM version

%
%       CREATED BY:     JASON CRAGGS
%       CREATED ON:     2017_09_29
%
%       USAGE:          GATHER INFORMATION ABOUT MOTION CORRECTION ESTIMATES
%                       AND WRITE OUT THE RESULTS

%		MODIFIED BY:	JASON CRAGGS
%		MODIFIED ON:	2017_10_04
%		MODIFIED ON:	2017_10_05
%       MODIFIED ON:    2017_10_07
%
%%       DETAILS OF WHAT THE SCRIPT DOES
%       SPECIFY LIBRARIES TO BE USED - IS THIS NECESSARY (?)
%       MOTION PARAMETERS ARE STORED IN THE rp*.txt files
%       SPECIFY THE BASE DIRECTORY OF DATA TO BE PROCESSED
%       LOAD THE TXT FILES CONTAINING THE 3DMC ESTIMATES
%           *   THE ESTIMATES IN THE LAST 3 COLUMNS ARE IN DEGREES
%           *   CONVERTING RADIANS INTO DEGREES
%       SPECIFY THE THRESHOLD FOR OUTLIER DETECTION
%       OPEN AN OUTPUT FILE (BEFORE THE LOOP)
%       LOOP OVER SUBJECTS
%           *   COMPUTE FRAME-WISE DISPLACEMENT
%              	*   A MORE SENSITIVE INDEX OF MOVEMENT
%            	*   IDENTIFY OUTLIERS
%           	*   WRITE EACH ITERATION TO THE OUTLIERS OUTPUT TXT FILE
%           *   COMPUTE MAXIMUM AMOUNT OF MOTION
%           	*   WRITE EACH ITERATION TO THE MAX MOTION OUTPUT TXT FILE
%           *   COMPUTE MAXIMUM AMOUT OF FRAME-WISE MOTION
%           	*   WRITE EACH ITERATION TO THE FRAMEWISE MOTION OUTLIERS OUTPUT TXT FILE
%       *   END THE LOOP
%       CLOST THE OUTPUT FILES
%       STORE THE INFORMATION AS A MAT FILE
%
%%


%       configure root path and subject pattern, as well as file patterns

rootpath = '/Volumes/Data/Imaging/R01/preprocessed/_Jason_0/';
subpattern = 'Sub*_v*';
rppattern = 'rp_aRSrun*.txt';

% get subjects folders
cd(rootpath);
subfolders = dir(subpattern);
subfolders(~cat(1, subfolders.isdir)) = [];

%       thresholds FOR OUTLIER DETECTION
transthresh = 1;
rotthresh = 1;

%       iterate over subject folders
rpfiles = {subfolders.name};
rpfolders = rpfiles;
for fc = 1:numel(rpfiles)

    %   locate all rp_*.txt files
    rpfiles{fc} = dir([subfolders(fc).name '/' rppattern]);
    rpfolders{fc} = repmat({subfolders(fc).name}, numel(rpfiles{fc}), 1);
end
rpfiles = cat(1, rpfiles{:});
rpfolders = cat(1, rpfolders{:});

%   DECLARE ADDITIONAL VARIABLES
rpfiles = {rpfiles.name};
rpoutliers = cell(numel(rpfiles), 1);
rpMax = cell(numel(rpfiles), 1);
maxFwMot = cell(numel(rpfiles), 1);

%   open an output file
file_id1 = fopen([rootpath 'motion_outLiers.txt'], 'w');
file_id2 = fopen([rootpath 'max_Motion.txt'], 'w');
file_id3 = fopen([rootpath 'maxFw_Motion.txt'], 'w');
file_id4 = fopen([rootpath 'Motion_estimates.txt'], 'w');
formatSpec = '%-32s: (maxMot) %.2f %.2f %.2f %.2f %.2f %.2f (maxFW) %.2f %.2f %.2f %.2f %.2f %.2f \n';

%   iterate over all files
for fc = 1:numel(rpfiles)
    
    %   generate filename with path
    rpfiles{fc} = [rpfolders{fc} '/' rpfiles{fc}];
    
    %   load RP file (into variable rp)
    rp = load(rpfiles{fc});
    
     %  convert radian into degrees (multiplying columns 4:6 with 180/pi)
    rp(:, 4:6) = (180/pi) .* rp(:, 4:6);
    
    %   compute frame-wise (volume-to-volume) displacement
    rpd = diff(rp);
    
    %   compute additional motion indicies
    %   *   MAX MOTION & MAX FRAMEWISE
    rpMax{fc} = (max(rp));
    maxFwMot{fc} = (max(rpd));
    
    %   find outliers (1+ as it is the first derivative) > thresholds
    rpoutliers{fc} = 1 + find(any(abs(rpd(:, 1:3)) > transthresh, 2) | ...
        any(abs(rpd(:, 4:6)) > rotthresh, 2));

    %   print out information
    olist = sprintf('%d, ', rpoutliers{fc});
    olist(end-1:end) = [];
    
    fprintf('%-32s: %d ([%s])\n', rpfiles{fc}, numel(rpoutliers{fc}), olist);
    fprintf(file_id1,'%-32s: %d ([%s])\n', rpfiles{fc}, numel(rpoutliers{fc}), olist);
    fprintf(file_id2,'%-32s(maxMot): %.2f %.2f %.2f %.2f %.2f %.2f \n', rpfiles{fc}, rpMax{fc});
    fprintf(file_id3,'%-32s(maxFwMot): %.2f %.2f %.2f %.2f %.2f %.2f \n', rpfiles{fc}, maxFwMot{fc});
    fprintf(file_id4,formatSpec,rpfiles{fc},rpMax{fc},maxFwMot{fc});
        
end

%fprintf(file_id4,'%.4f %.4f %.4f %.4f %.4f %.4f \n', rpd);

%   close file
fclose(file_id1);
fclose(file_id2);
fclose(file_id3);
fclose(file_id4);

%   TRANSPOSE THE SUBJECT LIST
rpSubs = rpfiles';

%   also store as a MAT file
save([rootpath 'motion_indices_v1.mat'], 'rpSubs', 'rpoutliers', 'rpMax', 'maxFwMot');


%%
%           END SCRIPT
%%