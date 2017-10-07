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
%		MODIFIED ON:	2017_10_07
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
%               	*   A MORE SENSITIVE INDEX OF MOVEMENT
%            	*   IDENTIFY OUTLIERS
%           		*   WRITE EACH ITERATION TO THE OUTLIERS OUTPUT FILE
%           *   COMPUTE MAXIMUM AMOUNT OF MOTION
%           		*   WRITE EACH ITERATION TO THE MAX MOTION OUTPUT FILE
%           *   COMPUTE MAXIMUM AMOUT OF FRAME-WISE MOTION
%           		*   WRITE EACH ITERATION TO THE FRAMEWISE MOTION OUTLIERS OUTPUT FILE
%           *   END THE LOOP
%       CLOST THE OUTPUT FILE
%
%%

%   MOTION PARAMETERS ARE STORED IN THE rp*.txt files
%   THE FIRST STEP IS TO LOAD THE TXT FILES
%   BE SURE TO CONVERTING RADIANS INTO DEGREES


% % libraries
% n = neuroelf;
% x = xff;
%           SPECIFY VARIABLES AND PARAMETERS
% configure root path and subject pattern, as well as file patterns
%rootpath = '/cluster/folder/craggs/study/preprocessed/';
%rootpath = '/Volumes/Data/Imaging/R01/preprocessed/';
%rootpath = '/Volumes/Data/Imaging/R01/preprocessed/_Jason_0/';
rootpath = '/Volumes/Data/Imaging/R01/preprocessed/_Jason_0/';
subpattern = 'Sub*_v*';
rppattern = 'rp_aRSrun*.txt';
%anatpattern = 'T1_*.nii';
%funcpattern = 'RSrun*.nii';
%numruns = 4;


% thresholds FOR OUTLIER DETECTION
transthresh = 1;
rotthresh = 1;


% find subjects in root folder
% dirinfo = dir([rootpath subpattern]);
% subjlist = {dirinfo.name};

% get subjects folders
cd(rootpath);
subfolders = dir(subpattern);
subfolders(~cat(1, subfolders.isdir)) = [];

% iterate over subject folders
rpfiles = {subfolders.name};
rpfolders = rpfiles;
for fc = 1:numel(rpfiles)

    % locate all rp_*.txt files
    rpfiles{fc} = dir([subfolders(fc).name '/' rppattern]);
    rpfolders{fc} = repmat({subfolders(fc).name}, numel(rpfiles{fc}), 1);
end
rpfiles = cat(1, rpfiles{:});
rpfolders = cat(1, rpfolders{:});

% open an output file
file_id1 = fopen([rootpath 'motion_outLiers.txt'], 'w');
file_id2 = fopen([rootpath 'max_Motion.txt'], 'w');
file_id3 = fopen([rootpath 'maxFw_Motion.txt'], 'w');
formatSpec = '%.2f %.2f %.2f %.2f %.2f %.2f \n';

%   DECLARE VARIABLES
rpfiles = {rpfiles.name};
rpoutliers = cell(numel(rpfiles), 1);
rpMax = cell(numel(rpfiles), 1);
maxFwMot = cell(numel(rpfiles), 1);
C = {};

%   iterate over all files
for fc = 1:numel(rpfiles)

    % generate filename with path
    rpfiles{fc} = [rpfolders{fc} '/' rpfiles{fc}];

    % load RP file (into variable rp)
    rp = load(rpfiles{fc});

    % convert radian into degrees (multiplying columns 4:6 with 180/pi)
    rp(:, 4:6) = (180/pi) .* rp(:, 4:6);

    % compute frame-wise (volume-to-volume) displacement
    rpd = diff(rp);

    %   compute additional motion indicies
    %   *   MAX MOTION & MAX FRAMEWISE
    rpMax{fc} = (max(rp));
    maxFwMot{fc} = (max(rpd));

    % find outliers (1+ as it is the first derivative) > thresholds
    rpoutliers{fc} = 1 + find(any(abs(rpd(:, 1:3)) > transthresh, 2) | ...
        any(abs(rpd(:, 4:6)) > rotthresh, 2));

%     % print out information
    %   olist = sprintf('%d, ', rpoutliers{fc});
    olist = sprintf('%d, ', rpoutliers{fc});
    olist(end-1:end) = [];
    %fprintf('%-32s: %d ([%s])\n', rpfiles{fc}, numel(rpoutliers{fc}), olist);

    %maxList = sprintf('%.2f %.2f %.2f %.2f %.2f %.2f\n', rpMax{fc});
    maxList = sprintf('%d, ', rpMax{fc});
    %fprintf('%.2f %.2f %.2f %.2f %.2f %.2f\n', rpMax{fc});
    fprintf(file_id1,'%-32s: %d ([%s])\n', rpfiles{fc}, numel(rpoutliers{fc}), olist);
    fprintf(file_id2,'%-32s(maxMot): %.2f %.2f %.2f %.2f %.2f %.2f \n', rpfiles{fc}, rpMax{fc});
    fprintf(file_id3,'%-32s(maxFwMot): %.2f %.2f %.2f %.2f %.2f %.2f \n', rpfiles{fc}, maxFwMot{fc});

    C = [rpMax];

    %fprintf('%-32s: %d ([%s])\n', rpfiles{fc}, numel(rpoutliers{fc}), olist);
   % fprintf(file_id,'%-32s: %d ([%s])\n', rpfiles{fc}, numel(rpoutliers{fc}), olist);
    %fprintf('%.2f %.2f %.2f %.2f %.2f %.2f\n', rpMax);
 %   fprintf('%.2f %.2f %.2f %.2f %.2f %.2f\n', maxMot);


end

%
% [nrows,ncols] = size(C);
% for row = 1:nrows
%     fprintf(file_id2,formatSpec,C{row,:});
% end



% T = cell2table(C);
%writetable(T,'tabledata.dat')


% close file
fclose(file_id1);
fclose(file_id2);
fclose(file_id3);

% also store as a MAT file
save([rootpath 'motion_indices_v1.mat'], 'rpfiles', 'rpoutliers');





%%
%       THIS SECTION IS FOR WRITING THE RESULTS TO A FILE
%
%save('test1.txt','rpfolders','rpoutliers', '-ascii');
% A = rpfiles;
% A = A';
%
% T = table(A);
% T2 = table(rpfolders);
% T3 = table(rpoutliers);
% T4 = cell2table(rpoutliers(:,:));

% fileID = fopen('test3.dat','w');
% formatSpec = '%s %d %2.1f %s\n';
% save('test3.dat','rpfolders','rpoutliers');
% fclose(fileID);


% fileID = fopen('testData.dat','w');
% formatSpec = '%s %d %2.1f %s\n';
% [nrows,ncols] = size(rpfiles);
% for row = 1:nrows
%     fprintf(fileID,formatSpec,rpfiles{row,:});
% end
% fclose(fileID);
%
%type testData.dat

%   fileID = fopen('text.txt','w');
%   fprintf(fileID,'%-32s: %d ([%s])\n', rpfiles{fc}, numel(rpoutliers{fc}), olist);
%   fclose(fileID);



%%





% pick subject according to job number
% for sc = 1:numel(subjlist)
%
%     % set primary path
%     primary_path = [rootpath subjlist{sc} filesep];
%     cd(primary_path);
% end

%rpfiles = n.findfiles([pwd '/Sub*'], 'rp*.txt', '-d1');
%rpfiles = n.findfiles(['rp*.txt'], '-d1');

% % cd(rootpath);
% %
% % rpfiles = n.findfiles(['rp*.txt']);
% %
% % %   MOTION PARAMETERS ARE STORED IN THE rp*.txt files
% % %   THE FIRST STEP IS TO LOAD THE TXT FILES
% % %   BE SURE TO CONVERTING RADIANS INTO DEGREES
% % %
% % %rp=load('rp_aRSrun1_v1.txt');
% %
% % rp=load(rpfiles{1,1});
% % rp2 = rp;
% %
% %
% % % convert radian into degrees VERSION 1
% % mot = (rp * diag([1, 1, 1, 180 ./ [pi, pi, pi]]));
% % % convert radian into degrees VERSION 2 (multiplying columns 4:6 with 180/pi)
% % rp2(:, 4:6) = (180/pi) .* rp2(:, 4:6);
% %
% %
% % %figure(1), plot(mot, 'DisplayName', 'mot');
% %
% % %   FRAME BY FRAME MOTION ESTIMATES VERSION 1
% % fbf_motion = abs(diff(rp * diag([1, 1, 1, 180 ./ [pi, pi, pi]])));
% %
% % % compute frame-wise (volume-to-volume) displacement VERSION 2
% % rpd = diff(rp2);
% %
% % % compute frame-wise (volume-to-volume) displacement VERSION 2
% % rpd2 = abs(diff(rp2));
% %
% % %total_motion = sum(fbf_motion);
% % max_motion = max(mot);
% % max_fbfmotion = max(fbf_motion);
% %
% % %%
% %
% % %   iterate over all files
% % %rpfiles = {rpfiles.name};
% % rpoutliers = cell(numel(rpfiles), 1);
% % %for fc = 1:numel(rpfiles)
% %
% %     % generate filename with path
% %  %   rpfiles{fc} = [rpfolders{fc} '/' rpfiles{fc}];
% %
% % % %     % load RP file (into variable rp)
% % % %     rp = load(rpfiles{fc});
% % % %
% % % %     % convert radian into degrees (multiplying columns 4:6 with 180/pi)
% % % %     rp(:, 4:6) = (180/pi) .* rp(:, 4:6);
% % % %
% % % %     % compute frame-wise (volume-to-volume) displacement
% % % %     rpd = diff(rp);
% % % %
% %     % find outliers (1+ as it is the first derivative) > thresholds
% % %    rpoutliers{fc} = 1 + find(any(abs(rpd(:, 1:3)) > transthresh, 2) | ...
% %     rpoutliers = 1 + find(any(abs(rpd(:, 1:3)) > transthresh, 2) | ...
% %         any(abs(rpd(:, 4:6)) > rotthresh, 2));

%     % print out information
%     olist = sprintf('%d, ', rpoutliers);
%     olist(end-1:end) = [];
%     fprintf('%-32s: %d ([%s])\n', rpfiles, numel(rpoutliers), olist);
% %
%
% %     % print out information
%     olist = sprintf('%d, ', rpoutliers{fc});
%     olist(end-1:end) = [];
%     fprintf('%-32s: %d ([%s])\n', rpfiles{fc}, numel(rpoutliers{fc}), olist);
%end

%%

%%

% find outliers (1+ as it is the first derivative) > thresholds
% rpoutliers{fc} = 1 + find(any(abs(rpd(:, 1:3)) > transthresh, 2) | ...
%     any(abs(rpd(:, 4:6)) > rotthresh, 2));

%     % print out information
% olist = sprintf('%d, ', rpoutliers{fc});
% olist(end-1:end) = [];
% fprintf('%-32s: %d ([%s])\n', rpfiles{fc}, numel(rpoutliers{fc}), olist);


% save('text2.txt','fbf_motion', '-ascii')
% save('text3.txt','fbf_motion')
% save('text3.txt','fbf_motion') -ascii
% save('text3.txt','fbf_motion', '-ascii')
% save('text3.txt','fbf_motion', 'maxFbf', '-ascii')
% type('text3.txt')
% save('text3.txt', 'maxFbf', '-ascii')

%%






%figure(2), plot(fbf_motion * diag([1, 1, 1, 180 ./ [pi, pi, pi]]));
%figure(2), plot(fbf_motion);

% cd ../

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

% title('frame by frame displacement', 'FontWeight','bold','FontSize',16);



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
