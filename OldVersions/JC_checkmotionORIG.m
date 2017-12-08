% read the motion parameter files of all subjects and print out
% (number of) outliers as well as overall motion estimates

% folders
rootpath = '/Volumes/Data/Imaging/R01/preprocessed/_Jason_0/';
subpattern = 'Sub*_v*';
rppattern = 'rp_aRSrun*.txt';

% thresholds
transthresh = 1;
rotthresh = 1;

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
file_id = fopen([rootpath 'motion_outliers.txt'], 'w');

%   iterate over all files
rpfiles = {rpfiles.name};
rpoutliers = cell(numel(rpfiles), 1);
for fc = 1:numel(rpfiles)
    
    % generate filename with path
    rpfiles{fc} = [rpfolders{fc} '/' rpfiles{fc}];
    
    % load RP file (into variable rp)
    rp = load(rpfiles{fc});
    
    % convert radian into degrees (multiplying columns 4:6 with 180/pi)
    rp(:, 4:6) = (180/pi) .* rp(:, 4:6);
    
    % compute frame-wise (volume-to-volume) displacement
    rpd = diff(rp);
    
    % find outliers (1+ as it is the first derivative) > thresholds
    rpoutliers{fc} = 1 + find(any(abs(rpd(:, 1:3)) > transthresh, 2) | ...
        any(abs(rpd(:, 4:6)) > rotthresh, 2));

%     % print out information
    olist = sprintf('%d, ', rpoutliers{fc});
    olist(end-1:end) = [];
        fprintf('%-32s: %d ([%s])\n', rpfiles{fc}, numel(rpoutliers{fc}), olist);
    fprintf(file_id, '%-32s: %d ([%s])\n', rpfiles{fc}, numel(rpoutliers{fc}), olist);
end

% close file
fclose(file_id);

% also store as a MAT file
save([rootpath 'motion_outliers.mat'], 'rpfiles', 'rpoutliers');

