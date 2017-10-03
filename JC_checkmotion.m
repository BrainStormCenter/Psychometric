% read the motion parameter files of all subjects and print out
% (number of) outliers as well as overall motion estimates

% folders
rootpath = '/Volumes/Data/Imaging/R01/preprocessed/_Jason/';
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
for fc = 1:numel(rpfiles)

    % locate all rp_*.txt files
    rpfiles{fc} = dir([subfolders(fc).name '/' rppattern]);
end
rpfiles = cat(1, rpfiles{:});

%   iterate over all files
%rpfolders = {rpfiles.folder};
rpfolders = {subfolders.name};
rpfiles = {rpfiles.name};
rpoutliers = cell(numel(rpfiles), 1);
for fc = 1:numel(rpfiles)
    rpfiles{fc} = [rpfolders{fc} '/' rpfiles{fc}];
    rp = load(rpfiles{fc});
    rp(:, 4:6) = (180/pi) .* rp(:, 4:6);
    rpd = diff(rp);
    rpoutliers{fc} = 1 + find(any(abs(rpd(:, 1:3)) > transthresh, 2) | ...
        any(abs(rpd(:, 4:6)) > rotthresh, 2));

%     % print out information
    olist = sprintf('%d, ', rpoutliers{fc});
    olist(end-1:end) = [];
    fprintf('%-88s: %d ([%s])\n', rpfiles{fc}, numel(rpoutliers{fc}), olist);
end
