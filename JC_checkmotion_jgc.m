% read the motion parameter files of all subjects and print out
% (number of) outliers as well as overall motion estimates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%		CREATED BY:		JOCHEN WEBER
%		CREATED ON:		2017-10-02
%		USAGE:			EXTRACT MOTION CORRECTION INFORMATION
%
%		MODIFIED BY:	JASON CRAGGS
%		MODIFIED ON:	2017_10_03
%		REASON:			FOLDER-FILE MISMATCH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%       CLEAR PREVIOUS ITERATION
clear all;


% folders
% configure root path and subject pattern, as well as file patterns
rootpath = '/Volumes/Data/Imaging/R01/preprocessed/_Jason/';
subpattern = 'Sub*_v*';
rppattern = 'rp_aRSrun*.txt';
numruns = 4;

% get subjects
dirinfo = dir([rootpath subpattern]);
sublist = {dirinfo.name};


% pick subject
% for sc = 1:numel(subjlist)
% %sc = 1:numel(subjlist);
% % set primary path
%   primary_path = {[rootpath subjlist{sc} filesep]};
%
% %  test_path = {primary_path,{sc}};
% end


% thresholds
transthresh = 1;
rotthresh = 1;

% get subjects folders
subfolders = dir(subpattern);
subfolders(~cat(1, subfolders.isdir)) = [];
rpfolders = {subfolders.name};
%


% iterate over subject folders to find the motion correction files
%   3DMC FILES ARE TEXT FILES WITH THE PREFIX rp (e.g., rp*.txt)
rpfiles = {subfolders.name};

for fc = 1:numel(rpfiles)

    % locate all rp_*.txt files, one for each run - typically 4 runs
    rpfiles{fc} = dir([subfolders(fc).name '/' rppattern]);
end
rpfiles = cat(1, rpfiles{:});
rpfiles = {rpfiles.name};


%   After a day working on this, I could not work out
%   how to get the nested loop to work properly.

for fc = 1:numel(rpfiles)

        for rc = 1:numel(rpfolders)
       %disp(rpfolders(rc));
        X = [rpfolders{rc},'/',rpfiles{fc}];
        disp(X);
%       i = 268 % the number of rpfiles
       %for run = 1:4
%           disp(rpfiles(run));
%           rpfiles{rc} = [rpfolders{rc} '/' rpfiles{rc}];
%           disp(rpfiles{rc});
%
%          %disp(rc);
%          % disp(rpfolders(run))


   end
end

% %


rpoutliers = cell(numel(rpfiles), 1);

%  for fc = 1:numel(rpfiles)
%      for rc = 1:numruns
%      rpfiles{fc} = [rpfolders{fc} '/' rpfiles{fc}];
%
%      test = (rpfiles{fc});
%
     %     rp = load(rpfiles{fc});
%     rp(:, 4:6) = (180/pi) .* rp(:, 4:6);
%     rpd = diff(rp);
%     rpoutliers{fc} = 1 + find(any(abs(rpd(:, 1:3)) > transthresh, 2) | ...
%         any(abs(rpd(:, 4:6)) > rotthresh, 2));
%
% %     % print out information
%     olist = sprintf('%d, ', rpoutliers{fc});
%     olist(end-1:end) = [];
%     fprintf('%-88s: %d ([%s])\n', rpfiles{fc}, numel(rpoutliers{fc}), olist);
%      end
% end
%
%
%       END SCRIPT
%



% mot = cat(1,[rpfolders(1),rpfiles(1)]);

% %   iterate over all files
% %rpfolders = {rpfiles.folder};
% rpfolders = {subfolders.name};
% rpfiles = {rpfiles.name};
