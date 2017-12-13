%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%		CREATED BY:		JOCHEN WEBER
%		CREATED ON:		2017-12-13
%
%		USAGE:			TESTING ROI CORRELATIONS ACROSS GROUPS
%
%    	MODIFIED ON:	 2017_12_13
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%

% neuroelf library
n = neuroelf;

%         SET PRIMARY PATH
rootpath = '/Volumes/Data/Imaging/R01/preprocessed/';

cd(rootpath);

% load variable
load FCvars.mat % contains slistd!

% load VOI
voi = xff('Craggs_VOIs.voi');
voinames = voi.VOINames;

% indices for pain
pain = find(~cellfun('isempty', regexpi(voinames, '^pain')));
dmn = find(~cellfun('isempty', regexpi(voinames, '^dmn')));
voiorder = [pain; dmn];
nvs = numel(voiorder);

% find subjects in three groups
g1 = find(slistd(:, 3) == 1 & ~any(isnan(slistd(:, 4:11)), 2));
g2 = find(slistd(:, 3) == 2 & ~any(isnan(slistd(:, 4:11)), 2));
g3 = find(slistd(:, 3) == 3 & ~any(isnan(slistd(:, 4:11)), 2));
g123 = [g1; g2; g3];
i1 = 1:numel(g1);
i2 = i1(end) + (1:numel(g2));
i3 = i2(end) + (1:numel(g3));
ns = numel(g123);
glistd = slistd(g123, :);
rlistd = glistd(:, 4:11);

% create cc arrays
afcccs = cat(3, fcccs{:});
gfcccs = reshape(afcccs(voiorder, voiorder, rlistd(:)), [nvs, nvs, ns, 2, 2, 2]);

% fisher transform
zgfcccs = n.fisherr2z(gfcccs);
zgfcccs(isinf(zgfcccs)) = 0;

% average over first and second run of each half-session
zgfcccs = squeeze(mean(zgfcccs, 4));

% this leaves dimensions
% 1, 2 (ROIs, in order of networks, 1-16 pain, 17-22 DMN)
% 3 subjects (in order of groups, 1-31 HC, 32-73 CLBP, 74-90 FM)
% 4 pre/post treatment (1 pre, 2 post)
% 5 neg/pos session (1 neg, 2 pos)

% to split into within network matrices
painzgfcccs = zgfcccs(1:16, 1:16, :, :, :);
dmnzgfcccs = zgfcccs(17:22, 17:22, :, :, :);

% average connectivity strengths
painnet = squeeze(mean(mean(painzgfcccs, 1), 2));
dmnnet = squeeze(mean(mean(dmnzgfcccs, 1), 2));

% to then run differences in group tests HC vs. CLBP
[h, p, ci, stats] = ttest2(painnet(i1, 1, 1), painnet(i2, 1, 1), ...
    'tail', 'both', 'vartype', 'unequal');
p
stats

% to unpack:
% - i1 and i2 are the indices for groups HC and CLBP
% - the next ", 1" is the "pre" (treatment) selection
% - the next ", 1" is the "neg session" selection

% to compute a compound "pre" score
[h, p, ci, stats] = ttest2(mean(painnet(i1, 1, :), 3), mean(painnet(i2, 1, :), 3), ...
    'tail', 'both', 'vartype', 'unequal');
p
stats

% to compute a "difference" for pre-post neg treatment (between groups)
%[h, p, ci, stats] = ttest2(painnet(i1, 2, 1) - painnet(i1, 1, 1), painnet(i2, 1, 2) - painnet(i2, 1, 1), ...
[h, p, ci, stats] = ttest2(painnet(i1, 2, 1) - painnet(i1, 1, 1), painnet(i2, 2, 1) - painnet(i2, 1, 1), ...
    'tail', 'both', 'vartype', 'unequal');
p
stats

% as a paired t-test (across all participants)
[h, p, ci, stats] = ttest(painnet(:, 1, 1), painnet(:, 1, 2), ...
    'tail', 'both');
p
stats

% as a one-sample t-test
