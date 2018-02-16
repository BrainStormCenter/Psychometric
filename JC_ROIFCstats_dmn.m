%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%		CREATED BY:		JOCHEN WEBER
%		CREATED ON:		2017-12-13
%
%		USAGE:			TESTING ROI CORRELATIONS ACROSS GROUPS
%
%    	MODIFIED ON:	 2017_12_13
%         MODIFIED ON:	  2018_02_01
%         MODIFIED ON:	  2018_02_02
%         MODIFIED ON:	  2018_02_12
%         MODIFIED ON:	 2018_02_13
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%       BASIC SETUP FOR ANALYSES
% neuroelf library
n = neuroelf;

%         SET PRIMARY PATH
rootpath = '/Volumes/Data/Imaging/R01/preprocessed/';
cd(rootpath);

% load variable
% contains slistd! WHICH CONATAINS THE SUBJECT LIST
load FCvars.mat

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
%   THESE ARE THE CROSS CORRELATIONS OF ALL THE BRAIN REGIONS LISTED IN THE VOI FILE
afcccs = cat(3, fcccs{:});
gfcccs = reshape(afcccs(voiorder, voiorder, rlistd(:)), [nvs, nvs, ns, 2, 2, 2]);

% fisher transform
zgfcccs = n.fisherr2z(gfcccs);
zgfcccs(isinf(zgfcccs)) = 0;

% average over first and second run of each half-session
%   THESE ARE THE 'PRE' MANIPULATION RESTING STATE SCANS
zgfcccs = squeeze(mean(zgfcccs, 4));

% this leaves 5-dimensions
%   THE 5 DIMENSIONS FOR THE zgfcccs ARRAY ARE:
% 1 (ROIs IN THE PAIN NETWORK (i.e., REGIONS [1-16], AS OF February 12, 2018)
% 2 (ROIs IN THE DMN NETWORK (i.e., REGIONS [17-22], AS OF February 12, 2018)
% 3 subjects (in order of groups, 1-31 HC, 32-73 CLBP, 74-90 FM)
% 4 pre/post treatment (1 pre, 2 post)
% 5 neg/pos session (1 neg, 2 pos)

% to split into within network matrices
%   THESE ARE ARRAYS OF CROSS CORRELATIONS AMONG REGIONS IN EACH NETWORK
%   THE ARRAYS ARE ORGANIZED AS (REGIONS^REGIONS, ALL SUBJECTS, PRE & POST, POS & NEG)
painzgfcccs = zgfcccs(1:16, 1:16, :, :, :);
dmnzgfcccs = zgfcccs(17:22, 17:22, :, :, :);

% average connectivity strengths
%   THESE ARE THE CCs FOR ALL 90 SUBJECTS ACROSS THE 16 PAIN REGIONS
painnet = squeeze(mean(mean(painzgfcccs, 1), 2));
%   THESE ARE THE CCs FOR ALL 90 SUBJECTS ACROSS THE 6 DMN REGIONS
dmnnet = squeeze(mean(mean(dmnzgfcccs, 1), 2));

% to unpack:
% - i1 and i2 are the indices for groups HC and CLBP
% - the next ", 1" is the "pre" (treatment) selection
% - the next ", 1" is the "neg session" selection
%
%           LIST OUT THE BRAIN REGIONS IN THE PAIN AND DMN NETWORKS
char(voinames(voiorder));
painnames = char(voinames(voiorder(1:16)));
dmnnames = char(voinames(voiorder(17:22)));
%%           ANALYSIS #0 (3 GROUP ANOVA FOR PRE)
%       COMPUTING 3-GROUP ANOVA FOR THE PRE-MANIPULATION RESTING STATE SCANS
%       STEP 1 = CREATE VARIABLES OF THE MEAN CORRELATION OF ALL PAIN REGIONS
%       FOR EACH GROUP OF THE PRE SCANS ACROSS BOTH VISITS
preDMNHC = mean(dmnnet(i1,1,:),3);        % MEAN OF HC
preDMNCLBP = mean(dmnnet(i2,1,:),3);      % MEAN OF CLBP
preDMNFM = mean(dmnnet(i3,1,:),3);        % MEAN OF FM
gpNames = {'HC','CLBP','FM'};               % VARIABLE OF GROUP NAMES
%       STEP 2 = CREATE AN ARRAY OF THE COMBINED VARIABLES FROM ABOVE
%       THE ARRAY NEEDS TO BE PADDED BECAUSE OF UNEVEN GROUP SIZES
%       IDENTIFY THE LARGEST GROUP
A = max([length(i1),length(i2),length(i3)]);
A = zeros(A,3);    % INITIALIZE ARRAY OF ALL ZEROS  FOR LARGEST GROUP
A(A == 0) = NaN;    % CONVERT ALL '0' TO 'NaN' (MISSING VALUES)
A(1:length(preDMNHC),1) = preDMNHC;       % HC TO COLUMN 1
A(1:length(preDMNCLBP),2) = preDMNCLBP;   % CLBP TO COLUMN 2
A(1:length(preDMNFM),3) = preDMNFM;       % FM TO COLUMN 3
%       STEP 3 = RUNNING THE ANOVA AND MULTIPLE COMPARISONS
%       CREATE A TABLE OF OVERALL F-TEST
[p,tbl,stats] = anova1(A,gpNames);      % TABLE OF OVERALL RESULTS
ftestNames = tbl(1,:);                  % VARIABLE NAMES FOR THE TABLE
ftestNames{1,6} = 'Prob_F';             % FIX THE SYMBOL ISSUE
tableFtest = array2table(tbl(2:4,:),'VariableNames',ftestNames);
figure;
%[~,~,stats] = anova1(A,gpNames);        % I AM NOT SURE WHAT THIS DOES ...
[c,~,~,gnames] = multcompare(stats);    % EVALUATE MULTIPLE COMPARISONS
%       STEP 4 = PREPARING STATS OUTPUT
%       CREATE AN ARRAY OF ANOVA OUTPUT
anovaPreOutput = [gnames(c(:,1)), gnames(c(:,2)), num2cell(c(:,3:6))];
%       INITIAL ORDER OF OUTPUT FROM THE MULTICOMPARISON STEP
%       COLUMNS 1-6 =  {'gp1','gp2','lCI','gpDiff','uCI','pval'}
%       CHANGING THE VARIABLE ORDER IN THE OUTPUT ARRAY TO
%       COLUMNS 1-6 = {'gp1','gp2', 'pval','gpDiff','lCI','uCI'}) AND THEN
%       CREATE A TABLE OF THE MULTICOMPARISON OUTPUT
anovaPreDMNOutput = anovaPreOutput(:,[1 2 6 4 3 5]);
tableAnovaDMNPre = array2table(anovaPreOutput, 'VariableNames',{'gp1','gp2', 'pval','gpDiff','lCI','uCI'});
