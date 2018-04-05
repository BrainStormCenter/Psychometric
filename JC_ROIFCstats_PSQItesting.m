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
%    	MODIFIED ON:	 2018_03_20
%    	MODIFIED ON:	 2018_04_05
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
%load FCvars_PSQI.mat

% load VOI
%voi = xff('Craggs_VOIs.voi');
voi = xff('/Users/jcraggs/Documents/GitHub/Psychometric/ROIs/AALmasks1.voi');
voinames = voi.VOINames;

% indices for pain
pain = find(~cellfun('isempty', regexpi(voinames, '^Pain')));
dmn = find(~cellfun('isempty', regexpi(voinames, '^DMN')));
both = find(~cellfun('isempty', regexpi(voinames, '^Both')));
%voiorder = [pain; dmn];
voiorder = [pain; dmn; both];
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

%painstart = bothend +1;
painstart = 1;
painend = length(pain);
dmnstart = painend + 1;
dmnend = length(pain) + length(dmn);
bothstart = dmnend + 1;
bothend = length(both) + length(pain) + length(dmn);

% painend = length(both)+length(pain);
% dmnstart = painend +1;
% dmnend = length(both) + length(pain) + length(dmn);
% bothstart = 1;
% bothend = length(both)+1;

bothzgfcccs = zgfcccs(bothstart:bothend, bothstart:bothend, :, :, :);
painzgfcccs = zgfcccs(painstart:painend, painstart:painend, :, :, :);
dmnzgfcccs = zgfcccs(dmnstart:dmnend, dmnstart:dmnend, :, :, :);


%painzgfcccs = zgfcccs(1:16, 1:16, :, :, :);
%dmnzgfcccs = zgfcccs(17:22, 17:22, :, :, :);

% average connectivity strengths
%   THESE ARE THE CCs FOR ALL 90 SUBJECTS ACROSS THE PAIN REGIONS
painnet = squeeze(sum(sum(painzgfcccs, 1), 2)) ./ (length(pain) * (length(pain) -1));
%   THESE ARE THE CCs FOR ALL 90 SUBJECTS ACROSS THE DMN REGIONS
dmnnet = squeeze(sum(sum(dmnzgfcccs, 1), 2)) ./ (length(dmn) * (length(dmn) -1));
%
%%
%   THESE ARE THE CCs FOR ALL 90 SUBJECTS ACROSS THE 16 PAIN REGIONS
%painnet = squeeze(sum(sum(painzgfcccs, 1), 2)) ./ (16 * 15);
%   THESE ARE THE CCs FOR ALL 90 SUBJECTS ACROSS THE 6 DMN REGIONS
%dmnnet = squeeze(sum(sum(dmnzgfcccs, 1), 2)) ./ (6 * 5);

% get left amygdala (region 2) to left anterior insula (region 8) from pain network
%pain_lamyg_2_linsula = squeeze(painzgfcccs(2, 8, :, :, :));

% to unpack:
% - i1 and i2 are the indices for groups HC and CLBP
% - the next ", 1" is the "pre" (treatment) selection
% - the next ", 1" is the "neg session" selection
%
%           LIST OUT THE BRAIN REGIONS IN THE PAIN AND DMN NETWORKS
char(voinames(voiorder));
painnames = char(voinames(voiorder(painstart:painend)));
dmnnames = char(voinames(voiorder(dmnstart:dmnend)));
bothnames = char(voinames(voiorder(bothstart:bothend)));
%painnames = char(voinames(voiorder(1:16)));
%dmnnames = char(voinames(voiorder(17:22)));

% computing the ANOVA for all pairs
pain_anovaresults_effect = zeros(length(pain),length(pain));
pain_anovaresults_pvalue = zeros(length(pain),length(pain));
%pain_anovaresults_effect = zeros(16, 16);
%pain_anovaresults_pvalue = zeros(16, 16);
% for node1 = 1:16
%     for node2 = 1:16
for node1 = 1:length(pain)
     for node2 = 1:length(pain)
        prePainHC = mean(squeeze(painzgfcccs(node1, node2, i1, 1, :)), 2);
        prePainCLBP = mean(squeeze(painzgfcccs(node1, node2, i2, 1, :)), 2);
        prePainFM = mean(squeeze(painzgfcccs(node1, node2, i3, 1, :)), 2);
%
%         % place the code between lines 116 and 131 here
        gpNames = {'HC','CLBP','FM'};               % VARIABLE OF GROUP NAMES
%       STEP 2 = CREATE AN ARRAY OF THE COMBINED VARIABLES FROM ABOVE
%       THE ARRAY NEEDS TO BE PADDED BECAUSE OF UNEVEN GROUP SIZES
%       IDENTIFY THE LARGEST GROUP
        A = max([length(i1),length(i2),length(i3)]);
        A = zeros(A,3);    % INITIALIZE ARRAY OF ALL ZEROS  FOR LARGEST GROUP
        A(A == 0) = NaN;    % CONVERT ALL '0' TO 'NaN' (MISSING VALUES)
        A(1:length(prePainHC),1) = prePainHC;       % HC TO COLUMN 1
        A(1:length(prePainCLBP),2) = prePainCLBP;   % CLBP TO COLUMN 2
        A(1:length(prePainFM),3) = prePainFM;       % FM TO COLUMN 3
%       STEP 3 = RUNNING THE ANOVA AND MULTIPLE COMPARISONS
%       CREATE A TABLE OF OVERALL F-TEST
        [p,tbl,stats] = anova1(A,gpNames, 'off');      % TABLE OF OVERALL RESULTS
%         ftestNames = tbl(1,:);                  % VARIABLE NAMES FOR THE TABLE
%         ftestNames{1,6} = 'Prob_F';             % FIX THE SYMBOL ISSUE
%         tableFtest = array2table(tbl(2:4,:),'VariableNames',ftestNames);
        % pain_anovaresults_effect(node1, node2) = SOME_VALUE;
        pain_anovaresults_pvalue(node1, node2) = p;
    end
end


%
%
%
%                   END SCRIPT
