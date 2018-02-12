%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%		CREATED BY:		JOCHEN WEBER
%		CREATED ON:		2017-12-13
%
%		USAGE:			TESTING ROI CORRELATIONS ACROSS GROUPS
%
%    	MODIFIED ON:	 2017_12_13
%   	MODIFIED ON:	  2018_02_01
%   	MODIFIED ON:	  2018_02_02
%       MODIFIED ON:	  2018_02_12
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
prePainHC = mean(painnet(i1,1,:),3);        % MEAN OF HC
prePainCLBP = mean(painnet(i2,1,:),3);      % MEAN OF CLBP
prePainFM = mean(painnet(i3,1,:),3);        % MEAN OF FM
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
anovaPreOutput = anovaPreOutput(:,[1 2 6 4 3 5]);
tableAnovaPre = array2table(anovaPreOutput, 'VariableNames',{'gp1','gp2', 'pval','gpDiff','lCI','uCI'});

%%           ANALYSIS #1 (3 GROUP ANOVA FOR POST COLLAPSED ACROSS CONDITIONS)
%       COMPUTING 3-GROUP ANOVA FOR THE PRE-MANIPULATION RESTING STATE SCANS
%       STEP 1 = CREATE VARIABLES OF THE MEAN CORRELATION OF ALL PAIN REGIONS
%       FOR EACH GROUP OF THE PRE SCANS ACROSS BOTH VISITS
postPainHC = mean(painnet(i1,2,:),3);        % MEAN OF HC
postPainCLBP = mean(painnet(i2,2,:),3);      % MEAN OF CLBP
postPainFM = mean(painnet(i3,2,:),3);        % MEAN OF FM
gpNames = {'HC post','CLBP post','FM post'};               % VARIABLE OF GROUP NAMES
%       STEP 2 = CREATE AN ARRAY OF THE COMBINED VARIABLES FROM ABOVE
%       THE ARRAY NEEDS TO BE PADDED BECAUSE OF UNEVEN GROUP SIZES
%       IDENTIFY THE LARGEST GROUP
A = max([length(i1),length(i2),length(i3)]);
A = zeros(A,3);    % INITIALIZE ARRAY OF ALL ZEROS  FOR LARGEST GROUP
A(A == 0) = NaN;    % CONVERT ALL '0' TO 'NaN' (MISSING VALUES)
A(1:length(postPainHC),1) = postPainHC;       % HC TO COLUMN 1
A(1:length(postPainCLBP),2) = postPainCLBP;   % CLBP TO COLUMN 2
A(1:length(postPainFM),3) = postPainFM;       % FM TO COLUMN 3
%       STEP 3 = RUNNING THE ANOVA AND MULTIPLE COMPARISONS
%       CREATE A TABLE OF OVERALL F-TEST
[p,tbl,stats] = anova1(A,gpNames,'off');      % TABLE OF OVERALL RESULTS
ftestNames = tbl(1,:);                  % VARIABLE NAMES FOR THE TABLE
ftestNames{1,6} = 'Prob_F';             % FIX THE SYMBOL ISSUE
tableFtest = array2table(tbl(2:4,:),'VariableNames',ftestNames);
figure;
%[~,~,stats] = anova1(A,gpNames);        % I AM NOT SURE WHAT THIS DOES ...
[c,~,~,gnames] = multcompare(stats);    % EVALUATE MULTIPLE COMPARISONS
%       STEP 4 = PREPARING STATS OUTPUT
%       CREATE AN ARRAY OF ANOVA OUTPUT
anovaPostOutput = [gnames(c(:,1)), gnames(c(:,2)), num2cell(c(:,3:6))];
%       INITIAL ORDER OF OUTPUT FROM THE MULTICOMPARISON STEP
%       COLUMNS 1-6 =  {'gp1','gp2','lCI','gpDiff','uCI','pval'}
%       CHANGING THE VARIABLE ORDER IN THE OUTPUT ARRAY TO
%       COLUMNS 1-6 = {'gp1','gp2', 'pval','gpDiff','lCI','uCI'}) AND THEN
%       CREATE A TABLE OF THE MULTICOMPARISON OUTPUT
anovaPostOutput = anovaPostOutput(:,[1 2 6 4 3 5]);
tableAnovaPost = array2table(anovaPostOutput, 'VariableNames',{'gp1','gp2', 'pval','gpDiff','lCI','uCI'});

%%           ANALYSIS #2 (3 GROUP ANOVA FOR POST NEGAVIVE MANIPULATION)
%       COMPUTING 3-GROUP ANOVA FOR THE PRE-MANIPULATION RESTING STATE SCANS
%       STEP 1 = CREATE VARIABLES OF THE MEAN CORRELATION OF ALL PAIN REGIONS
%       FOR EACH GROUP OF THE PRE SCANS ACROSS BOTH VISITS
postNegPainHC = mean(painnet(i1,2,1),3);        % MEAN OF HC
postNegPainCLBP = mean(painnet(i2,2,1),3);      % MEAN OF CLBP
postNegPainFM = mean(painnet(i3,2,1),3);        % MEAN OF FM
gpNames = {'HC neg','CLBP neg','FM neg'};               % VARIABLE OF GROUP NAMES
%       STEP 2 = CREATE AN ARRAY OF THE COMBINED VARIABLES FROM ABOVE
%       THE ARRAY NEEDS TO BE PADDED BECAUSE OF UNEVEN GROUP SIZES
%       IDENTIFY THE LARGEST GROUP
A = max([length(i1),length(i2),length(i3)]);
A = zeros(A,3);    % INITIALIZE ARRAY OF ALL ZEROS  FOR LARGEST GROUP
A(A == 0) = NaN;    % CONVERT ALL '0' TO 'NaN' (MISSING VALUES)
A(1:length(postNegPainHC),1) = postNegPainHC;       % HC TO COLUMN 1
A(1:length(postNegPainCLBP),2) = postNegPainCLBP;   % CLBP TO COLUMN 2
A(1:length(postNegPainFM),3) = postNegPainFM;       % FM TO COLUMN 3
%       STEP 3 = RUNNING THE ANOVA AND MULTIPLE COMPARISONS
%       CREATE A TABLE OF OVERALL F-TEST
[p,tbl,stats] = anova1(A,gpNames,'off');      % TABLE OF OVERALL RESULTS
ftestNames = tbl(1,:);                  % VARIABLE NAMES FOR THE TABLE
ftestNames{1,6} = 'Prob_F';             % FIX THE SYMBOL ISSUE
tableFtest = array2table(tbl(2:4,:),'VariableNames',ftestNames);
figure;
%[~,~,stats] = anova1(A,gpNames);        % I AM NOT SURE WHAT THIS DOES ...
[c,~,~,gnames] = multcompare(stats);    % EVALUATE MULTIPLE COMPARISONS
%       STEP 4 = PREPARING STATS OUTPUT
%       CREATE AN ARRAY OF ANOVA OUTPUT
anovaPostNegOutput = [gnames(c(:,1)), gnames(c(:,2)), num2cell(c(:,3:6))];
%       INITIAL ORDER OF OUTPUT FROM THE MULTICOMPARISON STEP
%       COLUMNS 1-6 =  {'gp1','gp2','lCI','gpDiff','uCI','pval'}
%       CHANGING THE VARIABLE ORDER IN THE OUTPUT ARRAY TO
%       COLUMNS 1-6 = {'gp1','gp2', 'pval','gpDiff','lCI','uCI'}) AND THEN
%       CREATE A TABLE OF THE MULTICOMPARISON OUTPUT
anovaPostNegOutput = anovaPostNegOutput(:,[1 2 6 4 3 5]);
tableAnovaNegPost = array2table(anovaPostNegOutput, 'VariableNames',{'gp1','gp2', 'pval','gpDiff','lCI','uCI'});

%%           ANALYSIS #3 (3 GROUP ANOVA FOR POST POSITIVE MANIPULATION)
%       COMPUTING 3-GROUP ANOVA FOR THE PRE-MANIPULATION RESTING STATE SCANS
%       STEP 1 = CREATE VARIABLES OF THE MEAN CORRELATION OF ALL PAIN REGIONS
%       FOR EACH GROUP OF THE PRE SCANS ACROSS BOTH VISITS
postPosPainHC = mean(painnet(i1,2,2),3);        % MEAN OF HC
postPosPainCLBP = mean(painnet(i2,2,2),3);      % MEAN OF CLBP
postPosPainFM = mean(painnet(i3,2,2),3);        % MEAN OF FM
gpNames = {'HC pos','CLBP pos','FM pos'};       % VARIABLE OF GROUP NAMES
%       STEP 2 = CREATE AN ARRAY OF THE COMBINED VARIABLES FROM ABOVE
%       THE ARRAY NEEDS TO BE PADDED BECAUSE OF UNEVEN GROUP SIZES
%       IDENTIFY THE LARGEST GROUP
A = max([length(i1),length(i2),length(i3)]);
A = zeros(A,3);    % INITIALIZE ARRAY OF ALL ZEROS  FOR LARGEST GROUP
A(A == 0) = NaN;    % CONVERT ALL '0' TO 'NaN' (MISSING VALUES)
A(1:length(postPosPainHC),1) = postPosPainHC;       % HC TO COLUMN 1
A(1:length(postPosPainCLBP),2) = postPosPainCLBP;   % CLBP TO COLUMN 2
A(1:length(postPosPainFM),3) = postPosPainFM;       % FM TO COLUMN 3
%       STEP 3 = RUNNING THE ANOVA AND MULTIPLE COMPARISONS
%       CREATE A TABLE OF OVERALL F-TEST
[p,tbl,stats] = anova1(A,gpNames,'off');      % TABLE OF OVERALL RESULTS
ftestNames = tbl(1,:);                  % VARIABLE NAMES FOR THE TABLE
ftestNames{1,6} = 'Prob_F';             % FIX THE SYMBOL ISSUE
tableFtest = array2table(tbl(2:4,:),'VariableNames',ftestNames);
figure;
%[~,~,stats] = anova1(A,gpNames);        % I AM NOT SURE WHAT THIS DOES ...
[c,~,~,gnames] = multcompare(stats);    % EVALUATE MULTIPLE COMPARISONS
%       STEP 4 = PREPARING STATS OUTPUT
%       CREATE AN ARRAY OF ANOVA OUTPUT
anovaPostPosOutput = [gnames(c(:,1)), gnames(c(:,2)), num2cell(c(:,3:6))];
%       INITIAL ORDER OF OUTPUT FROM THE MULTICOMPARISON STEP
%       COLUMNS 1-6 =  {'gp1','gp2','lCI','gpDiff','uCI','pval'}
%       CHANGING THE VARIABLE ORDER IN THE OUTPUT ARRAY TO
%       COLUMNS 1-6 = {'gp1','gp2', 'pval','gpDiff','lCI','uCI'}) AND THEN
%       CREATE A TABLE OF THE MULTICOMPARISON OUTPUT
anovaPostPosOutput = anovaPostPosOutput(:,[1 2 6 4 3 5]);
tableAnovaPosPost = array2table(anovaPostPosOutput, 'VariableNames',{'gp1','gp2', 'pval','gpDiff','lCI','uCI'});



%%          ANALYSIS #4 (T-TEST, HC vs. CLBP [pre, both visits])
% to compute a compound "pre" score
%   ttest2(group-1, pre/post, neg/post, group-2, pre/post, neg/post)
%   Below = (HC, pre, both visits vs. CLBP, pre, both visits)
[h, p, ci, stats] = ttest2(mean(painnet(i1, 1, :), 3), mean(painnet(i2, 1, :), 3), ...
    'tail', 'both', 'vartype', 'unequal');
% p
% stats
    Gps = (['HC v CLBP']);
    pval = ([p]);
    tVal = ([stats.tstat]);
    scanSet = (['combined pre']);
    cond = (['both']);
    Hy = ([h]);
%%          ANALYSIS #5 (T-TEST [HC v CLBP], COMPARING PRE FROM THE NEG VISITS, NOT OVERLY USEFUL)
% to then run differences in group tests HC vs. CLBP
%   ttest2(group-1, pre/post, neg/post, group-2, pre/post, neg/post)
%   Below = (HC, pre, neg vs. CLBP, pre, neg)
[h, p, ci, stats] = ttest2(painnet(i1, 1, 1), painnet(i2, 1, 1), ...
    'tail', 'both', 'vartype', 'unequal');
% p
% stats
    Gps = ([Gps;{'HC v CLBP'}]);
    pval = ([pval;p]);
    tVal = ([tVal;stats.tstat]);
    scanSet = ([scanSet;{'pre'}]);
    cond = ([cond;{'neg'}]);
    Hy = ([Hy;h]);

%%           ANALYSIS #6 (T-TEST [HC v CLBP], COMPARING PRE FROM THE POS VISITS, NOT OVERLY USEFUL)
%   ttest2(group-1, pre/post, neg/post, group-2, pre/post, neg/post)
%   Below = (HC, pre, pos vs. CLBP, pre, pos)
[h, p, ci, stats] = ttest2(painnet(i1, 1, 2), painnet(i2, 1, 2), ...
    'tail', 'both', 'vartype', 'unequal');
    Gps = ([Gps;{'HC v CLBP'}]);
    pval = ([pval;p]);
    tVal = ([tVal;stats.tstat]);
    scanSet = ([scanSet;{'pre'}]);
    cond = ([cond;{'pos'}]);
    Hy = ([Hy;h]);

%%           ANALYSIS #7 (T-TEST [HC v CLBP], POST > PRE (NEG), WHAT DID THE NEG MANIPULATION INCREASE)
%   TO TEST GROUP DIFFERENCES IN POST-PRE IN THE NEGATIVE CONDITION
%   THAT IS, GROUP DIFFERENCES IN HOW MUCH THE NEGATIVE MOOD MANIPULATION
%   INCREASED ACTIVITY IN THE ROIs
[h, p, ci, stats] = ttest2(painnet(i1, 2, 1) - painnet(i1, 1, 1), painnet(i2, 2, 1) - painnet(i2, 1, 1), ...
    'tail', 'both', 'vartype', 'unequal');
    Gps = ([Gps;{'HC v CLBP'}]);
    pval = ([pval;p]);
    tVal = ([tVal;stats.tstat]);
    scanSet = ([scanSet;{'post > pre'}]);
    cond = ([cond;{'neg'}]);
    Hy = ([Hy;h]);

%%           ANALYSIS #8 (T-TEST [HC v CLBP], POST > PRE (POS), WHAT DID THE POS MANIPULATION INCREASE)
%   TO TEST GROUP DIFFERENCES IN POST-PRE IN THE POSITIVE CONDITION
%   THAT IS, GROUP DIFFERENCES IN HOW MUCH THE POSITIVE MOOD MANIPULATION
%   INCREASED ACTIVITY IN THE ROIs
[h, p, ci, stats] = ttest2(painnet(i1, 2, 2) - painnet(i1, 1, 2), painnet(i2, 2, 2) - painnet(i2, 1, 2), ...
    'tail', 'both', 'vartype', 'unequal');
    Gps = ([Gps;{'HC v CLBP'}]);
    pval = ([pval;p]);
    tVal = ([tVal;stats.tstat]);
    scanSet = ([scanSet;{'post > pre'}]);
    cond = ([cond;{'pos'}]);
    Hy = ([Hy;h]);

    %%           ANALYSIS #9 (T-TEST [HC v FM], POST > PRE (NEG), WHAT DID THE NEG MANIPULATION INCREASE)
%   TO TEST GROUP DIFFERENCES IN POST-PRE IN THE NEGATIVE CONDITION
%   THAT IS, GROUP DIFFERENCES IN HOW MUCH THE NEGATIVE MOOD MANIPULATION
%   INCREASED ACTIVITY IN THE ROIs
[h, p, ci, stats] = ttest2(painnet(i1, 2, 1) - painnet(i1, 1, 1), painnet(i3, 2, 1) - painnet(i3, 1, 1), ...
    'tail', 'both', 'vartype', 'unequal');
    Gps = ([Gps;{'HC v FM'}]);
    pval = ([pval;p]);
    tVal = ([tVal;stats.tstat]);
    scanSet = ([scanSet;{'post > pre'}]);
    cond = ([cond;{'neg'}]);
    Hy = ([Hy;h]);

%%           ANALYSIS #10 (T-TEST [HC v FM], POST > PRE (POS), WHAT DID THE POS MANIPULATION INCREASE)
%   TO TEST GROUP DIFFERENCES IN POST-PRE IN THE POSITIVE CONDITION
%   THAT IS, GROUP DIFFERENCES IN HOW MUCH THE POSITIVE MOOD MANIPULATION
%   INCREASED ACTIVITY IN THE ROIs
[h, p, ci, stats] = ttest2(painnet(i1, 2, 2) - painnet(i1, 1, 2), painnet(i3, 2, 2) - painnet(i3, 1, 2), ...
    'tail', 'both', 'vartype', 'unequal');
    Gps = ([Gps;{'HC v FM'}]);
    pval = ([pval;p]);
    tVal = ([tVal;stats.tstat]);
    scanSet = ([scanSet;{'post > pre'}]);
    cond = ([cond;{'pos'}]);
    Hy = ([Hy;h]);

%%       CREATING A TABLE OF OUTPUT VARIABLES
tableTtest = table(Gps,scanSet,cond,Hy,pval,tVal,'VariableNames',{'group','scanSet','Condition','Sig','pvalue','tvalue'});


%%         SAVE WORKSPACE
ROIFCstatsOutput = ['ROIFCstats_',datestr(now, 'dd-mm-yyyy'),'.mat']
save(ROIFCstatsOutput);






% %   Below = (HC, pre, neg vs. CLBP, pre, pos)
% [h, p, ci, stats] = ttest2(painnet(i1, 1, 1), painnet(i2, 1, 1), ...
%     'tail', 'both', 'vartype', 'unequal');
%



%%              EXTRA STUFF

% pval2 = p;

%Qs = ([Qs;Q2]);
%Qs = ([Qs;{'HC ver CLBP'}]);


% % to compute a "difference" for pre-post neg treatment (between groups)
% %[h, p, ci, stats] = ttest2(painnet(i1, 2, 1) - painnet(i1, 1, 1), painnet(i2, 1, 2) - painnet(i2, 1, 1), ...
% [h, p, ci, stats] = ttest2(painnet(i1, 2, 1) - painnet(i1, 1, 1), painnet(i2, 2, 1) - painnet(i2, 1, 1), ...
%     'tail', 'both', 'vartype', 'unequal');
% p
% stats
%pval = {pval;p};

% as a paired t-test (across all participants)
% [h, p, ci, stats] = ttest(painnet(:, 1, 1), painnet(:, 1, 2), ...
%     'tail', 'both');
% p
% stats
% pval = ([pval;p]);
% tVal = ([tVal;stats.tstat]);
% pval3 = p;

% as a one-sample t-test
%
%
%
%
%
%
%
%                   END SCRIPT
