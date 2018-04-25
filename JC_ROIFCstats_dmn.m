%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%		CREATED BY:		JOCHEN WEBER
%		CREATED ON:		2017-12-13
%
%		USAGE:			TESTING DMN ROI CORRELATIONS ACROSS GROUPS
%
%    	MODIFIED ON:	 2017_12_13
%         MODIFIED ON:	  2018_02_01
%         MODIFIED ON:	  2018_02_02
%         MODIFIED ON:	  2018_02_12
%         MODIFIED ON:	 2018_02_13
%    	MODIFIED ON:	 2018_04_24
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

painzgfcccs = zgfcccs(painstart:painend, painstart:painend, :, :, :);
dmnzgfcccs = zgfcccs(dmnstart:dmnend, dmnstart:dmnend, :, :, :);
bothzgfcccs = zgfcccs(bothstart:bothend, bothstart:bothend, :, :, :);

% average connectivity strengths
%   THESE ARE THE CCs FOR ALL 90 SUBJECTS ACROSS THE PAIN REGIONS
painnet = squeeze(sum(sum(painzgfcccs, 1), 2)) ./ (length(pain) * (length(pain) -1));
%   THESE ARE THE CCs FOR ALL 90 SUBJECTS ACROSS THE DMN REGIONS
dmnnet = squeeze(sum(sum(dmnzgfcccs, 1), 2)) ./ (length(dmn) * (length(dmn) -1));
%

% to unpack:
% - i1 and i2 are the indices for groups HC and CLBP
% - the next ", 1" is the "pre" (treatment) selection
% - the next ", 1" is the "neg session" selection
%
%           LIST OUT THE BRAIN REGIONS IN THE PAIN AND DMN NETWORKS
%           LIST OUT THE BRAIN REGIONS IN THE PAIN AND DMN NETWORKS
char(voinames(voiorder));
painnames = char(voinames(voiorder(painstart:painend)));
dmnnames = char(voinames(voiorder(dmnstart:dmnend)));
bothnames = char(voinames(voiorder(bothstart:bothend)));

% computing the ANOVA for all pairs
% computing the ANOVA for all pairs
Dmn_anovaresults_effect = zeros(length(pain),length(pain));
Dmn_anovaresults_pvalue = zeros(length(pain),length(pain));

for node1 = 1:length(dmn)
     for node2 = 1:length(dmn)
        preDmnHC = mean(squeeze(dmnzgfcccs(node1, node2, i1, 1, :)), 2);
        preDmnCLBP = mean(squeeze(dmnzgfcccs(node1, node2, i2, 1, :)), 2);
        preDmnFM = mean(squeeze(dmnzgfcccs(node1, node2, i3, 1, :)), 2);
%
%         % place the code between lines 116 and 131 here
        gpNames = {'HC','CLBP','FM'};               % VARIABLE OF GROUP NAMES
%       STEP 2 = CREATE AN ARRAY OF THE COMBINED VARIABLES FROM ABOVE
%       THE ARRAY NEEDS TO BE PADDED BECAUSE OF UNEVEN GROUP SIZES
%       IDENTIFY THE LARGEST GROUP
        A = max([length(i1),length(i2),length(i3)]);
        A = zeros(A,3);    % INITIALIZE ARRAY OF ALL ZEROS  FOR LARGEST GROUP
        A(A == 0) = NaN;    % CONVERT ALL '0' TO 'NaN' (MISSING VALUES)
        A(1:length(preDmnHC),1) = preDmnHC;       % HC TO COLUMN 1
        A(1:length(preDmnCLBP),2) = preDmnCLBP;   % CLBP TO COLUMN 2
        A(1:length(preDmnFM),3) = preDmnFM;       % FM TO COLUMN 3
%       STEP 3 = RUNNING THE ANOVA AND MULTIPLE COMPARISONS
%       CREATE A TABLE OF OVERALL F-TEST
        [p,tbl,stats] = anova1(A,gpNames, 'off');      % TABLE OF OVERALL RESULTS
%         ftestNames = tbl(1,:);                  % VARIABLE NAMES FOR THE TABLE
%         ftestNames{1,6} = 'Prob_F';             % FIX THE SYMBOL ISSUE
%         tableFtest = array2table(tbl(2:4,:),'VariableNames',ftestNames);
        % pain_anovaresults_effect(node1, node2) = SOME_VALUE;
        dmn_anovaresults_pvalue(node1, node2) = p;
    end
end

%%           ANALYSIS #0 (3 GROUP ANOVA FOR PRE - DMN REGIONS)
%       COMPUTING 3-GROUP ANOVA FOR THE PRE-MANIPULATION RESTING STATE SCANS
%       STEP 1 = CREATE VARIABLES OF THE MEAN CORRELATION OF ALL DMN REGIONS
%       FOR EACH GROUP OF THE PRE SCANS ACROSS BOTH VISITS
preDmnHC = mean(dmnnet(i1,1,:),3);        % MEAN OF HC
preDmnCLBP = mean(dmnnet(i2,1,:),3);      % MEAN OF CLBP
preDmnFM = mean(dmnnet(i3,1,:),3);        % MEAN OF FM
gpNames = {'HC','CLBP','FM'};               % VARIABLE OF GROUP NAMES
%       STEP 2 = CREATE AN ARRAY OF THE COMBINED VARIABLES FROM ABOVE
%       THE ARRAY NEEDS TO BE PADDED BECAUSE OF UNEVEN GROUP SIZES
%       IDENTIFY THE LARGEST GROUP
A = max([length(i1),length(i2),length(i3)]);
A = zeros(A,3);    % INITIALIZE ARRAY OF ALL ZEROS  FOR LARGEST GROUP
A(A == 0) = NaN;    % CONVERT ALL '0' TO 'NaN' (MISSING VALUES)
A(1:length(preDmnHC),1) = preDmnHC;       % HC TO COLUMN 1
A(1:length(preDmnCLBP),2) = preDmnCLBP;   % CLBP TO COLUMN 2
A(1:length(preDmnFM),3) = preDmnFM;       % FM TO COLUMN 3
%       STEP 3 = RUNNING THE ANOVA AND MULTIPLE COMPARISONS
%       CREATE A TABLE OF OVERALL F-TEST
[p,tbl,stats] = anova1(A,gpNames);      % TABLE OF OVERALL RESULTS
ftestNamesPreDmn = tbl(1,:);                  % VARIABLE NAMES FOR THE TABLE
ftestNamesPreDmn{1,6} = 'Prob_F';             % FIX THE SYMBOL ISSUE
tableFtestPreDmn = array2table(tbl(2:4,:),'VariableNames',ftestNamesPreDmn);
figure;
%[~,~,stats] = anova1(A,gpNames);        % I AM NOT SURE WHAT THIS DOES ...
[c,~,~,gnames] = multcompare(stats);    % EVALUATE MULTIPLE COMPARISONS
%       STEP 4 = PREPARING STATS OUTPUT
%       CREATE AN ARRAY OF ANOVA OUTPUT
anovaOutputPreDmn = [gnames(c(:,1)), gnames(c(:,2)), num2cell(c(:,3:6))];
%       INITIAL ORDER OF OUTPUT FROM THE MULTICOMPARISON STEP
%       COLUMNS 1-6 =  {'gp1','gp2','lCI','gpDiff','uCI','pval'}
%       CHANGING THE VARIABLE ORDER IN THE OUTPUT ARRAY TO
%       COLUMNS 1-6 = {'gp1','gp2', 'pval','gpDiff','lCI','uCI'}) AND THEN
%       CREATE A TABLE OF THE MULTICOMPARISON OUTPUT
anovaOutputPreDmn = anovaOutputPreDmn(:,[1 2 6 4 3 5]);
tableAnovaPreDmn = array2table(anovaOutputPreDmn, 'VariableNames',{'gp1','gp2', 'pval','gpDiff','lCI','uCI'});


%%             SYNTAX FOR PERFORMING A MANOVA
%              START BY CREATING A VECTOR REPRESENTING ALL THE GROUPS
gp1 = ones(length(i1),1);
gp2 = 2*ones(length(i2),1);
gp3 = 3*ones(length(i3),1);
gps123 = cat(1,gp1,gp2,gp3);
groups = nominal(gps123);               % SPECIFY THIS AS AN ORDINAL VARIABLE
groups = setlabels(groups,{'HC','CLBP','FM'});    % SET THE VARIABLE LABELS
preDmn123 = cat(1,preDmnHC, preDmnCLBP,preDmnFM);  % CREATE ANOTHER VECTOR TO INCLUDE
%              RUN THE MANOVA
[d,p,stats] = manova1(preDmn123,groups)

%           ANALYSIS #1 (3 GROUP ANOVA FOR POST COLLAPSED ACROSS CONDITIONS)
%              COMPUTING 3-GROUP ANOVA FOR THE POST-MANIPULATION RESTING STATE SCANS
%       STEP 1 = CREATE VARIABLES OF THE MEAN CORRELATION OF ALL DMN REGIONS
%         FOR EACH GROUP OF THE POST SCANS ACROSS BOTH VISITS
postDmnHC = mean(dmnnet(i1,2,:),3);        % MEAN OF HC
postDmnCLBP = mean(dmnnet(i2,2,:),3);      % MEAN OF CLBP
postDmnFM = mean(dmnnet(i3,2,:),3);        % MEAN OF FM
gpNames = {'HC post','CLBP post','FM post'};               % VARIABLE OF GROUP NAMES
%       STEP 2 = CREATE AN ARRAY OF THE COMBINED VARIABLES FROM ABOVE
%         THE ARRAY NEEDS TO BE PADDED BECAUSE OF UNEVEN GROUP SIZES
%         IDENTIFY THE LARGEST GROUP
A = max([length(i1),length(i2),length(i3)]);
A = zeros(A,3);    % INITIALIZE ARRAY OF ALL ZEROS  FOR LARGEST GROUP
A(A == 0) = NaN;    % CONVERT ALL '0' TO 'NaN' (MISSING VALUES)
A(1:length(postDmnHC),1) = postDmnHC;       % HC TO COLUMN 1
A(1:length(postDmnCLBP),2) = postDmnCLBP;   % CLBP TO COLUMN 2
A(1:length(postDmnFM),3) = postDmnFM;       % FM TO COLUMN 3
%       STEP 3 = RUNNING THE ANOVA AND MULTIPLE COMPARISONS
%         CREATE A TABLE OF OVERALL F-TEST
[p,tbl,stats] = anova1(A,gpNames,'off');      % TABLE OF OVERALL RESULTS
ftestNamesPostDmn = tbl(1,:);                  % VARIABLE NAMES FOR THE TABLE
ftestNamesPostDmn{1,6} = 'Prob_F';             % FIX THE SYMBOL ISSUE
tableFtestPostDmn = array2table(tbl(2:4,:),'VariableNames',ftestNamesPostDmn);
figure;
%[~,~,stats] = anova1(A,gpNames);        % I AM NOT SURE WHAT THIS DOES ...
[c,~,~,gnames] = multcompare(stats);    % EVALUATE MULTIPLE COMPARISONS
%       STEP 4 = PREPARING STATS OUTPUT
%       CREATE AN ARRAY OF ANOVA OUTPUT
anovaOutputPostDmn = [gnames(c(:,1)), gnames(c(:,2)), num2cell(c(:,3:6))];
%       INITIAL ORDER OF OUTPUT FROM THE MULTICOMPARISON STEP
%       COLUMNS 1-6 =  {'gp1','gp2','lCI','gpDiff','uCI','pval'}
%       CHANGING THE VARIABLE ORDER IN THE OUTPUT ARRAY TO
%       COLUMNS 1-6 = {'gp1','gp2', 'pval','gpDiff','lCI','uCI'}) AND THEN
%       CREATE A TABLE OF THE MULTICOMPARISON OUTPUT
anovaOutputPostDmn = anovaOutputPostDmn(:,[1 2 6 4 3 5]);
tableAnovaPostDmn = array2table(anovaOutputPostDmn, 'VariableNames',{'gp1','gp2', 'pval','gpDiff','lCI','uCI'});







%%         SAVE WORKSPACE
DmnROIFCstatsOutput = ['DmnROIFCstats_',datestr(now, 'yyyy-mm-dd'),'.mat']
save(DmnROIFCstatsOutput);


%
% % painzgfcccs = zgfcccs(1:16, 1:16, :, :, :);
% % dmnzgfcccs = zgfcccs(17:22, 17:22, :, :, :);
%
% % average connectivity strengths
% %   THESE ARE THE CCs FOR ALL 90 SUBJECTS ACROSS THE 16 PAIN REGIONS
% painnet = squeeze(sum(sum(painzgfcccs, 1), 2)) ./ (16 * 15);
% %   THESE ARE THE CCs FOR ALL 90 SUBJECTS ACROSS THE 6 DMN REGIONS
% dmnnet = squeeze(sum(sum(dmnzgfcccs, 1), 2)) ./ (6 * 5);
%
% % get left amygdala (region 2) to left anterior insula (region 8) from pain network
% pain_lamyg_2_linsula = squeeze(painzgfcccs(2, 8, :, :, :));

%
% char(voinames(voiorder));
% painnames = char(voinames(voiorder(1:16)));
% dmnnames = char(voinames(voiorder(17:22)));

%
% pain_anovaresults_effect = zeros(16, 16);
% pain_anovaresults_pvalue = zeros(16, 16);
%
% for node1 = 1:16
%     for node2 = 1:16
%         preDmnHC = mean(squeeze(painzgfcccs(node1, node2, i1, 1, :)), 2);
%         preDmnCLBP = mean(squeeze(painzgfcccs(node1, node2, i2, 1, :)), 2);
%         preDmnFM = mean(squeeze(painzgfcccs(node1, node2, i3, 1, :)), 2);
