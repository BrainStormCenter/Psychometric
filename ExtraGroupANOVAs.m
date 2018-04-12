
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
writetable(tableTtest,'t-tests.txt','Delimiter',' ');

%%         SAVE WORKSPACE
ROIFCstatsOutput = ['ROIFCstats_',datestr(now, 'yyyy-mm-dd'),'.mat']
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
