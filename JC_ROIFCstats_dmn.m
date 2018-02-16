
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
