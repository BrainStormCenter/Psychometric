%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%		CREATED BY:		JOCHEN WEBER
%		CREATED ON:		2017-12-13
%
%		USAGE:			TESTING ROI CORRELATIONS ACROSS GROUPS
%
%         LATEST MODIFICATION:     2018_04_26
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
voi = xff('/Users/jcraggs/Documents/GitHub/Psychometric/ROIs/AALmasks1.voi');
voinames = voi.VOINames;

% indices for pain
pain = find(~cellfun('isempty', regexpi(voinames, '^Pain')));
dmn = find(~cellfun('isempty', regexpi(voinames, '^DMN')));
both = find(~cellfun('isempty', regexpi(voinames, '^Both')));
%voiorder = [pain; dmn];
voiorder = [pain; dmn; both];
nvs = numel(voiorder);

%    DETERMINING THE NUMBER OF ROIS FOR EACH NETWORK
painstart = 1;
painend = length(pain);
dmnstart = painend + 1;
dmnend = length(pain) + length(dmn);
bothstart = dmnend + 1;
bothend = length(both) + length(pain) + length(dmn);

%         LIST OUT THE BRAIN REGIONS IN THE PAIN AND DMN NETWORKS
char(voinames(voiorder));
painnames = char(voinames(voiorder(painstart:painend)));
dmnnames = char(voinames(voiorder(dmnstart:dmnend)));
bothnames = char(voinames(voiorder(bothstart:bothend)));

%         CLEANING UP ROI NAMES
painnames2 = cellstr(painnames);  % CONVERT FROM CHAR TO CELL
dmnnames2 = cellstr(dmnnames);  % CONVERT FROM CHAR TO CELL
painExpr = '([A-Z][a-z].+_)(rMNI_)([A-Za-z].+)(_)(roi.nii)';     % REGEX EXPRESSION
dmnExpr = '([A-Z].+_)(rMNI_)([A-Za-z].+)(_)(roi.nii)';           % REGEX EXPRESSION
newROI = '$1$3'; % NEW ROI NAME BASED ON REGEX EXPRESSION ABOVE
painnames2 = regexprep(painnames2, painExpr, newROI);
dmnnames2 = regexprep(dmnnames2, dmnExpr, newROI);
painnames3 = char(painnames2);
dmnnames3 = char(dmnnames2);

%         COMBINE BEHAVIORAL AND DEMOGRAPHIC DATA
%         THIS NEEDS TO BE DONE BEFORE IDENTIFYING SUBJECTS BELOW
%         TO AVOID MISSING DATA
slistdORIG = slistd;                              % PRESERVE ORIGINAL DATA
slistd = [slistd,struct2array(psqiStruct)];       % ADD BEHAVIORAL DATA

%         find subjects in three groups
g1 = find(slistd(:, 3) == 1 & ~any(isnan(slistd(:, 4:24)), 2));
g2 = find(slistd(:, 3) == 2 & ~any(isnan(slistd(:, 4:24)), 2));
g3 = find(slistd(:, 3) == 3 & ~any(isnan(slistd(:, 4:24)), 2));
g123 = [g1; g2; g3];
g23 = [g2;g3];      % BOTH PAIN GROUPS
i1 = 1:numel(g1);
i2 = i1(end) + (1:numel(g2));
i3 = i2(end) + (1:numel(g3));
i4 = [i2, i3];
ns = numel(g123);
glistd = slistd(g123, :);
rlistd = glistd(:, 4:11);
rlistd2 = glistd(:,[1 4:11]); % KEEP SUBJECT NUMBERS

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

%    this leaves 5-dimensions
%    THE 5 DIMENSIONS FOR THE zgfcccs ARRAY ARE:
%    1 (ROIs IN THE PAIN NETWORK (i.e., REGIONS [1-16], AS OF February 12, 2018)
%    2 (ROIs IN THE DMN NETWORK (i.e., REGIONS [17-22], AS OF February 12, 2018)
%    3 subjects (in order of groups, 1-31 HC, 32-73 CLBP, 74-90 FM)
%    4 pre/post treatment (1 pre, 2 post)
%    5 neg/pos session (1 neg, 2 pos)
%         to split into within network matrices
%         THESE ARE ARRAYS OF CROSS CORRELATIONS AMONG REGIONS IN EACH NETWORK
%         THE ARRAYS ARE ORGANIZED AS (REGIONS^REGIONS, ALL SUBJECTS, PRE & POST, POS & NEG)

%         THESE ARE THE CCs MATRICIES OF BRAIN REGIONS FOR EACH SUBJECT IN EACH ROI SET
bothzgfcccs = zgfcccs(bothstart:bothend, bothstart:bothend, :, :, :);
painzgfcccs = zgfcccs(painstart:painend, painstart:painend, :, :, :);
dmnzgfcccs = zgfcccs(dmnstart:dmnend, dmnstart:dmnend, :, :, :);

%              average connectivity strengths
%              THESE ARE THE AVERAGE CCs FOR ALL THE SUBJECTS ACROSS THE PAIN REGIONS
painnet = squeeze(sum(sum(painzgfcccs, 1), 2)) ./ (length(pain) * (length(pain) -1));
%              THESE ARE THE CCs FOR ALL THE SUBJECTS ACROSS THE DMN REGIONS
dmnnet = squeeze(sum(sum(dmnzgfcccs, 1), 2)) ./ (length(dmn) * (length(dmn) -1));
%
%
%              CREATE MATRIX OF PAIN CC'S AND BEHAVIORAL DATA
%              RESHAPE THE PAINNET CC'S, COLUMNS 1&2=PRE/POST; COLUMNS 3&4=NEG/POS
sub_by_painCCs = [painnet(:,1:2,1),painnet(:,1:2,2)];
psqiData = glistd(:,[1 12:24]);
psqiNames = {'ID','TiB_hrs', 'SoL_min','WASO_min','TST_hrs','SleepEfficiency','psqi_Durat','psqi_Distb', ...
               'psqi_Latency','psqi_DayDys','psqi_SE','psqi_BadSQual','psqi_Meds','PSQI_total'};

PSQIandPainCCs = [psqiData,sub_by_painCCs];

%         to unpack:
%         - i1 and i2 and i3 are the indices for groups HC and CLBP and FM
%         - the next ", 1" is the "pre" (treatment) selection
%         - the next ", 1" is the "neg session" selection
%
%         COMPUTING THE ANOVA FOR ALL PAIRS OF PAIN ROIS
pain_anovaresults_effect = zeros(length(pain),length(pain));
pain_anovaresults_pvalue = zeros(length(pain),length(pain));

%%        ANALYSIS #1 (3 GROUP ANOVA FOR PRE - PAIN REGIONS)
%         COMPUTING 3-GROUP ANOVA FOR THE PRE-MANIPULATION RESTING STATE SCANS
%              STEP 1 = CREATE VARIABLES OF THE MEAN CORRELATION OF ALL PAIN REGIONS
%                   FOR EACH GROUP OF THE PRE SCANS ACROSS BOTH VISITS
HC_PainPre = mean(painnet(i1,1,:),3);        % MEAN OF HC
CLBP_PainPre = mean(painnet(i2,1,:),3);      % MEAN OF CLBP
FM_PainPre = mean(painnet(i3,1,:),3);        % MEAN OF FM
gpNames = {'HC','CLBP','FM'};               % VARIABLE OF GROUP NAMES
%              STEP 2 = CREATE AN ARRAY OF THE COMBINED VARIABLES FROM ABOVE
%                   THE ARRAY NEEDS TO BE PADDED BECAUSE OF UNEVEN GROUP SIZES
%                   IDENTIFY THE LARGEST GROUP
A = max([length(i1),length(i2),length(i3)]);
A = zeros(A,3);    % INITIALIZE ARRAY OF ALL ZEROS  FOR LARGEST GROUP
A(A == 0) = NaN;    % CONVERT ALL '0' TO 'NaN' (MISSING VALUES)
A(1:length(HC_PainPre),1) = HC_PainPre;       % HC TO COLUMN 1
A(1:length(CLBP_PainPre),2) = CLBP_PainPre;   % CLBP TO COLUMN 2
A(1:length(FM_PainPre),3) = FM_PainPre;       % FM TO COLUMN 3
%              STEP 3 = RUNNING THE ANOVA AND MULTIPLE COMPARISONS
%                   CREATE A TABLE OF OVERALL F-TEST
[p_PainPre,tbl_PainPre,stats_PainPre] = anova1(A,gpNames);      % TABLE OF OVERALL RESULTS
ftestNamesPain = tbl_PainPre(1,:);                  % VARIABLE NAMES FOR THE TABLE
ftestNamesPain{1,6} = 'Prob_F';             % FIX THE SYMBOL ISSUE
tableFtestPainPre = array2table(tbl_PainPre(2:4,:),'VariableNames',ftestNamesPain);
figure;
%[~,~,stats] = anova1(A,gpNames);        % I AM NOT SURE WHAT THIS DOES ...
[c,~,~,gnames] = multcompare(stats_PainPre);    % EVALUATE MULTIPLE COMPARISONS
%       STEP 4 = PREPARING STATS OUTPUT
%       CREATE AN ARRAY OF ANOVA OUTPUT
anovaOutputPainPre = [gnames(c(:,1)), gnames(c(:,2)), num2cell(c(:,3:6))];
%       INITIAL ORDER OF OUTPUT FROM THE MULTICOMPARISON STEP
%       COLUMNS 1-6 =  {'gp1','gp2','lCI','gpDiff','uCI','pval'}
%       CHANGING THE VARIABLE ORDER IN THE OUTPUT ARRAY TO
%       COLUMNS 1-6 = {'gp1','gp2', 'pval','gpDiff','lCI','uCI'}) AND THEN
%       CREATE A TABLE OF THE MULTICOMPARISON OUTPUT
anovaOutputPainPre = anovaOutputPainPre(:,[1 2 6 4 3 5]);
tableAnovaPainPre = array2table(anovaOutputPainPre, 'VariableNames',{'gp1','gp2', 'pval','gpDiff','lCI','uCI'});
%
%%              T-TEST COMPARING GROUPS ON THE REGION-TO-REGION CROSS-CORRELATIONS
for node1 = 1:length(pain)         % PAIN REGION #1
     for node2 = 1:length(pain)    % PAIN REGION #1
        prePainHC = mean(squeeze(painzgfcccs(node1, node2, i1, 1, :)), 2);
        % prePainCLBP = mean(squeeze(painzgfcccs(node1, node2, i2, 1, :)), 2);
        % prePainFM = mean(squeeze(painzgfcccs(node1, node2, i3, 1, :)), 2);
        prePainGps = mean(squeeze(painzgfcccs(node1, node2, i4, 1, :)), 2);     % BOTH CP GROUPS
        gpNames = {'HC','CLBP','FM'};           % VARIABLE: GROUP NAMES
        gpNames2 = {'HC','Pain'};               % VARIABLE: GROUP NAMES COLLAPSED ACROSS CP GROUPS
%       STEP 2 = CREATE AN ARRAY OF THE COMBINED VARIABLES FROM ABOVE
%       THE ARRAY NEEDS TO BE PADDED BECAUSE OF UNEVEN GROUP SIZES
%       IDENTIFY THE LARGEST GROUP
        % A = max([length(i1),length(i2),length(i3)]);
        % A = zeros(A,3);    % INITIALIZE ARRAY OF ALL ZEROS  FOR LARGEST GROUP
        % A(A == 0) = NaN;    % CONVERT ALL '0' TO 'NaN' (MISSING VALUES)
        % A(1:length(prePainHC),1) = prePainHC;       % HC TO COLUMN 1
        % A(1:length(prePainCLBP),2) = prePainCLBP;   % CLBP TO COLUMN 2
        % A(1:length(prePainFM),3) = prePainFM;       % FM TO COLUMN 3
%       STEP 3 = RUNNING THE ANOVA AND MULTIPLE COMPARISONS
%       CREATE A TABLE OF OVERALL F-TEST
        %[p,tbl,stats] = anova1(A,gpNames, 'off');      % TABLE OF OVERALL RESULTS
        [H_PainPre,P_PainPre,CI_PainPre,STATS_PainPre] = ttest2(prePainHC,prePainGps);           % 2 SAMPLE T-TEST
%         ftestNames = tbl(1,:);                  % VARIABLE NAMES FOR THE TABLE
%         ftestNames{1,6} = 'Prob_F';             % FIX THE SYMBOL ISSUE
%         tableFtest = array2table(tbl(2:4,:),'VariableNames',ftestNames);
        % pain_anovaresults_effect(node1, node2) = cell2mat(tbl(2,5)); % THESE ARE F-VALUES
        % pain_anovaresults_pvalue(node1, node2) = p;
        ttest_tval_PainPre(node1, node2) = STATS_PainPre.tstat; % tstat
        ttest_pval_PainPre(node1, node2) = P_PainPre; % pvalue

    end
end

for i = 1:24
     for j=i+1:24
          ttest_pval_PainPre(i,j)=NaN;
          ttest_tval_PainPre(i,j)=NaN;
     end
end

%         THESE IDENTIFY GROUP DIFFERENCES IN ROI-TO-ROI FUNCTIONAL CONNECTIVITY
%         FDR CORRECTION
pid = FDR(ttest_pval_PainPre,.05);
[I,J] = find(ttest_pval_PainPre <= pid);
[I J];

%% Fancy code to print-out RESULTS

% first, create an empty array of character strings
OUT_Text_PainPre = [];
% Now, add strings to that OUT_TEXT
% and you can use the "sprintf" command for syntax like tabs, line breaks

OUT_Text_PainPre = [OUT_Text_PainPre 'For the pain ROIs, pre mood manipulation, there are sig ' ...
                    'group diffs in the CCs of these ROI-pairs:' sprintf('\t') ...
                    sprintf('\n')];
for i=1:numel(I)
     pairNum = num2str(i);
     roi1num = num2str(I(i));
     roi1str = painnames3(I(i),:);
     roi2num = num2str(J(i));
     roi2str = painnames3(J(i),:);
     PainROI_tval = ttest_tval_PainPre(I(i),J(i));
     PainROI_pval = ttest_pval_PainPre(I(i),J(i));
     % OUT_Text_PainPre = [OUT_Text_PainPre roi1str ' (#' roi1num ') with ' roi2str ' (#' roi2num ')' sprintf('\t') ...
     %  't-score: ' num2str(PainROI_tval) sprintf('\t') ' p-value: ' sprintf('%0.05f',PainROI_pval) ...
     %  sprintf('\n') ];
      % OUT_Text_PainPre = [OUT_Text_PainPre roi1str '(#',roi1num,') with ',roi2str '(#', roi2num ')' ...
      % 't-score: ' num2str(PainROI_tval) ' p-value: ' sprintf('%0.05f',PainROI_pval) ...
      % sprintf('\n') ];
      OUT_Text_PainPre = [OUT_Text_PainPre pairNum, '. ' roi1str, '(#',roi1num,') with ',roi2str, ...
      '(#', roi2num,')', sprintf('\t'), 't-val: ', num2str(PainROI_tval), sprintf('\t'), 'p-val: ' ...
      sprintf('%0.05f',PainROI_pval) sprintf('\n')];
end

OUT_Text_PainPre
%         NEGATIVE T-VALUES = CONTROL LESS THAN PAIN GROUP
%         WRITE OUT THE SIGNIFICANT ROI-TO-ROI CORRELATION BETWEEN GROUPS
%dlmwrite('PainRoiTtest_PainPre.txt',OUT_Text_PainPre,'')
whenRun = datestr(now, 'yyyy-mm-dd_HHMM');
% file_id1 = fopen([rootpath 'GpFCdiffs_',datestr(now, 'yyyy-mm-dd_HH:MM'),'.txt'], 'w');
file_id1 = fopen([rootpath 'ROIttest_PainPre', whenRun,'.txt'], 'w');
fprintf(file_id1, OUT_Text_PainPre,'');
fclose(file_id1);

figure;
imagesc(ttest_tval_PainPre);colorbar;colormap('jet');

%% Now that we have identified ROI-ROI correlations that differ between groups
%  We want to test what behavioral variables (if any) contribute to those differences in fxnl connectivity



psqiPlusRoiNames = [psqiNames, 'Y'];

%             RUNNING MULTIPLE LINERAR REGRESSION USING fitlm (WITH A TABLE)
%
%             DECALRE A STRUCT TO HOLD ALL THE RESULTS FROM THE fitlm LOOP
LM_PainPre = struct();
for i=1:numel(I);
     node1 = I(i);
     node2 = J(i);
     Y = mean(squeeze(painzgfcccs(node1, node2, :, 1, :)), 2);
     psqiPlusROI = [psqiData, Y];
     tablePsqiPlusROI = array2table(psqiPlusROI, 'VariableNames',psqiPlusRoiNames);
     LM_PainPre.(strcat('lm', num2str(i))) = fitlm(tablePsqiPlusROI, 'Y~TiB_hrs+SoL_min+WASO_min');
end
% %         OUTPUT THE OVERALL F AND P-VALUES FOR EACH MODEL
for i=1:numel(I);
     LM_PainPre.Model = anova(LM_PainPre.(strcat('lm', num2str(i))), 'summary');
     LM_PainPre.modelSummary(i, :) = LM_PainPre.Model(2,4:5);
end
% %         IDENTIFY SIGNIGCANT MODELS
 LM_PainPre.modelSummary.sig = [LM_PainPre.modelSummary.pValue < 0.05];
LM_PainPre.modelSummary.node1 = I;
LM_PainPre.modelSummary.node2 = J;
for i=1:numel(I);
  node1 = I(i);
  node2 = J(i);
  roi1num = num2str(I(i));
  roi1str = painnames2(I(i),:);
  roi2num = num2str(J(i));
  roi2str = painnames2(J(i),:);
  DUMMY_ROI1{i} = roi1str;
  DUMMY_ROI2{i} = roi2str;
end
LM_PainPre.modelSummary.roi1 = DUMMY_ROI1';
LM_PainPre.modelSummary.roi2 = DUMMY_ROI2';
LM_PainPre.modelSummary.model = (1:numel(I))';
LM_PainPre.modelSummary = LM_PainPre.modelSummary(:, [8 3 2 1 4 5 6 7]);
%
%
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%   END SCRIPT     %%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%   OTHER ANALYSES    %%%%%%%%%%%%%%%%%%%%%%
% Method A:  Simple linear regression
% NOTE: VERY Simple
% Instead, we are going to use multiple regression
% (but keeping this code in case need for later)
% for each ROI-ROI correlation (Y) from above,
%    Test if behavioral variable (X, i.e. a sleep variable) predicts Y

% for i=1:numel(I)
%      node1=I(i);
%      node2=J(i);
%      for behavior_index=3:5 %14 %2:5%numel(psqiNames) % skip first column, which is ID
%           BEHAVIOR = psqiData(:,behavior_index);
%           BEHAVIOR_label=psqiNames{behavior_index};
%           X = [ones(size(BEHAVIOR)) BEHAVIOR];
%           Y = mean(squeeze(painzgfcccs(node1, node2, :, 1, :)), 2);
%           [B,BINT,R,RINT,STATS] = regress(Y,X);
%           % B(1) is beta for constant term; B(2) is beta for behavioral
%           % STATS lists (1) R2, (2) F stat, (3) p-value, (4) error variance
%           % can estimate t-value as square-root of F-stat times sign of B(2) (+ or -)
%
%           % scatterplot
%           % if P<0.05 (uncorrected)
%           %   scatterplot of relationship
%           if STATS(3) < 0.05
%                figure
%                plot(X(:,2),Y,'k+');
%                hold on
%                l = lsline;
%                set(l,'LineWidth',2)
%                xlabel(psqiNames{behavior_index})
%                ylabel(['RS-FC of ROI' num2str(node1) '-ROI' num2str(node2)])
%                R2=STATS(1);
%                xposition = max(X(:,2)) - 0.1*range(X(:,2));
%                yposition = max(Y) - 0.1*range(Y);
%                text(xposition,yposition,['R^2 = ' sprintf('%0.3f',R2)])
%           end
%
%           Regress_P_out(i,behavior_index) = STATS(3);
%           Regress_T_out(i,behavior_index) = sqrt(STATS(2))*sign(B(2));
%           Regress_B_out(i,behavior_index) = B(2);
%           % Note: in line 287, we are looping through variables 3-5
%           % in line 317, we are assigning to columns 3-5
%           % therefore, columns 1 and 2 will be empty (all zeros)
%           %
%           % You can change code so that line 316 assigns to (behavior_index-2), but this could cause problems later
%           % Recommend keep as is
%           %
%           % So each row of Regress_B_out is a different ROI-to-ROI pairs
%           % and each column corresponds to a variable in psqiData (if tested)
%           % For example: we did not test variable #1 (subject ID), so this column is all zeros
%
%      end
% end
%
%
%
%%%%%%%%%%%%%%%%%%%%%%   OTHER ANALYSES    %%%%%%%%%%%%%%%%%%%%%%







%              RUNNING A STEPWISE REGRESSION USING stepwiselm (WITH A TABLE)

%              RUNNING MULTIPLE LINERAR REGRESSION USING fitlm (WITHOUT A TABLE)
%
%              DECALRE A STRUCT TO HOLD ALL THE RESULTS FROM THE fitlm LOOP
% LM = struct();
% for i=1:numel(I);
%      node1 = I(i);
%      node2 = J(i);
%      roi1num = num2str(I(i));
%      roi1str = painnames(I(i),:);
%      roi2num = num2str(J(i));
%      roi2str = painnames(J(i),:);
%      xVars = [psqiData(:, 3:5)];
%      % nodePairs = [nodePairs, roi1str, roi2str, sprintf('\n')];
%      Y = mean(squeeze(painzgfcccs(node1, node2, :, 1, :)), 2);
%      %
%      % LM.(strcat('rois')) = [nodePairs, node1, node2];
%      % % {roi1str, roi2str};
%      LM.(strcat('lm', num2str(i))) = fitlm(xVars,Y);
%
% end
%
% %         OUTPUT THE OVERALL F AND P-VALUES FOR EACH MODEL
% for i=1:numel(I);
% ztbl2 = anova(LM.(strcat('lm', num2str(i))), 'summary');
% modelSummary(i, :) = ztbl2(2,4:5);
% end
%

% zlm = fitlm(tablePsqiPlusROI, 'Y~TiB_hrs+SoL_min+WASO_min')

% tablePSQI = array2table(psqiData, 'VariableNames',psqiNames)





% X2 = [psqiData(:, 3:5)];
% lm= fitlm(X2,Y);
% lm

% Linear regression model:
%     y ~ 1 + x1 + x2 + x3
%
% Estimated Coefficients:
%                     Estimate          SE         tStat     pValue
%                    ___________    __________    _______    _______
%
%     (Intercept)       -0.09716       0.12003    -0.8095    0.42049
%     x1              0.00039373    0.00085143    0.46244    0.64495
%     x2             -0.00031876    0.00029569     -1.078    0.28407
%     x3                0.013969      0.015571    0.89715    0.37217
%
%
% Number of observations: 89, Error degrees of freedom: 85
% Root Mean Squared Error: 0.145
% R-squared: 0.0378,  Adjusted R-Squared 0.00388
% F-statistic vs. constant model: 1.11, p-value = 0.348

% AJ: Consider running stepwise REGRESSION
% [BB,SE,PVAL,INCLUDE,NEXTSTEP,HISTORY]=stepwisefit(X,Y);
% NOTE:  stepwise fit automatically adds a constant column of ones;
%  (so you don't have to...)


% If you want to test relationship only in group 1
%     i1 = group 1 SUBJECTS
% X = [ones(size(BEHAVIOR)) BEHAVIOR];
% Y = mean(squeeze(painzgfcccs(node1, node2, :, 1, :)), 2);
% X = X(i1,:)
% Y = Y(i1,:)

% for i=1:numel(I)
%      node1=I(i);
%      node2=J(i);
%      for behavior_index=3:5 %14 %2:5%numel(psqiNames) % skip first column, which is ID
%           BEHAVIOR = psqiData(:,behavior_index);
%           BEHAVIOR_label=psqiNames{behavior_index};
%           X = [ones(size(BEHAVIOR)) BEHAVIOR];
%           Y = mean(squeeze(painzgfcccs(node1, node2, :, 1, :)), 2);
%            X3 = X(i4,:);        % REGRESSION FOR ONLY GROUP 2
%            Y3 = Y(i4,:);        % REGRESSION FOR ONLY GROUP 2
%           % X3 = X;
%           % Y3 = Y;
%           [B3,BINT3,R3,RINT3,STATS3] = regress(Y3,X3);
%           % B(1) is beta for constant term; B(2) is beta for behavioral
%           % STATS lists (1) R2, (2) F stat, (3) p-value, (4) error variance
%           % can estimate t-value as square-root of F-stat times sign of B(2) (+ or -)
%
%           % scatterplot
%           % if P<0.05 (uncorrected)
%           %   scatterplot of relationship
%           if STATS3(3) < 0.05
%                figure
%                plot(X3(:,2),Y3,'k+');
%                hold on
%                l = lsline;
%                set(l,'LineWidth',2)
%                xlabel(psqiNames{behavior_index})
%                ylabel(['RS-FC of ROI' num2str(node1) '-ROI' num2str(node2)])
%                R2=STATS3(1);
%                xposition = max(X3(:,2)) - 0.1*range(X3(:,2));
%                yposition = max(Y3) - 0.1*range(Y3);
%                text(xposition,yposition,['R^2 = ' sprintf('%0.3f',R2)])
%           end
%
%           Regress_P_out3(i,behavior_index) = STATS3(3);
%           Regress_T_out3(i,behavior_index) = sqrt(STATS3(2))*sign(B3(2));
%           Regress_B_out3(i,behavior_index) = B3(2);
%           % Note: in line 287, we are looping through variables 3-5
%           % in line 317, we are assigning to columns 3-5
%           % therefore, columns 1 and 2 will be empty (all zeros)
%           %
%           % You can change code so that line 316 assigns to (behavior_index-2), but this could cause problems later
%           % Recommend keep as is
%           %
%           % So each row of Regress_B_out is a different ROI-to-ROI pairs
%           % and each column corresponds to a variable in psqiData (if tested)
%           % For example: we did not test variable #1 (subject ID), so this column is all zeros
%
%      end
% end


% Method B:  multiple linear regression
% As simple, but reading in all behaviors

% psqiNames = {'ID','TiB_hrs', 'SoL_min','WASO_min','TST_hrs','SleepEfficiency','psqi_Durat','psqi_Distb', ...
%               'psqi_Latency','psqi_DayDys','psqi_SE','psqi_BadSQual','psqi_Meds','PSQI_total'}
% for i=1:numel(I)
%      node1=I(i);
%      node2=J(i);
%      Y = mean(squeeze(painzgfcccs(node1, node2, :, 1, :)), 2);
%      % X = [psqiData(:,2:end)]; % note: only doing 2-5 due to multicollinarity - address later
%      X = [psqiData(:,2:5)];
%      lm= fitlm(X,Y);
%      % B(1) is beta for constant term; B(2) is beta for behavioral
%      % STATS lists (1) R2, (2) F stat, (3) p-value, (4) error variance
%      % can estimate t-value as square-root of F-stat times sign of B(2) (+ or -)
%
%      % scatterplot
%      % if P<0.05 (uncorrected)
%      %   scatterplot of relationship
% end



%%
%         THERE IS A SIGNIGCANT MAIN EFFECT FOR GROUP
%         BOTH PAIN GROUPS DIFFER FROM THE HC GROUP
%         NEITHER PAIN GROUP DIFFERES FROM EACH OTHER
%         GOING FORWARD, COLLAPSING ACROSS PAIN GROUPS
%%











%
%
% %%             SYNTAX FOR PERFORMING A MANOVA
% %              START BY CREATING A VECTOR REPRESENTING ALL THE GROUPS
% gp1 = ones(length(i1),1);
% gp2 = 2*ones(length(i2),1);
% gp3 = 3*ones(length(i3),1);
% gps123 = cat(1,gp1,gp2,gp3);
% groups = nominal(gps123);               % SPECIFY THIS AS AN ORDINAL VARIABLE
% groups = setlabels(groups,{'HC','CLBP','FM'});    % SET THE VARIABLE LABELS
% prePain123 = cat(1,prePainHC, prePainCLBP,prePainFM);  % CREATE ANOTHER VECTOR TO INCLUDE
% %              RUN THE MANOVA
% [d,p,stats] = manova1(prePain123,groups)
%

% as a one-sample t-test
%
%
%
%
%
%
%
%                   END SCRIPT
% for node1 = 1:16
%     for node2 = 1:16

%         THIS CREATES A MEAN CORRELATION MATRIX FOR THE PAIN ROIS OF EACH GROUP
% for node1 = 1:length(pain)
%      for node2 = 1:length(pain)
%         prePainHC = mean(squeeze(painzgfcccs(node1, node2, i1, 1, :)), 2);
%         prePainCLBP = mean(squeeze(painzgfcccs(node1, node2, i2, 1, :)), 2);
%         prePainFM = mean(squeeze(painzgfcccs(node1, node2, i3, 1, :)), 2);
% %
% %         % place the code between lines 116 and 131 here
%         gpNames = {'HC','CLBP','FM'};               % VARIABLE OF GROUP NAMES
% %       STEP 2 = CREATE AN ARRAY OF THE COMBINED VARIABLES FROM ABOVE
% %       THE ARRAY NEEDS TO BE PADDED BECAUSE OF UNEVEN GROUP SIZES
% %       IDENTIFY THE LARGEST GROUP
%         A = max([length(i1),length(i2),length(i3)]);
%         A = zeros(A,3);    % INITIALIZE ARRAY OF ALL ZEROS  FOR LARGEST GROUP
%         A(A == 0) = NaN;    % CONVERT ALL '0' TO 'NaN' (MISSING VALUES)
%         A(1:length(prePainHC),1) = prePainHC;       % HC TO COLUMN 1
%         A(1:length(prePainCLBP),2) = prePainCLBP;   % CLBP TO COLUMN 2
%         A(1:length(prePainFM),3) = prePainFM;       % FM TO COLUMN 3
% %       STEP 3 = RUNNING THE ANOVA AND MULTIPLE COMPARISONS
% %       CREATE A TABLE OF OVERALL F-TEST
%         [p,tbl,stats] = anova1(A,gpNames, 'off');      % TABLE OF OVERALL RESULTS
% %         ftestNames = tbl(1,:);                  % VARIABLE NAMES FOR THE TABLE
% %         ftestNames{1,6} = 'Prob_F';             % FIX THE SYMBOL ISSUE
% %         tableFtest = array2table(tbl(2:4,:),'VariableNames',ftestNames);
%         pain_anovaresults_effect(node1, node2) = cell2mat(tbl(2,5)); % THESE ARE F-VALUES
%         pain_anovaresults_pvalue(node1, node2) = p;
%     end
% end


%painzgfcccs = zgfcccs(1:16, 1:16, :, :, :);
%dmnzgfcccs = zgfcccs(17:22, 17:22, :, :, :);


%%
%    THESE ARE THE CCs FOR ALL 90 SUBJECTS ACROSS THE 16 PAIN REGIONS
%    painnet = squeeze(sum(sum(painzgfcccs, 1), 2)) ./ (16 * 15);
%    THESE ARE THE CCs FOR ALL 90 SUBJECTS ACROSS THE 6 DMN REGIONS
%    dmnnet = squeeze(sum(sum(dmnzgfcccs, 1), 2)) ./ (6 * 5);


%get left amygdala (region 2) to left anterior insula (region 8) from pain network
%pain_lamyg_2_linsula = squeeze(painzgfcccs(2, 8, :, :, :));


%pain_anovaresults_effect = zeros(16, 16);
%pain_anovaresults_pvalue = zeros(16, 16);
