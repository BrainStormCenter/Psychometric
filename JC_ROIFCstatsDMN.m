%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%		CREATED BY:		JOCHEN WEBER
%		CREATED ON:		2017-12-13
%         MODIFIED WITH ANDY  2018_04_21
%
%		USAGE:			TESTING ROI CORRELATIONS ACROSS GROUPS
%
%         LATEST MODIFICATION:     2018_05_08
%
%         AS OF May 4, 2018, THERE ARE ISSUES WITH THE PSQI DATA
%              SUB 161 = WASO SCORE OF -30
%              SUB 172 = SLEEP EFFICIENCY SCORE GREATER THAN 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
painvoi = find(~cellfun('isempty', regexpi(voinames, '^Pain')));
dmnvoi = find(~cellfun('isempty', regexpi(voinames, '^DMN')));
bothvoi = find(~cellfun('isempty', regexpi(voinames, '^Both')));
%voiorder = [pain; dmn];
voiorder = [painvoi; dmnvoi; bothvoi];
nvs = numel(voiorder);

%    DETERMINING THE NUMBER OF ROIS FOR EACH NETWORK
painstart = 1;
painend = length(painvoi);
dmnstart = painend + 1;
dmnend = length(painvoi) + length(dmnvoi);
bothstart = dmnend + 1;
bothend = length(bothvoi) + length(painvoi) + length(dmnvoi);

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

%              find subjects in three groups
%              COLUMN 3 IDENTIFIES SUBJECT GROUP
%              1=HC, 2=CLBP, 3=FM
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
rlistd2 = glistd(:,[1 3:11]); % KEEP SUBJECT NUMBERS

%         GROUP NAMES
gpNames = {'HC','CLBP','FM'}; % ALL THREE GROUPS
gpNames2 = {'HC','Pain'};     % GROUP NAMES COLLAPSED ACROSS CP GROUPS

%                             CREATE CC ARRAYS
%         THESE ARE THE CROSS CORRELATIONS OF ALL THE BRAIN REGIONS LISTED IN THE VOI FILE
afcccs = cat(3, fcccs{:});
gfcccs = reshape(afcccs(voiorder, voiorder, rlistd(:)), [nvs, nvs, ns, 2, 2, 2]);
%         fisher transform
zgfcccs = n.fisherr2z(gfcccs);
zgfcccs(isinf(zgfcccs)) = 0;
%         average over first and second run of each half-session
%         THESE ARE THE 'PRE' MANIPULATION RESTING STATE SCANS
zgfcccs = squeeze(mean(zgfcccs, 4));

%    this leaves 5-dimensions
%    THE 5 DIMENSIONS FOR THE zgfcccs ARRAY ARE:
%    1 (ROIs IN THE PAIN NETWORK (i.e., REGIONS [1-16], AS OF February 12, 2018)
%    2 (ROIs IN THE DMN NETWORK (i.e., REGIONS [17-22], AS OF February 12, 2018)
%    3 subjects (in order of groups, 1-31 HC, 32-73 CLBP, 74-90 FM)
%    4 pre/post treatment (1 pre, 2 post)
%    5 neg/pos session (1 neg, 2 pos)
%         to unpack:
%         - i1 and i2 and i3 are the indices for groups HC and CLBP and FM
%         - the next ", 1" is the "pre" (treatment) selection
%         - the next ", 1" is the "neg session" selection
%
%         to split into within network matrices
%         THESE ARE ARRAYS OF CROSS CORRELATIONS AMONG REGIONS IN EACH NETWORK
%         THE ARRAYS ARE ORGANIZED AS (REGIONS^REGIONS, ALL SUBJECTS, PRE & POST, POS & NEG)
%
%         THESE ARE THE CCs MATRICIES OF BRAIN REGIONS FOR EACH SUBJECT IN EACH ROI SET
bothzgfcccs = zgfcccs(bothstart:bothend, bothstart:bothend, :, :, :);
painzgfcccs = zgfcccs(painstart:painend, painstart:painend, :, :, :);
dmnzgfcccs = zgfcccs(dmnstart:dmnend, dmnstart:dmnend, :, :, :);

%              average connectivity strengths
%              THESE ARE THE AVERAGE CCs FOR ALL THE SUBJECTS ACROSS THE PAIN REGIONS
painnet = squeeze(sum(sum(painzgfcccs, 1), 2)) ./ (length(painvoi) * (length(painvoi) -1));
%              THESE ARE THE CCs FOR ALL THE SUBJECTS ACROSS THE DMN REGIONS
dmnnet = squeeze(sum(sum(dmnzgfcccs, 1), 2)) ./ (length(dmnvoi) * (length(dmnvoi) -1));
%
%              CREATE MATRIX OF PAIN CC'S AND BEHAVIORAL DATA
%              RESHAPE THE PAINNET CC'S, COLUMNS 1&2=PRE/POST; COLUMNS 3&4=NEG/POS
sub_by_painCCs = [painnet(:,1:2,1),painnet(:,1:2,2)];
psqiData = glistd(:,[1 3 12:24]);
psqiNames = {'ID','GP','TiB_hrs', 'SoL_min','WASO_min','TST_hrs','SleepEfficiency','psqi_Durat','psqi_Distb', ...
               'psqi_Latency','psqi_DayDys','psqi_SE','psqi_BadSQual','psqi_Meds','PSQI_total'};

PSQIandPainCCs = [psqiData,sub_by_painCCs];  % WHY DID I MAKE THIS?
%%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   THE ANOVA LOOKING AT THE OVERALL MAIN EFFECT OF GROUP
%                   CCS AVERAGED ACROSS ALL PAIN ROIS; FOR THE PRE-MOOD MANIPULATION RESTING STATE SCANS
%                   FINISHED May 1, 2018; CODE BLOCK = LINES 145-188
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              DECLARE ARRAY TO STORE ANOVA RESULTS
Pain.Pre.n1_Gp_AvgPain = struct();
%              STEP 1 = CREATE VARIABLES OF THE MEAN CORRELATION OF ALL PAIN REGIONS
%                   FOR EACH GROUP OF THE PRE SCANS ACROSS BOTH VISITS
fcccs_PainPreAvg_HC = mean(painnet(i1,1,:),3);        % MEAN OF HC
fcccs_PainPreAvg_CLBP = mean(painnet(i2,1,:),3);      % MEAN OF CLBP
fcccs_PainPreAvg_FM = mean(painnet(i3,1,:),3);        % MEAN OF FM
%              STEP 2 = CREATE AN ARRAY OF THE COMBINED VARIABLES FROM ABOVE
%                   THE ARRAY NEEDS TO BE PADDED BECAUSE OF UNEVEN GROUP SIZES
%                   IDENTIFY THE LARGEST GROUP
A = max([numel(i1),numel(i2),numel(i3)]);
A = zeros(A,3);    % INITIALIZE ARRAY OF ALL ZEROS  FOR LARGEST GROUP
A(A == 0) = NaN;    % CONVERT ALL '0' TO 'NaN' (MISSING VALUES)
A(1:numel(fcccs_PainPreAvg_HC),1) = fcccs_PainPreAvg_HC;       % HC TO COLUMN 1
A(1:numel(fcccs_PainPreAvg_CLBP),2) = fcccs_PainPreAvg_CLBP;   % CLBP TO COLUMN 2
A(1:numel(fcccs_PainPreAvg_FM),3) = fcccs_PainPreAvg_FM;       % FM TO COLUMN 3
%              STEP 3 = RUNNING THE ANOVA AND MULTIPLE COMPARISONS
%                   CREATE A TABLE OF OVERALL F-TEST
[Pain.Pre.n1_Gp_AvgPain.p_PainPre,Pain.Pre.n1_Gp_AvgPain.modelSummary_PainPre,Pain.Pre.n1_Gp_AvgPain.stats_PainPre] ...
     = anova1(A,gpNames);
%              THIS CONVERTS THE ABOVE RESULTS TO TABLE FORMAT
Pain.Pre.n1_Gp_AvgPain.ftest_tblHdr = Pain.Pre.n1_Gp_AvgPain.modelSummary_PainPre(1,:);    % VARIABLE NAMES FOR THE TABLE
Pain.Pre.n1_Gp_AvgPain.ftest_tblHdr{1,6} = 'Prob_F';             % FIX THE SYMBOL ISSUE
Pain.Pre.n1_Gp_AvgPain.overallFtest = array2table(Pain.Pre.n1_Gp_AvgPain.modelSummary_PainPre(2:4,:), ...
     'VariableNames',Pain.Pre.n1_Gp_AvgPain.ftest_tblHdr);
figure;
%              POST-HOC GROUP COMPARISONS
[Pain.Pre.n1_Gp_AvgPain.multcompare.c,~,~,Pain.Pre.n1_Gp_AvgPain.multcompare.gnames] = ...
     multcompare(Pain.Pre.n1_Gp_AvgPain.stats_PainPre);    % EVALUATE MULTIPLE COMPARISONS
%              STEP 4 = PREPARING POST-HOC STATS OUTPUT
%                   CREATE AN ARRAY OF ANOVA STATISTICAL OUTPUT
Pain.Pre.n1_Gp_AvgPain.multcompare.multout_PainPre = ...
     [Pain.Pre.n1_Gp_AvgPain.multcompare.gnames(Pain.Pre.n1_Gp_AvgPain.multcompare.c(:,1)), ...
     Pain.Pre.n1_Gp_AvgPain.multcompare.gnames(Pain.Pre.n1_Gp_AvgPain.multcompare.c(:,2)), ...
     num2cell(Pain.Pre.n1_Gp_AvgPain.multcompare.c(:,3:6))];
%              INITIAL ORDER OF OUTPUT FROM THE MULTICOMPARISON STEP
%                   COLUMNS 1-6 =  {'gp1','gp2','lCI','gpDiff','uCI','pval'}
%              CHANGING THE VARIABLE ORDER IN THE OUTPUT ARRAY TO
%                   COLUMNS 1-6 = {'gp1','gp2', 'pval','gpDiff','lCI','uCI'}) AND THEN
%              CREATE A TABLE OF THE MULTICOMPARISON OUTPUT
Pain.Pre.n1_Gp_AvgPain.multcompare.multout_PainPre = Pain.Pre.n1_Gp_AvgPain.multcompare.multout_PainPre(:,[1 2 6 4 3 5]);
Pain.Pre.n1_Gp_AvgPain.multcompare.multoutTbl_PainPre = ...
     array2table(Pain.Pre.n1_Gp_AvgPain.multcompare.multout_PainPre, ...
     'VariableNames',{'gp1','gp2', 'pval','gpDiff','lCI','uCI'});
%
%%
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   ANOVAS TESTING FOR A MAIN EFFECT OF GROUP ON ALL PAIN ROI PAIRS
%                   FOR THE PRE CONDITION ONLY
%                   FINISHED May 2, 2018; CODE BLOCK = LINES 196-273
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              DECLARE A STRUCT TO HOLD THE RESULTS
Pain.Pre.n2_Gp_RoiPairPain = struct();
%              INITIALIZE ARRAY OF ALL ZEROS TO HOLD STAT VALUES
Pain.Pre.n2_Gp_RoiPairPain.anova_fvals = zeros(length(painvoi),length(painvoi));
Pain.Pre.n2_Gp_RoiPairPain.anova_pvals = zeros(length(painvoi),length(painvoi));
%              THIS LOOP COMPUTES AN ANOVA FOR ALL ROI-TO-ROI PAIRS
%              THIS IS FOR THE PRE-MOOD MANIPULATION SCANS ONLY
for node1 = 1:length(painvoi)         % PAIN REGION #1
     for node2 = 1:length(painvoi)    % PAIN REGION #2
        fcccs_PainPreRoi_HC = mean(squeeze(painzgfcccs(node1, node2, i1, 1, :)), 2);
        fcccs_PainPreRoi_CLBP = mean(squeeze(painzgfcccs(node1, node2, i2, 1, :)), 2);
        fcccs_PainPreRoi_FM = mean(squeeze(painzgfcccs(node1, node2, i3, 1, :)), 2);
        fcccs_PainPreRoi_PainGps = mean(squeeze(painzgfcccs(node1, node2, i4, 1, :)), 2);     % BOTH CP GROUPS
        gpNames = {'HC','CLBP','FM'};           % VARIABLE: GROUP NAMES
        roi1str = painnames3(node1,:);
        roi2str = painnames3(node2,:);
%              STEP 2 = CREATE AN ARRAY OF THE COMBINED VARIABLES FROM ABOVE
%                   THE ARRAY NEEDS TO BE PADDED BECAUSE OF UNEVEN GROUP SIZES
%                   IDENTIFY THE LARGEST GROUP
        A = max([numel(i1),numel(i2),numel(i3)]);
        A = zeros(A,3);       % INITIALIZE ARRAY OF ALL ZEROS  FOR LARGEST GROUP
        A(A == 0) = NaN;      % CONVERT ALL '0' TO 'NaN' (MISSING VALUES)
        A(1:numel(fcccs_PainPreRoi_HC),1) = fcccs_PainPreRoi_HC;       % HC TO COLUMN 1
        A(1:numel(fcccs_PainPreRoi_CLBP),2) = fcccs_PainPreRoi_CLBP;   % CLBP TO COLUMN 2
        A(1:numel(fcccs_PainPreRoi_FM),3) = fcccs_PainPreRoi_FM;       % FM TO COLUMN 3
%              STEP 3 = RUNNING THE ANOVA AND MULTIPLE COMPARISONS
%                   CREATE A TABLE OF OVERALL F-TEST
        [Pain.Pre.n2_Gp_RoiPairPain.p,Pain.Pre.n2_Gp_RoiPairPain.tbl,Pain.Pre.n2_Gp_RoiPairPain.stats] = anova1(A,gpNames, 'off');      % TABLE OF OVERALL RESULTS
        Pain.Pre.n2_Gp_RoiPairPain.ftestHdr = Pain.Pre.n2_Gp_RoiPairPain.tbl(1,:);        % VARIABLE NAMES FOR THE TABLE
        Pain.Pre.n2_Gp_RoiPairPain.ftestHdr{1,6} = 'Prob_F';                    % FIX THE SYMBOL ISSUE
        Pain.Pre.n2_Gp_RoiPairPain.ftestTable = array2table(Pain.Pre.n2_Gp_RoiPairPain.tbl(2:4,:),'VariableNames',Pain.Pre.n2_Gp_RoiPairPain.ftestHdr);
        Pain.Pre.n2_Gp_RoiPairPain.anova_fvals(node1, node2) = cell2mat(Pain.Pre.n2_Gp_RoiPairPain.tbl(2,5)); % THESE ARE F-VALUES
        Pain.Pre.n2_Gp_RoiPairPain.anova_pvals(node1, node2) = Pain.Pre.n2_Gp_RoiPairPain.p;
    end
end
%              THIS LOOP REMOVES THE UPPER DIAGONAL PVALS AND FVALS MATRICES
for i = 1:length(painvoi)
     for j=i+1:length(painvoi)
          Pain.Pre.n2_Gp_RoiPairPain.anova_pvals(i,j)=NaN;
          Pain.Pre.n2_Gp_RoiPairPain.anova_pvals(i,j)=NaN;
          Pain.Pre.n2_Gp_RoiPairPain.anova_fvals(i,j)=NaN;
          Pain.Pre.n2_Gp_RoiPairPain.anova_fvals(i,j)=NaN;
     end
end
%              IDENTIFY THE ROI-TO-ROI PAIRS WITH SIGNIFICANT GROUP DIFFERENCES IN CROSS CORRECTIONS
%              SIGNIFICANT PAIRS ARE IDENTIFIED USING THE P-VALUE SPECIFIED BELOW
sigpval = 0.005;
[I,J] = find(Pain.Pre.n2_Gp_RoiPairPain.anova_pvals < sigpval);
[I J]          % THE SIGNIFICANT ROI PAIRS IS SENT TO THE SCREEN
Pain.Pre.n2_Gp_sigRoiPairPain = [I,J];
%
%              INITIALIZE A VARIABLE TO STORE THE SIGNIFICANT PAIRS
Pain.Pre.n2_Gp_RoiPairPain.n2_Gp_sigRoiPairPain = [];
%              Now, add strings to that OUT_TEXT
%              and you can use the "sprintf" command for syntax like tabs, line breaks
%              ADD A HEADER TO EXPLAIN WHAT INFORMATION IS BEING STORED
Pain.Pre.n2_Gp_RoiPairPain.n2_Gp_sigRoiPairPain = [Pain.Pre.n2_Gp_RoiPairPain.n2_Gp_sigRoiPairPain ...
     'There are significant group difference in these ROI-to_ROI CCs:' sprintf('\t') sprintf('\n')];
%              THIS LOOP WRITES OUT THE SIGNIFICANT ROI PAIRS IDENTIFIED ABOVE TO THE NEW VARIABLE
for i=1:numel(I)
     pairNum = num2str(i);
     roi1num = num2str(I(i));
     roi1str = painnames3(I(i),:);
     roi2num = num2str(J(i));
     roi2str = painnames3(J(i),:);
     this_pval1 = Pain.Pre.n2_Gp_RoiPairPain.anova_pvals(I(i),J(i));
     this_fval1 = Pain.Pre.n2_Gp_RoiPairPain.anova_fvals(I(i),J(i));
     Pain.Pre.n2_Gp_RoiPairPain.n2_Gp_sigRoiPairPain = [Pain.Pre.n2_Gp_RoiPairPain.n2_Gp_sigRoiPairPain pairNum, '. ' roi1str 'with ' roi2str ...
          sprintf('\t') '(#' roi1num ') <-->' ' (#' roi2num ')  ' sprintf('\t') ...
          'f-val: ' num2str(this_fval1) sprintf('\t') 'p-val: ' sprintf('%0.04f',this_pval1) sprintf('\n')];
end
%              THIS PRINTS THE OUTPUT FROM THE LOOP ABOVE TO THE SCREEN
Pain.Pre.n2_Gp_RoiPairPain.n2_Gp_sigRoiPairPain
%              THIS CREATES A TEXT FILE WITH THE SAME OUTPUT AS ABOVE
whenRun = datestr(now, 'yyyy-mm-dd_HHMM');
file_id1 = fopen([rootpath 'Sig_Pain_Pre_RoiPairs', whenRun,'.txt'], 'w');
fprintf(file_id1, Pain.Pre.n2_Gp_RoiPairPain.n2_Gp_sigRoiPairPain,'');
fclose(file_id1);
%
%%
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   FOLLOW UP POST HOC TESTING FOR THE SIGNIFICANT ROI PAIRS IDENTIFIED ABOVE
%                   FOR THE PRE CONDITION ONLY
%                   FINISHED May 2, 2018; CODE BLOCK = LINES 280-327
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              DECLARE A STRUCT TO HOLD THE RESULTS
Pain.Pre.n3_Gp_pHocRoiPairPain = struct();
for i=1:numel(I);
     node1 = I(i);
     node2 = J(i);
     fcccs_PainPreRoi_HC = mean(squeeze(painzgfcccs(node1, node2, i1, 1, :)), 2);
     fcccs_PainPreRoi_CLBP = mean(squeeze(painzgfcccs(node1, node2, i2, 1, :)), 2);
     fcccs_PainPreRoi_FM = mean(squeeze(painzgfcccs(node1, node2, i3, 1, :)), 2);
     gpNames = {'HC','CLBP','FM'};           % VARIABLE: GROUP NAMES
%              STEP 2 = CREATE AN ARRAY OF THE COMBINED VARIABLES FROM ABOVE
%                   THE ARRAY NEEDS TO BE PADDED BECAUSE OF UNEVEN GROUP SIZES
%                   IDENTIFY THE LARGEST GROUP
     A = max([numel(i1),numel(i2),numel(i3)]);
     A = zeros(A,3);       % INITIALIZE ARRAY OF ALL ZEROS  FOR LARGEST GROUP
     A(A == 0) = NaN;      % CONVERT ALL '0' TO 'NaN' (MISSING VALUES)
     A(1:numel(fcccs_PainPreRoi_HC),1) = fcccs_PainPreRoi_HC;       % HC TO COLUMN 1
     A(1:numel(fcccs_PainPreRoi_CLBP),2) = fcccs_PainPreRoi_CLBP;   % CLBP TO COLUMN 2
     A(1:numel(fcccs_PainPreRoi_FM),3) = fcccs_PainPreRoi_FM;       % FM TO COLUMN 3
%                   STEP 3 = RUNNING THE ANOVA AND MULTIPLE COMPARISONS
%                        CREATE A TABLE OF OVERALL F-TEST
     [Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).model_p,Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).model_tbl, ...
          Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).model_stats] = anova1(A,gpNames, 'off');  % TABLE OF OVERALL RESULTS
     %              THIS CONVERTS THE ABOVE RESULTS TO TABLE FORMAT
     Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).ftest_tblHdr = Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).model_tbl(1,:);  % VARIABLE NAMES FOR THE TABLE
     Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).ftest_tblHdr{1,6} = 'Prob_F';    % FIX THE SYMBOL ISSUE
     Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).Ftable = ...
          array2table(Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).model_tbl(2:4,:), ...
          'VariableNames', Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).ftest_tblHdr);
%                   THESE TWO LINES PROVIDE SOMEWHAT REDUNDANT INFORMATION
     Pain.Pre.n3_Gp_pHocRoiPairPain.ph_Fvals2(node1, node2) = cell2mat(Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).model_tbl(2,5)); % THESE ARE F-VALUES
     Pain.Pre.n3_Gp_pHocRoiPairPain.ph_Pvals2(node1, node2) = cell2mat(Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).model_tbl(2,6)); % THESE ARE P-VALUES
%                   THE NEXT 3 LINES PERFORM THE POST-HOC GROUP COMPARISONS USING: multcompare - ON THE 3RD LINE
%                   OUTPUT: c=MATRIX OF RESULTS, h=FIGURE HANDLE, m=GROUP MEAN and STD ERROR
     [Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).c,Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).m, ...
          Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).h,Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).gpNames] ...
          = multcompare(Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).model_stats);%, 'CType','bonferroni'); % POTENTIAL CORRECTION
%                   THE NEXT 3 SEGMENTS CLEAN UP THE OUTPUT AND IMPROVES READIBILITY
     Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).gpPostHoc = ...
          [Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).gpNames(Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).c(:,1)), ...
          Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).gpNames(Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).c(:,2)), ...
          num2cell(Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).c(:,3:6))];

     Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).gpPostHoc = ...
          Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).gpPostHoc(:,[1 2 6 4 3 5]);

     Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).gpPostHocTable = ...
          array2table(Pain.Pre.n3_Gp_pHocRoiPairPain.(strcat('posthoc', num2str(i))).gpPostHoc, ...
          'VariableNames',{'gp1','gp2', 'pval','gpDiff','lCI','uCI'});
end
%
%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   THE FIRST ANOVA (PAIN-PRE) FOUND NO SIGNIFICANT DIFFERENCES BETWEEN THE PAIN GROUPS
%                   THE PAIN GROUPS WERE COLLAPSED FOR THE FOLLOWING ANALYSES (PAIN, PRE-MOOD MANIPULATION)
%                   T-TEST COMPARING GROUPS ON THE REGION-TO-REGION CROSS-CORRELATIONS
%                   FINISHED May 2, 2018; CODE BLOCK = LINES 335-390
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              DECLARE STRUCT TO HOLD T-TEST RESULTS
Pain.Pre.n4_Gp_ttestRoiPairPain = struct();
for node1 = 1:numel(painvoi)         % PAIN REGION #1
     for node2 = 1:numel(painvoi)    % PAIN REGION #2
          fcccs_PainPreRoi_HC = mean(squeeze(painzgfcccs(node1, node2, i1, 1, :)), 2);         % HC GROUP
          fcccs_PainPreRoi_PainGps = mean(squeeze(painzgfcccs(node1, node2, i4, 1, :)), 2);     % BOTH CP GROUPS
          gpNames2 = {'HC','Pain'};               % VARIABLE: GROUP NAMES COLLAPSED ACROSS CP GROUPS
          [Pain.Pre.n4_Gp_ttestRoiPairPain.h,Pain.Pre.n4_Gp_ttestRoiPairPain.p,Pain.Pre.n4_Gp_ttestRoiPairPain.CI,Pain.Pre.n4_Gp_ttestRoiPairPain.stats] ...
               = ttest2(fcccs_PainPreRoi_HC,fcccs_PainPreRoi_PainGps);               % 2 SAMPLE T-TEST
          Pain.Pre.n4_Gp_ttestRoiPairPain.ttest_tval(node1, node2) = Pain.Pre.n4_Gp_ttestRoiPairPain.stats.tstat;      % tstat
          Pain.Pre.n4_Gp_ttestRoiPairPain.ttest_pval(node1, node2) = Pain.Pre.n4_Gp_ttestRoiPairPain.p;                % pvalue
    end
end
%              REMOVE THE UPPER DIAGONAL FROM THE MATRICES BELOW
for i = 1:24
     for j=i+1:24
          Pain.Pre.n4_Gp_ttestRoiPairPain.ttest_tval(i,j)=NaN;
          Pain.Pre.n4_Gp_ttestRoiPairPain.ttest_pval(i,j)=NaN;
     end
end
%              THESE IDENTIFY GROUP DIFFERENCES IN ROI-TO-ROI FUNCTIONAL CONNECTIVITY
%              FDR CORRECTION
pid = FDR(Pain.Pre.n4_Gp_ttestRoiPairPain.ttest_pval,.05);
[I,J] = find(Pain.Pre.n4_Gp_ttestRoiPairPain.ttest_pval <= pid);
[I J]
Pain.Pre.n4_Gp_sigRoiPairPain = [I,J];
%%             FANCY CODE TO PRINT-OUT RESULTS
%              FIRST, CREATE AN EMPTY ARRAY OF CHARACTER STRINGS
Pain.Pre.n4_Gp_ttestRoiPairPain.ttest_results = [];
%              NOW, ADD STRINGS TO THAT ARRAY
%              AND YOU CAN USE THE "SPRINTF" COMMAND FOR SYNTAX LIKE TABS, LINE BREAKS
Pain.Pre.n4_Gp_ttestRoiPairPain.ttest_results = [Pain.Pre.n4_Gp_ttestRoiPairPain.ttest_results 'Significant t-test of pain ROIs, pre mood change, ' ...
     'in the CCs of these ROI-pairs:' sprintf('\t') sprintf('\n')];
%              THIS LOOP WRITES OUT THE SIGNIFICANT ROI PAIRS IDENTIFIED BY THE T-TEST ABOVE
for i=1:numel(I)
     pairNum = num2str(i);
     roi1num = num2str(I(i));
     roi1str = painnames3(I(i),:);
     roi2num = num2str(J(i));
     roi2str = painnames3(J(i),:);
     Pain.Pre.n4_Gp_ttestRoiPairPain.ttest_FDRtval = Pain.Pre.n4_Gp_ttestRoiPairPain.ttest_tval(I(i),J(i));
     Pain.Pre.n4_Gp_ttestRoiPairPain.ttest_FDRpval = Pain.Pre.n4_Gp_ttestRoiPairPain.ttest_pval(I(i),J(i));
     Pain.Pre.n4_Gp_ttestRoiPairPain.ttest_results = [Pain.Pre.n4_Gp_ttestRoiPairPain.ttest_results pairNum, '. ' roi1str, 'with ' roi2str ...
      sprintf('\t') '(#' roi1num ') <-->' ' (#' roi2num ')  ' sprintf('\t') ...
      't-val: ', num2str(Pain.Pre.n4_Gp_ttestRoiPairPain.ttest_FDRtval), sprintf('\t'), ...
      'p-val: ' sprintf('%0.04f',Pain.Pre.n4_Gp_ttestRoiPairPain.ttest_FDRpval) sprintf('\n')];
end
%              SEND THE OUTPUT TO THE SCREEN
Pain.Pre.n4_Gp_ttestRoiPairPain.ttest_results
%              NEGATIVE T-VALUES = HEALTHY CONTROL LESS THAN THE COMBINED PAIN GROUPS
figure;
imagesc(Pain.Pre.n4_Gp_ttestRoiPairPain.ttest_tval);colorbar;colormap(parula);

%              WRITE OUT THE SIGNIFICANT ROI-TO-ROI CORRELATION BETWEEN GROUPS
whenRun = datestr(now, 'yyyy-mm-dd_HHMM');
file_id2 = fopen([rootpath 'Sig_Pain_Pre_ROIttest', whenRun,'.txt'], 'w');
fprintf(file_id2, Pain.Pre.n4_Gp_ttestRoiPairPain.ttest_results,'');
fclose(file_id2);
%
%%
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   NOW THAT WE HAVE IDENTIFIED ROI-ROI CORRELATIONS THAT DIFFER BETWEEN GROUPS VIA T-TEST
%                   WE WILL USE REGRESSION ANALYSES TO SEE WHICH, IF ANY, BEHAVIORAL VARIABLES SIGNIFICANLY
%                   CONTRIBUTE TO GROUP DIFFERENCES IN THE FUNCTIONAL CONNECTIVITY BETWEEN THE ROI PAIRS
%                   FINISHED May 3, 2018; CODE BLOCK = LINES 400-446
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              CREATE A HEADER FOR THE PSQI DATA MATRIX AND 'Y'
%              'Y' WILL BE THE ROI-TO-ROI CC BEING PREDICTED IN THE ANALYSES BELOW
%              THE FOLLOWING ARE MULTIPLE LINEAR REGRESSION USING fitlm (WITH A TABLE)
%              CREATE A HEADER FOR THE PSQI DATA MATRIX AND 'Y'
%              'Y' WILL BE THE ROI-TO-ROI CC BEING PREDICTED IN THE ANALYSES BELOW
%              THE FOLLOWING ARE MULTIPLE LINEAR REGRESSION USING fitlm (WITH A TABLE)
psqiPlusRoiNames = [psqiNames, 'Y'];
%              DECALRE A STRUCT TO HOLD ALL THE RESULTS FROM THE REGRESSIONS BELOW
Pain.Pre.n4_Gp_regrssRoiPairPain = struct();
for i=1:numel(I);
     node1 = I(i);
     node2 = J(i);
     roi1str = painnames2(I(i),:);
     roi2str = painnames2(J(i),:);
     Ypair = strcat(roi1str,'_with_',roi2str);
     Y = mean(squeeze(painzgfcccs(node1, node2, :, 1, :)), 2);
     psqiPlusROI = [psqiData, Y];
     tblPsqi_ROI = array2table(psqiPlusROI, 'VariableNames',psqiPlusRoiNames);
     tblPsqi_ROI.GP = nominal(tblPsqi_ROI.GP);
     tblPsqi_ROIreduced = tblPsqi_ROI(:,{'TiB_hrs','SoL_min','WASO_min','PSQI_total','GP','Y'});
%              LINEAR REGRESSION (fitlm) WITH SPECIFIC INTERACTION TERMS
     Pain.Pre.n4_Gp_regrssRoiPairPain.(strcat('LMpair_', num2str(i))) = fitlm(tblPsqi_ROIreduced, ...
          'Y~TiB_hrs+SoL_min+WASO_min+PSQI_total+TiB_hrs*GP+SoL_min*GP+WASO_min*GP+PSQI_total*GP');
          %    STEPWISE LINEAR REGRESSION (stepwiselm) WITH SPECIFIC INTERACTION TERMS
     Pain.Pre.n4_Gp_regrssRoiPairPain.(strcat('SLMpair_', num2str(i))) = stepwiselm(tblPsqi_ROIreduced, ...
          'Y~TiB_hrs+SoL_min+WASO_min+PSQI_total+TiB_hrs*GP+SoL_min*GP+WASO_min*GP+PSQI_total*GP');
%              INTERCEPT ONLY MODEL (stepwiselm)
     Pain.Pre.n4_Gp_regrssRoiPairPain.(strcat('SLMintercept_', num2str(i))) = stepwiselm(tblPsqi_ROIreduced, ...
          'constant', 'ResponseVar', 'Y');
%              INTERACTION MODEL (stepwiselm)
     Pain.Pre.n4_Gp_regrssRoiPairPain.(strcat('SLMinteraction_', num2str(i))) = stepwiselm(tblPsqi_ROIreduced, ...
          'interactions', 'ResponseVar', 'Y');
%              QUADRATIC MODEL (stepwiselm)
     Pain.Pre.n4_Gp_regrssRoiPairPain.(strcat('SLMquadratic_', num2str(i))) = stepwiselm(tblPsqi_ROIreduced, ...
          'quadratic','ResponseVar','Y','Upper','quadratic');
end
%         EVALUATE THE RSQUARED VALUES OF THE REGRESSION MODELS (T-TEST ROI PAIRS)
for i=1:numel(I)
     Pain.Pre.n4_Gp_regrssRoiPairPain.(strcat('RSq_', num2str(i))) = ...
          [Pain.Pre.n4_Gp_regrssRoiPairPain.(strcat('LMpair_', num2str(i))).Rsquared.Adjusted ...
          Pain.Pre.n4_Gp_regrssRoiPairPain.(strcat('SLMpair_', num2str(i))).Rsquared.Adjusted ...
          Pain.Pre.n4_Gp_regrssRoiPairPain.(strcat('SLMintercept_', num2str(i))).Rsquared.Adjusted ...
          Pain.Pre.n4_Gp_regrssRoiPairPain.(strcat('SLMinteraction_', num2str(i))).Rsquared.Adjusted ...
          Pain.Pre.n4_Gp_regrssRoiPairPain.(strcat('SLMquadratic_', num2str(i))).Rsquared.Adjusted];
end
%         OUTPUT THE OVERALL F AND P-VALUES FOR EACH MODEL
for i=1:numel(I);
     Pain.Pre.n4_Gp_regrssRoiPairPain.lmModel = anova(Pain.Pre.n4_Gp_regrssRoiPairPain.(strcat('LMpair_', num2str(i))), 'summary');
     Pain.Pre.n4_Gp_regrssRoiPairPain.lmModelSummary(i, :) = Pain.Pre.n4_Gp_regrssRoiPairPain.lmModel(2,4:5);
end
%         IDENTIFY SIGNIGCANT MODELS BASED ON THE PVALUES STORED IN lmModelSummary
%         THE PVALUE CRITERIA IS SPECIFIED BELOW
%         SIG = '1'; NON-SIG = '0'
%         THE MODEL NUMBER IS ADDED
%         THE NUMBERS AND NAMES OF EACH ROI PAIR ARE ADDED
%         THE NEW COLUMN ORDER IS 'model','sig','pvalue','F','node1','node2','roi1','roi2'
Pain.Pre.n4_Gp_regrssRoiPairPain.lmModelSummary.sig = [Pain.Pre.n4_Gp_regrssRoiPairPain.lmModelSummary.pValue < 0.05];
Pain.Pre.n4_Gp_regrssRoiPairPain.lmModelSummary.node1 = I;
Pain.Pre.n4_Gp_regrssRoiPairPain.lmModelSummary.node2 = J;
%
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
%
Pain.Pre.n4_Gp_regrssRoiPairPain.lmModelSummary.roi1 = DUMMY_ROI1';
Pain.Pre.n4_Gp_regrssRoiPairPain.lmModelSummary.roi2 = DUMMY_ROI2';
Pain.Pre.n4_Gp_regrssRoiPairPain.lmModelSummary.model = (1:numel(I))';
Pain.Pre.n4_Gp_regrssRoiPairPain.lmModelSummary = Pain.Pre.n4_Gp_regrssRoiPairPain.lmModelSummary(:, [8 3 2 1 4 5 6 7]);

%{
LM_PainPre.modelSummary.roi1 = DUMMY_ROI1';
LM_PainPre.modelSummary.roi2 = DUMMY_ROI2';
LM_PainPre.modelSummary.model = (1:numel(I))';
LM_PainPre.modelSummary = LM_PainPre.modelSummary(:, [8 3 2 1 4 5 6 7]);
%
%
%{
%
%         STEPWISE LINEAR REGRESSION
% tablePsqiPlusROIstep = tablePsqiPlusROI(:,3:7); % SUBTABLE FOR STEPWISE REGRESSION
% Pain.Pre.LM.(strcat('pair_'+1, num2str(i))) = stepwisefit(tablePsqiPlusROIstep,Y,'penter',0.05,'premove',0.10);
%    TO INCLUDE COVARIATES, OR INTERACTION TERMS, CONSIDER MODIFYING THE CODE BELOW
%         zmdl2 = fitlm(tablePsqiPlusROI, 'Y~TiB_hrs+TiB_hrs*GP')
%    FOR DOING SCATTERPLOTS OF THE DATA CONSIDER USING THE CODE BELOW
%         gscatter(tablePsqiPlusROI.Y,tablePsqiPlusROI.SleepEfficiency,tablePsqiPlusROI.GP,'bgr','x.o')


%
%
%}




%%%%%%%%%%%%%%%%%%%%%%   END SCRIPT     %%%%%%%%%%%%%%%%%%%%%%
% [H_PainPre,P_PainPre,CI_PainPre,STATS_PainPre] = ttest2(prePainHC,prePainGps);           % 2 SAMPLE T-TEST
% [H_PainPre,P_PainPre,CI_PainPre,STATS_PainPre] = ttest2(prePainHC,prePainGps);           % 2 SAMPLE T-TEST

% ttest_tval_PainPre(node1, node2) = STATS_PainPre.tstat; % tstat
% ttest_pval_PainPre(node1, node2) = P_PainPre; % pvalue
% [c,~,~,gnames] = multcompare(stats);    % EVALUATE MULTIPLE COMPARISONS

%    THESE NEXT 4 LINES MIGHT NOT BE RELATED TO THE T-TEST
% multcomparePreOutput = [gnames(c(:,1)), gnames(c(:,2)), num2cell(c(:,3:6))];
% multcomparePreOutput = multcomparePreOutput(:,[1 2 6 4 3 5]);
% tableAnovaPre = array2table(multcomparePreOutput, 'VariableNames',{'gp1','gp2', 'pval','gpDiff','lCI','uCI'});
% multiTest.OUT_TEXT2 = [multiTest.OUT_TEXT2 roi1str 'with ' roi2str 'new line' sprintf('\n')];
%


%{
% pid = FDR(pain_anovaresults_pvals,.05);
% [x,y] = find(pain_anovaresults_pvals <= pid);
% [x y];
%
%}


%{
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
%}
%}
