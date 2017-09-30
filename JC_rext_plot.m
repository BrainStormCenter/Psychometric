% go to correct folder
cd /Volumes/Data/Imaging/R01/preprocessed

% libraries
n = neuroelf;
x = xff;

% load GLM
try
    glm = x.Document('ChronicPain_R01_zFC_maps.glm');
catch
    glm = xff('ChronicPain_R01_zFC_maps.glm');
end

% load VOI
try
    voi = x.Document('DMN-10_BA13-2_PAIN-12_clusters.voi');
catch
    voi = xff('DMN-10_BA13-2_PAIN-12_clusters.voi');
end

% load LUT
try
    lut = x.Document('fullspec_extended.olt');
catch
    lut = xff('~/Documents/MATLAB/NeuroElf_v11_6501/_core/lut/fullspec_extended.olt');
end
llut = lut.Colors;

% store indices to HC and CLBP groups (excluding 42)
hc = [1,2,3,4,6,8:18,21,22,23,24,26,27,28,30,31,33,35,36];
clbp = [5,7,19,20,25,29,32,34,37,38,39,43,46,47,49,51:56,58:61,63];

% extract VOI betas
vb = glm.VOIBetas(voi);

% DMN pre conditions (2, 4, 10, and 12) and BA13 pre (1, 3, 9, and 11)
hcbpre = [squeeze(mean(vb(hc, [2,4,10,12], :), 2)), squeeze(mean(vb(hc, [1,3,9,11], :), 2))];
cpbpre = [squeeze(mean(vb(clbp, [2,4,10,12], :), 2)), squeeze(mean(vb(clbp, [1,3,9,11], :), 2))];
hccpre = corrcoef(hcbpre);
cpcpre = corrcoef(cpbpre);

% same for pos post and neg post
hcbpsp = [squeeze(mean(vb(hc, [6,8], :), 2)), squeeze(mean(vb(hc, [5,7], :), 2))];
hcbpsn = [squeeze(mean(vb(hc, [14,16], :), 2)), squeeze(mean(vb(hc, [13,15], :), 2))];
cpbpsp = [squeeze(mean(vb(clbp, [6,8], :), 2)), squeeze(mean(vb(clbp, [5,7], :), 2))];
cpbpsn = [squeeze(mean(vb(clbp, [14,16], :), 2)), squeeze(mean(vb(clbp, [13,15], :), 2))];
hccpsp = corrcoef(hcbpsp);
hccpsn = corrcoef(hcbpsn);
cpcpsp = corrcoef(cpbpsp);
cpcpsn = corrcoef(cpbpsn);

% combine into image
hcallc = [hccpre, zeros(48,8), hccpsp, zeros(48,8), hccpsn];
cpallc = [cpcpre, zeros(48,8), cpcpsp, zeros(48,8), cpcpsn];
allc = [hcallc; zeros(8, 160); cpallc];
allmap = n.threshlutc(allc, llut);

% show image
figure; image(allmap)
ylabel('HC (top) vs. CLBP (bottom)');
xlabel('pre (left) vs. post-pos (middle) vs. post-neg (right)');
title('cross-correlation of ROIs across subjects');
set(gca, 'FontSize', 18)

% extract DMN portion
[ltrx, ltry] = ndgrid(1:10, 1:10);
ltry(ltrx <= ltry) = 0;
ltrx(ltry == 0) = 0;
ltrx = ltrx(:) + 10 .* max(0, (ltry(:) - 1));
ltrx(ltrx == 0) = [];
hcdmnpre = hccpre(1:10, 1:10);
cpdmnpre = cpcpre(1:10, 1:10);
hcdmnppost = hccpsp(1:10, 1:10);
cpdmnppost = cpcpsp(1:10, 1:10);
hcdmnnpost = hccpsn(1:10, 1:10);
cpdmnnpost = cpcpsn(1:10, 1:10);

% extract area-to-area correlations (across subjects) = network flex
hcdmnpreflex = hcdmnpre(ltrx);
cpdmnpreflex = cpdmnpre(ltrx);
hcdmnppostflex = hcdmnppost(ltrx);
cpdmnppostflex = cpdmnppost(ltrx);
hcdmnnpostflex = hcdmnnpost(ltrx);
cpdmnnpostflex = cpdmnnpost(ltrx);

% compute means and approx. SEs
mhcdmnpreflex = mean(hcdmnpreflex);
shcdmnpreflex = std(hcdmnpreflex) / 3;
mcpdmnpreflex = mean(cpdmnpreflex);
scpdmnpreflex = std(cpdmnpreflex) / 3;
mhcdmnppostflex = mean(hcdmnppostflex);
shcdmnppostflex = std(hcdmnppostflex) / 3;
mcpdmnppostflex = mean(cpdmnppostflex);
scpdmnppostflex = std(cpdmnppostflex) / 3;
mhcdmnnpostflex = mean(hcdmnnpostflex);
shcdmnnpostflex = std(hcdmnnpostflex) / 3;
mcpdmnnpostflex = mean(cpdmnnpostflex);
scpdmnnpostflex = std(cpdmnnpostflex) / 3;

% error bars
figure; errorbar( ...
    [nan, mhcdmnpreflex, mcpdmnpreflex, nan, mhcdmnppostflex, mcpdmnppostflex, nan, mhcdmnnpostflex, mcpdmnnpostflex, nan], ...
    [0, shcdmnpreflex, scpdmnpreflex, 0, shcdmnppostflex, scpdmnppostflex, 0, shcdmnnpostflex, scpdmnnpostflex, 0]);
set(gca, 'XLim', [1, 10]);

% scatter plot
r1 = 1; r2 = 10; figure; subplot(1,2,1); scatter(hcbpre(:, r1), hcbpre(:, r2)); subplot(1,2,2); scatter(cpbpre(:, r1), cpbpre(:, r2));
