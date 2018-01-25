%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%		CREATED BY:		JASON CRAGGS
%		CREATED ON:		2018-01-25
%
%		USAGE:			TO CREATE GRAPHS OF RESULTS
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% libraries
n = neuroelf;
x = xff;

% load LUT
try
    lut = x.Document('fullspec_extended.olt');
catch
     lut = xff('~/Documents/MATLAB/neuroelf-matlab/_core/lut/fullspec_extended.olt');
end

llut = lut.Colors;

%   GRAPHS FOR THE PAIN NETWORK
allmap1 = n.threshlutc(mean(mean(painzgfcccs( : , : , i1 , : , 1), 4), 3) + eps, llut);
allmap2 = n.threshlutc(mean(mean(painzgfcccs( : , : , i2 , : , 1), 4), 3) + eps, llut);
[h, p, ci, stats] = ttest2(mean(painzgfcccs(: , : , i1, 1, :), 5), mean(painzgfcccs(: , : , i2, 1, :), 5), 'dim',3, 'tail', 'both', 'vartype', 'unequal');



%     GRAPHS FOR THE DMN NETWORK
allmap3 = n.threshlutc(mean(mean(dmnzgfcccs( : , : , i1 , : , 1), 4), 3) + eps, llut);
allmap4 = n.threshlutc(mean(mean(dmnzgfcccs( : , : , i2 , : , 1), 4), 3) + eps, llut);
[h, p, ci, stats2] = ttest2(mean(dmnzgfcccs(: , : , i1, 1, :), 5), mean(dmnzgfcccs(: , : , i2, 1, :), 5), 'dim',3, 'tail', 'both', 'vartype', 'unequal');


figure; subplot(2, 2, 1); image(allmap1); subplot(2, 2, 2); image(allmap2); subplot(2, 2, 3); imagesc(stats.tstat);

figure; subplot(2, 2, 1); image(allmap3); subplot(2, 2, 2); image(allmap4); subplot(2, 2, 3); imagesc(stats2.tstat);
