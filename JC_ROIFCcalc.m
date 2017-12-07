% check SPM version
if ~strcmpi(spm('ver'), 'spm12')
    error('spm:version:wrongSPMVersion', 'This script requires SPM12.');
end

% initialize SPM defaults and job manager
% spm12path = '/cluster/folder/craggs/software/spm12';
spm12path = '/Users/jcraggs/Applications/spm12';
spm('defaults','FMRI');
spm_jobman('initcfg');
clear matlabbatch;

% configure root path and subject pattern, as well as file patterns
%rootpath = '/cluster/folder/craggs/study/preprocessed/';
%rootpath = '/Volumes/Data/Imaging/R01/preprocessed/';
rootpath = '/Volumes/Data/Imaging/R01/preprocessed/_Jason_0/';
cd(rootpath);

% load VOI
voi = xff('*.voi');
nvoi = numel(voi.VOI);
voigroups = {'DMN', 'BA13', 'PAIN'};
voigi = voigroups;
for vc = 1:numel(voigi)
    voigi{vc} = find(~cellfun('isempty', regexpi(voi.VOINames, ['^' voigroups{vc}])));
end
voigii = cell(numel(voigi) + round(0.5 * numel(voigi) * (numel(voigi) - 1)), 1);
voigin = voigii;
gi = 1;
for vc = 1:numel(voigi)
    for vc2 = vc:numel(voigi)
        voigin{gi} = sprintf('%s <-> %s', voigroups{vc}, voigroups{vc2});
        [gix, giy] = ndgrid(voigi{vc}(:)', voigi{vc2}(:)');
        gii = [gix(:), giy(:)];
        gii(gii(:, 1) == gii(:, 2), :) = [];
        voigii{gi} = sub2ind([nvoi, nvoi], gii(:, 1), gii(:, 2));
        gi = gi + 1;
    end
end

% subject pattern
subpattern = 'Sub*_v*';

% requires neuroelf
n = neuroelf;

% locate FC prepared files
fcswa = n.findfiles([pwd '/' subpattern], 'fcsw*.nii', '-d1');
if isempty(fcswa)
    error('JCfunc:filematch:notFilesFound', 'No FC-prepared files found.');
end

% load first volume
vol = spm_vol([fcswa{1} ',1']);
voldim = vol.dim;
volvox = prod(voldim(1:3));
volmat = vol.mat;
volimat = inv(volmat);

% get indices
voiidx = voi.VOI;
voiidx = {voiidx.Voxels};
voi.ClearObject;
for vc = 1:numel(voiidx)

    % transform to voxel coordinates
    vox = round(volimat * ([voiidx{vc}, ones(size(voiidx{vc}, 1), 1)])')';
    vox(:, 4) = [];
    vox(any(vox < 1, 2) | vox(:, 1) > voldim(1) | vox(:, 2) > voldim(2) | vox(:, 3) > voldim(3), :) = [];
    if isempty(vox)
        voiidx{vc} = zeros(0, 1);
        continue;
    end
    voiidx{vc} = unique(sub2ind(voldim, vox(:, 1), vox(:, 2), vox(:, 3)));
end

% generate cell array for time courses, CC matrices, and averages
fctcs = cell(numel(fcswa), numel(voiidx));
fcccs = cell(numel(fcswa), 1);
fccas = zeros(numel(fcswa), numel(voigii));

% loop over data
for fc = 1:numel(fcswa)

    % load data
    fprintf('Processing %s...\n', fcswa{fc});
    fcvol = spm_vol(fcswa{fc});
    nvol = numel(fcvol);
    fcdata = spm_read_vols(fcvol);
    fcdata = reshape(fcdata, volvox, nvol);

    % extract for each voi
    for vc = 1:numel(voiidx)

        % empty?
        if isempty(voiidx{vc})
            fctcs{fc, vc} = zeros(nvol, 1);
            continue;
        end

        % time courses
        voxtcs = fcdata(voiidx{vc}, :);

        % remove bad time courses
        badtcs = (any(isnan(voxtcs), 2) | all(voxtcs == 0, 2));
        voxtcs(badtcs, :) = [];
        if isempty(voxtcs)
            fctcs{fc, vc} = zeros(nvol, 1);
        else
            fctcs{fc, vc} = mean(voxtcs, 1)';
        end
    end

    % compute correlation matrix
    fcccs{fc} = corrcoef(cat(2, fctcs{fc, :}));

    % extract averages
    for vc = 1:numel(voigii)
        fccas(fc, vc) = mean(atanh(fcccs{fc}(voigii{vc})));
    end
end
