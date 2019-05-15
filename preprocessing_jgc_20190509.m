% check SPM version

%%
%       CREATED BY:     JOCHEN WEBER
%       CREATED ON:     11/22/16
%
%       MODIFIED BY:    JASON CRAGGS
%      LATEST MODIFICATION:     2019_05_09 (modifying for McCrae study)
%
%       USAGE:          PREPROCESS CHRISTINA'S R01 DATA
%       MODIFIED TO:    JASON IS PREPROCESSING ADDITIONAL SUBJECTS
%                       WHICH REQUIRES THAT THEY BE IN A DIFFERENT FOLDER
%
%       2017_12_06
%       TESTING THE ORDER THE SCRIPTS SHOULD BE RUN
%       THIS IS THE FIRST ONE
%%

if ~strcmpi(spm('ver'), 'spm12')
    error('spm:version:wrongSPMVersion', 'This script requires SPM12.');
end

% initialize SPM defaults and job manager
% spm12path = '/cluster/folder/craggs/software/spm12';
%spm12path = '/Users/jcraggs/Applications/spm12';
spm12path = '/storage/hpc/group/sleeplab/software/spm12';
spm('defaults','FMRI');
spm_jobman('initcfg');
clear matlabbatch;

% configure root path and subject pattern, as well as file patterns
%rootpath = '/cluster/folder/craggs/study/preprocessed/';
%rootpath = '/Volumes/Data/Imaging/R01/preprocessed/';
%rootpath = '/Volumes/Data/Imaging/R01/preprocessed/_Jason/';
%rootpath = '/Volumes/Data/Imaging/R01/preprocessed/_Jason_0/';
%rootpath = '/Volumes/Data/Imaging/R01/preprocessed/_Jason_step1/';
rootpath = '/storage/hpc/group/sleeplab/preprocessed/';
subpattern = 'SPIN2_013_tes*';
anatpattern = 'SPIN2_*_t1.nii';
funcpattern = 'SPIN2_*_run*_*.nii';
%anatpattern = 'T1_*.nii';
%funcpattern = 'RSrun*.nii';
numruns = 5;

% set variables, number of volumes and functional slices, TR
nvols = 120;
nslices = 36;
TR = 2.46;

%   SLICE ACQUISITION ORDER IS INTERLEAVED BUT THE
%   FIRST SLICE ACQUIRED IS SLICE #2; THEREFORE THE
%   CORRECT SPM SPECIFICATION IS [2:2:42 1:2:42]
% slice aquisition order (interleaved bottom) and reference slice
if mod(nslices, 2) == 1
    sliceorder = [1:2:nslices, 2:2:nslices];
else
    sliceorder = [2:2:nslices, 1:2:nslices];
end
refslice = sliceorder(1);

% reslice voxel size and bounding box
wvox = 2;
wbbox = [-78, -112, -70; 78, 76, 85];

% smoothing kernel (EPI)
epismk = 6;

% compute TA
TA = TR - (TR / nslices);

% find subjects in root folder
dirinfo = dir([rootpath subpattern]);
subjlist = {dirinfo.name};

% pick subject according to job number
for sc = 1:numel(subjlist)

    % set primary path
    primary_path = [rootpath subjlist{sc} filesep];
    cd(primary_path);

    % gzip -d everything
    gzipfiles = dir('*.gz');
    if ~isempty(gzipfiles)
        system('gzip -d *.gz');
    end

    % locate files
    anatfile = dir([primary_path anatpattern]);
    funcfiles = dir([primary_path funcpattern]);
    if numel(anatfile) ~= 1 || numel(funcfiles) ~= numruns
        warning('spm:prepro:invalidNumberOfFiles', ...
            'Number of files incorrect for %s.', subjlist{sc});
        continue;
    end
    if exist(['swa' funcfiles(end).name], 'file') == 2
        fprintf('%s already preprocessed.\n', subjlist{sc});
        continue;
    end

    % check number of slices
    vol = spm_vol(funcfiles(1).name);
    if vol(1).dim(3) ~= nslices
        warning('spm:prepro:invalidNumberOfSlices', ...
            'Number of slices mismatch for %s.', subjlist{sc});
        continue;
    end

    % change filetype of functional data to float32 (reduce precision loss)
    for fc = 1:numruns
        vol = spm_vol(funcfiles(fc).name);
        if vol(1).dt(1) == 16
            continue;
        end
        fprintf('Converting %s in %s to float32...\n', funcfiles(fc).name, subjlist{sc});
        voldata = spm_read_vols(vol);
        [vol.dt] = deal([16, 0]);
        for vc = 1:numel(vol)
            vol(vc).pinfo(3) = vol(1).pinfo(3) + (vc - 1) * prod([4, vol(1).dim]);
        end
        delete(funcfiles(fc).name);
        for vc = 1:numel(vol)
            spm_write_vol(vol(vc), voldata(:, :, :, vc));
        end
    end

    % smooth the structural
    anat = [primary_path anatfile.name];
    matlabbatch{1}.spm.spatial.smooth.data = {anat};
    matlabbatch{1}.spm.spatial.smooth.fwhm = [12 12 12];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    fprintf('Smoothing %s in %s...\n', anatfile.name, subjlist{sc});
    spm_jobman('run', matlabbatch);
    clear matlabbatch;

    % coregister smoothed struct to old T1 template
    sanat = [primary_path 's' anatfile.name];
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[spm12path '/toolbox/OldNorm/T1.nii,1']};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {sanat};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {anat};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4, 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 ...
        0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7, 7];
    fprintf('Coregistering (roughly) s%s in %s to T1 template...\n', anatfile.name, subjlist{sc});
    spm_jobman('run', matlabbatch);
    clear matlabbatch;
    %delete(sanat);

    % segment the structural
    DefField = [primary_path 'y_' anatfile.name];
    manat = [primary_path 'm' anatfile.name];
    tpms = repmat({[spm12path filesep 'tpm' filesep 'TPM.nii,']}, 1, 6);
    for tc = 1:6
        tpms{tc} = {[tpms{tc} char(48+tc)]};
    end
    matlabbatch{1}.spm.spatial.preproc.channel.vols = {anat};
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel.write = [0, 1];
    matlabbatch{1}.spm.spatial.preproc.tissue = struct('tpm', tpms, ...
        'ngaus', {1, 1, 2, 3, 4, 2}, ...
        'native', {[1, 1], [1, 1], [1, 1], [0, 0], [0, 0], [0, 0]}, ...
        'warped', {[0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]});
    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0, 0.001, 0.5, 0.05, 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [1, 1];
    if exist(DefField, 'file') ~= 2 || exist(manat, 'file') ~= 2
        spm_jobman('run', matlabbatch);
    end
    clear matlabbatch;

    % normalize the segmented structural to T1 template
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {DefField};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {manat};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = ...
        [-78, -112, -70; 78, 76, 85];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    fprintf('Warping anatomical %s in %s to MNI space...\n', anatfile.name, subjlist{sc});
    spm_jobman('run', matlabbatch);
    clear matlabbatch;

    % skull strip the segmented structural
    c1 = [primary_path 'c1' anatfile.name];
    c2 = [primary_path 'c2' anatfile.name];
    c3 = [primary_path 'c3' anatfile.name];
    sx = ['xm' anatfile.name];
    matlabbatch{1}.spm.util.imcalc.input = {manat; c1; c2; c3};
    matlabbatch{1}.spm.util.imcalc.output = sx;
    matlabbatch{1}.spm.util.imcalc.outdir = {''};
    matlabbatch{1}.spm.util.imcalc.expression = 'i1.*((i2+i3+i4)>=.5)';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    fprintf('Skull stripping m%s in %s...\n', anatfile.name, subjlist{sc});
    spm_jobman('run', matlabbatch);
    clear matlabbatch;

    % generate func file names
    funcfiles = {funcfiles.name};
    meanfuncfile = sprintf('%smeana%s,1', primary_path, funcfiles{1});
    smeanfuncfile = sprintf('%ssmeana%s,1', primary_path, funcfiles{1});
    afuncfiles = funcfiles;
    wafuncfiles = funcfiles;
    for fc = 1:numruns
        funcfiles{fc} = repmat(funcfiles(fc), nvols, 1);
        afuncfiles{fc} = funcfiles{fc};
        wafuncfiles{fc} = funcfiles{fc};
        for vc = 1:nvols
            funcfiles{fc}{vc} = sprintf('%s%s,%d', primary_path, funcfiles{fc}{vc}, vc);
            afuncfiles{fc}{vc} = sprintf('%sa%s,%d', primary_path, afuncfiles{fc}{vc}, vc);
            wafuncfiles{fc}{vc} = sprintf('%swa%s,%d', primary_path, wafuncfiles{fc}{vc}, vc);
        end
    end

    % slice-timing of functional data
    matlabbatch{1}.spm.temporal.st.scans = funcfiles;
    matlabbatch{1}.spm.temporal.st.nslices = nslices;
    matlabbatch{1}.spm.temporal.st.tr = TR;
    matlabbatch{1}.spm.temporal.st.ta = TA;
    matlabbatch{1}.spm.temporal.st.so = sliceorder;
    matlabbatch{1}.spm.temporal.st.refslice = refslice;
    matlabbatch{1}.spm.temporal.st.prefix = 'a';
    fprintf('Slice-time correcting data in %s...\n', subjlist{sc});
    spm_jobman('run', matlabbatch);
    clear matlabbatch;

    % motion correction/realignment (+ mean image)
    matlabbatch{1}.spm.spatial.realign.estwrite.data = afuncfiles;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions = struct( ...
        'quality', 0.9, 'sep', 4, 'fwhm', 5, 'rtm', 1, 'interp', 2, 'wrap', [0, 0, 0], 'weight', '');
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions = struct( ...
        'which', [0, 1], 'interp', 4, 'wrap', [0, 0, 0,], 'mask', 1, 'prefix', 'r');
    fprintf('Realignment of data in %s...\n', subjlist{sc});
    spm_jobman('run', matlabbatch);
    clear matlabbatch;

    % smooth mean functional (for first step of coreg)
    matlabbatch{1}.spm.spatial.smooth.data = {meanfuncfile};
    matlabbatch{1}.spm.spatial.smooth.fwhm = [12, 12, 12];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    fprintf('Smoothing %s in %s (for coreg)...\n', meanfuncfile, subjlist{sc});
    spm_jobman('run', matlabbatch);
    clear matlabbatch;

    % coregister funcs to EPI template
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = ...
        {[spm12path filesep 'toolbox' filesep 'OldNorm' filesep 'EPI.nii,1']};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {smeanfuncfile};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = [{meanfuncfile}; cat(1, afuncfiles{:})];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4, 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 ...
        0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7, 7];
    fprintf('Coregistering smoothed mean-func in %s to EPI template...\n', subjlist{sc});
    spm_jobman('run', matlabbatch);
    clear matlabbatch;

    % coregister funcs to T1
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[primary_path sx]};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {meanfuncfile};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = cat(1, afuncfiles{:});
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4, 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 ...
        0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7, 7];
    fprintf('Coregistering mean-func in %s to skull-stripped anatomical...\n', subjlist{sc});
    spm_jobman('run', matlabbatch);
    clear matlabbatch;

    % MNI-normalize EPI data
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {DefField};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cat(1, afuncfiles{:});
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = wbbox;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = repmat(wvox, 1, 3);
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    fprintf('MNI-warping EPI data in %s using %gmm voxel size...\n', subjlist{sc}, wvox);
    spm_jobman('run', matlabbatch);
    clear matlabbatch

    % smooth EPI data
    matlabbatch{1}.spm.spatial.smooth.data = cat(1, wafuncfiles{:});
    matlabbatch{1}.spm.spatial.smooth.fwhm = repmat(epismk, 1, 3);
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    fprintf('Smoothing warped EPI data in %s with %gmm kernel...\n', subjlist{sc}, epismk);
    spm_jobman('run', matlabbatch);
    clear matlabbatch;
end
