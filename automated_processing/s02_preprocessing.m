%
%       CREATED BY:	JOCHEN WEBER
%       CREATED ON:	11/22/16
%
%       MODIFIED BY:	JASON CRAGGS
%       MODIFIED ON:	2019-05-09
%
%	MODIFIED BY:	JACOB GOTBERG
%	MODIFIED ON:	2019-03-02
%
%       PURPOSE:	PREPROCESS FMRI SCANS FOR ANALYSIS
%
%	USAGE:		matlab -r "prog 2_preprocessing.m dirname" 
%

%'/storage/hpc/group/sleeplab/preprocessed/SPIN2_013_V1_unknown'

function []=s02_preprocessing(dirname)
    disp("*** Starting 2_preprocessing ***")
    fprintf("Pre Processing file %s\n", dirname)
   
    % initialize SPM defaults and job manager
    disp("Loading SPM12")
    spm12path = '/storage/hpc/group/sleeplab/software/spm12';
    addpath(spm12path)
    
    spm('defaults','FMRI');
    spm_jobman('initcfg');
    clear matlabbatch;
    disp("SPM12 Loaded")

    % check SPM version
    if ~strcmpi(spm('ver'), 'spm12')
        error('spm:version:wrongSPMVersion', 'This script requires SPM12.');
    end

    % configure root path and subject pattern, as well as file patterns
    anatpattern = 'SPIN2_*_t1_*.nii';
    funcpattern = 'SPIN2_*_fMRI_*.nii';
    
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

    % set primary path
    cd(dirname);
    % hard code to 5 per protocol at MRI machine
    %runz = dir('*fMRI*');
    numruns = 5; %size(runz, 1)

    % gzip -d everything
    gzipfiles = dir('*.gz');
    if ~isempty(gzipfiles)
        system('gzip -d *.gz');
    end

    % locate files
    anatfile = dir(anatpattern);
    funcfiles = dir(funcpattern);
    
     if numel(anatfile) ~= 1 || numel(funcfiles) ~= numruns
         warning('spm:prepro:invalidNumberOfFiles', ...
             'Number of files incorrect for %s.', dirname);
         return
     end

     if exist(['swa' funcfiles(end).name], 'file') == 2
         fprintf('%s already preprocessed.\n', dirname);
         return
     end

     % check number of slices
     vol = spm_vol(funcfiles(1).name);
     if vol(1).dim(3) ~= nslices
         warning('spm:prepro:invalidNumberOfSlices', ...
             'Number of slices mismatch for %s.', dirname);
         return
     end

    % change filetype of functional data to float32 (reduce precision loss)
    for fc = 1:numruns
        vol = spm_vol(funcfiles(fc).name);
        if vol(1).dt(1) == 16
            continue;
        end
        fprintf('Converting %s in %s to float32...\n', funcfiles(fc).name, dirname);
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
    %matlabbatch{1}.spm.spatial.smooth.data = {'/storage/hpc/group/sleeplab/preprocessed/SPIN2_013_V1_unknown/SPIN2_013_t1_space_sag_p2_iso_s02.nii,1'};
    
    anat_file = fullfile(dirname, anatfile.name);
    file_to_smooth = strcat(anat_file, ",1");
    anat = cellstr(file_to_smooth);

    matlabbatch{1}.spm.spatial.smooth.data = anat;
    matlabbatch{1}.spm.spatial.smooth.fwhm = [12 12 12];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    fprintf('Smoothing %s in %s...\n', anatfile.name, dirname);
    spm_jobman('run', matlabbatch);
    %spm_jobman('interactive',matlabbatch);
    clear matlabbatch;

    % coregister smoothed struct to old T1 template
    sanat_file = fullfile(dirname, strcat("s",anatfile.name));
    sanat = cellstr(sanat_file);

    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[spm12path '/toolbox/OldNorm/T1.nii,1']};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = sanat;
    matlabbatch{1}.spm.spatial.coreg.estimate.other = anat;
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4, 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 ...
        0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7, 7];
    fprintf('Coregistering (roughly) s%s in %s to T1 template...\n', anatfile.name, dirname);
    spm_jobman('run', matlabbatch);
    clear matlabbatch;
    %delete(sanat);

    % segment the structural
    def_file = fullfile(dirname, strcat("y_", anatfile.name));
    DefField = cellstr(def_file);
    
    manat_file = fullfile(dirname, strcat("m", anatfile.name));
    manat = cellstr(fullfile(dirname, strcat("m", anatfile.name)));
    
    tpms = repmat({[spm12path filesep 'tpm' filesep 'TPM.nii,']}, 1, 6);
    for tc = 1:6
        tpms{tc} = {[tpms{tc} char(48+tc)]};
    end
    matlabbatch{1}.spm.spatial.preproc.channel.vols = anat;
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
    if exist(def_file, 'file') ~= 2 || exist(manat_file, 'file') ~= 2
        spm_jobman('run', matlabbatch);
    end
    clear matlabbatch;

    % normalize the segmented structural to T1 template
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = DefField;
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = manat;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = ...
        [-78, -112, -70; 78, 76, 85];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    fprintf('Warping anatomical %s in %s to MNI space...\n', anatfile.name, dirname);
    spm_jobman('run', matlabbatch);
    clear matlabbatch;

     % skull strip the segmented structural
    
    c1 = fullfile(dirname, strcat("c1", anatfile.name));
    c2 = fullfile(dirname, strcat("c2", anatfile.name));
    c3 = fullfile(dirname, strcat("c3", anatfile.name));
    sx = ['xm' anatfile.name];
    matlabbatch{1}.spm.util.imcalc.input = cellstr([manat; c1; c2; c3]);
    matlabbatch{1}.spm.util.imcalc.output = sx;
    matlabbatch{1}.spm.util.imcalc.outdir = {''};
    matlabbatch{1}.spm.util.imcalc.expression = 'i1.*((i2+i3+i4)>=.5)';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    fprintf('Skull stripping m%s in %s...\n', anatfile.name, dirname);
    spm_jobman('run', matlabbatch);
    clear matlabbatch;
  
    % generate func file names
    fprintf('Starting func file names');
    funcfiles = {strcat('/', funcfiles.name)};
    meanfuncfile = sprintf('%smeana%s,1', dirname, funcfiles{1});
    smeanfuncfile = sprintf('%ssmeana%s,1', dirname, funcfiles{1});
    afuncfiles = funcfiles;
    wafuncfiles = funcfiles;
    for fc = 1:numruns
        disp(fc);
        
        funcfiles{fc} = repmat(funcfiles(fc), nvols, 1);
        afuncfiles{fc} = funcfiles{fc};
        wafuncfiles{fc} = funcfiles{fc};
        for vc = 1:nvols
            funcfiles{fc}{vc} = sprintf('%s%s,%d', dirname, funcfiles{fc}{vc}, vc);
            afuncfiles{fc}{vc} = sprintf('%sa%s,%d', dirname, afuncfiles{fc}{vc}, vc);
            wafuncfiles{fc}{vc} = sprintf('%swa%s,%d', dirname, wafuncfiles{fc}{vc}, vc);
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
    fprintf('Slice-time correcting data in %s...\n', dirname);
    spm_jobman('run', matlabbatch);
    clear matlabbatch;

    % motion correction/realignment (+ mean image)
    matlabbatch{1}.spm.spatial.realign.estwrite.data = afuncfiles;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions = struct( ...
        'quality', 0.9, 'sep', 4, 'fwhm', 5, 'rtm', 1, 'interp', 2, 'wrap', [0, 0, 0], 'weight', '');
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions = struct( ...
        'which', [0, 1], 'interp', 4, 'wrap', [0, 0, 0,], 'mask', 1, 'prefix', 'r');
    fprintf('Realignment of data in %s...\n', dirname);
    spm_jobman('run', matlabbatch);
    clear matlabbatch;

    % smooth mean functional (for first step of coreg)
    matlabbatch{1}.spm.spatial.smooth.data = {meanfuncfile};
    matlabbatch{1}.spm.spatial.smooth.fwhm = [12, 12, 12];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    fprintf('Smoothing %s in %s (for coreg)...\n', meanfuncfile, dirname);
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
    fprintf('Coregistering smoothed mean-func in %s to EPI template...\n', dirname);
    spm_jobman('run', matlabbatch);
    clear matlabbatch;

    % coregister funcs to T1
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[dirname sx]};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {meanfuncfile};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = cat(1, afuncfiles{:});
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4, 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 ...
        0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7, 7];
    fprintf('Coregistering mean-func in %s to skull-stripped anatomical...\n', dirname);
    spm_jobman('run', matlabbatch);
    clear matlabbatch;

    % MNI-normalize EPI data
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {DefField};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cat(1, afuncfiles{:});
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = wbbox;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = repmat(wvox, 1, 3);
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    fprintf('MNI-warping EPI data in %s using %gmm voxel size...\n', dirname, wvox);
    spm_jobman('run', matlabbatch);
    clear matlabbatch

    % smooth EPI data
    matlabbatch{1}.spm.spatial.smooth.data = cat(1, wafuncfiles{:});
    matlabbatch{1}.spm.spatial.smooth.fwhm = repmat(epismk, 1, 3);
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    fprintf('Smoothing warped EPI data in %s with %gmm kernel...\n', dirname, epismk);
    spm_jobman('run', matlabbatch);
    clear matlabbatch;
    
end
