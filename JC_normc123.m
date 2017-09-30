% check SPM version

%%
%       CREATED BY:     JOCHEN WEBER
%       CREATED ON:     11/22/16
%       MODIFIED BY:	JASON CRAGGS
%		MODIFIED ON:	2017_09_27
%       
%       USAGE:          PREPROCESS JASON'S R01 DATA
%       MODIFIED TO:    JASON IS PREPROCESSING ADDITIONAL SUBJECTS
%                       WHICH REQUIRES THAT THEY BE IN A DIFFERENT FOLDER
%%

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
rootpath = '/Volumes/Data/Imaging/R01/preprocessed/_Jason/';
subpattern = 'Sub*_v*';

% get subjects
dirinfo = dir([rootpath subpattern]);
subjlist = {dirinfo.name};

% iterate over subjects
for sc = 1:numel(subjlist)
    
    % find warpfield
    cd([rootpath filesep subjlist{sc}]);
    deffield = dir('y_*.nii');
    wcfiles = dir('wc*.nii');
    if isempty(deffield) || ~isempty(wcfiles)
        continue;
    end
    deffield = [pwd filesep deffield.name];
    
    % find c1/c2/c3 files
    c123 = [dir('c1T*.nii'); dir('c2T*.nii'); dir('c3T*.nii')];
    if numel(c123) ~= 3
        continue;
    end
    c123 = {c123.name};
    for c = 1:3
        c123{c} = [pwd filesep c123{c} ',1'];
    end
    
    % resample job
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {deffield};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = c123(:);
    matlabbatch{1}.spm.spatial.normalise.write.woptions = ...
        struct('bb', [-78, -112, -70; 78, 76, 85], 'vox', [2, 2, 2], 'interp', 1, 'prefix', 'w');
    spm_jobman('run', matlabbatch);
    clear matlabbatch;
end
