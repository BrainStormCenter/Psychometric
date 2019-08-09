% List of open inputs
% Slice Timing: Session - cfg_files
nrun = X; % enter the number of runs here
jobfile = {'/Users/jcraggs/Library/Mobile Documents/com~apple~CloudDocs/Documents/git/McCrae/Psychometric/automated_processing/McC_preprocessing_jgc_20190709_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(1, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Slice Timing: Session - cfg_files
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
