function JC_FCprepro(ppfile, rpfile, trfilt, gsfiles)
%   JC_FCPREPRO  Preprocess functional data file(s) for functional connectivity
%   JC_FCPREPRO(PPFILE, RPFILE) removes the variance associated with
%   columns in RPFILE (text file or matfile) from PPFILE.
%
%   JC_FCPREPRO(PPFILE, RPFILE, TRFILT) also removes low-frequency drifts
%   present in the data up to TRFILT wavelength (using DCT filtering).
%
%   JC_FCPREPRO(PPFILE, RPFILE, TRFILT, GSFILES) also removes the mean
%   covariate from each file in GSFILES (resampled to the space in PPFILE)
%
%   This function uses SPM functions, and thus is dependent on SPM8 or
%   higher. Functions used include spm_vol, spm_read_vols, spm_write_vol,
%   spm_filter, and others.
%
%   The function will make a copy of (each of) the input file(s) with a
%   prefix of "fc" (e.g. "fcswraRUN.nii").



%%
%       December 7, 2017
%       TESTING THE SCRIPT
% THESE COMMANDS WORKED FOR A SINGLE SUBJECT ?
%rpfiles = n.findfiles([pwd '/Sub*'], 'rp*.txt', '-d1');
% wc1 = n.findfiles([pwd '/Sub*'], 'wc1*.nii', '-d1');
% wc2 = n.findfiles([pwd '/Sub*'], 'wc2*.nii', '-d1');
% wc3 = n.findfiles([pwd '/Sub*'], 'wc3*.nii', '-d1');
%
%JC_FCprepro(swa,rpfiles,120/2.8, [wc1; wc2; wc3])
%
%
%
%%

% requires CHAR or CELL input
if ischar(ppfile)
    ppfile = cellstr(ppfile);
elseif ~iscell(ppfile)
    error('JCfunc:prepro:invalidArgument', 'Invalid PPFILE argument.');
end
ppfile = ppfile(:);
if ~all(cellfun(@ischar, ppfile)) || any(cellfun('isempty', ppfile)) || ...
    any(cellfun('isempty', regexpi(ppfile, '\.nii$'))) || ...
    sum(cellfun(@exist, ppfile)) ~= (2 * numel(ppfile))
    error('JCfunc:prepro:invalidArgument', 'Invalid PPFILE argument.');
end
if ischar(rpfile)
    rpfile = cellstr(rpfile);
elseif ~iscell(rpfile)
    error('JCfunc:prepro:invalidArgument', 'Invalid RPFILE argument.');
end
if numel(rpfile) ~= numel(ppfile)
    error('JCfunc:prepro:invalidArgument', 'Invalid RPFILE argument.');
end
rpfile = rpfile(:);
if ~all(cellfun(@ischar, rpfile)) || any(cellfun('isempty', rpfile)) || ...
    any(cellfun('isempty', regexpi(rpfile, '(rp.*\.txt|\.mat)$'))) || ...
    sum(cellfun(@exist, rpfile)) ~= (2 * numel(rpfile))
    error('JCfunc:prepro:invalidArgument', 'Invalid RPFILE argument.');
end
if nargin < 3 || ~isa(trfilt, 'double') || numel(trfilt) ~= 1 || ...
    isnan(trfilt) || trfilt < 2
    trfilt = Inf;
end
if nargin < 4
    gsfiles = {};
elseif ischar(gsfiles)
    gsfiles = cellstr(gsfiles);
elseif ~iscell(gsfiles)
    error('JCfunc:prepro:invalidArgument', 'Invalid GSFILES argument.');
end
gsfiles = gsfiles(:);
if ~isempty(gsfiles) && ...
   (~all(cellfun(@ischar, gsfiles)) || any(cellfun('isempty', gsfiles)) || ...
    any(cellfun('isempty', regexpi(gsfiles, '\.nii$'))) || ...
    sum(cellfun(@exist, gsfiles)) ~= (2 * numel(gsfiles)))
    error('JCfunc:prepro:invalidArgument', 'Invalid GSFILES argument.');
end

% iterate over files
for fc = 1:numel(ppfile)

    % try
    try

        % load RPFILE
        rpdata = load(rpfile{fc});
        while isstruct(rpdata)
            rpfields = fieldnames(rpdata);
            rpdata = rpdata.(rpfields{1});
        end
        if ~isa(rpdata, 'double') || size(rpdata, 2) ~= 6
            warning('JCfunc:prepro:invalidFileContent', ...
                'Invalid RPFILE content in %s.', rpfile{fc});
            continue;
        end

        % MAT file exists as well
        niifile = ppfile{fc};
        matfile = [niifile(1:end-2) 'mat'];
        matexists = (exist(matfile, 'file') == 2);

        % new filename
        [fcpath, fcname] = fileparts(niifile);
        if isempty(fcpath)
            fcpath = pwd;
        end
        fcmatfile = [fcpath filesep 'fc' fcname '.mat'];
        fcfile = [fcpath filesep 'fc' fcname '.nii'];

        % load file
        fcvols = spm_vol(niifile);
        nvols = numel(fcvols);
        if nvols ~= size(rpdata, 1)
            warning('JCfunc:prepro:invalidFileContent', ...
                'Invalid RPFILE content in %s.', rpfile{fc});
            continue;
        end

        % copy file, then re-load with those headers
        copyfile(niifile, fcfile);
        if matexists
            copyfile(matfile, fcmatfile);
        end
        fcvols = spm_vol(fcfile);

        % load data
        fcdata = spm_read_vols(fcvols);

        % transpose
        szdata = size(fcdata);
        fcdata = reshape(fcdata, prod(szdata(1:3)), szdata(4));
        udata = ~any(isinf(fcdata) | isnan(fcdata), 2) & any(fcdata ~= 0, 2);

        % filters
        if trfilt > (2 * nvols)
            dcts = zeros(nvols, 0);
        else
            dcts = spm_filter(struct('RT', 1, 'row', (1:nvols)', 'HParam', trfilt));
            dcts = dcts.X0;
        end

        % reslice gsimage in space of first volume
        gmatrix = zeros(nvols, numel(gsfiles));
        if ~isempty(gsfiles)
            spm_reslice([{[niifile ',1']}; gsfiles(:)], struct('mean', false, 'which', 1));
        end

        % global signal files
        for gc = numel(gsfiles):-1:1
            [gspath, gsfile, gsext] = fileparts(gsfiles{gc});
            if isempty(gspath)
                gspath = pwd;
            end
            gsvol = spm_vol([gspath filesep 'r' gsfile gsext]);
            gsidx = spm_read_vols(gsvol);
            gsidx = gsidx(:, :, :, 1);
            gsidx = find(gsidx(:) >= 0.5);
            if isempty(gsidx)
                gmatrix(:, gc) = [];
                continue;
            end
            gsdata = fcdata(gsidx, :);
            gsdata(~udata(gsidx), :) = [];
            if isempty(gsdata)
                gmatrix(:, gc) = [];
            else
                gmatrix(:, gc) = ztrans(mean(gsdata)', 1);
            end
        end

        % create complete filter matrix
        fmatrix = [ones(nvols, 1), ztrans(rpdata, 1), dcts, gmatrix];

        % filter data
        fltbetas = fcdata(udata, :) * ((fmatrix' * fmatrix) \ fmatrix')';

        % remove variance
        fcdata(udata, :) = ztrans(fcdata(udata, :) - fltbetas * fmatrix', 2);

        % write back into volumes
        fcdata = reshape(fcdata, szdata);
        for vc = 1:nvols
            spm_write_vol(fcvols(vc), fcdata(:, :, :, vc));
        end

    % error handling
    catch eobj
        warning(eobj.identifier, eobj.message);
    end
end



%% sub-functions

% ztrans (see @neuroelf/private/ztrans for comments)
function ztc = ztrans(tc, dim)
ts = size(tc);
td = ts(dim);
if dim == 1
    tc = reshape(tc, ts(1), prod(ts(2:end)));
elseif dim > 2 || numel(ts) > 2
    if dim == numel(ts)
        tc = reshape(tc, prod(ts(1:dim-1)), td);
    else
        tc = reshape(tc, prod(ts(1:dim-1)), td, prod(ts(dim+1:end)));
    end
    dim = 2;
end
if dim == 1
    sref = {ones(1, td), ':'};
elseif numel(ts) == 2
    sref = {':', ones(1, td)};
else
    sref = {':', ones(1, td), ':'};
end
zsh = sum(tc, dim) ./ td;
ztc = tc - zsh(sref{:});
zf = 1 ./ sqrt(sum((1 / (td - 1)) .* (ztc .* ztc), dim));
zf(isinf(zf) | isnan(zf)) = 0;
ztc = reshape(ztc .* zf(sref{:}), ts);
