%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%		CREATED BY:		JOCHEN WEBER
%
%		CREATED ON:		2016-11-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function JC_FCcalc(fcfiles, roifiles)
%   JC_FCCALC  Calculate functional connectivity maps
%   JC_FCCALC(FCFILES, ROIFILES) calculates the functional connectivity
%   between the mean of each ROI mask file (averaged across voxels) and
%   the functional connectivity files in FCFILES.
%
%   Note: the files MUST already be in the same space!

% requires CHAR or CELL input
if ischar(fcfiles)
    fcfiles = cellstr(fcfiles);
elseif ~iscell(fcfiles)
    error('JCfunc:calc:invalidArgument', 'Invalid FCFILES argument.');
end
fcfiles = fcfiles(:);
if ~all(cellfun(@ischar, fcfiles)) || any(cellfun('isempty', fcfiles)) || ...
    any(cellfun('isempty', regexpi(fcfiles, '\.nii$'))) || ...
    sum(cellfun(@exist, fcfiles)) ~= (2 * numel(fcfiles))
    error('JCfunc:calc:invalidArgument', 'Invalid FCFILES argument.');
end
if ischar(roifiles)
    roifiles = cellstr(roifiles);
elseif ~iscell(roifiles)
    error('JCfunc:calc:invalidArgument', 'Invalid ROIFILES argument.');
end
roifiles = roifiles(:);
if ~all(cellfun(@ischar, roifiles)) || any(cellfun('isempty', roifiles)) || ...
    any(cellfun('isempty', regexpi(roifiles, '\.nii$'))) || ...
    sum(cellfun(@exist, roifiles)) ~= (2 * numel(roifiles))
    error('JCfunc:calc:invalidArgument', 'Invalid FCFILES argument.');
end

% iterate over ROI files
roiinfo = cell(numel(roifiles), 4);
for fc = 1:numel(roifiles)
    roivoi = spm_vol([roifiles{fc} ',1']);
    roidata = spm_read_vols(roivoi);
    roiinfo{fc, 1} = find(roidata(:) >= 0.5);
    roiinfo{fc, 2} = roivoi.dim;
    roiinfo{fc, 3} = roivoi.mat;
    [roipath, roiinfo{fc, 4}] = fileparts(roifiles{fc});
end
if fc > 1 && (any(any(diff(cat(3, roiinfo{:, 2}), 1, 3), 3), 2) || ...
    any(any(any(diff(cat(3, roiinfo{:, 3}), 1, 3), 1), 2)))
    error('JCfunc:calc:invalidArgument', 'ROI files must match in space.');
end

% iterate over FC files
for fc = 1:numel(fcfiles)

    % try
    try

        % for new filename
        [fcpath, fcname] = fileparts(fcfiles{fc});
        if isempty(fcpath)
            fcpath = pwd;
        end

        % load file
        fcvols = spm_vol(fcfiles{fc});
        nvols = numel(fcvols);
        if nvols < 20 || ~isequal(fcvols(1).dim, roiinfo{1, 2}) || ...
           ~isequal(fcvols(1).mat, roiinfo{1, 3})
            warning('Jcfunc:calc:invalidFile', 'Bad FC file: %s', fcfiles{fc});
            continue;
        end

        % load data
        fcdata = spm_read_vols(fcvols);

        % transpose
        szdata = size(fcdata);
        fcdata = reshape(fcdata, prod(szdata(1:3)), szdata(4));
        udata = ~any(isinf(fcdata) | isnan(fcdata), 2) & any(fcdata ~= 0, 2);

        % for each ROI
        for rc = 1:numel(roifiles)

            % update udata
            mdata = false(size(udata));
            mdata(roiinfo{rc, 1}) = true;
            mdata = mdata & udata;

            % don't compute if no good voxels for this ROI
            if ~any(mdata(:))
                rdata = reshape(NaN .* ones(size(mdata)), szdata(1:3));

            % compute
            else
                rdata = ztrans(mean(fcdata(mdata(:), :))', 1);
                rdata = reshape(atanh((1 / (nvols - 1)) .* ...
                    sum(fcdata * sparse(1:nvols, 1:nvols, rdata, nvols, nvols, nvols), 2)), ...
                    szdata(1:3));
            end

            % create volume
            rvolname = [fcpath filesep 'z' fcname '_' roiinfo{rc, 4} '.nii'];
            fprintf('Writing %s...\n', rvolname);
            rvol = spm_create_vol( ...
                struct('fname', rvolname, 'dim', szdata(1:3), 'dt', [64, 0], ...
                'pinfo', [1; 0; 352], 'mat', roiinfo{1, 3}, 'n', [1, 1], ...
                'descrip', sprintf('FC (%s:%s)', fcname, roiinfo{rc, 4})));
            spm_write_vol(rvol, rdata);
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
