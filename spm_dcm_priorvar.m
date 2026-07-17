function [Cp,Ep,Cg] = spm_dcm_priorvar(priorvar,pV,pE,gV)
% spm_dcm_priorvar  Adjust prior variances for selected parameters
% FORMAT [pC,Ep,gC] = spm_dcm_priorvar(priorvar,pV,[pE],[gV])
%
% Each element of priorvar must have:
%   .struct       'pC' or 'gC'
%   .fieldindices 1x1 cellstr (parameter selector, spm_fieldindices format)
%   .scale        scalar > 0 (multiply variance by this factor)
%   or
%   .set          scalar >= 0 (assign this value directly)
%   .expand       logical (optional) -- only valid with .set, only for G
%                 fields. Expands pE.G and pV.G column-wise before
%                 applying the variance. For int{k}.G, also expands all
%                 corresponding int{k}.B{j} entries to match, preserving
%                 existing values. Freedom check is skipped for expanded
%                 elements.
%
% Parameters whose prior variance is zero in pV (or gV) are structurally
% fixed and are silently skipped, regardless of the requested operation.
%__________________________________________________________________________

if nargin < 3, pE = []; end
if nargin < 4, gV = []; end

Ep = pE;

% =========================================================================
% Pre-pass: expand G fields where .expand = true
% =========================================================================
for i = 1:numel(priorvar)

    if ~isfield(priorvar(i),'expand') || isempty(priorvar(i).expand) ...
            || ~priorvar(i).expand
        continue
    end

    pname = priorvar(i).fieldindices{1};
    tgt   = lower(priorvar(i).struct);

    switch tgt
        case 'pc'
            eval(['pV.' pname ' = 0;']);
            if isstruct(Ep)
                eval(['Ep.' pname ' = 0;']);
            end

        case 'gc'
            if isstruct(gV)
                eval(['gV.' pname ' = 0;']);
            end
    end

end

if isfield(pV,'int')
    for i = 1:numel(pV.int)

        if isfield(pV.int{i},'B')

            if iscell(pV.int{i}.G) && iscell(pV.int{i}.B)

                for b = 1:numel(pV.int{i}.B)
                    if length(pV.int{i}.G{b}) ~= length(pV.int{i}.B{b})
                        n = numel(pV.int{i}.G{b});
                        pV.int{i}.B{b}(n) = 0;
                        Ep.int{i}.B{b}(n) = 0;
                    end
                end

            elseif ~iscell(pV.int{i}.G) && iscell(pV.int{i}.B)

                for b = 1:numel(pV.int{i}.B)
                    if length(pV.int{i}.G) ~= length(pV.int{i}.B{b})
                        n = numel(pV.int{i}.G);
                        pV.int{i}.B{b}(n) = 0;
                        Ep.int{i}.B{b}(n) = 0;
                    end
                end

            else

                if length(pV.int{i}.G) ~= length(pV.int{i}.B)
                    n = numel(pV.int{i}.G);
                    pV.int{i}.B(n) = 0;
                    Ep.int{i}.B(n) = 0;
                end

            end
        end
    end
end

pv   = spm_vec(pV);
pStr = spm_fieldindices(pV,1:length(pv));

if isstruct(gV)
    gv   = spm_vec(gV);
    gStr = spm_fieldindices(gV,1:length(gv));
else
    gv   = spm_vec(0);
    gStr = {};
end

pv_ref = pv;
gv_ref = gv;

% =========================================================================
% Apply scalings / settings
% =========================================================================
for i = 1:numel(priorvar)

    tgt    = lower(priorvar(i).struct);
    pname  = priorvar(i).fieldindices{1};
    is_exp = isfield(priorvar(i),'expand') && ~isempty(priorvar(i).expand) ...
             && priorvar(i).expand;

    % Resolve index and check structural freedom
    %----------------------------------------------------------------------
    switch tgt
        case 'pc'
            idx = find(strcmp(pStr, pname));
            if isempty(idx)
                continue
            end
            if ~is_exp && ~pv_ref(idx)
                continue   
            end
        case 'gc'
            idx = find(strcmp(gStr, pname));
            if isempty(idx)
                continue
            end
            if ~is_exp && ~gv_ref(idx)
                continue  
            end
        otherwise
            error('priorvar(%d).struct must be ''pC'' or ''gC''.', i);
    end

    % Apply operation
    %----------------------------------------------------------------------
    if isfield(priorvar(i),'scale') && ~isempty(priorvar(i).scale)

        scale = priorvar(i).scale;
        switch tgt
            case 'pc'
                pv(idx) = pv(idx) * scale;
                fprintf('Scaled state-space prior variance: %s x %.4f\n', pname, scale);
            case 'gc'
                gv(idx) = gv(idx) * scale;
                fprintf('Scaled observer prior variance:    %s x %.4f\n', pname, scale);
        end

    elseif isfield(priorvar(i),'set') && ~isempty(priorvar(i).set)

        val = priorvar(i).set;
        switch tgt
            case 'pc'
                pv(idx) = val;
                if is_exp
                    fprintf('Expanded state-space prior variance: %s = %.4f\n', pname, val);
                else
                    fprintf('Set state-space prior variance:      %s = %.4f\n', pname, val);
                end
            case 'gc'
                gv(idx) = val;
                if is_exp
                    fprintf('Expanded observer prior variance:    %s = %.4f\n', pname, val);
                else
                    fprintf('Set observer prior variance:         %s = %.4f\n', pname, val);
                end
        end

    else
        warning('priorvar(%d) has neither .scale nor .set -- skipped.', i);
    end

end

% =========================================================================
% Reconstruct structured outputs
% =========================================================================
Cp = spm_unvec(pv,pV);
if isstruct(gV)
    Cg = spm_unvec(gv,gV);
else
    Cg = 0;
end

end