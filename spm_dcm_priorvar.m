function [Cp,Cg] = spm_dcm_priorvar(priorvar,pV,gV)
% spm_dcm_priorvar  Adjust prior variances for selected parameters
% FORMAT [pC,gC] = spm_dcm_priorvar(priorvar,pV,[gV])
%
% Each element of priorvar must have:
%   .struct       'pC' or 'gC'
%   .fieldindices 1x1 cellstr (parameter selector, spm_fieldindices format)
%   .scale        scalar > 0 (multiply variance)
%   or
%   .set          

if nargin < 3, gV = []; end

% Vectorise latent
pv   = spm_vec(pV);
pStr = spm_fieldindices(pV,1:length(pv));

% Vectorise observer if provided
if isstruct(gV)
    gv   = spm_vec(gV);
    gStr = spm_fieldindices(gV,1:length(gv));
else
    gv   = spm_vec(0);
    gStr = {};
end

% Apply scalings
for i = 1:numel(priorvar)
    tgt   = lower(priorvar(i).struct);
    pname = priorvar(i).fieldindices{1};

    if isfield(priorvar(i),'scale') && ~isempty(priorvar(i).scale)
    
        scale = priorvar(i).scale;

        if strcmp(tgt,'pc')
            idx = find(strcmp(pStr, pname));
            if ~isempty(idx)
                pv(idx) = pv(idx) * scale;
                fprintf('Scaled state-space prior variance: %s × %.3f\n', pname, scale);
            end
        elseif strcmp(tgt,'gc')
            if ~isempty(gStr)
                idx = find(strcmp(gStr, pname));
                if ~isempty(idx)
                    gv(idx) = gv(idx) * scale;
                    fprintf('Scaled observer prior variance: %s × %.3f\n', pname, scale);
                end
            end
        else
            error('priorvar(%d).struct must be ''pC'' or ''gC''.', i);
        end


    elseif isfield(priorvar(i),'set') && ~isempty(priorvar(i).set)

        set = priorvar(i).set;

        if strcmp(tgt,'pc')
            idx = find(strcmp(pStr, pname));
            if ~isempty(idx)
                pv(idx) = set;
                fprintf('Set state-space prior variance: %s - %.3f\n', pname, set);
            end
        elseif strcmp(tgt,'gc')
            if ~isempty(gStr)
                idx = find(strcmp(gStr, pname));
                if ~isempty(idx)
                    gv(idx) = gv(idx) * set;
                    fprintf('Set observer prior variance: %s - %.3f\n', pname, set);
                end
            end
        else
            error('priorvar(%d).struct must be ''pC'' or ''gC''.', i);
        end


    end

    
end

% Reconstruct outputs
Cp = spm_unvec(pv,pV);
if isstruct(gV)
    Cg = spm_unvec(gv,gV);
else
    Cg = 0;
end
end
