function [pE,pC] = spm_L_seq_priors(dipfit,model,pE,pC)


if nargin < 3
    pE = [];
end

if nargin < 4
    pC = [];
end

if iscell(dipfit)

    for d = 1:numel(dipfit)

        [lE,lC] = spm_L_seq_priors(dipfit{d},model,pE,pC);

        oE.obs{d} = lE;
        oC.obs{d} = lC;

    end


    if numel(oE.obs) == 1
        pE = oE.obs{1};
        pC = oC.obs{1};
    else
        pE = oE;
        pC = oC;
    end
    return
end

% get priors on observation model
%--------------------------------------------------------------------------
dipfit.model = model;
switch upper(dipfit.type)

    case {'ECD','IMG','LFP'}
        [pE,pC] = spm_L_erp_priors(dipfit,pE,pC);

    case {'COM'}
        [pE,pC] = spm_L_com_priors(dipfit,pE,pC);
end

