function [L] = spm_lx_seq(P,dipfit)


% extract dipfit from model if necessary
%--------------------------------------------------------------------------
if  isfield(dipfit,'dipfit'), dipfit = dipfit.dipfit; end
% 
% if ~iscell(dipfit), dipfit = {dipfit}; end

for d = 1:numel(dipfit)

    switch dipfit{d}.type

        case {'ECD','IMG','LFP'}

            l =   spm_lx_erp(P.obs{d},dipfit{d});

        case 'COM'
            l =  spm_lx_com(P.obs{d},dipfit{d});


            % l =   spm_vec(P.obs{d}.J)';
    end

    if d == 1
        L = l;
    else
        L = [L; l]; 
    end

end


