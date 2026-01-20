function [U] = spm_dcm_seq_channelmodes(dipfit,Nm,xY,CVA)


% number of channels and modes
%--------------------------------------------------------------------------
if nargin < 2, Nm = 8; end
if ~iscell(dipfit), dipfit = {dipfit}; end

% Determine maximum number of channels across all modalities
Nc_max = 0;
for d = 1:numel(dipfit)
    Nc_max = max(Nc_max, max(dipfit{d}.Iu));
end

% Build U with proper indexing
uc = {};
for d = 1:numel(dipfit)
    switch dipfit{d}.modality
        case {'EEG','MEG','MEGGRAD'}
            if CVA
                u = spm_dcm_eeg_channelmodes(dipfit{d}, Nm, xY);
            else
                u = spm_dcm_eeg_channelmodes(dipfit{d}, Nm);
            end
        case 'HR'
            % Nc = length(xY.2 Ic);
            % u = 1/(Nc_max-1);
              u = 1/sqrt(Nc_max);
    end
    
    % Create full-size matrix with u placed at correct rows
    [~, Nc_u]              = size(u);
    uc                  = sparse(Nc_max, Nc_u);
    uc(dipfit{d}.Iu, :) = u;
    
    if d == 1
        U = uc;
    else
        U = cat(2, U, uc);
    end
end
