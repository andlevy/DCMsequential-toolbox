function u = spm_seq_u(t, P, M)
% spm_seq_u generates a time series from 1 to 0 with a sigmoid transition
% Inputs:
%   t - time vector in seconds
%   P - parameter structure for scaling of sigmoid parameters
%       P should contain:
%           P.K - slope parameter for the sigmoid function
%   M - model structure for model specification
%       M should contain:
%           M.transons - location of the transition in milliseconds
% Output:
%   u - time series vector
%
% See spm_seq_u, spm_fx_seq.m and spm_erp_priors.m
%__________________________________________________________________________
% Copyright (C) 2024 Wellcome Trust Centre for Neuroimaging
% Andrew Levy

if isfield(M,'su') && ~isempty(M.su)

    for s = 1:numel(M.su)
        Q      = P;
        Q.R    = P.R(s,:);

        N     = M;
        N.ons = M.ons(s);
        N.dur = M.dur(s);
        u(:,s) = feval(M.su{s},t,Q,N);
    end

else
    u  = spm_erp_u(t,P,M);
end

% Parameter trajectory
%--------------------------------------------------------------------------

% Nb turn sigmoid into a function so that this part can be funcitonally
% specified

if ~isempty(P.K)

    t     = t*1000;
    v     = [];
    for i = 1:length(P.K)

        % delay  = M.ons(i+1) + 128*P.R(i+1,1);
        delay  = M.ons(i+1) + 128* (P.R(i+1,1) + P.K(i,1)); %Nb for future release to allow for variation in transition mid-point
        if isfield(M.dur,'K') && ~isempty(M.dur.K)
            scale  = M.durK(i+1) * exp(P.R(i+1,2));
        else
            scale  = M.dur(i+1) * exp(P.R(i+1,2));
        end

        % sigma = (scale/(sqrt(2)*2)) * exp(P.K(i));
        sigma  = (scale/(sqrt(2)*2)) * exp(P.K(i,2)); %Nb for future release to allow for variation in transition mid-point
        v(:,i) = 0.5*(1 + erf((t - delay)/(sqrt(2)*sigma)));

    end

    u = [u v];
    
end