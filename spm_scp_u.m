function [u] = spm_scp_u(t,P,M)
% [scalar] input for EEG models (Skew-Normal function)
% FORMAT [u] = spm_scp_u(t,P,M)
% t      - PST (seconds)
% P      - parameter structure
%   P.R  - scaling of [Gaussian] parameters
%
% u   - stimulus-related (subcortical) input
%
% See spm_fx_erp.m and spm_erp_priors.m
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________

% Andrew Levy



% Preliminaries - check durations (ms)
%--------------------------------------------------------------------------
try
    if length(M.dur) ~= length(M.ons)
        M.dur = M.dur(1) + M.ons - M.ons;
    end
catch
    M.dur = 32 + M.ons - M.ons;
end

% Stimulus - Skew-Normal (subcortical) impulse
%--------------------------------------------------------------------------
nu = length(M.ons);
u  = sparse(length(t),nu);
t  = t*1000;                       % ms

try
    alpha = M.alpha;
catch
    alpha = -1.25;                        % fixed skew: long left shoulder, faster right fall
end

for i = 1:nu

    % Bump function
    %----------------------------------------------------------------------
    delay   = M.ons(i) + 32*(P.R(i,1));                  % peak near onset (ms)
    f       = sqrt(exp(1) ./ (1 + exp(-(P.R(i,2)-1))));  % sqrt(e * sigmoid(z-1))
    scale   = M.dur(i) * f;                              % width scale (ms)

    Z       = (t - delay) ./ scale;
    phi     = exp(-0.5*Z.^2) / sqrt(2*pi);           % standard normal pdf
    Phi     = 0.5*(1 + erf((alpha*Z)/sqrt(2)));      % standard normal cdf
    U       = 2 * (phi./scale) .* Phi;               % skew-normal pdf

    U       = U / max(U);                        % peak normalise
    u(:,i)  = 16 * U;    


end


return

figure;
hold on
dur_span = linspace(0,2,5);
alp_span = [0 -1.5 -2.5];
str = {};
s = 1;

% Define a nice color scheme
colors = lines(length(dur_span));

for b = alp_span
    subplot(3,1,s)
    hold on
    
    for a = 1:length(dur_span)
        M.alpha = b;
        M.ons   = 300;
        M.dur   = 64;
        P.R     = zeros(1,2);
        P.R(2)  = dur_span(a);
        t       = linspace(0,.6,2000);
        [u]     = spm_scp_u(t,P,M);
        
        plot(t, u, 'LineWidth', 2, 'Color', colors(a,:), 'DisplayName', num2str(dur_span(a)))
    end
    
    % Formatting
    axis tight
    yticks([])
    % grid on
    % ylabel('Response', 'FontSize', 11, 'FontWeight', 'bold')
    title(['\alpha = ' num2str(b)], 'FontSize', 12, 'FontWeight', 'bold')
    set(gca, 'FontSize', 10, 'Box', 'on', 'LineWidth', 1.2)
    
    if s == 3
        xlabel('PST (s)', 'FontSize', 11, 'FontWeight', 'bold')
        legend('show', 'Location', 'best', 'FontSize', 9)
    end
    
    s = s + 1;
end

% Overall figure title
sgtitle('Effect of Duration Parameter across \alpha values', 'FontSize', 14, 'FontWeight', 'bold')

