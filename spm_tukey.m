
function H = spm_tukey(n,alpha)
% Return the n-point Tukey window in a column vector
% FORMAT H = spm_tukey(n,alpha)
% n      - length of Tukey window
% alpha  - taper proportion (0 = rectangular, 1 = Hanning)
% H      - Tukey window
%__________________________________________________________________________
%
% Andrew Levy

if nargin < 2
    alpha = 0.5; % default taper proportion
end

k = (0:n-1)';   % column vector
H = ones(n,1);

if alpha <= 0
    % rectangular window
    H(:) = 1;

elseif alpha >= 1
    % Hanning window
    H = (1 - cos(2*pi*(k+1)/(n+1)))/2;

else
    % Tukey definition
    per = alpha*(n-1)/2;

    % rising cosine
    idx = k < per;
    H(idx) = 0.5*(1 + cos(pi*((2*k(idx))/(alpha*(n-1)) - 1)));

    % falling cosine
    idx = k > (n-1)*(1 - alpha/2);
    H(idx) = 0.5*(1 + cos(pi*((2*k(idx))/(alpha*(n-1)) - 2/alpha + 1)));
end
