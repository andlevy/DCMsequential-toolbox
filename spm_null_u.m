function [u] = spm_null_u(t,P,M)

% Initialize the output variable based on the input parameters
nu    = length(M.ons);
u     = sparse(length(t),nu);