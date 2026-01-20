function [y] = spm_int_S(P,M,U)
% Integrate a MIMO nonlinear system (strategic updates at input events)
% FORMAT [y] = spm_int_S(P,M,U)
% P   - model parameters
% M   - model structure
% U   - input structure or matrix
%
% y   - (v x l)  response y = g(x,u,P)
%__________________________________________________________________________
% Integrates the MIMO system described by
%  
%    dx/dt = f(x,u,P,M)
%    y     = g(x,u,P,M)
% 
% using the update scheme:
%
%    x(t + dt) = x(t) + U*dx(t)/dt
%
%            U = (expm(dt*J) - I)*inv(J)
%            J = df/dx
%
% at strategically chosen input times. This integration scheme pre-computes
% update points based on input structure (Gaussian bumps, sigmoids) and
% evaluates the Jacobian at the onset, offset, and intermediate points
% during input transitions. This provides a balance between spm_int_U
% (updates when input changes) and spm_int_L (single Jacobian evaluation).
%
% If f returns the Jacobian (i.e. [fx J] = feval(f,M.x,u,P,M) it will
% be used. Otherwise it is evaluated numerically.
%
% spm_int will also handle static observation models by evaluating
% g(x,u,P,M)
%
% Configuration (optional fields in M):
%   M.nJ  - number of updates per input event (default: 6)
%   M.th  - threshold for detecting input events (default: 1e-6)
%
%--------------------------------------------------------------------------
%
% SPM solvers or integrators
%
% spm_int_ode:  uses ode45 (or ode113) which are one and multi-step solvers
% respectively. They can be used for any ODEs, where the Jacobian is
% unknown or difficult to compute; however, they may be slow.
%
% spm_int_J: uses an explicit Jacobian-based update scheme that preserves
% nonlinearities in the ODE: dx = (expm(dt*J) - I)*inv(J)*f. If the
% equations of motion return J = df/dx, it will be used; otherwise it is
% evaluated numerically, using spm_diff at each time point. This scheme is
% infallible but potentially slow, if the Jacobian is not available (calls
% spm_dx).
%
% spm_int_E: As for spm_int_J but uses the eigensystem of J(x(0)) to eschew 
% matrix exponentials and inversion during the integration. It is probably
% the best compromise, if the Jacobian is not available explicitly.
%
% spm_int_B: As for spm_int_J but uses a first-order approximation to J
% based on J(x(t)) = J(x(0)) + dJdx*x(t).
%
% spm_int_L: As for spm_int_B but uses J(x(0)).
%
% spm_int_U: like spm_int_J but only evaluates J when the input changes.
% This can be useful if input changes are sparse (e.g., boxcar functions).
% It is used primarily for integrating EEG models
%
% spm_int_S: like spm_int_U but pre-computes strategic update points based
% on input structure. Updates the Jacobian at onset, offset, and 
% intermediate points of input events (Gaussian bumps, sigmoids). Provides
% efficiency gains for structured inputs whilst maintaining accuracy during
% transitions.
%
% spm_int:   Fast integrator that uses a bilinear approximation to the
% Jacobian evaluated using spm_bireduce. This routine will also allow for
% sparse sampling of the solution and delays in observing outputs. It is
% used primarily for integrating fMRI models
%__________________________________________________________________________
 
% Andrew Quinn
% Copyright (C) 2024 Wellcome Centre for Human Neuroimaging


% convert U to U.u if necessary
%--------------------------------------------------------------------------
if ~isstruct(U), U.u = U; end
try, dt = U.dt; catch, dt = 1; end
try, ns = M.ns; catch, ns = length(U.u); end
try, nJ = M.nJ; catch, nJ = 5;    end
 
if nargin < 4; N = 1;              end

% add [0] states if not specified
%--------------------------------------------------------------------------
try
    f   = spm_funcheck(M.f);
catch
    f   = @(x,u,P,M) sparse(0,1);
    M.n = 0;
    M.x = sparse(0,0);
end
M.f = f;

% output nonlinearity, if specified
%--------------------------------------------------------------------------
try
    g = spm_funcheck(M.g);
    if isempty(g)
        g  = @(x,u,P,M) x;
    end
catch
    g = @(x,u,P,M) x;
end
M.g = g;

% Initial states and inputs
%--------------------------------------------------------------------------
try
    u   = U.u;
catch
    u   = sparse(1,M.m);
end

try
    x   = M.x;
catch
    x   = sparse(0,1);
    M.x = x;
end

% check for delay operator
%--------------------------------------------------------------------------
try
    [fx,dfdx,D] = f(x,u,P,M);
catch
    D = 1;
end

% Compute strategic update indices
%--------------------------------------------------------------------------
if size(u,2) ~= size(P.C,2)

    if size(u,2) < size(P.C,2)
        % detect bireduce
        u = u;
        v = 0;
    else
        Nu = size(P.C,2);
        v  = u(:,Nu+1:end);
        u  = u(:,1:Nu);
    end

    % updates triggered by parameter transitions
    %-------------------------------------------
    [~, on]   = max(v > 1e-6, [], 1);
    [~, mid]  = min(abs(v-.5),[],1);
    [~, off]  = max(v > 1-1e-6, [], 1);
    % uJ       = [on,mid,off];
    uJ       = [off];

    % % updates triggered by input onset and offset
    % %--------------------------------------------
    % du       = sqrt(diff([u(1,:);u]).^2);
    % du       = du > 1e-6;
    % [on,~]   = find(diff([zeros(1,Nu); du]) == 1);
    % [off,~]  = find(diff([du; zeros(1,Nu)]) == -1);
    % uJ       = [uJ,on',off'];

    uJ       = sort(unique(uJ));

else
    uJ = [];

end


cJ = false; w = true;
% integrate
%--------------------------------------------------------------------------
v     = spm_vec(x);
for i = 1:ns
    
    % input
    %----------------------------------------------------------------------
    try
        u = U.u(i,:);
    end

    % re-compute Jacobian at strategic points or when poorly conditioned
    %----------------------------------------------------------------------
    if i == 1 || ismember(i, uJ)  || cJ


        % dx(t)/dt and Jacobian df/dx (and check for delay operator)
        %------------------------------------------------------------------
        D   = 1;
        n   = spm_length(x);
        if nargout(f) >= 3
            [fx,dfdx,D] = f(x,u,P,M);
        elseif nargout(f) == 2
            [fx,dfdx]   = f(x,u,P,M);
        else
            dfdx = spm_cat(spm_diff(f,x,u,P,M,1));
        end

        % local linear update operator Q = (expm(dt*J) - I)*inv(J)
        %------------------------------------------------------------------
        dfdx  = dfdx - speye(n,n)*exp(-16);
        Q     = (spm_expm(dt*D*dfdx/N) - speye(n,n))/dfdx;
        k     = 0;

        % try
        %     [fx,dfdx] = feval(f,x,u,P,M);
        % catch
        %     dfdx = spm_diff(f,x,u,P,M,1);
        % end

        % % Comprehensive Jacobian health checks
        % %------------------------------------------------------------------
        % kappa = cond(full(dfdx));
        % nJ    = norm(dfdx, 'fro');
        % 
        % % Check 1: Condition number (catches singularity, near-singularity, rank deficiency)
        % bad_cond = kappa > 1e10 || isnan(kappa) || isinf(kappa);
        % 
        % % Check 2: Contains NaN or Inf (catches computational failures)
        % bad_vals = any(isnan(dfdx(:))) || any(isinf(dfdx(:)));
        % 
        % % Check 3: Norm too large or too small (catches extreme scaling)
        % bad_norm = nJ > 1e10 || nJ < 1e-10;
        % 
        % % Check 4: Eigenvalue analysis (catches stiffness - different timescales)
        lambda      = eig(full(dfdx));
        lambda      = lambda(abs(lambda) > 1e-12);
        if ~isempty(lambda)
            stiff     = max(abs(lambda)) / min(abs(lambda));
            bad_stiff = stiff > 1e4;
        else
            bad_stiff = false;
        end

        % Set flag
        if bad_stiff %||bad_cond || bad_vals || bad_norm 
            cJ = true;
            if  w
                warning('Jacobian stiff');
                w = false;
            end
        else
            cJ = false;
        end
    end



     try
        % update dx = (expm(dt*J) - I)*inv(J)*f(x,u)
        %------------------------------------------------------------------
        for j = 1:N
            v = v + Q*f(v,u,P,M);
        end
        
        % output - implement g(x)
        %------------------------------------------------------------------
        y(:,i) = g(v,u,P,M);

    catch
        
        % update dx = (expm(dt*J) - I)*inv(J)*f(x,u)
        %------------------------------------------------------------------
        for j = 1:N
            x = spm_vec(x) + Q*spm_vec(f(x,u,P,M));
            x = spm_unvec(x,M.x);
        end
        
        % output - implement g(x)
        %------------------------------------------------------------------
        y(:,i) = spm_vec(g(x,u,P,M));
 
        
     end
    % 
    % 
    % % dx(t)/dt
    % %----------------------------------------------------------------------
    % fx = f(x,u,P,M);
    % 
    % % update dx = (expm(dt*J) - I)*inv(J)*fx = Q*fx;
    % %----------------------------------------------------------------------
    % x  = spm_unvec(spm_vec(x) + spm_dx(D*dfdx,D*fx,dt),x);
    % 
    % % output - implement g(x)
    % %----------------------------------------------------------------------
    % if ~isempty(g)
    %     y(:,i) = g(x,u,P,M);
    % else
    %     y(:,i) = spm_vec(x);
    % end
    % 
end
 
% transpose
%--------------------------------------------------------------------------
y = real(y');