function [Cp,Cg] = spm_dcm_priorcorr(priorcorr,pV,pE,gV,gE)
% FORMAT [pC,gC] = spm_dcm_priorcorr(corr,pV,[gV])
%
% Purpose
% Inputs
%   corr    : struct array specifying pairwise correlations. Each element must have:
%             .struct         'pC' or 'gC'  (target prior: latent or observer)
%             .fieldindices   1x2 cellstr of field selectors, each selecting a
%                              SINGLE scalar parameter (as listed by spm_fieldindices)
%             .coeff          scalar ρ in [-1,1] (correlation coefficient)
%             Notes:
%               - Provide one corr(i) per pair you want to set.
%               - Pairs may overlap (i.e., share a parameter); feasibility is
%                 enforced by SPD projection, which may slightly adjust targets.
%
%   pV      : PRIOR VARIANCE structure for latent parameters (same shape as pE).
%             This is a variance *structure*, not a covariance matrix.
%
%   gV      : (optional) PRIOR VARIANCE structure for observer/noise parameters
%             (same shape as gE). If omitted or empty, gC is returned as zeros.
%
% Outputs
%   pC      : adjusted prior covariance matrix for latent parameters
%   gC      : adjsuted prior covariance matrix for observer parameters
%
% See also: spm_dcm_priorvar
%__________________________________________________________________________
% Andrew Levy

%{
% Example
%--------------------------------------------------------------------------
pC.A                      = 1./randi(32,3,3);
priorcorr                 = struct();
priorcorr(1).struct       = 'pC';
priorcorr(1).fieldindices = {'A(2,2)','A(2,1)'};
priorcorr(1).coeff        = .3;
priorcorr(2).struct       = 'pC';
priorcorr(2).fieldindices = {'A(2,1)','A(2,3)'};
priorcorr(2).coeff        = .6;
Cp                        = spm_dcm_priorcorr(priorcorr, pC);
%}

% Internal choice of SPD repair method
%   'shrinkage'  : default. Preserves diagonal exactly, shrinks off diagonals.
%   'higham'     : original SVD based projection.
method   = 'modified Cholesky';


% Argument checks
%--------------------------------------------------------------------------

if nargin < 2
    error('At least priorcorr and pV must be provided.');
end

% If pV is provided, enforce that pE is too
if nargin < 3 || isempty(pE)
    if ~isempty(pV) && (isstruct(pV) || (isnumeric(pV) && ~isscalar(pV)))
        error('pE must be supplied whenever pV is provided.');
    else
        pE = [];
    end
end

if nargin < 4
    gV = [];
    Cg = [];
end

% If gV is provided, enforce that gE is too
if nargin < 5 || isempty(gE)
    if ~isempty(gV) && (isstruct(gV) || (isnumeric(gV) && ~isscalar(gV)))
        error('gE must be supplied whenever gV is provided.');
    else
        gE = [];
    end
end


% Prepare latent model
%--------------------------------------------------------------------------
if isstruct(pV)

    pv    = spm_vec(pV);
    pC    = spm_diag(pv);
    Ip    = double(~~pC);
    Rp    = zeros(length(pv));
    pStr  = spm_fieldindices(pV,1:length(pv));
    pEig  = eig(full((pC+pC')/2));

end

% Prepare observer model
%--------------------------------------------------------------------------
Rg0 = 0;
if isnumeric(gV) && ismatrix(gV) && nnz(gV) >1
    Rg0     = spm_cov2corr(gV);
    Rg0     = Rg0-diag(diag(Rg0));
    gC      = (gV + gV')/2;
    gV      = spm_unvec(diag(gC),gE);
end

if isstruct(gV)

    gv      = spm_vec(gV);
    gC      = spm_diag(gv);
    Ig      = double(~~gC);
    Rg      = zeros(length(gv));
    gStr    = spm_fieldindices(gE,1:length(gv));
    gEig    = eig(full((gC+gC')/2));

else
    gv      = spm_vec(gV);
    gC      = spm_diag(0*gv);
    gEig    = 0;
    Eigg    = 0;
end

% Induce correlations
%--------------------------------------------------------------------------
for i = 1:numel(priorcorr)

    % Target prior
    %-----------------------------------------
    is_p =  strcmpi(priorcorr(i).struct,'pc');
    is_g =  strcmpi(priorcorr(i).struct,'gc');
    if (is_p || is_g)
        if is_p
            R   = Rp;
            Str = pStr;
            v   = pv;
        else
            R   = Rg;
            Str = gStr;
            v   = gv;
        end
    else
        error('corr(%d).struct must be ''pC'' or ''gC''.',i);
    end

    % Assign coefficient (R)
    %------------------------------------------
    r = sort(find(ismember(Str,priorcorr(i).fieldindices)));

    if ~isempty(r) && v(r(1)) && v(r(2)) && R(r(1),r(2)) == 0
        R(r(1),r(2)) = priorcorr(i).coeff;
        R(r(2),r(1)) = priorcorr(i).coeff;
    else
        % warning('%s to %s misspecified - skipping',Str{r(1)},Str{r(2)})
        continue
    end

    if is_p
        Rp = R;
    else
        Rg = R;
    end

end

Rp        = Rp + Ip;
[Cp,Eigp] = prepare_Cp(Rp,pC,pStr);

if isstruct(gV)
    Rg          = Rg0+Rg + Ig;
    [Cg,Eigg]   = prepare_Cp(Rg,gC,gStr);
end
diagnostic = true;

if diagnostic
    diagnostic_plot(pEig,Eigp,Cp,pStr)
    if isstruct(gV)
        diagnostic_plot(gEig,Eigg,Cg,gStr)
    end
end


function diagnostic_plot(pEig,Eigp,Cp,Str)

% Plot eigenvalues, which reveal SPD issues and condition
figure;
t = tiledlayout(5,1,"TileSpacing","tight","Padding","tight");
nexttile(t,[1 1])
plot(sort(pEig,'descend'),'o-'); hold on;
plot(sort(Eigp,'descend'),'.-');
xlabel('index');
ylabel('eigenvalue');
legend({'pC before','pC after'},'Location','best');
title(sprintf('pC: min-eig before %.3g, after %.3g', min(pEig), min(Eigp)));
axis square

nexttile(t,[4 1])
spy(Cp-diag(diag(Cp)),".",30)
axis square
set(gca,'FontSize',8)

drawnow

function [Cp,Eig] = prepare_Cp(R,pC,Str)

% Ensure SPD
%----------------------------------------
% R   = nearestSPD(R);

% Corr to cov
%----------------------------------------
D    = sqrt(pC);
R    = (R + R.')/2;
Cp   = (D*R*D);

% Diagnotics
%----------------------------------------
Eig = eig(full((Cp+Cp')/2));
if any(Eig < -1e-10)
    warning('pC is not SPD. Minimum eigenvalue %.3g corrected covariance to nearest SPD', min(Eig));
    R = modcholSPD(R);
    % R = shrinkageSPD(R);
    % R    = nearestSPD(R);
    D    = sqrt(pC);
    Cp   = (D*R*D);
end



Rp      = spm_cov2corr(Cp);
Rp      = tril(Rp,-1);
[rp,cp] = find(abs(Rp) > 1e-6);

fprintf('Pairs with nonzero prior correlation in pC:\n');
for k = 1:numel(rp)
    fprintf('%s with %s, rho = %+0.4f\n', ...
        Str{rp(k)}, Str{cp(k)}, Rp(rp(k),cp(k)));
end

function Cfix = modcholSPD(C, tau)
% Modified Cholesky style repair:
% 1) try Cholesky
% 2) if fails, use LDL' with permutation, lift 1x1 and 2x2 blocks minimally
% 3) reconstruct. This is local and usually small.
tau= 1e-6;
C   = (C + C.')/2;
[R,p] = chol(C);
if p == 0
    Cfix = C;
    return
end

% LDL' factorisation with permutation
[L,D,perm] = ldl(C,'vector');  % C(perm,perm) = L*D*L'
n = size(C,1);
Dplus = D;

k = 1;
while k <= n
    if k < n && D(k,k+1) ~= 0
        % 2x2 block
        B = D(k:k+1,k:k+1);
        [Q,lam] = eig((B+B')/2);
        lam = max(diag(lam), tau);
        Bp  = Q*diag(lam)*Q';
        % symmetrise to clean round-off
        Bp = (Bp + Bp.')/2;
        Dplus(k:k+1,k:k+1) = Bp;
        k = k + 2;
    else
        % 1x1 block
        Dplus(k,k) = max(D(k,k), tau);
        k = k + 1;
    end
end

% Reconstruct and unpermute
Cperm = L * Dplus * L.';
Cperm = (Cperm + Cperm.')/2;
invp  = zeros(n,1); invp(perm) = 1:n;
Cfix  = Cperm(invp,invp);

% Final symmetry and small tidy-up
Cfix = (Cfix + Cfix.')/2;

% As a belt-and-braces, if Cholesky still fails due to round-off, add tiny load
[~,pfinal] = chol(Cfix);
if pfinal ~= 0
    mineig = min(eig((Cfix+Cfix')/2));
    epsload = max(0, tau - mineig) + 1e-12*norm(Cfix,2);
    Cfix = Cfix + epsload*eye(n);
    Cfix = (Cfix + Cfix.')/2;
end


function Ahat = shrinkageSPD(A)

% Shrink towards identity with minimal alpha to reach eig_floor
% eig(R_alpha) = (1 - alpha)*eig(R) + alpha
eig_floor= 1e-10;     % small floor for strict positivity
e = eig(A);
lam_min = min(e);
if lam_min >= eig_floor
    alpha = 0;
else
    alpha = (eig_floor - lam_min) / (1 - lam_min);
    alpha = min(max(alpha,0),1);
end
Ahat = (1 - alpha)*A + alpha*eye(size(A));
% enforce symmetry and unit diagonal again
Ahat = (Ahat + Ahat.')/2;
Ahat(1:size(Ahat,1)+1:end) = 1;


function Ahat = nearestSPD(A)

% Higham method
B       = (A + A')/2;
[~,S,V] = svd(B);
H       = V*S*V';
Ahat    = (B + H)/2;
Ahat    = (Ahat + Ahat')/2;

k = 0;
while true
    [~,p] = chol(Ahat);  if p == 0, break; end

    mineig  = min(eig(Ahat));
    Ahat    = Ahat + (-mineig*eye(size(Ahat)))*(10^k);
    k       = k +1;
end

% normalise to unit diagonal
d       = 1./sqrt(diag(Ahat) + eps);
Ahat    = diag(d)*Ahat*diag(d);
Ahat(Ahat<1e-12) = 0;

%======================= FINAL DIAGNOSTIC SUMMARY =========================%
% Inputs expected in workspace:
%   Cp      final latent prior covariance
%   Cold_p  target latent covariance before repair (rename if needed)
%   tau     eigen floor used during repair
%   Optionally repeat for Cg with Cold_g

function diagnostic_summary(C, Cold, tau, name)

if nargin < 4, name = 'pC'; end

% Symmetry and SPD
Csym = (C + C')/2;
sym_err = norm(C - Csym,'fro')/max(1,norm(C,'fro'));
[~,pchol] = chol(Csym);
e = eig(Csym);
lambda_min = min(e);
lambda_max = max(e);
floor_mass = sum(e >= tau & e <= 10*tau);

% Conditioning
kappaC = condest(Csym);
R = spm_cov2corr(Csym);
kappaR = condest(R);

% Precision and partial correlations
Q = inv(Csym);
qdiag = diag(Q);
prec_min = min(qdiag);
prec_max = max(qdiag);
prec_cv  = std(qdiag)/max(mean(qdiag),eps);
Dq = diag(sqrt(max(qdiag, eps)));
Ppartial = - Dq \ Q / Dq;  % ones on diagonal
max_abs_partial = max(max(abs(tril(Ppartial,-1))));

% Correlations
Rlt = tril(R,-1);
max_abs_corr = full(max(abs(Rlt(:))));

% Deviation from target
if nargin >= 2 && ~isempty(Cold)
    relF = norm(Csym - Cold,'fro')/max(1,norm(Cold,'fro'));
else
    relF = NaN;
end

% Thresholds
thr_sym    = 1e-12;
thr_kappa  = 1e9;
thr_part   = 0.9;       % relax if scientifically intended
thr_corr   = 0.99;      % flag very high unintended correlations
thr_relF   = 0.25;      % example tolerance for global change
tau_lo     = tau;       % lambda_min should be at least tau
tau_hi     = 10*tau;    % many eigenvalues in [tau,10*tau] implies heavy regularisation

% Print summary
fprintf('\n---------------- PRIOR COVARIANCE DIAGNOSTICS [%s] ----------------\n', name);
fprintf('Cholesky              : %s (must pass)\n', ternary(pchol==0,'PASS','FAIL'));
fprintf('lambda_min            : % .3e  (>= %.1e)\n', lambda_min, tau_lo);
fprintf('lambda_max            : % .3e\n', lambda_max);
fprintf('eigs near floor       : %d      (count in [%.1e, %.1e])\n', floor_mass, tau_lo, tau_hi);
fprintf('cond2(C)              : % .3e  (<= %.1e)\n', kappaC, thr_kappa);
fprintf('cond2(R)              : % .3e  (smaller is better)\n', kappaR);
fprintf('max |corr|            : % .3f   (<= %.2f unless intended)\n', max_abs_corr, thr_corr);
fprintf('max |partial corr|    : % .3f   (<= %.2f unless intended)\n', max_abs_partial, thr_part);
fprintf('precision diag range  : [% .3e, % .3e]\n', prec_min, prec_max);
fprintf('precision diag CV     : % .3f   (smaller is better)\n', prec_cv);
if ~isnan(relF)
    fprintf('relative Fro norm dev : % .3f   (<= %.2f preferred)\n', relF, thr_relF);
end

% Decision hint
ok = (pchol==0) && (lambda_min >= tau_lo) && (kappaC <= thr_kappa) ...
    && (max_abs_partial <= thr_part || isnan(thr_part));
fprintf('Decision              : %s  (aim is stable Laplace inversion)\n', ternary(ok,'OK','CHECK'));
fprintf('-------------------------------------------------------------------\n');



function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end

return

% demo for correlated bivaraite log-normal pdf
%==========================================================================
vset    = 1./[2 16 32 128];             % log-variances (v)
rset    = [0.0 0.3 0.6 0.9];            % correlations in log space
mu      = [0;0];                        % log-means
x1      = linspace(0.1, 2, 100);
x2      = linspace(0.1, 2, 200);
[X1,X2] = meshgrid(x1,x2);

figure('Name','Bivariate log-normal: analytic PDF contours')
p=1;
for k = 1:numel(rset)
    for s = 1:numel(vset)
        v   = vset(s);
        rho = rset(k);
        Sig = [v, rho*v; rho*v, v];
        Sgi = inv(Sig);
        detS= det(Sig);
        L1  = log(X1) - mu(1);
        L2  = log(X2) - mu(2);
        % Quadratic form (vectorised): [L1 L2] * Sgi * [L1; L2]
        Q = Sgi(1,1).*L1.^2 + 2*Sgi(1,2).*L1.*L2 + Sgi(2,2).*L2.^2;
        f = exp(-0.5*Q) ./ (2*pi*sqrt(detS) .* X1 .* X2);

        subplot(numel(rset),numel(vset),p)
        contourf(X1, X2, f, 12, 'LineStyle','none')
        axis tight; axis square
        title(sprintf('v=1/%d, \\rho=%.1f', round(1/v), rho))
        xlabel('X_1'); ylabel('X_2')
        p=p+1;
    end
end
sgtitle('Bivariate log-normal density')

