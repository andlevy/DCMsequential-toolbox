function [f,J,Q] = spm_fx_seq(x,u,P,M)


% Sequential scheme
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
    seq = true;
else
    seq = false;

end

% model-specific parameters
%==========================================================================
% model or node-specific state equations of motions
%--------------------------------------------------------------------------
fx{1} = @spm_fx_erp;                     % ERP model
fx{2} = @spm_fx_cmc_2017;                     % CMC model
fx{3} = @spm_fx_bgt;                     % basal ganglia circuit
fx{4} = @spm_fx_mmc;                     % motor micro circuit
fx{5} = @spm_fx_ban3;                     % heart-rate micro circuit
fx{6} = @spm_fx_cvp;                     % heart-rate micro circuit

% indices of extrinsically coupled hidden states
%--------------------------------------------------------------------------
efferent(1,:) = [9 9 9 9];               % sources of ERP connections
afferent(1,:) = [4 8 5 8];               % targets of ERP connections

efferent(2,:) = [3 3 7 7];               % sources of CMC connections
afferent(2,:) = [2 8 4 6];               % targets of CMC connections

efferent(3,:) = [9 9 9 9];               % sources of BGT connections (thalamus)
afferent(3,:) = [2 6 2 6];               % targets of BGT connections (striatum & STN)

efferent(4,:) = [3 3 7 7];               % sources of MMC connections
afferent(4,:) = [2 4 8 0];               % targets of MMC connections

efferent(5,:) = [13 13 13 13];           % sources of BAN connections
afferent(5,:) = [2 4 2 4];               % targets of BAN connections

efferent(6,:) = [3 3 3 3];               % sources of CVP connections
afferent(6,:) = [0 0 2 2];               % targets of CVP connections


% sign of extrinsically coupled states
%--------------------------------------------------------------------------
polarity(1,:) = [1  1  1  1]; 
polarity(2,:) = [1  1 -1 -1]; 
% Sign(2,:) = [1  1  1  1];
polarity(3,:) = [1  1  1  1];
polarity(4,:) = [1  1  1  1];
polarity(5,:) = [1  1  1  1];
polarity(6,:) = [1  1  1 -1];

% scaling of afferent extrinsic connectivity (Hz)
%--------------------------------------------------------------------------
E(1,:) = [64 2 32 8]*1000;               % ERP connections
E(2,:) = [64 4 64 8]*1000;               % CMC connections
E(3,:) = [1.8 1.2 1.8 1.2]*100000;       % BGC connections
E(4,:) = [.9 .9 .11 0]*100000;           % MMC connections
E(5,:) = [.9 .9 1 0]*100000;                 % HMC connections
E(6,:) = [1 1 1 1]*10000;                  % CVP connections

if isfield(M,'ERP_E'); E(1,:)= M.ERP_E; end
if isfield(M,'CMC_E'); E(2,:)= M.CMC_E; end
if isfield(M,'BGT_E'); E(3,:)= M.BGT_E; end
if isfield(M,'MMC_E'); E(4,:)= M.MMC_E; end
if isfield(M,'BAN_E'); E(5,:)= M.BAN_E; end
if isfield(M,'CVP_E'); E(6,:)= M.CVP_E; end

% intrinsic delays ms (scaling)
%--------------------------------------------------------------------------
D(1) = 2;                                % ERP connections
D(2) = 1;                                % CMC connections
D(3) = 4;                                % BGT connections
D(4) = 1;                                % MMC connections
D(5) = 2;                                % HMC connections
D(6) = 4;                                % CVP connections

% get model specific operators
%==========================================================================
if isvector(x)
    x = spm_unvec(x,M.x);                               % neuronal states
end

% get the neural mass models {'ERP','CMC'}
%--------------------------------------------------------------------------
n     = numel(x);

model = M.model;
for i = 1:n
    if     strcmp(model(i).source,'ERP')
        nmm(i) = 1;
    elseif strcmp(model(i).source,'CMC')
        nmm(i) = 2;
    elseif strcmp(model(i).source,'BGT')
        nmm(i) = 3;
    elseif strcmp(model(i).source,'MMC')
        nmm(i) = 4;
    elseif strcmp(model(i).source,'BAN')
        nmm(i) = 5;
    elseif strcmp(model(i).source,'CVP')
        nmm(i) = 6;
    end
end

% synaptic activation function
%--------------------------------------------------------------------------
R     = 2/3;                      % gain of sigmoid activation function
B     = 0;                        % bias or background (sigmoid)
S     = @(x,R,B)1./(1 + exp(-R*x(:) + B)) - 1/(1 + exp(B));
dSdx  = @(x,R,B)(R*exp(B - R*x(:)))./(exp(B - R*x(:)) + 1).^2;

if isfield(P,'S') && iscell(P.S) && seq
    QS = spm_trajectory(P.S,v,true);
    R   = R.*QS;
else
    R   = R.*exp(P.S);              % posterior gain
end

for i = 1:n
    Sx{i} = S(x{i},R,B);          % sigmod firing rate function
end


% Extrinsic connections
%--------------------------------------------------------------------------
if iscell(P.A{1})  && seq
    QA   = vertcat(P.A{:});
    QA   = arrayfun(@(a) QA(:,a), 1:size(QA,2), 'UniformOutput', false);
    for a = 1:numel(QA)
        A{a} = spm_dcm_trajectory(QA{a},v,true);
    end
else
    for i = 1:numel(P.A)
        A{i} = exp(P.A{i});
    end
end
% detect and reduce the strength of reciprocal (lateral) connections
%--------------------------------------------------------------------------
TOL   = exp(-8);
for i = 1:numel(A)
    L    = (A{i} > TOL) & (A{i}' > TOL);
    A{i} = A{i}./(1 + 4*L);
end

% and scale of extrinsic connectivity (Hz)
%--------------------------------------------------------------------------
for i = 1:n
    for k = 1:numel(A)
        A{k}(i,:) = E(nmm(i),k)*polarity(nmm(i),k)*A{k}(i,:);
    end
end



% assemble flow
%==========================================================================
N     = M;
for i = 1:n

    % intrinsic flow
    %----------------------------------------------------------------------
    N.x  = M.x{i};
    if isfield(M,'u')
        ui = u(i,:);
    else
        ui = u;
    end
    Q    = P.int{i};

    if isfield(Q,'T') && iscell(Q.T)  && seq
        Q.T = log(spm_dcm_trajectory(Q.T,v,true));
    end

    if isfield(Q,'S') && iscell(Q.S)  && seq
        Q.S = log(spm_dcm_trajectory(Q.S,v,true));
    end

    if isfield(Q,'H') && iscell(Q.H)  && seq
        Q.H = log(spm_dcm_trajectory(Q.H,v,true));
    end

    if isfield(Q,'G') && iscell(Q.G)  && seq

        Q.G = log(spm_dcm_trajectory(Q.G,v,true));


    end

 
    Q.C  = P.C(i,:);
    f{i} = fx{nmm(i)}(x{i},ui,Q,N);


    % extrinsic flow
    %----------------------------------------------------------------------
    for j = 1:n
        for k = 1:numel(A)
            if abs(A{k}(i,j)) > TOL
                ik       = afferent(nmm(i),k);
                jk       = efferent(nmm(j),k);
                f{i}(ik) = f{i}(ik) + A{k}(i,j)*Sx{j}(jk);
            end
        end
    end
end


% concatenate flow
%--------------------------------------------------------------------------
f  = spm_vec(f);

if nargout < 2; return, end

% Jacobian
%==========================================================================
for i = 1:n
    for j = 1:n

        J{i,j} = zeros(numel(x{i}),numel(x{j}));
        if i == j

            % intrinsic flow
            %--------------------------------------------------------------
            N.x    = M.x{i};
            if isfield(M,'u')
                ui = u(i,:);
            else
                ui = u;
            end
            Q      = P.int{i};
            if isfield(Q,'T') && iscell(Q.T)  && seq
                Q.T = log(spm_dcm_trajectory(Q.T,v,true));
            end

            if isfield(Q,'S') && iscell(Q.S)  && seq
                Q.S = log(spm_dcm_trajectory(Q.S,v,true));
            end

            if isfield(Q,'H') && iscell(Q.H)  && seq
                Q.H = log(spm_dcm_trajectory(Q.H,v,true));
            end

            if isfield(Q,'G') && iscell(Q.G)  && seq
                Q.G = log(spm_dcm_trajectory(Q.G,v,true));

            end

            Q.C    = P.C(i,:);
            J{i,i} = spm_diff(fx{nmm(i)},x{i},ui,Q,N,1);

        else

            % extrinsic flow
            %--------------------------------------------------------------
            for k = 1:numel(A)
                if A{k}(i,j) > TOL
                    ik  = afferent(nmm(i),k);
                    jk  = efferent(nmm(j),k);
                    dS  = dSdx(x{j}(jk),R,B);
                    J{i,j}(ik,jk) = J{i,j}(ik,jk) + A{k}(i,j)*dS;
                end
            end
        end
    end
end

J  = spm_cat(J);


if nargout < 3; return, end


% delays
%==========================================================================
% Delay differential equations can be integrated efficiently (but
% approximately) by absorbing the delay operator into the Jacobian
%
%    dx(t)/dt     = f(x(t - d))
%                 = Q(d)f(x(t))
%
%    J(d)         = Q(d)df/dx
%--------------------------------------------------------------------------
% Implement: dx(t)/dt = f(x(t - d)) = inv(1 + D.*dfdx)*f(x(t))
%                     = Q*f = Q*J*x(t)
%--------------------------------------------------------------------------

% source specific (intrinsic) delays
%--------------------------------------------------------------------------
if iscell(P.D)
    P.D = log(spm_dcm_trajectory(P.D,v,true));
end

for i = 1:n
    P.D(i,i) = P.D(i,i) + log(D(nmm(i)));
end

% N-th order Taylor approximation to delay
%--------------------------------------------------------------------------
N     = 256;
for i = 1:8
    Q = spm_dcm_delay(P,M,J,N);
    if isfinite(norm(Q,'inf')) && (norm(Q*J,'inf') < norm(J,'inf')*4)
        break
    else
        N = N/2;
    end
end