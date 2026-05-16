function [QE, QC] = spm_dcm_posterior_trajectory(DCM)


M            = DCM.M;
U            = DCM.xU;
Ep           = DCM.Ep;
Cp           = DCM.Cp;

pst = (1:M.ns)*U(1).dt;
u   = spm_seq_u(pst, Ep, M);
if size(u,2) ~= size(Ep.C,2)
    Nu = size(Ep.C,2);
    v  = u(:,Nu+1:end);
end

% Extract covariances 
%==========================================================================
nevent       = numel(Ep.A);
nA           = numel(Ep.A{1});
nB           = numel(Ep.B{1});
nS           = size(Ep.A{1}{1},1);
q            = spm_find_pC(DCM);
Pstr         = spm_fieldindices(DCM.M.pE,q);

% for A and B
%--------------------------------------------------------------------------
AB = {};
for iE = 1:nevent
    for iA = 1:nA
        for iB = 1:nB
            AB{iA}{iB}{iE,1} = zeros(nS,nS);
            for iS1 = 1:nS
                for iS2 = 1:nS

                    Astr = sprintf('A{%i}{%i}(%i,%i)',iE,iA,iS1,iS2);
                    Bstr = sprintf('B{%i}{%i}(%i,%i)',iE,iB,iS1,iS2);

                    Aidx = q(find(strcmp(Pstr,Astr)));
                    Bidx = q(find(strcmp(Pstr,Bstr)));

                    if numel(Aidx) && numel(Bidx)
                        AB{iA}{iB}{iE,1}(iS1,iS2) = Cp(Aidx,Bidx);
                    end

                end
            end
        end
    end
end


% For G and B
%--------------------------------------------------------------------------
% this currently only handles the first intrinsic connection being
% modulated
if isfield(Ep,'int') && isfield(Ep.int{1},'G')

    GB = {};
    for iE = 1:nevent

        for iS = 1:nS

            % GB.int{iS}.G{iE} = spm_zeros(Ep.int{iS}.G{iE});
            for iG = 1:length(Ep.int{iS}.G{iE})
                for iB = 1:nB

                    Gstr = sprintf('int{%i}.G{%i}(1,%i)',iS,iE,iG);

                    if isfield(Ep.int{1},'B')
                        Bstr = sprintf('int{%i}.B{%i}(1,%i)',iS,iE,iG);

                    else
                        Bstr = sprintf('B{%i}{%i}(%i,%i)',iE,iB,iS,iS);
                    end

                    Gidx = q(find(strcmp(Pstr,Gstr)));
                    Bidx = q(find(strcmp(Pstr,Bstr)));

                    if numel(Gidx) && numel(Bidx)
                        GB.int{iS}.GB{iE}(1,iG) = Cp(Gidx,Bidx);
                    else
                        GB.int{iS}.GB{iE}(1,iG) = 0;
                    end

                end
            end
        end
    end


end


% Convert variances into structure and append covariances
%--------------------------------------------------------------------------
Cp         = spm_unvec(diag(Cp),Ep);

% Arrange efficacies
%--------------------------------------------------------------------------
Q.A      = vertcat(Ep.A{:});
Q.A      = arrayfun(@(a) Q.A(:,a), 1:size(Q.A,2), 'UniformOutput', false);
C.A      = vertcat(Cp.A{:});
C.A      = arrayfun(@(a) C.A(:,a), 1:size(C.A,2), 'UniformOutput', false);

Q.B      = vertcat(Ep.B{:});
Q.B      = arrayfun(@(a) Q.B(:,a), 1:size(Q.B,2), 'UniformOutput', false);
C.B      = vertcat(Cp.B{:});
C.B      = arrayfun(@(a) C.B(:,a), 1:size(C.B,2), 'UniformOutput', false);  % Fixed typo here

C.AB    = AB;


if isfield(Ep,'int') && isfield(Ep.int{1},'G')
    nI                = numel(Ep.int);
    for i = 1:nI
        Q.int{i}.G    = Ep.int{i}.G;
        C.int{i}.G    = Cp.int{i}.G;

        C.int{i}.GB   = GB.int{i}.GB;

        if isfield(Ep.int{1},'B')
            Q.int{i}.B = Ep.int{i}.B;
            C.int{i}.B = Cp.int{i}.B;
        end
    end

    C.GB    = GB;
end

% Synaptic trajectories
%--------------------------------------------------------------------------
QE = struct();
QC = struct();
nA = numel(Q.A);
nB = numel(Q.B);

for t = 1:length(v)

    for a = 1:nA
        QE.A{a}(t,:,:)        = spm_dcm_trajectory(Q.A{a},v(t,:),false);
        QC.A{a}(t,:,:)        = spm_dcm_trajectory(C.A{a},v(t,:),false);
    end

    for b = 1:nB
        QE.B{b}(t,:,:)        = spm_dcm_trajectory(Q.B{b},v(t,:),false);
        QC.B{b}(t,:,:)        = spm_dcm_trajectory(C.B{b},v(t,:),false);
    end

    for a = 1:nA
        for b = 1:nB
            Cov.AB{a}{b}(t,:,:)   = spm_dcm_trajectory(C.AB{a}{b},v(t,:),false);
        end
    end

    if isfield(Ep,'int') && isfield(Ep.int{1},'G')

        for i = 1:nI
            QE.int{i}.G(t,:)  = spm_dcm_trajectory(Q.int{i}.G,v(t,:),false);
            QC.int{i}.G(t,:)  = spm_dcm_trajectory(C.int{i}.G,v(t,:),false);
        end


        if isfield(Ep.int{1},'B')
            for i = 1:nI
                QE.int{i}.B(t,:)  = spm_dcm_trajectory(Q.int{i}.B,v(t,:),false);
                QC.int{i}.B(t,:)  = spm_dcm_trajectory(C.int{i}.B,v(t,:),false);
            end
        end

        for i = 1:nI
            for b = 1:nB
                Cov.int{i}.GB{b}(t,:) = spm_dcm_trajectory(C.int{i}.GB,v(t,:),false);
            end
        end

    end
end

if isfield(Ep,'int') && isfield(Ep.int{1},'T')
    for t = 1:length(v)
        for i = 1:nI
            QE.int{i}.T(t,:)  = spm_dcm_trajectory(Ep.int{i}.T,v(t,:),false);
            QC.int{i}.T(t,:)  = spm_dcm_trajectory(Cp.int{i}.T,v(t,:),false);
        end
    end
end


% Modulated trajectories
%--------------------------------------------------------------------------
X = [0;1];
for x = 1:length(X)
    Xx = X(x,:);
    for b = 1:nB

        for a = 1:nA
            QE.AB{x}{a} = QE.A{a} + Xx(b)*QE.B{b};
            QC.AB{x}{a} = QC.A{a} + Xx(b)^2*QC.B{b} + 2*Xx(b)*Cov.AB{a}{b}; 
        end

        if isfield(Ep,'int')
            for i = 1:nI
                if isfield(Ep.int{1},'B')
                    QE.int{i}.GB{x} = QE.int{i}.G + Xx(b)*QE.int{i}.B;
                    QC.int{i}.GB{x} = QC.int{i}.G + Xx(b)*QC.int{i}.B + 2*Xx(b)*Cov.int{i}.GB{b};
                else
                    QE.int{i}.GB{x} = QE.int{i}.G(:,1) + Xx(b)*QE.B{b}(:,i,i);
                    QC.int{i}.GB{x} = QC.int{i}.G(:,1) + Xx(b)^2*QC.B{b}(:,i,i) + 2*Xx(b)*Cov.int{i}.GB{b}(:,1);
                end

            end
        end

    end
end

