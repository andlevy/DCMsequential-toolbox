function QT = spm_dcm_synaptic_trajectory(P,M,U)




pst = (1:M.ns)*U(1).dt;
u   = spm_seq_u(pst, P, M);

if size(u,2) ~= size(P.C,2)
    Nu = size(P.C,2);
    v  = u(:,Nu+1:end);
end


% Arrange efficacies
%--------------------------------------------------------------------------
Q.A      = vertcat(P.A{:});
Q.A      = arrayfun(@(a) Q.A(:,a), 1:size(Q.A,2), 'UniformOutput', false);

Q.B      = vertcat(P.B{:});
Q.B      = arrayfun(@(a) Q.B(:,a), 1:size(Q.B,2), 'UniformOutput', false);

if isfield(P,'int')
    nI                = numel(P.int);
    for i = 1:nI
        Q.int{i}.G    = P.int{i}.G;
    end
end

% Synaptic trajectories
%--------------------------------------------------------------------------
QT = struct();
nA = numel(Q.A);
nB = numel(Q.B);

for t = 1:length(v)

    for a = 1:nA
        QT.A{a}(t,:,:)        = spm_dcm_trajectory(Q.A{a},v(t,:),false);
    end

    for b = 1:nB
        QT.B{b}(t,:,:)        = spm_dcm_trajectory(Q.B{b},v(t,:),false);
    end

    if isfield(P,'int')
        for i = 1:nI
            QT.int{i}.G(t,:)  = spm_dcm_trajectory(Q.int{i}.G,v(t,:),false);
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

            QT.AB{x}{a} = QT.A{a} + Xx(b)*QT.B{b};

        end

        if isfield(P,'int')
            for i = 1:nI
                
                QT.int{i}.GB{x} = QT.int{i}.G(:,1) + Xx(b)*QT.B{b}(:,i,i);

            end
        end
    end
end


if isfield(P,'int') && isfield(P.int{1},'T')
    for t = 1:length(v)
        for i = 1:nI
            QT.int{i}.T(t,:)  = spm_dcm_trajectory(P.int{i}.T,v(t,:),false);
        end
    end
end






