function [Q] = spm_seq_Q(P,X)
% Helper routine for spm_gen routines
% FORMAT [Q] = spm_gen_Q(P,X)
%
% P - parameters
% X - vector of between trial effects
% c - trial in question
%
% Q - trial or condition-specific parameters
%
% This routine computes the parameters of a DCM for a given trial, where
% trial-specific effects are deployed according to a design vector X. The
% parameterisation follows a standard naming protocol where, for example,
% X(1)*P.B{1} + X(2)*P.B{2}... adjusts P.A for all (input) effects encoded
% in P.B.
% P.BN and P.AN operate at NMDA receptors along extrinsic connections
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


% condition or trial specific parameters
%==========================================================================
if isfield(P,'B')
    Q = rmfield(P,'B');
else
    Q = P;
end


% trial-specific effects on C (first effect only)
%--------------------------------------------------------------------------
try
    Q.C = Q.C(:,:,1) + X(1)*P.C(:,:,2);
end

% trial-specific effects on A (connections)
%--------------------------------------------------------------------------
for i = 1:length(X)

    % extrinsic (driving) connections
    %----------------------------------------------------------------------
    if iscell(Q.A{1})

        for k = 1:numel(Q.A)
            for j = 1:length(Q.A{k})

                Q.A{k}{j} = Q.A{k}{j} + X(i)*P.B{k}{i};

                % CMM-NMDA specific modulation on extrinsic NMDA connections
                %------------------------------------------------------------------
                if isfield(P,'AN')
                    Q.AN{k}{j} = Q.AN{k}{j} + X(i)*P.BN{k}{i};
                end

            end
        end

    else

        for j = 1:length(Q.A)

            Q.A{j} = Q.A{j} + X(i)*P.B{i};

            % CMM-NMDA specific modulation on extrinsic NMDA connections
            %------------------------------------------------------------------
            if isfield(P,'AN')
                Q.AN{j} = Q.AN{j} + X(i)*P.BN{i};
            end
        end

    end

    % modulatory connections
    %----------------------------------------------------------------------
    if isfield(P,'M')
        if iscell(Q.M)
            for k = 1:numel(Q.M)
                Q.M{k} = Q.M{k} + X(i)*P.N{k}{i};
            end
        else
            Q.M  = Q.M + X(i)*P.N{i};
        end
    end

    % intrinsic connections
    %----------------------------------------------------------------------
    if isfield(Q,'G')

        if iscell(Q.G)
            for k = 1:numel(Q.G)
                Q.G{k}(:,1) = Q.G{k}(:,1) + X(i)*diag(P.B{k}{i});
            end
        else
            Q.G(:,1) = Q.G(:,1) + X(i)*diag(P.B{i});
        end

    end

    % intrinsic connections
    %----------------------------------------------------------------------
    if isfield(Q,'int')
        for j = 1:numel(Q.int)

            if isfield(Q.int{j},'B') && ~isempty(Q.int{j}.B) && iscell(Q.int{j}.B)

                for k = 1:numel(Q.int{j}.B)
                    % if isempty(Q.int{j}.B{k})
                    %     Q.int{j}.B{k} = spm_zeros(Q.int{j}.G{k});
                    % end
                    % 
                    % if numel(Q.int{j}.B{k}) < numel(Q.int{j}.G{k})
                    %     haha = 1;
                    % end
                    Q.int{j}.G{k} = Q.int{j}.G{k} + X(i)*Q.int{j}.B{k};
                end

            elseif isfield(Q.int{j},'G') && iscell(Q.int{j}.G)

                for k = 1:numel(Q.int{j}.G)
                    Q.int{j}.G{k}(:,1) = Q.int{j}.G{k}(:,1) + X(i)*P.B{k}{i}(j,j);
                end
                
            else

                if isfield(Q.int{j},'B') && ~isempty(Q.int{j}.B)
                    Q.int{j}.G = Q.int{j}.G + X(i)*Q.int{j}.B;
                else
                    if iscell(P.B{i})
                        G = {};
                        for k = 1:numel(P.B)
                            G{k}(:,1) = Q.int{j}.G(:,1) + X(i)*P.B{k}{i}(j,j);
                        end
                        Q.int{j}.G = G;
                    else
                        Q.int{j}.G(:,1) = Q.int{j}.G(:,1) + X(i)*P.B{i}(j,j);
                    end
                end
            end

        end
    end

end
