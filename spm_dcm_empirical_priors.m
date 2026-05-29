% function DCMp = spm_dcm_empirical_priors(DCMq, DCMp, n_free)
% % Update DCM priors from posterior estimates across a set of DCMs
% %
% % Usage: DCMp = spm_dcm_empirical_priors(DCMq, DCMp, n_free)
% %
% % Inputs:
% %   DCMq   - cell array of DCMs whose posteriors will be used to estimate
% %            updated priors (e.g. from previous inversions)
% %   DCMp   - DCM to be initialised with the updated priors
% %   n_free - number of G parameters to leave free per NMM type (those with
% %            the largest total variance across sources and DCMs)
% %
% % Outputs:
% %   DCMp   - updated DCM with new prior fields:
% %            DCMp.<NMM>_G  : updated G prior values
% %            DCMp.<NMM>_T  : updated T prior values
% %            DCMp.<NMM>_j  : indices of free G parameters (highest total variance
% %                            plus any mandatory free parameters)
% %
% % Notes on free parameter selection (G parameters only):
% %   Parameters are ranked by total variance, decomposed as:
% %     Sigma_total   = Sigma_between + Sigma_within
% %     Sigma_between = empirical scatter of collapsed means across DCMs
% %     Sigma_within  = mean scatter of sources around their collapsed mean
% %                     within each DCM, averaged across DCMs
% %   This ensures parameters are left free if they vary either across
% %   inversions (between-DCM) or across NMM sources within an inversion
% %   (within-DCM), or both. Certain G parameters are additionally forced
% %   free regardless of their variance rank (see g_forced below).
% %   T parameters are all passed through unconditionally.
% %
% %   Note on source ordering: within each DCM, parameters for sources of a
% %   given NMM type are assumed to be stacked in contiguous blocks of n_p
% %   within the vectorised posterior, in the same order as inmm. This
% %   follows directly from spm_fieldindices and the str_nmm construction
% %   below, and should be verified if the SPM field-ordering convention
% %   changes.
% 
% % =========================================================================
% % SETUP
% % =========================================================================
% 
% % Equal model weights (reserved for future BMA weighting by free energy)
% w = ones(1, length(DCMq)) / length(DCMq); %#ok<NASGU>
% 
% % Default prior values for G and T per NMM type
% priors.CMC.G = [4 4 4 4 4 2 4 4 2 1 2 4 4] * 200;
% priors.CMC.T = [2 2 16 28];
% priors.MMC.G = [2 4 2 2 2 2 2 2 2 2 4 2 2 2] * 200;
% priors.MMC.T = [3 2 12 18];
% priors.BCN.G = [4 6 2 2 2 4 1 2 1 2 1 2 6 2 4 0.01 2 4 6] * 200;
% priors.BCN.T = [4 4 4 128 32 256 4 0.1];
% priors.BGT.G = [2 2 2 2 2 2 2 2 2] * 200;
% priors.BGT.T = [8 8 4 8 8];
% 
% % G parameters that are always free regardless of variance rank.
% % Forced indices are inclusive within n_free: they claim slots first and
% % the remaining slots are filled by variance ranking.
% g_forced.CMC = [7];
% g_forced.MMC = [1];
% g_forced.BCN = [1];
% g_forced.BGT = [1];
% 
% % Parameter types
% params = {'G', 'T'};
% 
% if ischar(DCMq{1})
%     DCMq = spm_dcm_load(DCMq);
% end
% 
% % =========================================================================
% % VALIDATE AND ESTABLISH GLOBAL NMM TYPE SET
% % =========================================================================
% % unmm must be fixed before the DCM loop so that the averaging and
% % plotting loops index consistent fields regardless of DCM ordering.
% % All DCMs are required to contain the same set of NMM types.
% 
% nmm_first = extractfield(DCMq{1}.options.model, 'source');
% unmm      = unique(nmm_first);
% 
% for d = 2:length(DCMq)
%     nmm_d   = extractfield(DCMq{d}.options.model, 'source');
%     unmm_d  = unique(nmm_d);
%     if ~isequal(sort(unmm), sort(unmm_d))
%         error('spm_dcm_empirical_priors: DCM %d contains a different NMM type set to DCM 1. All DCMs must share the same NMM types.', d);
%     end
% end
% 
% % =========================================================================
% % OUTER LOOP OVER DCMs: COLLAPSE WITHIN MODEL
% % =========================================================================
% for d = 1:length(DCMq)
% 
%     nmm  = extractfield(DCMq{d}.options.model, 'source'); % NMM type per source
%     Ep   = spm_vec(DCMq{d}.Ep);                           % vectorised posterior means
%     Cp   = DCMq{d}.Cp;                                    % posterior covariance
%     Pstr = spm_fieldindices(DCMq{d}.M.pE, 1:spm_length(DCMq{d}.Ep));
% 
%     % Inner loop over NMM types
%     for i = 1:length(unmm)
%         inmm = find(strcmpi(nmm, unmm{i}));
% 
%         % Loop over parameter types G and T
%         for p = 1:length(params)
%             par = params{p};
% 
%             % Build parameter name stubs for all sources of this NMM type.
%             % Sources are appended in the order of inmm, which determines
%             % the contiguous block structure assumed by the Level-1 scatter
%             % computation below.
%             str_nmm = {};
%             for j = 1:length(inmm)
%                 str_nmm = [str_nmm sprintf('int{%i}.%s', inmm(j), par)];
%             end
% 
%             % Extract parameters and joint covariance block
%             ip    = contains(Pstr, str_nmm);
%             ep    = Ep(ip);
%             cp    = Cp(ip, ip);
%             pstr  = Pstr(ip);
%             n_p   = sum(contains(pstr, sprintf('int{%i}.%s', inmm(1), par)));
%             n_s   = length(inmm);
% 
%             % GLS collapse across sources to single representative vector.
%             % A_mat encodes the assumption that all sources share a common
%             % underlying parameter vector; the collapse is precision-weighted.
%             A_mat   = repmat(eye(n_p), n_s, 1);
%             Lp      = inv(cp);
%             Lp_coll = A_mat' * Lp * A_mat;
% 
%             Ap.(par)(d).(unmm{i}) = Lp_coll \ (A_mat' * Lp * ep);
%             Bp.(par)(d).(unmm{i}) = inv(Lp_coll);
% 
%             % -----------------------------------------------------------------
%             % Level-1: between-source scatter around the collapsed mean.
%             % Captures disagreement among NMM sources of the same type within
%             % this DCM. Parameters for source j occupy the contiguous block
%             % ep((j-1)*n_p + (1:n_p)), consistent with the str_nmm ordering.
%             % If only one source exists the scatter is zero and contributes
%             % nothing to Sigma_within at the group level.
%             % -----------------------------------------------------------------
%             mu_d           = Ap.(par)(d).(unmm{i});
%             Sigma_within_d = zeros(n_p, n_p);
%             for j = 1:n_s
%                 ep_j           = ep((j - 1) * n_p + (1:n_p));
%                 resid          = ep_j - mu_d;
%                 Sigma_within_d = Sigma_within_d + resid * resid';
%             end
%             if n_s > 1
%                 Sigma_within_d = Sigma_within_d / (n_s - 1);
%             end
%             Cp_within.(par)(d).(unmm{i}) = Sigma_within_d;
%         end
%     end
% 
%     % Store free energy for this DCM (reserved for future BMA weighting)
%     Fp(d) = DCMq{d}.F;
% end
% 
% % =========================================================================
% % AVERAGE ACROSS DCMs AND UPDATE PRIORS
% % =========================================================================
% 
% % Warn if only one DCM is provided: Sigma_between cannot be estimated and
% % the ranking will be based on within-DCM scatter alone.
% if length(DCMq) == 1
%     warning('spm_dcm_empirical_priors: only one DCM provided. Between-DCM variance cannot be estimated; free parameter ranking will be based on within-DCM between-source scatter only.');
% end
% 
% for i = 1:length(unmm)
%     for p = 1:length(params)
%         par = params{p};
% 
%         % -----------------------------------------------------------------
%         % Precision-weighted average across DCMs (fixed-effects pooling)
%         % -----------------------------------------------------------------
%         Lp_total = zeros(size(Bp.(par)(1).(unmm{i})));
%         Lp_mu    = zeros(size(Ap.(par)(1).(unmm{i})));
% 
%         for d = 1:length(DCMq)
%             Lp_d     = inv(Bp.(par)(d).(unmm{i}));
%             Lp_total = Lp_total + Lp_d;
%             Lp_mu    = Lp_mu + Lp_d * Ap.(par)(d).(unmm{i});
%         end
% 
%         mu_bma = Lp_total \ Lp_mu;
%         Cp_bma = inv(Lp_total);
% 
%         % Store BMA estimates
%         Ap_bma.(par).(unmm{i}) = mu_bma;
%         Bp_bma.(par).(unmm{i}) = Cp_bma;
% 
%         % Within-model variance (retained for plotting of T parameters)
%         Var_bma.(par).(unmm{i}) = diag(Cp_bma);
% 
%         % Apply BMA average to default priors via exponential weighting
%         updated_priors.(unmm{i}).(par) = priors.(unmm{i}).(par)(:) .* exp(mu_bma);
% 
%         % -----------------------------------------------------------------
%         % Free parameter selection (G only)
%         % Rank by total variance = between-DCM scatter + mean within-DCM
%         % between-source scatter. Forced indices claim slots first; the
%         % remaining slots up to n_free are filled by variance ranking.
%         % -----------------------------------------------------------------
%         if strcmpi(par, 'G')
% 
%             % Level-2: empirical between-DCM scatter of collapsed means.
%             % With a single DCM this remains an unnormalised outer product
%             % and is set to zero so ranking falls back to Sigma_within alone.
%             Sigma_between = zeros(size(Cp_bma));
%             if length(DCMq) > 1
%                 for d = 1:length(DCMq)
%                     resid         = Ap.(par)(d).(unmm{i}) - mu_bma;
%                     Sigma_between = Sigma_between + resid * resid';
%                 end
%                 Sigma_between = Sigma_between / (length(DCMq) - 1);
%             end
% 
%             % Level-1: mean between-source scatter across DCMs.
%             % DCMs where n_s == 1 contribute a zero matrix, which correctly
%             % reflects the absence of within-DCM source disagreement.
%             Sigma_within = zeros(size(Cp_bma));
%             for d = 1:length(DCMq)
%                 Sigma_within = Sigma_within + Cp_within.(par)(d).(unmm{i});
%             end
%             Sigma_within = Sigma_within / length(DCMq);
% 
%             % Total variance and ranking
%             Sigma_total = Sigma_between + Sigma_within;
% 
%             % Store components for plotting
%             Var_between.(par).(unmm{i}) = diag(Sigma_between);
%             Var_within.(par).(unmm{i})  = diag(Sigma_within);
%             Var_total.(par).(unmm{i})   = diag(Sigma_total);
% 
%             % Forced indices claim their slots first
%             idx_forced  = g_forced.(unmm{i})(:);
%             n_remaining = max(0, n_free - length(idx_forced));
% 
%             % Fill remaining slots from variance ranking, excluding forced
%             [~, idx_sorted]         = sort(diag(Sigma_total), 'descend');
%             idx_ranked_pool         = idx_sorted(~ismember(idx_sorted, idx_forced));
%             idx_ranked              = idx_ranked_pool(1:min(n_remaining, length(idx_ranked_pool)));
% 
%             % Final free set: total count is always n_free (forced + ranked)
%             idx_ranked_new               = idx_ranked(~ismember(idx_ranked, idx_forced));
%             Idx_free.(par).(unmm{i})     = [idx_forced(:); idx_ranked_new(:)];
% 
%         else
%             % All T parameters passed through unconditionally
%             Var_between.(par).(unmm{i}) = NaN(size(Var_bma.(par).(unmm{i})));
%             Var_within.(par).(unmm{i})  = NaN(size(Var_bma.(par).(unmm{i})));
%             Var_total.(par).(unmm{i})   = NaN(size(Var_bma.(par).(unmm{i})));
%             Idx_free.(par).(unmm{i})    = (1:length(diag(Cp_bma)))';
%         end
% 
%         % Write updated priors and free parameter indices into DCMp
%         DCMp.M.([unmm{i} '_' par])    = spm_vec(updated_priors.(unmm{i}).(par))';
%         if strcmp(par,'G')
%             DCMp.M.([lower(unmm{i}) 'j']) = Idx_free.(par).(unmm{i});
%         end
%     end
% end
% 
% % =========================================================================
% % PLOT: PRIORS + VARIANCE DECOMPOSITION
% % =========================================================================
% % Layout: 4 columns per NMM type row
% %   Col 2k-1 : prior comparison for params{k}
% %   Col 2k   : variance decomposition for params{k}
% %              G: stacked between + within bars, markers on free parameters
% %                 (* = variance-ranked free; filled circle = forced free)
% %              T: within-model (Cp_bma) variance for reference only
% 
% figure;
% n_rows = length(unmm);
% n_cols = length(params) * 2;   % 4 columns
% 
% for i = 1:length(unmm)
%     for p = 1:length(params)
%         par      = params{p};
%         col_pri  = (p - 1) * 2 + 1;   % columns 1, 3
%         col_var  = (p - 1) * 2 + 2;   % columns 2, 4
%         free_idx = Idx_free.(par).(unmm{i});
% 
%         % --- Prior comparison subplot ---
%         subplot(n_rows, n_cols, (i - 1) * n_cols + col_pri);
%         hold on;
%         bar(priors.(unmm{i}).(par),          'FaceAlpha', 0.4);
%         bar(updated_priors.(unmm{i}).(par)', 'FaceAlpha', 0.4);
%         hold off;
%         legend('Original', 'Updated', 'Location', 'best');
%         ylabel(unmm{i}, 'Interpreter', 'none');
%         if i == 1
%             title(sprintf('Priors: %s', par));
%         end
%         xlabel(sprintf('%s parameter index', par));
%         grid on;
% 
%         % --- Variance subplot ---
%         subplot(n_rows, n_cols, (i - 1) * n_cols + col_var);
% 
%         if strcmpi(par, 'G')
%             % Stacked bar: within (bottom) + between (top)
%             v_within  = Var_within.(par).(unmm{i});
%             v_between = Var_between.(par).(unmm{i});
%             b = bar([v_within, v_between], 'stacked', 'FaceAlpha', 0.6);
%             b(1).FaceColor = [0.2 0.6 0.4];
%             b(2).FaceColor = [0.6 0.2 0.8];
%             hold on;
% 
%             % Build legend entries and marker plots conditionally to avoid
%             % phantom legend entries when either set is empty
%             legend_handles = [b(1) b(2)];
%             legend_labels  = {'Within-DCM', 'Between-DCM'};
% 
%             v_total    = Var_total.(par).(unmm{i});
%             idx_forced = g_forced.(unmm{i})(:);
%             idx_ranked = setdiff(free_idx, idx_forced);
% 
%             if ~isempty(idx_ranked)
%                 h_ranked = plot(idx_ranked, v_total(idx_ranked) * 1.05, ...
%                                 'k*', 'MarkerSize', 8, 'LineWidth', 1.2);
%                 legend_handles(end + 1) = h_ranked;
%                 legend_labels{end + 1}  = 'Variance-ranked';
%             end
%             if ~isempty(idx_forced)
%                 h_forced = plot(idx_forced, v_total(idx_forced) * 1.05, ...
%                                 'ko', 'MarkerSize', 8, 'LineWidth', 1.2, ...
%                                 'MarkerFaceColor', 'k');
%                 legend_handles(end + 1) = h_forced;
%                 legend_labels{end + 1}  = 'Forced free';
%             end
% 
%             hold off;
%             legend(legend_handles, legend_labels, 'Location', 'best');
%             if i == 1
%                 title(sprintf('Total variance: %s  (* ranked, o forced)', par));
%             end
% 
%         else
%             % T parameters: show within-model posterior variance for reference
%             bar(Var_bma.(par).(unmm{i}), ...
%                 'FaceColor', [0.3 0.5 0.8], 'FaceAlpha', 0.5);
%             if i == 1
%                 title(sprintf('Within-model variance: %s  (all free)', par));
%             end
%         end
% 
%         xlabel(sprintf('%s parameter index', par));
%         grid on;
%     end
% end
% 
% sgtitle('Prior update and variance decomposition by NMM type');



function DCMp = spm_dcm_empirical_priors(DCMp)
% Estimate empirical priors from a set of previously inverted DCMs and
% update DCMp accordingly.
%
% Usage: DCMp = spm_dcm_empirical_priors(DCMp)
%
% DCMp must contain DCMp.options.empirical with the following fields
% (all optional; defaults are applied where missing):
%
%   empirical.dcmpths        - cell array of DCM paths or structs to use
%   empirical.params         - parameter types (default: {'G','T'})
%   empirical.cmc.n          - number of free G params for CMC (default: 4)
%   empirical.cmc.j          - forced free G indices for CMC (default: [7])
%   empirical.mmc.n          - number of free G params for MMC (default: 4)
%   empirical.mmc.j          - forced free G indices for MMC (default: [1])
%   empirical.bcn.n          - number of free G params for BCN (default: 4)
%   empirical.bcn.j          - forced free G indices for BCN (default: [1])
%   empirical.bgt.n          - number of free G params for BGT (default: 4)
%   empirical.bgt.j          - forced free G indices for BGT (default: [1])
%
% Outputs:
%   DCMp   - updated DCM with new prior fields in DCMp.M:
%            DCMp.M.<NMM>_G  : updated G prior values
%            DCMp.M.<NMM>_T  : updated T prior values
%            DCMp.M.<nmm>j   : indices of free G parameters

% =========================================================================
% DEFAULTS
% =========================================================================

emp = DCMp.options.empirical;

if ~isfield(emp, 'params'); emp.params = {'G','T'}; end

% Per-NMM defaults: n=4, j=forced indices as defined below
nmm_defaults.cmc_n = 4;  nmm_defaults.cmc_j = [7];
nmm_defaults.mmc_n = 4;  nmm_defaults.mmc_j = [1];
nmm_defaults.bcn_n = 4;  nmm_defaults.bcn_j = [1];
nmm_defaults.bgt_n = 4;  nmm_defaults.bgt_j = [1];

nmm_list = fieldnames(nmm_defaults);
for k = 1:length(nmm_list)
    nm = nmm_list{k};
    % if ~isfield(emp, nm) || ~isfield(emp.(nm), 'n')
        % emp.(nm) = nmm_defaults.(nm).n;
    % end
    if ~isfield(emp, nm)
        emp.(nm) = nmm_defaults.(nm);
    end
end

params = emp.params;

% =========================================================================
% HARDCODED PRIOR VALUES FOR G AND T PER NMM TYPE
% =========================================================================
priors.CMC.G = [4 4 4 4 4 2 4 4 2 1 2 4 4] * 200;
priors.CMC.T = [2 2 16 28];
priors.MMC.G = [2 4 2 2 2 2 2 2 2 2 4 2 2 2] * 200;
priors.MMC.T = [3 2 12 18];
priors.BCN.G = [4 6 2 2 2 4 1 2 1 2 1 2 6 2 4 0.01 2 4 6] * 200;
priors.BCN.T = [4 4 4 128 32 256 4 0.1];
priors.BGT.G = [2 2 2 2 2 2 2 2 2] * 200;
priors.BGT.T = [8 8 4 8 8];

% =========================================================================
% LOAD DCMs
% =========================================================================
DCMq = emp.dcmpths;
if ischar(DCMq{1})
    DCMq = spm_dcm_load(DCMq);
end

% =========================================================================
% VALIDATE AND ESTABLISH GLOBAL NMM TYPE SET
% =========================================================================
nmm_first = extractfield(DCMq{1}.options.model, 'source');
unmm      = unique(nmm_first);

for d = 2:length(DCMq)
    nmm_d  = extractfield(DCMq{d}.options.model, 'source');
    unmm_d = unique(nmm_d);
    if ~isequal(sort(unmm), sort(unmm_d))
        error(['spm_dcm_empirical_priors: DCM %d contains a different ' ...
               'NMM type set to DCM 1. All DCMs must share the same ' ...
               'NMM types.'], d);
    end
end

% =========================================================================
% OVERRIDE DEFAULTS WITH EXISTING EMPIRICAL PRIORS IF PRESENT IN DCMq
% =========================================================================
for i = 1:length(unmm)
    for p = 1:length(params)
        par   = params{p};
        field = [unmm{i} '_' par];
        if isfield(DCMq{1}.M, field)
            ref = DCMq{1}.M.(field)(:)';
            for d = 2:length(DCMq)
                if isfield(DCMq{d}.M, field) && ~isequal(ref, DCMq{d}.M.(field)(:)')
                    warning('spm_dcm_empirical_priors: %s differs across DCMs; using DCM 1.', field);
                end
            end
            priors.(unmm{i}).(par) = ref;
        end
    end
end

% =========================================================================
% OUTER LOOP OVER DCMs: COLLAPSE WITHIN MODEL
% =========================================================================
for d = 1:length(DCMq)

    nmm  = extractfield(DCMq{d}.options.model, 'source');
    Ep   = spm_vec(DCMq{d}.Ep);
    Cp   = DCMq{d}.Cp;
    Pstr = spm_fieldindices(DCMq{d}.M.pE, 1:spm_length(DCMq{d}.Ep));

    for i = 1:length(unmm)
        inmm = find(strcmpi(nmm, unmm{i}));

        for p = 1:length(params)
            par = params{p};

            str_nmm = {};
            for j = 1:length(inmm)
                str_nmm = [str_nmm sprintf('int{%i}.%s', inmm(j), par)];
            end

            ip    = contains(Pstr, str_nmm);
            ep    = Ep(ip);
            cp    = Cp(ip, ip);
            pstr  = Pstr(ip);
            n_p   = sum(contains(pstr, sprintf('int{%i}.%s', inmm(1), par)));
            n_s   = length(inmm);

            % GLS collapse across sources
            A_mat   = repmat(eye(n_p), n_s, 1);
            Lp      = inv(cp);
            Lp_coll = A_mat' * Lp * A_mat;

            Ap.(par)(d).(unmm{i}) = Lp_coll \ (A_mat' * Lp * ep);
            Bp.(par)(d).(unmm{i}) = inv(Lp_coll);

            % Level-1: between-source scatter around the collapsed mean
            mu_d           = Ap.(par)(d).(unmm{i});
            Sigma_within_d = zeros(n_p, n_p);
            for j = 1:n_s
                ep_j           = ep((j - 1) * n_p + (1:n_p));
                resid          = ep_j - mu_d;
                Sigma_within_d = Sigma_within_d + resid * resid';
            end
            if n_s > 1
                Sigma_within_d = Sigma_within_d / (n_s - 1);
            end
            Cp_within.(par)(d).(unmm{i}) = Sigma_within_d;
        end
    end

    Fp(d) = DCMq{d}.F;
end

% =========================================================================
% AVERAGE ACROSS DCMs AND UPDATE PRIORS
% =========================================================================
if length(DCMq) == 1
    warning(['spm_dcm_empirical_priors: only one DCM provided. ' ...
             'Between-DCM variance cannot be estimated; ranking will ' ...
             'be based on within-DCM between-source scatter only.']);
end
for i = 1:length(unmm)

    % Retrieve per-NMM options using lowercase NMM name as field
    nmm_opts.n = emp.([lower(unmm{i}) '_n']);
    nmm_opts.j = emp.([lower(unmm{i}) '_j']);

    for p = 1:length(params)
        par = params{p};

        % Precision-weighted average across DCMs
        Lp_total = zeros(size(Bp.(par)(1).(unmm{i})));
        Lp_mu    = zeros(size(Ap.(par)(1).(unmm{i})));

        for d = 1:length(DCMq)
            Lp_d     = inv(Bp.(par)(d).(unmm{i}));
            Lp_total = Lp_total + Lp_d;
            Lp_mu    = Lp_mu + Lp_d * Ap.(par)(d).(unmm{i});
        end

        mu_bma = Lp_total \ Lp_mu;
        Cp_bma = inv(Lp_total);

        Ap_bma.(par).(unmm{i})  = mu_bma;
        Bp_bma.(par).(unmm{i})  = Cp_bma;
        Var_bma.(par).(unmm{i}) = diag(Cp_bma);

        updated_priors.(unmm{i}).(par) = priors.(unmm{i}).(par)(:) .* exp(mu_bma);

        if strcmpi(par, 'G')

            % Level-2: between-DCM scatter
            Sigma_between = zeros(size(Cp_bma));
            if length(DCMq) > 1
                for d = 1:length(DCMq)
                    resid         = Ap.(par)(d).(unmm{i}) - mu_bma;
                    Sigma_between = Sigma_between + resid * resid';
                end
                Sigma_between = Sigma_between / (length(DCMq) - 1);
            end

            % Level-1: mean between-source scatter across DCMs
            Sigma_within = zeros(size(Cp_bma));
            for d = 1:length(DCMq)
                Sigma_within = Sigma_within + Cp_within.(par)(d).(unmm{i});
            end
            Sigma_within = Sigma_within / length(DCMq);

            Sigma_total = Sigma_between + Sigma_within;

            Var_between.(par).(unmm{i}) = diag(Sigma_between);
            Var_within.(par).(unmm{i})  = diag(Sigma_within);
            Var_total.(par).(unmm{i})   = diag(Sigma_total);

            % Use existing j-index if present, otherwise rank by total variance
            jfield = [lower(unmm{i}) 'j'];
            % if isfield(DCMq{1}.M, jfield)
            %     Idx_free.(par).(unmm{i}) = DCMq{1}.M.(jfield)(:);
            % else
                idx_forced  = nmm_opts.j(:);
                n_free      = nmm_opts.n;
                n_remaining = max(0, n_free - length(idx_forced));

                [~, idx_sorted]     = sort(diag(Sigma_total), 'descend');
                idx_ranked_pool     = idx_sorted(~ismember(idx_sorted, idx_forced));
                idx_ranked          = idx_ranked_pool(1:min(n_remaining, length(idx_ranked_pool)));

                Idx_free.(par).(unmm{i}) = [idx_forced(:); idx_ranked(:)];
            % end

        else
            Var_between.(par).(unmm{i}) = NaN(size(Var_bma.(par).(unmm{i})));
            Var_within.(par).(unmm{i})  = NaN(size(Var_bma.(par).(unmm{i})));
            Var_total.(par).(unmm{i})   = NaN(size(Var_bma.(par).(unmm{i})));
            Idx_free.(par).(unmm{i})    = (1:length(diag(Cp_bma)))';
        end

        DCMp.M.([unmm{i} '_' par]) = spm_vec(updated_priors.(unmm{i}).(par))';
        if strcmp(par, 'G')
            DCMp.M.([lower(unmm{i}) 'j']) = Idx_free.(par).(unmm{i});
        end
    end
end


% =========================================================================
% PLOT: PRIORS + VARIANCE DECOMPOSITION
% =========================================================================
figure;
n_rows = length(unmm);
n_cols = length(params) * 2;

for i = 1:length(unmm)
    nmm_opts.n = emp.([lower(unmm{i}) '_n']);
    nmm_opts.j = emp.([lower(unmm{i}) '_j']);

    for p = 1:length(params)
        par      = params{p};
        col_pri  = (p - 1) * 2 + 1;
        col_var  = (p - 1) * 2 + 2;
        free_idx = Idx_free.(par).(unmm{i});

        % Prior comparison
        subplot(n_rows, n_cols, (i - 1) * n_cols + col_pri);
        hold on;
        bar(priors.(unmm{i}).(par),          'FaceAlpha', 0.4);
        bar(updated_priors.(unmm{i}).(par)', 'FaceAlpha', 0.4);
        hold off;
        legend('Original', 'Updated', 'Location', 'best');
        ylabel(unmm{i}, 'Interpreter', 'none');
        if i == 1; title(sprintf('Priors: %s', par)); end
        xlabel(sprintf('%s parameter index', par));
        grid on;

        % Variance decomposition
        subplot(n_rows, n_cols, (i - 1) * n_cols + col_var);

        if strcmpi(par, 'G')
            b = bar([Var_within.(par).(unmm{i}), Var_between.(par).(unmm{i})], ...
                    'stacked', 'FaceAlpha', 0.6);
            b(1).FaceColor = [0.2 0.6 0.4];
            b(2).FaceColor = [0.6 0.2 0.8];
            hold on;

            legend_handles = [b(1) b(2)];
            legend_labels  = {'Within-DCM', 'Between-DCM'};

            v_total    = Var_total.(par).(unmm{i});
            idx_forced = nmm_opts.j(:);
            idx_ranked = setdiff(free_idx, idx_forced);

            if ~isempty(idx_ranked)
                h = plot(idx_ranked, v_total(idx_ranked) * 1.05, ...
                         'k*', 'MarkerSize', 8, 'LineWidth', 1.2);
                legend_handles(end + 1) = h;
                legend_labels{end + 1}  = 'Variance-ranked';
            end
            if ~isempty(idx_forced)
                h = plot(idx_forced, v_total(idx_forced) * 1.05, ...
                         'ko', 'MarkerSize', 8, 'LineWidth', 1.2, ...
                         'MarkerFaceColor', 'k');
                legend_handles(end + 1) = h;
                legend_labels{end + 1}  = 'Forced free';
            end

            hold off;
            legend(legend_handles, legend_labels, 'Location', 'best');
            if i == 1; title(sprintf('Total variance: %s  (* ranked, o forced)', par)); end

        else
            bar(Var_bma.(par).(unmm{i}), 'FaceColor', [0.3 0.5 0.8], 'FaceAlpha', 0.5);
            if i == 1; title(sprintf('Within-model variance: %s  (all free)', par)); end
        end

        xlabel(sprintf('%s parameter index', par));
        grid on;
    end
end

sgtitle('Prior update and variance decomposition by NMM type');