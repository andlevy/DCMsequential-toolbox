function DCMp = spm_dcm_peb_nmm(DCMp)
% *** WORK IN PROGRESS -- INCOMPLETE / NOT FULLY FUNCTIONAL ***
% Not currently used in any pipeline. Known issues before this can be
% relied on:
%   - M.pC handling does not branch on the per-condition (cell) case the
%     way M.pE/Ep do; when M.pC is a struct and the parameter is not
%     stored per-condition, the {c} indexing below will error.
%   - The second-level design matrix X is intercept-only. Xmodel is
%     computed but not included, and the region-coding path (Xregion,
%     spm_dcm_peb_nmm_region_name) is fully commented out / dead. The
%     "controlling for region and originating DCM" description below is
%     therefore aspirational, not yet implemented.
% Do not use for real inference until these are resolved.
% *******************************************************************
%
% Estimate empirical priors from a set of previously inverted DCMs using PEB
% run separately
% for each neural mass model (NMM) type, with the second-level design
% controlling for region and originating DCM ("model").
%
% Usage: DCMp = spm_dcm_peb_nmm(DCMp)
%
% DCMp must contain DCMp.options.empirical with the following fields
% (all optional; defaults are applied where missing):
%
%   empirical.dcmpths  - cell array of DCM paths or structs to use
%   empirical.params   - parameter types to pool, default {'G','T'}
%   empirical.cmc_n     - number of free G params for CMC (default: 4)
%   empirical.cmc_j     - forced free G indices for CMC (default: [7])
%   empirical.mmc_n     - number of free G params for MMC (default: 4)
%   empirical.mmc_j     - forced free G indices for MMC (default: [1])
%   empirical.bcn_n     - number of free G params for BCN (default: 4)
%   empirical.bcn_j     - forced free G indices for BCN (default: [1])
%   empirical.bgt_n     - number of free G params for BGT (default: 4)
%   empirical.bgt_j     - forced free G indices for BGT (default: [1])
%
% Outputs:
%   DCMp   - updated DCM with new prior fields in DCMp.M:
%            DCMp.M.<NMM>_G  : updated G prior values
%            DCMp.M.<NMM>_T  : updated T prior values
%            DCMp.M.<nmm>j   : forced-plus-ranked free G indices
%   DCMp.PEB.<NMM>  - the full PEB structure for that NMM type (region
%                      and model effect estimates, second-level
%                      covariance, free energy), for inspection.


% =========================================================================
% DEFAULTS
% =========================================================================
emp = DCMp.options.empirical;
if ~isfield(emp,'params'); emp.params = {'G','T'}; end
params = emp.params;

% Per-NMM defaults for the forced-plus-ranked free G parameter selection,
% matching the GLS script's configuration interface exactly.
nmm_defaults.cmc_n = 4;  nmm_defaults.cmc_j = [7];
nmm_defaults.mmc_n = 4;  nmm_defaults.mmc_j = [1];
nmm_defaults.bcn_n = 4;  nmm_defaults.bcn_j = [1];
nmm_defaults.bgt_n = 4;  nmm_defaults.bgt_j = [1];

nmm_list = fieldnames(nmm_defaults);
for k = 1:length(nmm_list)
    nm = nmm_list{k};
    if ~isfield(emp, nm)
        emp.(nm) = nmm_defaults.(nm);
    end
end

% =========================================================================
% HARDCODED BASELINE PRIOR VALUES FOR G AND T PER NMM TYPE
% =========================================================================
priors.CMC.G = [4 4 4 4 4 2 4 4 2 1 2 4 4] * 200;
priors.CMC.j = [7 4 12 9 11 6 5 8 1 3 10 13 2];
priors.CMC.T = [2 2 16 28];

priors.MMC.G = [2 4 2 2 2 2 2 2 2 2 4 2 2 2] * 200;
priors.MMC.j = 1:14;
priors.MMC.T = [3 2 12 18];

priors.BCN.G = [4 4 2 2 2 4 .125 2 1 2 1 2 4 2 4 0.01 2 4 6] * 200;
priors.BCN.j = 1:19;
priors.BCN.T = [4 4 4 128 12 256 4 .1];

priors.BGT.G = [2 2 2 2 2 2 2 2 2] * 200;
priors.BGT.j = 1:9;
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
        error(['spm_dcm_peb_nmm: DCM %d contains a different NMM type ' ...
            'set to DCM 1. All DCMs must share the same NMM types.'], d);
    end
end

% =========================================================================
% OVERRIDE BASELINE PRIORS WITH EXISTING EMPIRICAL PRIORS IF PRESENT
% =========================================================================
for i = 1:length(unmm)
    for p = 1:length(params)

        par   = params{p};
        field = [unmm{i} '_' par];

        if isfield(DCMq{1}.M, field)

            ref = DCMq{1}.M.(field)(:)';
            for d = 2:length(DCMq)
                if isfield(DCMq{d}.M, field) && ~isequal(ref, DCMq{d}.M.(field)(:)')
                    warning('spm_dcm_peb_nmm: %s differs across DCMs; using DCM 1.', field);
                end
            end

            priors.(unmm{i}).(par) = ref;

        end
    end
end

% =========================================================================
% PEB PER NMM TYPE
% =========================================================================
updated_priors = struct();
Idx_free       = struct();
Var_total      = struct();

for d = 1:length(DCMq)
    Pstr{d} = spm_fieldindices(DCMq{d}.M.pE, 1:spm_length(DCMq{d}.Ep));
end

for i = 1:length(unmm)

    nmtype     = unmm{i};
    nmm_opts.n = emp.([lower(nmtype) '_n']);
    nmm_opts.j = emp.([lower(nmtype) '_j']);

    % ---------------------------------------------------------------
    % Build the first-level array: one synthetic single-source "DCM"
    % per instance of this NMM type, across every DCM in the set.
    % ---------------------------------------------------------------
    for p = 1:length(params)

        P_DCM = {}; region_label = {}; model_id = [];

        for d = 1:length(DCMq)

            nmm_d = extractfield(DCMq{d}.options.model, 'source');
            nmm_i = find(strcmpi(nmm_d, nmtype));

            if isempty(nmm_i); continue; end

            pstr   = Pstr{d};
            jfield = [lower(nmtype) 'j'];
            if isfield(DCMq{d}.M, jfield)
                j_d = DCMq{d}.M.(jfield)(:)';
            else
                j_d = priors.(nmtype).j;
            end

            for kk = 1:numel(nmm_i)
                k = nmm_i(kk);

                if ~isfield(DCMq{d}.M,'pE') || ~isfield(DCMq{d}.M.pE,'int') || numel(DCMq{d}.M.pE.int) < k
                    warning(['spm_dcm_peb_nmm: DCM %d source %d (%s) is missing M.pE.int{%d}; ' ...
                        'skipping this source.'], d, k, nmtype, k);
                    continue
                end

                ip_full = {};
                if iscell(DCMq{d}.M.pE.int{k}.(params{p}))
                    n_cell = numel(DCMq{d}.M.pE.int{k}.(params{p}));
                    for c = 1:n_cell
                        ip_full{c} = find(contains(pstr, sprintf('int{%d}.%s{%i}', k,params{p},c)));
                    end
                else
                    n_cell     = 0;
                    ip_full{1} = find(contains(pstr, sprintf('int{%d}.%s', k,params{p})));
                end

                for c = 1:numel(ip_full)

                    DCM_inst      = struct();
                    if n_cell
                        DCM_inst.M.pE = DCMq{d}.M.pE.int{k}.(params{p}){c};
                        DCM_inst.Ep   = DCMq{d}.Ep.int{k}.(params{p}){c};
                        DCM_inst.Cp   = DCMq{d}.Cp(ip_full{c}, ip_full{c});
                    else
                        DCM_inst.M.pE = DCMq{d}.M.pE.int{k}.(params{p});
                        DCM_inst.Ep   = DCMq{d}.Ep.int{k}.(params{p});
                        DCM_inst.Cp   = DCMq{d}.Cp(ip_full{c}, ip_full{c});
                    end

                    if isstruct(DCMq{d}.M.pC)
                        DCM_inst.M.pC = DCMq{d}.M.pC.int{k}.(params{p}){c};
                    else
                        DCM_inst.M.pC = DCMq{d}.M.pC(ip_full{c}, ip_full{c});
                    end

                    DCM_inst.F              = 1;%DCMq{d}.F;
                    DCM_inst.name           = sprintf('DCM%d_%s_src%d', d, nmtype, k);
                    P_DCM{end+1,1}          = DCM_inst;
                    % region_label{end+1,1} = spm_dcm_peb_nmm_region_name(DCMq{d}, k);
                    model_id(end+1,1)     = d;
                end
            end
        end

        if isempty(P_DCM)
            warning('spm_dcm_peb_nmm: no sources of type %s found; skipping.', nmtype);
            continue
        end


        % ---------------------------------------------------------------
        % Second-level design: intercept + effects-coded region + model
        % ---------------------------------------------------------------
        % u_region = unique(region_label, 'stable');
        % Xregion = spm_dcm_peb_nmm_effects_code(region_label, u_region);
        u_model  = unique(model_id, 'stable');
        Xmodel   = spm_dcm_peb_nmm_effects_code(model_id, u_model);
        X        = [ones(numel(P_DCM),1)];% Xmodel];%, Xregion];
        %
        % u_region_row = reshape(u_region(1:end-1), 1, []);
        % u_model_row  = reshape(u_model(1:end-1), 1, []);
        %
        % if ~isempty(u_region_row)
        %     region_names = strcat({'Region: '}, u_region_row);
        % else
        %     region_names = {};
        % end
        %
        % if ~isempty(u_model_row)
        %     model_names = strcat({'Model: DCM'}, cellfun(@num2str, num2cell(u_model_row), 'UniformOutput', false));
        % else
        %     model_names = {};
        % end

        Xnames      = {'Mean'};%,'Model'}];%, region_names, model_names];
        M_T        = struct();
        M_T.X      = X;
        M_T.Xnames = Xnames;
        M_T.Q      = 'single';

        % ---------------------------------------------------------------
        % Run PEB
        % ---------------------------------------------------------------
        [PEB_T, ~]          = spm_dcm_peb(P_DCM, M_T);
        mu      = PEB_T.Ep(:,1)';
        Np      = numel(PEB_T.Pnames);

        if strcmpi(params{p},'G')
            n = length(mu);
            if length(j_d) < n
                missing   = setdiff(1:n, j_d);
                j_d       = [j_d, missing];
            end
            inv_j       = zeros(1, n);
            inv_j(j_d)  = 1:n;
            G_new       = priors.(nmtype).G .* exp(mu(inv_j));
            Var_G       = diag(PEB_T.Cp(1:Np,1:Np));

            updated_priors.(nmtype).G = G_new;
            idx_forced                = nmm_opts.j(:)';
            n_free                    = nmm_opts.n;
            n_remaining               = max(0, n_free - numel(idx_forced));
            rankable                  = setdiff(1:numel(G_new), idx_forced);
            rank_var                  = Var_G(rankable);
            rank_var(isnan(rank_var)) = -Inf;
            [~, order]                = sort(rank_var, 'descend');
            idx_ranked                = rankable(order(1:min(n_remaining, numel(rankable))));
            Idx_free.(nmtype)         = [idx_forced(:); idx_ranked(:)];

            % stash for the prior/variance plot below
            Var_total.G.(nmtype)      = Var_G;

        elseif strcmpi(params{p},'T')
            T_new = priors.(nmtype).T .* exp(mu);
            updated_priors.(nmtype).T = T_new;

            % stash for the prior/variance plot below
            Var_total.T.(nmtype)      = diag(PEB_T.Cp(1:Np,1:Np));

        end
    end
end

% =========================================================================
% WRITE OUTPUT
% =========================================================================
for i = 1:length(unmm)
    if ~isfield(updated_priors, unmm{i}); continue; end
    DCMp.M.([unmm{i} '_G']) = full(updated_priors.(unmm{i}).G);
    DCMp.M.([unmm{i} '_T']) = full(updated_priors.(unmm{i}).T);
    if isfield(Idx_free, unmm{i})
        DCMp.M.([lower(unmm{i}) 'j']) = Idx_free.(unmm{i});
    end
end

% =========================================================================
% PLOT: PRIORS + VARIANCE
% =========================================================================
figure;
n_rows = length(unmm);
n_cols = length(params) * 2;

for i = 1:length(unmm)
    nmtype = unmm{i};
    if ~isfield(updated_priors, nmtype); continue; end

    nmm_opts.n = emp.([lower(nmtype) '_n']);
    nmm_opts.j = emp.([lower(nmtype) '_j']);

    for p = 1:length(params)
        par     = params{p};
        col_pri = (p - 1) * 2 + 1;
        col_var = (p - 1) * 2 + 2;

        if ~isfield(updated_priors.(nmtype), par); continue; end

        % Prior comparison
        subplot(n_rows, n_cols, (i - 1) * n_cols + col_pri);
        hold on;
        bar(priors.(nmtype).(par),          'FaceAlpha', 0.4);
        bar(updated_priors.(nmtype).(par)', 'FaceAlpha', 0.4);
        hold off;
        legend('Original', 'Updated', 'Location', 'best');
        ylabel(nmtype, 'Interpreter', 'none');
        if i == 1; title(sprintf('Priors: %s', par)); end
        xlabel(sprintf('%s parameter index', par));
        grid on;

        % Variance
        subplot(n_rows, n_cols, (i - 1) * n_cols + col_var);

        if strcmpi(par, 'G') && isfield(Idx_free, nmtype)

            v_total    = Var_total.G.(nmtype);
            free_idx   = Idx_free.(nmtype);
            idx_forced = nmm_opts.j(:);
            idx_ranked = setdiff(free_idx, idx_forced);

            bar(v_total, 'FaceColor', [0.2 0.6 0.4], 'FaceAlpha', 0.6);
            hold on;

            legend_handles = [];
            legend_labels  = {};

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
            if ~isempty(legend_handles)
                legend(legend_handles, legend_labels, 'Location', 'best');
            end
            if i == 1; title(sprintf('Second-level variance: %s  (* ranked, o forced)', par)); end

        elseif isfield(Var_total, par) && isfield(Var_total.(par), nmtype)
            bar(Var_total.(par).(nmtype), 'FaceColor', [0.3 0.5 0.8], 'FaceAlpha', 0.5);
            if i == 1; title(sprintf('Second-level variance: %s  (all free)', par)); end
        end

        xlabel(sprintf('%s parameter index', par));
        grid on;
    end
end

sgtitle('Prior update and second-level variance by NMM type');


function name = spm_dcm_peb_nmm_region_name(DCM, k)
% Get a consistent region label for source k of DCM, trying DCM.Sname
% then DCM.xY(k).name, with a positional fallback and a warning if
% neither is present.
name = '';
if isfield(DCM,'Sname') && numel(DCM.Sname) >= k && ~isempty(DCM.Sname{k})
    name = DCM.Sname{k};
elseif isfield(DCM,'xY') && numel(DCM.xY) >= k && isfield(DCM.xY(k),'name') && ~isempty(DCM.xY(k).name)
    name = DCM.xY(k).name;
end
if isempty(name)
    name = sprintf('source%d', k);
    warning(['spm_dcm_peb_nmm: no region name (Sname/xY.name) found for source %d; ' ...
        'using positional fallback ''%s''.'], k, name);
end


function E = spm_dcm_peb_nmm_effects_code(labels, levels)
% Sum-to-zero effects coding of a categorical label vector against a set
% of levels, dropping the last level as the (implicit) reference. Returns
% an N x (numel(levels)-1) matrix, empty in the second dimension if there
% is only one level.
N = numel(labels);
L = numel(levels);
E = zeros(N, max(L-1,0));
for n = 1:N
    if iscell(labels)
        li = find(strcmp(labels{n}, levels));
    else
        li = find(labels(n) == levels);
    end
    if isempty(li); continue; end
    if li == L
        E(n,:) = -1;
    else
        E(n,li) = 1;
    end
end