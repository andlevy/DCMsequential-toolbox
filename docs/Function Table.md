# DCM for Sequential Responses Toolbox: Function Table

| Function | Description |
|---|---|
| `spm_cell` | Splits a matrix into cells like `M.x` |
| `spm_dcm_empirical_priors` | Empirical priors from inverted DCMs |
| `spm_dcm_peb_nmm` | PEB empirical priors (work in progress) |
| `spm_dcm_posterior_trajectory` | Posterior trajectory moments |
| `spm_dcm_priorcorr` | Sets prior correlations |
| `spm_dcm_priorvar` | Sets prior variances |
| `spm_dcm_review_network_graph` | Network graph plot per event |
| `spm_dcm_seq` | Model specification and inversion |
| `spm_dcm_seq_channelmodes` | Channel modes, including HR channels |
| `spm_dcm_seq_collapse` | Collapses event-specific A, B, C into standard DCM A, B, C |
| `spm_dcm_seq_data` | Data preparation |
| `spm_dcm_seq_dipfit` | Forward model specification |
| `spm_dcm_seq_neural_priors` | Neural priors |
| `spm_dcm_seq_sensitivity` | Source sensitivity to parameter contrasts |
| `spm_dcm_seq_x_neural` | Initial states and equations of motion |
| `spm_dcm_synaptic_trajectory` | Synaptic parameter trajectories |
| `spm_dcm_trajectory` | Parameter trajectory across events |
| `spm_diff_par` | Parallel `spm_diff` |
| `spm_eeg_erp_tukey` | ERP preprocessing with a Tukey window |
| `spm_fx_seq` | Equations of motion with event-varying parameters |
| `spm_gen_seq` | Generative model |
| `spm_int_S` | Integrator, updates at input events |
| `spm_int_U` (in `docs/`) | Reference copy of SPM's `spm_int_U` |
| `spm_L_com_priors` | COM observation priors |
| `spm_L_erp_priors` | ECD/IMG/LFP observation priors |
| `spm_L_seq_priors` | Observation priors, multimodal |
| `spm_lx_com` | COM observation model |
| `spm_lx_seq` | Observation model, multimodal |
| `spm_nlsi_N_par` | Parallel `spm_nlsi_N` |
| `spm_null_u` | Null input |
| `spm_scp_u` | Skew-normal input |
| `spm_seq_Q` | Trial-specific parameters |
| `spm_seq_u` | Inputs with sigmoid event transitions |
| `spm_tukey` | Tukey window |
