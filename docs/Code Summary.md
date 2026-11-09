# Code Summary: DCM for Sequential Responses Toolbox

## Main Sequential DCM functions

## spm_dcm_seq

Main estimation function that performs Variational Laplace inversion of DCM models for sequential responses using expectation-maximization

## spm_fx_seq

State equations for sequential designs that combine multiple neural mass models with time-varying connectivity and intrinsic parameters. Mention this is an extension of generic DCM spm_fx_gen

## spm_seq_u

Generates augemnted inputs timeseries u including standard driving bump functions with sigmoid transitions %v_i(t)$.

### spm_dcm_trajectory

Computes parameter trajectories by interpolating between event-specific parameter values using a weighted combination piecewise functions.

### spm_seq_Q

Helper routine that computes trial or condition-specific parameters for DCM by deploying trial-specific effects according to a design vector, adjusting connectivity matrices A, B, C and intrinsic parameters.

## dcm_seq_collapse

Collapses time-varying connectivity matrices (A, B, C) from sequential DCM specifications into standard DCM format by accumulating connections across events to enable use of general DCM plotting routines

## spm_dcm_seq_data

Prepares data structures for forward modeling in EEG, MEG and LFP including time windowing, down-sampling, detrending, and windowing operations.

## spm_gen_seq

Generates predictions of trial-specific source activity by integrating the neural mass model for each trial with condition-specific parameters.

## spm_L_priors_seq

Prepares prior moments for lead-field parameters of ERP models including spatial locations, orientations, and contributing states.

## spm_dcm_neural_priors_seq

Prepares priors on parameters of neural mass models with support for sequential designs where parameters can vary across events.

## spm_dcm_x_neural_seq

Returns the initial states and equations of motion for various neural mass models (ERP, CMC, LFP, NMM, etc.).

## Helper functions

## spm_tukey

Returns an n-point Tukey window with specified taper proportion for signal windowing applications. Rather than using the current handing window, two key windows are preferred for longer durations, as handing windows tend to penalize, tends to suppress things happening towards the boundary a bit too severely for these types of experimental designs.

## spm_eeg_erp_tukey

Applies corrections to ERPs or single trials including detrending and Tukey windowing to make sensor-level analysis consistent with DCM-ERP.

## spm_dcm_priorcorr

Induces pairwise correlations between prior parameters by constructing correlation matrices and ensuring positive semi-definite covariance structures.

## spm_dcm_priorvar

Adjusts prior variances for selected parameters by applying scaling factors to specified fields in the prior covariance structures.

## spm_nlsi_N_par

Bayesian inversion of a linear-nonlinear model using Fisher scoring and Laplace approximation with parallel numerical differentiation.

## spm_dcm_seq_review_network_graph

Visualization function that draws network graphs showing connectivity structure with nodes, edges, and modulations for sequential DCM models.

## spm_scp_u

Generates skew-normal input functions for EEG models providing stimulus-related subcortical inputs with asymmetric temporal profiles.

## spm_null_u

Initializes a sparse matrix for null input functions with dimensions based on the number of onsets specified in the model.

## spm_diff_par

Matrix high-order numerical differentiation using parallel computation for efficient gradient calculations with finite differences.

## spm_dcm_seq_synaptic_trajectory

Computes posterior expectation and covariance synaptic parameter trajectories over time by evaluating connectivipty matrices at each time point during sigmoid transitions.
