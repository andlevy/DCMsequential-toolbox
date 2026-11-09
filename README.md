# Sequential DCM Toolbox

An SPM toolbox for dynamic causal modelling (DCM) of sequential responses in EEG, MEG and LFP data.

It extends DCM so that model parameters can change across a sequence of events within a trial, rather than staying fixed. It works with the neural mass models already in SPM and follows the standard SPM DCM workflow.

## Requirements

- MATLAB
- SPM
- Parallel Computing Toolbox (optional)

## Installation

Add the toolbox folder to the MATLAB path:

```matlab
addpath('path/to/dcm_seq')
```

## Usage

Set up a `DCM` structure as for `spm_dcm_erp`, then invert it with:

```matlab
DCM = spm_dcm_seq(DCM);
```

For a list of functions, see [docs/Function Table.md](docs/Function%20Table.md).

## Reference

Levy AD, Zeidman P, Friston K (2026). Time-Varying Dynamic Causal Modelling for Sequential Responses: Neural Mechanisms of Slow Cortical Potentials, Preparation, Planning and Beyond. *bioRxiv* 2026.03.24.714008. https://doi.org/10.64898/2026.03.24.714008

## Author

Andrew Levy
