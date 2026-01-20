# Methods: Modeling Synaptic Transients with Trajectory Functions

## Overview

Sequential Dynamic Causal Modeling (DCM) extends conventional DCM by allowing synaptic parameters to evolve smoothly over time, thereby capturing plasticity and non-stationary dynamics in neural and physiological systems. The fundamental mechanism enabling this temporal evolution is the **trajectory function** (`spm_dcm_trajectory`), which models smooth transitions between discrete synaptic states using sigmoid-weighted interpolation.

This document describes the mathematical formulation and computational implementation of parameter trajectories for modeling time-varying synaptic efficacy, intrinsic connectivity, and other neurophysiological quantities.

---

## Conceptual Framework

### The Problem: Time-Varying Parameters

In conventional DCM, synaptic parameters (e.g., connection strengths, time constants, gain parameters) are assumed constant throughout the observation period. This assumption is violated when:

1. **Plasticity occurs** - synaptic weights change due to learning, adaptation, or state transitions
2. **Experimental conditions change** - tasks involve sequential phases with different connectivity patterns
3. **Physiological states transition** - e.g., shifts from rest to activity, or between cardiac phases

Rather than estimating separate models for each discrete epoch (which is statistically inefficient and introduces discontinuities), sequential DCM estimates a **single generative model** with parameters that evolve continuously according to smooth trajectory functions.

### The Solution: Sigmoid-Weighted Parameter Interpolation

The core idea is to represent time-varying parameters as **weighted sums of discrete parameter states**, where the weights are smooth sigmoid functions that transition from 0 to 1 at specified times. This allows parameters to "flow" from one value to another in a differentiable, physiologically plausible manner.

---

## Mathematical Formulation

### Trajectory Function Definition

Let **θ(t)** denote a time-varying parameter (e.g., a synaptic connection strength). In sequential DCM, θ(t) is modeled as:

```
θ(t) = θ₁ + Σᵢ₌₁ᴷ [θᵢ₊₁ - θᵢ] · vᵢ(t)
```

where:
- **K** = number of transitions
- **θᵢ** = parameter value in state *i* (estimated from data)
- **vᵢ(t)** = sigmoid transition function for transition *i*
- **t** = time

Each sigmoid function vᵢ(t) transitions from 0 to 1, effectively "blending" between consecutive parameter states.

### Sigmoid Transition Functions

The transition functions are cumulative error functions (sigmoid-shaped):

```
vᵢ(t) = 0.5 · [1 + erf((t - τᵢ) / (√2 · σᵢ))]
```

where:
- **τᵢ** = temporal location (center) of transition *i*
- **σᵢ** = temporal scale (width) of transition *i*
- **erf()** = Gaussian error function

This formulation ensures:
- **Smooth transitions**: No discontinuities in parameters or their derivatives
- **Bounded range**: vᵢ(t) ∈ [0, 1]
- **Interpretable parameters**: τᵢ controls *when* transitions occur; σᵢ controls *how fast*

The transition width σᵢ is parameterized as:

```
σᵢ = (dᵢ / 2.563103) · exp(Kᵢ)
```

where:
- **dᵢ** = nominal duration of transition (milliseconds)
- **Kᵢ** = log-scale parameter (estimated from data)
- **2.563103** = scaling constant (≈ 2σ for 98% of sigmoid transition)

### Log-Scale Parameterization

For parameters that must remain positive (e.g., connection strengths, time constants), the trajectory operates in log-space:

```
θ(t) = exp[log(θ₁) + Σᵢ [log(θᵢ₊₁) - log(θᵢ)] · vᵢ(t)]
     = exp[Σᵢ log(θᵢ₊₁) · vᵢ(t) + log(θ₁) · (1 - Σᵢ vᵢ(t))]
```

This ensures θ(t) > 0 for all t, and the estimation operates in an unconstrained (log) space.

---

## Computational Implementation

### Function Signature

```matlab
function Q = spm_dcm_trajectory(P, v, logscale)
```

**Inputs:**
- `P` - Cell array of parameter states {θ₁, θ₂, ..., θₖ₊₁}
- `v` - Matrix of sigmoid values [v₁(t), v₂(t), ..., vₖ(t)] evaluated at time t
- `logscale` - Boolean flag:
  - `true` (default): parameters in log-space (exponentiated)
  - `false`: parameters in linear space

**Output:**
- `Q` - Interpolated parameter value at time t

### Algorithm

The function implements the trajectory equation efficiently:

```matlab
function Q = spm_dcm_trajectory(P, v, logscale)

    % Default: log-scale parameterization
    if nargin < 3
        logscale = true;
    end

    % Initialize with first state
    if logscale
        Q = exp(P{1});  % θ₁ in natural space
    else
        Q = P{1};
    end

    % Accumulate weighted transitions
    for k = 1:size(v, 2)

        if logscale
            q1 = exp(P{k+1});  % θₖ₊₁
            q2 = exp(P{k});    % θₖ
        else
            q1 = P{k+1};
            q2 = P{k};
        end

        % Add weighted difference: (θₖ₊₁ - θₖ) · vₖ(t)
        Q = Q + (q1 - q2) * v(:, k);

    end
end
```

**Key Features:**
1. **Vectorized**: Supports vector/matrix parameters (e.g., connectivity matrices)
2. **Flexible scaling**: Handles both log-scale (positive parameters) and linear-scale parameters
3. **Efficient**: Single-pass accumulation of weighted differences

### Computational Complexity

- **Time**: O(K) where K = number of transitions
- **Space**: O(dim(θ)) where dim(θ) = dimensionality of parameter

---

## Application to DCM Parameters

### Intrinsic Connectivity (A matrices)

Time-varying intrinsic connectivity between neural populations:

```matlab
A(t) = spm_dcm_trajectory(Q.A, v(t,:), true)
```

where `Q.A = {A₁, A₂, ..., Aₖ₊₁}` are connectivity matrices for each state.

**Example**: Transition from resting-state to task-active connectivity.

### Modulatory Connectivity (B matrices)

Time-varying modulation of connections by experimental inputs:

```matlab
B(t) = spm_dcm_trajectory(Q.B, v(t,:), false)
```

Note: B parameters are typically *not* log-scaled (can be positive or negative).

**Example**: Changing attentional modulation across task phases.

### Neural Mass Parameters

Time-varying intrinsic parameters of neural populations:

```matlab
G(t) = spm_dcm_trajectory(Q.int.G, v(t,:), true)  % Gain parameters
T(t) = spm_dcm_trajectory(Q.int.T, v(t,:), true)  % Time constants
H(t) = spm_dcm_trajectory(Q.int.H, v(t,:), true)  % Synaptic efficacies
```

**Example**: Cardiac cycle modulation of neural excitability.

### Conduction Delays

Time-varying axonal/conduction delays:

```matlab
D(t) = log(spm_dcm_trajectory(P.D, v(t,:), true))
```

Note: Result is log-transformed for use in delay operators.

---

## Integration with Sequential DCM

### Workflow Overview

1. **Prior specification**: Define discrete parameter states {θ₁, ..., θₖ₊₁} with Gaussian priors
2. **Sigmoid generation**: Compute sigmoid functions v(t) based on transition times/widths (see `spm_seq_u.m`)
3. **Forward model**: Evaluate trajectory at each integration time step within `spm_fx_seq.m`
4. **Variational inversion**: Estimate posterior distributions over discrete states using `spm_nlsi_N_par.m`
5. **Post-hoc reconstruction**: Compute full time-varying trajectories using `spm_dcm_synaptic_trajectory.m`

### Example: Two-State Transition Model

Consider a simple model with one transition at t = 500ms:

```matlab
% Define discrete states
P.A = {A₁, A₂};  % Pre- and post-transition connectivity

% Define sigmoid transition
τ = 500;         % Transition center (ms)
σ = 50;          % Transition width (ms)
t = 0:1:1000;    % Time vector (ms)
v = 0.5 * (1 + erf((t - τ) / (sqrt(2) * σ)));

% Compute trajectory
for i = 1:length(t)
    A(i,:,:) = spm_dcm_trajectory(P.A, v(i), true);
end
```

**Result**: Smooth connectivity matrix A(t) that transitions from A₁ to A₂ around 500ms.

---

## Statistical Estimation

### Parameter Space

The trajectory function separates parameters into:

1. **Structural parameters** (estimated): {θ₁, ..., θₖ₊₁}
   - Prior: Multivariate Gaussian with mean θ₀ and covariance Σθ
   - Estimated via variational Bayes (EM algorithm)

2. **Temporal parameters** (fixed or estimated):
   - Transition times {τ₁, ..., τₖ}: Often fixed based on experimental design
   - Transition widths {K₁, ..., Kₖ}: Can be estimated or fixed

### Priors on Transition Width

The log-scale width parameter K typically has a weakly informative prior:

```
K ~ N(0, 1)  ⟹  σ ∈ [0.14d, 7.4d] (95% interval)
```

This allows transitions to range from very sharp to very gradual, while regularizing extreme values.

### Advantages for Inference

1. **Parameter efficiency**: Estimate K+1 states instead of continuous function
2. **Regularization**: Gaussian priors on discrete states provide smoothing
3. **Interpretability**: Discrete states correspond to meaningful physiological epochs
4. **Model comparison**: Different numbers of states can be compared using Bayesian model evidence

---

## Physiological Interpretation

### Temporal Resolution of Plasticity

The sigmoid width σ determines the temporal scale of synaptic change:

- **σ ~ 10-50 ms**: Fast plasticity (e.g., short-term depression/facilitation)
- **σ ~ 100-500 ms**: Intermediate dynamics (e.g., attention shifts, state transitions)
- **σ ~ 1-5 s**: Slow adaptation (e.g., learning, habituation)

### Biological Plausibility

The sigmoid functional form approximates several biological processes:

1. **Hebbian plasticity**: Gradual strengthening/weakening of synapses
2. **Neuromodulation**: Slow accumulation/clearance of modulatory substances
3. **Network reconfiguration**: Distributed synaptic changes across populations
4. **Cardiac-phase transitions**: Smooth changes across cardiac cycle phases

### Contrast with Step Functions

Unlike discrete switching models, trajectories:
- Avoid unphysiological discontinuities
- Enable gradient-based optimization
- Better fit gradual empirical transitions
- Reduce model comparison artifacts

---

## Practical Considerations

### Number of Transitions

**Model selection**: Use Bayesian model evidence (free energy) to compare models with different K.

**Guidelines**:
- Start with K based on experimental design (e.g., number of task phases)
- Test K ± 1 for robustness
- Penalize excessive transitions (prefer parsimony)

### Transition Timing

**Fixed vs. estimated**:
- **Fixed**: When timing is known (e.g., stimulus onsets) → more power for structural parameters
- **Estimated**: When timing is uncertain → requires informative priors, reduces power

### Numerical Stability

**Log-scale parameters**:
- Always use log-scale for positive parameters (prevents negative values)
- Check for overflow: exp(P) > 10¹⁵ may indicate numerical issues

**Sigmoid evaluation**:
- Error function `erf()` is numerically stable for |z| < 6
- For extreme values, sigmoid saturates to {0, 1} (no numerical issues)

---

## Extensions and Variations

### Non-Monotonic Trajectories

Current formulation assumes monotonic transitions (always 0 → 1). For non-monotonic changes:

- Use multiple transitions with opposing directions
- Example: θ₁ → θ₂ → θ₁ requires K=2 transitions

### Input-Dependent Trajectories

Transitions can be triggered by experimental inputs:

```
vᵢ(t) = 0.5 · [1 + erf((t - τᵢ(u)) / (√2 · σᵢ))]
```

where τᵢ(u) depends on input u(t).

### Hierarchical Trajectories

Different parameters can have different transition timings:

```matlab
A(t) = spm_dcm_trajectory(Q.A, vₐ(t,:), true)   % Connectivity transitions
G(t) = spm_dcm_trajectory(Q.G, vᵍ(t,:), true)   % Gain transitions
```

where vₐ ≠ vᵍ (different sigmoid sets).

---

## Summary

The `spm_dcm_trajectory` function provides a principled, efficient mechanism for modeling smooth parameter evolution in dynamic systems. Key advantages include:

1. **Biological realism**: Continuous, differentiable parameter changes
2. **Statistical efficiency**: Estimate discrete states rather than continuous functions
3. **Computational simplicity**: Single-pass weighted interpolation
4. **Flexibility**: Applicable to any DCM parameter (connectivity, gains, delays, etc.)
5. **Identifiability**: Regularized by Gaussian priors on discrete states

This framework enables sequential DCM to capture non-stationary dynamics, plasticity, and state transitions in neural and physiological systems, extending DCM's applicability to cardiac neurophysiology, learning paradigms, and other time-varying processes.

---

## References

**Implementation Files**:
- [`spm_dcm_trajectory.m`](../spm_dcm_trajectory.m) - Core trajectory interpolation function
- [`spm_dcm_synaptic_trajectory.m`](../spm_dcm_synaptic_trajectory.m) - Full trajectory reconstruction
- [`spm_seq_u.m`](../spm_seq_u.m) - Sigmoid generation and input specification
- [`spm_fx_seq.m`](../spm_fx_seq.m) - Sequential neural dynamics (uses trajectories)

**Related Theory**:
- Friston, K. J., et al. (2003). "Dynamic causal modelling." *NeuroImage*, 19(4), 1273-1302.
- Friston, K. J., et al. (2016). "Variational Laplace and AIC." *NeuroImage*, 151, 293-305.
