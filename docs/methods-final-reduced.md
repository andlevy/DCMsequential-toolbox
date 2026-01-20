# Methods

## Time-Varying DCM

Standard dynamic causal modelling (DCM) assumes time-invariant parameters throughout observation, constraining flexibility when neural dynamics exhibit non-stationarity. Medrano et al. (2024) proposed a framework using discrete cosine transform (DCT) basis sets to model time-varying modulatory connectivity (B parameters) whilst maintaining fixed baseline connectivity (A parameters). This approach segments timeseries into consecutive epochs, treating each as a separate condition with effective connectivity expressed as weighted DCT basis functions evaluated at epoch midpoints, demonstrating improved model fit for non-stationary data with regular temporal structure.

Sequential DCM (seq-DCM) extends this foundation by expanding parameter variation beyond modulatory connections to encompass all model parameters. Whilst basis-DCM segments timeseries into discrete epochs, seq-DCM models continuous observations where parameters transition smoothly between discrete physiological states via sigmoid functions. These transitions can be fixed to experimental events, estimated from data, or locked to measured physiological markers, affording advantages for paradigms with irregular transition timing, mixed exogenous and endogenous state changes, and multi-scale dynamics.

### Basis DCM

Basis DCM achieves time-varying connectivity through projection onto temporal basis functions. The observation period is partitioned into $K$ consecutive epochs, with effective connectivity of type $c$ at epoch $k$ expressed as:

$$E_{ck} = A_c + \sum_{m=1}^{M} B_m X_m(t_k)$$

where $A_c$ represents baseline connectivity, $B_m$ are modulatory effect matrices, $X_m(t_k)$ are DCT basis functions evaluated at epoch midpoint $t_k$, and $M$ is the basis set order.

Computational efficiency derives from parameter reduction: rather than estimating separate parameters for each of $K$ epochs, the method estimates only $M$ basis weights where $M \ll K$, exploiting trajectory smoothness. The approach accommodates irregular sampling through basis function evaluation at arbitrary time points and integrates seamlessly with existing DCM inversion schemes. However, restricting temporal variation to modulatory connections may limit applicability when time-variation across all parameters, irregular state transitions, or discontinuous connectivity changes are required.

### Sequential DCM

Seq-DCM permits time-varying trajectories across all neural parameters including baseline connectivity (A), intrinsic synaptic efficacies (G), neural mass parameters (T, H), and transmission delays (D), enabling representation of comprehensive physiological state changes. The framework models continuous timeseries where parameters transition smoothly between discrete physiological states via sigmoid functions, naturally accommodating carryover effects between sequential events such as nonlinear hysteresis, refractoriness, and state-dependent responses.

Transition timing flexibly accommodates both exogenous and endogenous causes of state changes. Transitions may be fixed to experimental event markers, estimated to discover latent state sequences, or locked to measured physiological signals such as heart rate or respiratory phase. The architecture readily accommodates irregular intervals between state transitions of arbitrary duration, from milliseconds to minutes, naturally handling paradigms with variable inter-stimulus intervals, blocked designs, or mixed timescales.

These extensions introduce increased model complexity. The expanded parameter space and explicit transition modelling increase computational demands and necessitate careful consideration of model identifiability and parsimony, addressed through structured priors, parameter constraints, and efficient inference schemes.

## Technical Details

### Trajectory Function Definition

Sequential DCM represents time-varying parameters through piecewise smooth trajectories transitioning between discrete physiological states. Let $\theta(t)$ denote a time-varying parameter. The parameter trajectory is constructed as:

$$\theta(t) = \theta_1 + \sum_{i=1}^{N} [\theta_{i+1} - \theta_i] \cdot v_i(t)$$

where $N$ denotes the number of transitions, $\theta_i$ represents the parameter value in physiological state $i$ (estimated from data), and $v_i(t)$ is a sigmoid transition function expressing $\theta(t)$ as an initial state $\theta_1$ plus cumulative increments weighted by transition functions.

Each transition function $v_i(t)$ smoothly interpolates from 0 to 1. Before transition $i$ occurs, $v_i(t) \approx 0$ and the parameter remains at $\theta_i$. After completion, $v_i(t) \approx 1$ and the parameter has incremented by $\theta_{i+1} - \theta_i$. The sigmoid shape ensures $C^\infty$ continuity, eliminating discontinuities in parameter trajectories and derivatives, essential for gradient-based optimisation and physiologically realistic dynamics.

Transition functions are implemented using cumulative Gaussian error functions:

$$v_i(t) = \frac{1}{2} \left[1 + \text{erf}\left(\frac{t - \tau_i}{\sqrt{2} \sigma_i}\right)\right]$$

where $\tau_i$ specifies the temporal centre of transition $i$ (milliseconds), $\sigma_i$ determines the transition width (standard deviation), and $\text{erf}(\cdot)$ is the Gaussian error function. The centre parameter $\tau_i$ controls when the transition occurs, whilst width parameter $\sigma_i$ controls transition speed.

Transition width is parameterised in log-space to facilitate unconstrained estimation:

$$\sigma_i = \frac{d_i}{2.563103} \exp(K_i)$$

where $d_i$ represents a nominal transition duration, $K_i$ is a log-scale parameter estimated via variational Bayes, and the normalisation constant ensures approximately 98% of the sigmoid transition occurs within duration $d_i$ ($\pm 2\sigma$ coverage). This log-space parameterisation permits $K_i$ to vary freely during optimisation whilst ensuring $\sigma_i > 0$.

Transition timing $\tau_i$ may be fixed to experimental event markers, estimated freely to discover latent state sequences, or locked to measured physiological signals, synchronising parameter dynamics to measured cardiac phase, respiratory cycle, or other biomarkers.

#### Log-Space Parameterisation for Positive Parameters

Neural parameters such as synaptic efficacies, membrane time constants, and connection strengths are inherently positive. To enforce positivity whilst permitting unconstrained optimisation, seq-DCM implements trajectory functions in log-space:

$$\theta(t) = \exp\left[\log(\theta_1) + \sum_{i=1}^{N} [\log(\theta_{i+1}) - \log(\theta_i)] \cdot v_i(t)\right]$$

This formulation guarantees $\theta(t) > 0$ for all $t$ regardless of estimated state values. The log-space representation improves numerical stability by avoiding parameter bounds and naturally accommodates multiplicative parameter changes, often reflecting physiological mechanisms more faithfully than additive changes.

#### Trajectory Covariance and Identifiability

Variational Bayesian inversion estimates posterior parameter covariance, quantifying uncertainty in parameter estimates. For seq-DCM, this uncertainty propagates through the trajectory function to characterise confidence in time-varying parameter estimates across the observation period.

The state-specific parameters $\{\theta_i\}_{i=1}^{N+1}$ are normally distributed random variables. The posterior mean evolves as:

$$\mathbb{E}[\theta(t)] = \mathbb{E}[\theta_1] + \sum_{i=1}^{N} [\mathbb{E}[\theta_{i+1}] - \mathbb{E}[\theta_i]] \cdot v_i(t)$$

The trajectory covariance between time points $t$ and $s$ captures how estimation uncertainty covaries:

$$\text{Cov}[\theta(t), \theta(s)] = \text{Cov}[\theta_1, \theta_1] + \sum_{i=1}^{N} \sum_{j=1}^{N} \text{Cov}[\theta_{i+1} - \theta_i, \theta_{j+1} - \theta_j] \cdot v_i(t) v_j(s)$$

The trajectory variance at time $t$ (obtained by setting $s = t$) provides time-resolved confidence bounds indicating when parameter estimates are well-constrained by data versus uncertain. High trajectory variance at specific times indicates limited data constraint on parameter values, suggesting potential identifiability issues. The temporal structure of covariance reveals whether uncertainty is localised to specific transitions or distributed across the trajectory, informing decisions about model refinement such as fixing poorly identified transition timings.

### Model Implementation

#### Generic DCM Framework and Extension

Sequential DCM builds upon the generic DCM framework (van Wijk et al., 2018), which enables coupling of heterogeneous neural mass models within a unified architecture through standardised interfaces for connecting arbitrary models via extrinsic connections. This modularity maintains mathematical consistency by ensuring extrinsic coupling operates through well-defined afferent and efferent connections between specified hidden states.

By building on the generic framework, seq-DCM accommodates parameter trajectories for any neural mass model available in SPM without requiring model-specific modifications. The generic framework provides a validated foundation, with seq-DCM introducing a single architectural extension: time-varying parameters.

Network dynamics computation proceeds in two stages. First, intrinsic flows for each neural mass model are evaluated through model-specific functions. Second, extrinsic contributions are added by summing weighted influences from connected nodes, with each connection contributing proportionally to source population activity and current connection strength $\theta(t)$. This separation preserves modularity whilst accommodating time-varying connectivity. Conduction delays are incorporated through a Taylor series approximation that absorbs the delay operator into the system Jacobian.

#### Trajectory Functions as Augmented Inputs

Sequential DCM incorporates parameter trajectories into numerical integration by absorbing transition functions into an augmented input matrix, ensuring endogenous parameter changes trigger appropriate Jacobian updates during integration whilst maintaining numerical accuracy.

The standard DCM input structure specifies exogenous driving inputs $u(t)$. In seq-DCM, the input matrix is augmented to include transition functions:

$$U_{\text{aug}}(t) = [u(t), v_1(t), v_2(t), \ldots, v_N(t)]$$

where $u(t)$ represents conventional experimental inputs and $\{v_i(t)\}$ are sigmoid transition functions.

The integration scheme `spm_int_U` re-evaluates the system Jacobian $J = \partial f/\partial x$ whenever inputs change substantially, detected by $\|u(t) - u(t-\Delta t)\|^2 > \epsilon$. By including transition functions $v_i(t)$ as input channels, changes in parameter trajectories are automatically detected as input changes. When a transition function evolves significantly, the Jacobian is recomputed to reflect the current system linearisation, maintaining integration accuracy during rapid parameter transitions whilst avoiding unnecessary computation during stability periods.

The trajectory-dependent parameters $\theta(t)$ are computed at each integration step using current values of $\{v_i(t)\\}$ extracted from $U_{\text{aug}}(t)$, ensuring parameter trajectories and exogenous inputs are evaluated synchronously at identical time points.

#### Extended Input Function Repertoire

Standard DCM employs Gaussian bump functions to model transient stimulus-driven inputs. Whilst appropriate for event-related paradigms with discrete stimuli, this representation proves restrictive for sequential paradigms exhibiting slower, asymmetric dynamics such as sustained anticipatory activity.

Sequential DCM introduces `spm_scp_u`, a skew-normal bump function modelling slow cortical potentials with asymmetric temporal structure:

$$u(t) = \frac{2}{\sigma} \phi\left(\frac{t - \mu}{\sigma}\right) \Phi\left(\alpha \frac{t - \mu}{\sigma}\right)$$

where $\mu$ represents input peak time, $\sigma$ controls temporal width, and $\alpha$ determines skewness. Negative $\alpha$ produces extended left tails (slow onset, rapid offset), positive $\alpha$ yields extended right tails (rapid onset, slow offset), and $\alpha = 0$ recovers the symmetric Gaussian case.

This parameterisation enables flexible representation of stimulus dynamics across multiple timescales. Short-duration symmetric inputs ($\alpha \approx 0$, $\sigma$ small) model phasic sensory stimuli, whilst longer asymmetric inputs ($|\alpha| > 1$, $\sigma$ large) capture sustained processes. The width parameter $\sigma$ is estimated via log-scaling factors permitting data-driven inference of input duration. By matching input function structure to temporal characteristics of driving perturbations, the model reduces parameter trade-offs between input shape and intrinsic dynamics.

#### Parallelized Gradient Computation

Sequential DCM substantially increases parameter space dimensionality, as each time-varying parameter requires estimation of multiple state-specific values $\{\theta_i\}$ plus transition parameters $\{\tau_i, K_i\}$. This expanded parameterisation necessitates computation of high-dimensional Jacobians during variational inversion, where each iteration evaluates partial derivatives $\partial f/\partial \theta$ across all free parameters.

The implementation addresses this through parallelised numerical differentiation via `spm_diff_par`, distributing Jacobian column calculations across multiple processor cores. The numerical derivative for parameter $j$ is computed via finite differences:

$$\frac{\partial f}{\partial \theta_j} \approx \frac{f(\theta + \delta_j \epsilon) - f(\theta)}{\epsilon}$$

where $\delta_j$ is the unit vector in the $j$-th parameter direction and $\epsilon = \exp(-8)$ is the step size. Each column is computed independently, making gradient computation embarrassingly parallel with minimal inter-process communication.

Parallelisation achieves near-linear speedup proportional to available processor cores. On systems with $N$ cores, computational time reduces by a factor approaching $N$, with parallel efficiency $E = S/N$ typically exceeding 0.9 where $S$ denotes observed speedup. This efficiency derives from the independence of finite difference evaluations and negligible overhead of distributing parameter perturbations.

#### Model Complexity and Scaling

The enhanced flexibility of sequential DCM relative to basis-DCM comes at the cost of increased model complexity, quantified by free parameters requiring estimation during inversion. Whilst basis-DCM constrains time-varying parameters to modulatory connectivity through a fixed number of basis function weights, seq-DCM permits arbitrary parameters to vary across multiple discrete states, with each state and transition contributing additional free parameters.

To illustrate scaling properties, consider a canonical two-region model with bidirectional extrinsic connections, yielding four directed connections. This minimal network architecture provides a concrete basis for comparing parameter counts between modelling frameworks as the number of temporal states increases, revealing fundamental trade-offs between temporal resolution and model parsimony that inform practical design choices.
