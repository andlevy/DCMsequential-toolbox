# Methods

## Time-Varying DCM

Standard dynamic causal modelling (DCM) assumes time-invariant parameters throughout the observation period, constraining model flexibility when neural dynamics exhibit non-stationarity. Recent extensions have introduced temporal variation into DCM parameters through basis function expansions. Medrano et al. (2024) proposed a framework using discrete cosine transform (DCT) basis sets to model time-varying modulatory connectivity (B parameters), whilst maintaining fixed baseline connectivity (A parameters). This approach segments timeseries into consecutive epochs, treating each as a separate condition with effective connectivity expressed as a weighted sum of DCT basis functions evaluated at epoch midpoints. The method demonstrated improved model fit for non-stationary data, particularly in paradigms with regular temporal structure and experimentally driven modulation.

The sequential DCM (seq-DCM) framework introduced here builds upon this foundation by extending the scope and flexibility of time-varying parameter modelling. Like basis-DCM, seq-DCM represents slowly evolving parameter trajectories, but expands parameter variation beyond modulatory connections to encompass all model parameters. Another key distinction is that whilst basis-DCM segments timeseries into discrete epochs and models connectivity within each using basis function weights, seq-DCM models the timeseries as a single continuous observation in which parameters transition smoothly between discrete physiological states via sigmoid functions. These transitions can be fixed to experimental events, estimated from data, or locked to measured physiological markers. This continuous modelling approach affords particular advantages for paradigms with irregular transition timing, mixed exogenous and endogenous state changes, and multi-scale dynamics, which are explored through direct comparison in subsequent sections.

The following sections describe the mathematical formulation of seq-DCM, its relationship to the generic DCM framework, and the specific implementation features that enable modelling of sequential physiological responses.

### Basis DCM

Basis DCM achieves time-varying connectivity through projection onto temporal basis functions. The observation period is partitioned into $K$ consecutive epochs, with effective connectivity of type $c$ at epoch $k$ expressed as:

$$E_{ck} = A_c + \sum_{m=1}^{M} B_m X_m(t_k)$$

where $A_c$ represents baseline connectivity, $B_m$ are modulatory effect matrices (free parameters), $X_m(t_k)$ are discrete cosine transform (DCT) basis functions evaluated at epoch midpoint $t_k$, and $M$ is the basis set order. The modulatory matrices $B_m$ weight each basis function, capturing the spectral expansion of time-varying connectivity.

The computational efficiency derives from parameter reduction: rather than estimating separate parameters for each of $K$ epochs, the method estimates only $M$ basis weights where $M \ll K$. This exploits trajectory smoothness, as slow temporal variations can be captured with few DCT basis functions whilst preserving model expressivity. The approach accommodates irregular sampling intervals through basis function evaluation at arbitrary time points and integrates seamlessly with existing DCM inversion schemes by treating epochs as conditions within the standard modulatory framework.

Basis DCM has demonstrated utility for paradigms with regular temporal structure and experimentally driven connectivity changes. However, restricting temporal variation to modulatory connections (with baseline connectivity A and neurobiological parameters fixed) may limit applicability when time-variation across all parameters, irregular state transitions, or discontinuous connectivity changes are required. These scenarios motivate the seq-DCM framework described in subsequent sections.

### Sequential DCM

seq-DCM extends time-varying parameter modelling to encompass the full neural parameter space whilst modelling state transitions explicitly. Where basis-DCM restricts temporal variation to modulatory connectivity, seq-DCM permits time-varying trajectories across all neural parameters including baseline connectivity (A), intrinsic synaptic efficacies (G), neural mass parameters (T, H), and transmission delays (D). This expanded parameterisation enables representation of comprehensive physiological state changes rather than isolated connectivity fluctuations.

The framework models continuous timeseries observations in which parameters transition smoothly between discrete physiological states via sigmoid functions. This continuous state-space formulation naturally accommodates carryover effects between sequential events, capturing phenomena such as nonlinear hysteresis, refractoriness, and state-dependent responses that emerge from incomplete relaxation to baseline between transitions. By embedding transition dynamics directly into the generative model, seq-DCM provides mechanistic representation of temporal dependencies across events.

Transition timing in seq-DCM is flexibly specified to accommodate both exogenous and endogenous causes of state changes. Transitions may be fixed to experimental event markers, approximate endogenously induced transitions, locked to measured physiological signals such as heart rate or respiratory phase. This versatility supports paradigms ranging from structured task designs to spontaneous fluctuations, enabling investigation of both stimulus-driven and intrinsic dynamics within a unified framework.

The architecture readily accommodates irregular intervals between state transitions of arbitrary duration, from milliseconds to minutes. Unlike epoch-based approaches requiring regular temporal segmentation, seq-DCM places no constraints on transition spacing or state duration, naturally handling paradigms with variable inter-stimulus intervals, blocked designs, or mixed timescales. This temporal flexibility proves particularly valuable for naturalistic experimental designs and clinical investigations where precise timing control is impractical.

These extensions afford substantial representational advantages but introduce increased model complexity. The expanded parameter space and explicit transition modelling increase computational demands and necessitate careful consideration of model identifiability and parsimony. The mathematical formulation and inversion procedures described in subsequent sections address these challenges through structured priors, parameter constraints, and efficient inference schemes.

## Technical Details

### Trajectory Function Definition

Sequential DCM represents time-varying parameters through piecewise smooth trajectories that transition between discrete physiological states. Let $\theta(t)$ denote a time-varying parameter (e.g., connectivity strength, synaptic gain, or neural mass time constant). The parameter trajectory is constructed as a weighted sum of discrete state values:

$$\theta(t) = \theta_1 + \sum_{i=1}^{N} [\theta_{i+1} - \theta_i] \cdot v_i(t)$$

where $N$ denotes the number of transitions, $\theta_i$ represents the parameter value in physiological state $i$ (estimated from data), and $v_i(t)$ is a sigmoid transition function. This formulation expresses $\theta(t)$ as an initial state $\theta_1$ plus cumulative increments weighted by transition functions.

Each transition function $v_i(t)$ smoothly interpolates from 0 to 1, effectuating a gradual shift from state $i$ to state $i+1$. Before transition $i$ occurs, $v_i(t) \approx 0$ and the parameter remains at $\theta_i$. After completion, $v_i(t) \approx 1$ and the parameter has incremented by $\theta_{i+1} - \theta_i$. The sigmoid shape ensures $C^\infty$ continuity (infinite differentiability), eliminating discontinuities in parameter trajectories and their derivatives, which is essential for gradient-based optimisation and physiologically realistic dynamics.

Transition functions are implemented using cumulative Gaussian error functions:

$$v_i(t) = \frac{1}{2} \left[1 + \text{erf}\left(\frac{t - \tau_i}{\sqrt{2} \sigma_i}\right)\right]$$

where $\tau_i$ specifies the temporal centre of transition $i$ (in milliseconds), $\sigma_i$ determines the transition width (standard deviation), and $\text{erf}(\cdot)$ is the Gaussian error function. The centre parameter $\tau_i$ controls when the transition occurs, whilst the width parameter $\sigma_i$ controls the transition speed; small $\sigma_i$ produces rapid state changes, whilst large $\sigma_i$ yields gradual evolution.

Transition width is parameterised in log-space to facilitate unconstrained estimation:

$$\sigma_i = \frac{d_i}{2.563103} \exp(K_i)$$

where $d_i$ represents a nominal transition duration (typically set to the inter-transition interval or experimentally determined time constant), $K_i$ is a log-scale parameter estimated via variational Bayes, and the normalisation constant 2.563103 ensures that approximately 98% of the sigmoid transition occurs within duration $d_i$ (corresponding to $\pm 2\sigma$ coverage). This log-space parameterisation permits $K_i$ to vary freely during optimisation whilst ensuring $\sigma_i > 0$.

The transition timing $\tau_i$ may be fixed to experimental event markers, estimated freely to discover latent state sequences, or locked to measured physiological signals. Fixed timing enforces transitions at known stimulus onsets or task events, estimated timing enables inference of endogenous state changes from data, and physiologically locked timing synchronises parameter dynamics to measured cardiac phase, respiratory cycle, or other biomarkers.

#### Log-Space Parameterisation for Positive Parameters

Neural parameters such as synaptic efficacies, membrane time constants, and connection strengths are inherently positive. To enforce positivity whilst permitting unconstrained optimisation, seq-DCM implements trajectory functions in log-space for such parameters:

$$\theta(t) = \exp\left[\log(\theta_1) + \sum_{i=1}^{N} [\log(\theta_{i+1}) - \log(\theta_i)] \cdot v_i(t)\right]$$

This formulation guarantees $\theta(t) > 0$ for all $t$ regardless of the estimated state values $\{\theta_i\}$. The log-space representation improves numerical stability during gradient descent by avoiding parameter bounds and naturally accommodates multiplicative parameter changes, which often reflect physiological mechanisms more faithfully than additive changes. For parameters without positivity constraints (e.g., connectivity that may be excitatory or inhibitory), the linear trajectory formulation is retained.

#### Trajectory Covariance and Identifiability

Variational Bayesian inversion estimates not only the posterior mean of model parameters but also their posterior covariance, quantifying uncertainty in parameter estimates. For seq-DCM, this uncertainty propagates through the trajectory function to characterise confidence in time-varying parameter estimates across the entire observation period. This trajectory covariance provides critical information about parameter identifiability and the reliability of inferred state sequences.

The state-specific parameters $\{\theta_i\}_{i=1}^{N+1}$ estimated during inversion are normally distributed random variables. For a given parameter trajectory $\theta(t)$, the posterior mean evolves as:

$$\mathbb{E}[\theta(t)] = \mathbb{E}[\theta_1] + \sum_{i=1}^{N} [\mathbb{E}[\theta_{i+1}] - \mathbb{E}[\theta_i]] \cdot v_i(t)$$

The trajectory covariance between time points $t$ and $s$ captures how estimation uncertainty at these times covaries:

$$\text{Cov}[\theta(t), \theta(s)] = \text{Cov}[\theta_1, \theta_1] + \sum_{i=1}^{N} \sum_{j=1}^{N} \text{Cov}[\theta_{i+1} - \theta_i, \theta_{j+1} - \theta_j] \cdot v_i(t) v_j(s)$$

The covariance terms $\text{Cov}[\theta_{i+1} - \theta_i, \theta_{j+1} - \theta_j]$ derive directly from the posterior parameter covariance matrix estimated during variational inversion. The trajectory variance at time $t$ is obtained by setting $s = t$, providing time-resolved confidence bounds that indicate when parameter estimates are well-constrained by data versus when they remain uncertain.

This covariance structure reveals several important features of parameter identifiability. High trajectory variance at specific times indicates that the data provide limited constraint on parameter values during those periods, suggesting potential issues with model identifiability or insufficient information in the observed signals. Conversely, low variance indicates strong constraint by data. The temporal structure of covariance also reveals whether uncertainty is localised to specific transitions or distributed across the entire trajectory, informing decisions about model refinement such as fixing poorly identified transition timings or simplifying state sequences.

For models with multiple time-varying parameters, the full trajectory covariance is a matrix-valued function capturing correlations between different parameters across time. This enables assessment of trade-offs between parameters and identification of redundant parameterisations. Such diagnostic information proves invaluable for model validation and comparison, as models with lower trajectory variance at critical time points demonstrate superior identifiability despite potentially similar model evidence values.


### Model Implementation

#### Generic DCM Framework and Extension

Sequential DCM builds upon the generic DCM framework (van Wijk et al., 2018), which enables coupling of heterogeneous neural mass models within a unified architecture. Generic DCM establishes a standardised interface for connecting arbitrary models through extrinsic connections, allowing cortical regions modelled with laminar microcircuits to interact with subcortical structures using anatomically specific formulations. This modularity maintains mathematical consistency by ensuring extrinsic coupling operates through well-defined afferent and efferent connections between specified hidden states.

Network dynamics computation proceeds in two stages. First, intrinsic flows for each neural mass model are evaluated through model-specific functions, capturing local population dynamics. Second, extrinsic contributions are added by summing weighted influences from connected nodes, with each connection contributing proportionally to source population activity and current connection strength $\theta(t)$. This separation of intrinsic and extrinsic dynamics preserves modularity whilst accommodating time-varying connectivity. Conduction delays are incorporated through a Taylor series approximation that absorbs the delay operator into the system Jacobian, maintaining numerical stability whilst capturing phase relationships, with each neural mass model assigned characteristic delays reflecting typical conduction times within its anatomical substrate.

This flexibility proves essential for sequential paradigms engaging multiple interacting systems with distinct intrinsic dynamics. By building on the generic framework, seq-DCM accommodates parameter trajectories for any neural mass model available in SPM without requiring model-specific modifications. The generic framework provides a validated foundation, with seq-DCM introducing a single architectural extension: time-varying parameters.


<!-- Claude do not address this comment - Discuss how we can therefore easily accommodate new models that can account for things outside of just neural mass models.  
We can actually start incorporating different physiological systems, which occur over different time scales.  -->


#### Trajectory Functions as Augmented Inputs

Sequential DCM incorporates parameter trajectories into the numerical integration scheme by absorbing transition functions into an augmented input matrix. This design ensures that endogenous parameter changes trigger appropriate updates to the system Jacobian during integration, maintaining numerical accuracy whilst leveraging efficient sparse-input integration schemes.

The standard DCM input structure specifies exogenous driving inputs $u(t)$ that enter the system dynamics through designated input nodes. In seq-DCM, the input matrix is augmented to include transition functions as additional input channels:

$$U_{\text{aug}}(t) = [u(t), v_1(t), v_2(t), \ldots, v_N(t)]$$

where $u(t)$ represents conventional experimental inputs (e.g., stimulus onsets) and $\{v_i(t)\}$ are the sigmoid transition functions defined previously. This augmented input structure is evaluated at each integration time point, providing current values for both exogenous stimuli and endogenous state transitions.

The integration scheme `spm_int_U` re-evaluates the system Jacobian $J = \partial f/\partial x$ whenever inputs change substantially, as detected by the condition $\|u(t) - u(t-\Delta t)\|^2 > \epsilon$ where $\epsilon$ is a numerical tolerance threshold. By including transition functions $v_i(t)$ as input channels, changes in parameter trajectories are automatically detected as input changes. When a transition function evolves significantly between integration steps—such as during the steep portion of a sigmoid—the Jacobian is recomputed to reflect the current system linearisation. This adaptive Jacobian evaluation maintains integration accuracy during rapid parameter transitions whilst avoiding unnecessary computation during periods of parameter stability.

The trajectory-dependent parameters $\theta(t)$ are computed at each integration step by evaluating:

$$\theta(t) = \theta_1 + \sum_{i=1}^{N} [\theta_{i+1} - \theta_i] \cdot v_i(t)$$

using the current values of $\{v_i(t)\}$ extracted from $U_{\text{aug}}(t)$. This ensures that parameter trajectories and exogenous inputs are evaluated synchronously at identical time points, maintaining temporal consistency between stimulus-driven dynamics and endogenous state evolution. The approach accommodates arbitrary combinations of fixed, estimated, and physiologically locked transition timings through appropriate specification of the transition function time series prior to integration.

#### Extended Input Function Repertoire

Standard DCM employs Gaussian bump functions to model transient stimulus-driven inputs, representing brief subcortical perturbations that evoke rapid cortical responses. Whilst appropriate for event-related paradigms with discrete stimuli, this representation proves restrictive for sequential paradigms exhibiting slower, asymmetric dynamics such as sustained anticipatory activity or gradual post-stimulus relaxation.

Sequential DCM extends the input function repertoire to accommodate diverse temporal profiles of exogenous perturbation. The implementation introduces `spm_scp_u`, a skew-normal bump function that models slow cortical potentials with asymmetric temporal structure. The skew-normal distribution combines a Gaussian probability density function $\phi(z)$ with a cumulative distribution function $\Phi(\alpha z)$ to produce asymmetric profiles:

$$u(t) = \frac{2}{\sigma} \phi\left(\frac{t - \mu}{\sigma}\right) \Phi\left(\alpha \frac{t - \mu}{\sigma}\right)$$

where $\mu$ represents the input peak time, $\sigma$ controls the temporal width, and $\alpha$ determines the skewness. Negative $\alpha$ produces distributions with extended left tails (slow onset, rapid offset), positive $\alpha$ yields extended right tails (rapid onset, slow offset), and $\alpha = 0$ recovers the symmetric Gaussian case.

This parameterisation enables flexible representation of stimulus dynamics across multiple timescales. Short-duration symmetric inputs ($\alpha \approx 0$, $\sigma$ small) model phasic sensory stimuli, whilst longer asymmetric inputs ($|\alpha| > 1$, $\sigma$ large) capture sustained processes such as preparatory states or prolonged modulatory influences. The width parameter $\sigma$ is estimated via log-scaling factors that permit data-driven inference of input duration, allowing the model to discover optimal input profiles without strong prior constraints.

The expanded input function repertoire accommodates sequential paradigms where stimulus effects unfold over extended periods or exhibit characteristic rise and decay profiles. By matching input function structure to the temporal characteristics of driving perturbations, the model reduces parameter trade-offs between input shape and intrinsic dynamics, improving identifiability and interpretability of estimated neural parameters.

#### Parallelized Gradient Computation

Sequential DCM substantially increases the dimensionality of parameter space relative to static DCM, as each time-varying parameter requires estimation of multiple state-specific values $\{\theta_i\}$ plus transition parameters $\{\tau_i, K_i\}$. This expanded parameterization necessitates computation of high-dimensional Jacobians during variational inversion, where each iteration evaluates partial derivatives $\partial f/\partial \theta$ across all free parameters. For models with long observation periods or numerous time-varying connections, gradient computation becomes the dominant computational bottleneck.

The implementation addresses this through parallelized numerical differentiation via `spm_diff_par`, which distributes Jacobian column calculations across multiple processor cores. The numerical derivative for parameter $j$ is computed via finite differences:

$$\frac{\partial f}{\partial \theta_j} \approx \frac{f(\theta + \delta_j \epsilon) - f(\theta)}{\epsilon}$$

where $\delta_j$ is the unit vector in the $j$-th parameter direction and $\epsilon = \exp(-8)$ is the step size. Each column of the Jacobian is computed independently by perturbing a single parameter and evaluating the model, making gradient computation embarrassingly parallel with minimal inter-process communication.

Parallelization achieves near-linear speedup proportional to the number of available processor cores. On systems with $N$ cores, computational time reduces by a factor approaching $N$, with parallel efficiency $E = S/N$ typically exceeding 0.9 where $S$ denotes observed speedup. For example, an 8-core system achieves approximately 8-fold acceleration of gradient computation, substantially reducing the time required for model inversion. This efficiency derives from the independence of finite difference evaluations and the negligible overhead of distributing parameter perturbations to parallel workers.

The parallelization strategy integrates seamlessly with existing SPM infrastructure, requiring no modifications to model specification or inversion procedures. Users benefit from automatic acceleration when MATLAB's parallel computing toolbox is available, with the implementation reverting to sequential computation otherwise. This computational optimization proves critical for applications with extensive timeseries or complex network architectures, enabling practical deployment of sequential DCM to datasets that would otherwise require prohibitive computation time.

#### Model Complexity and Scaling

The enhanced flexibility of sequential DCM relative to basis-DCM comes at the cost of increased model complexity, quantified by the number of free parameters requiring estimation during inversion. Whilst basis-DCM constrains time-varying parameters to modulatory connectivity through a fixed number of basis function weights, seq-DCM permits arbitrary parameters to vary across multiple discrete states, with each state and transition contributing additional free parameters. Understanding this complexity scaling proves essential for practical model specification, as excessive parameterisation degrades identifiability and computational efficiency whilst insufficient parameterisation limits representational capacity.

To illustrate the scaling properties, consider a canonical two-region model with bidirectional extrinsic connections, yielding four directed connections (forward and backward projections in each direction). This minimal network architecture provides a concrete basis for comparing parameter counts between modelling frameworks as the number of temporal states increases. The analysis examines how parameter dimensionality grows with the number of transitions $N$ for seq-DCM versus the basis set order $M$ for basis-DCM, revealing fundamental trade-offs between temporal resolution and model parsimony that inform practical design choices for sequential paradigms.