# Continuous Size-Spectrum (`mizer`) and State-Space Food-Web Models (`EwE`, `Rpath`, `EcoState`) {#sec-size-spectrum}

## Continuous Multispecies Size-Spectrum Models (`mizer`)

### Theoretical Lineage and Core Philosophy
Size-spectrum theory originated from the empirical discovery by [@sheldon1972] and [@sheldon1977] that pelagic marine biomass is approximately constant across logarithmic body-mass bins ($B(w) dw \propto w^0 d(\ln w) \implies N(w) \propto w^{-2}$), and from the continuous mass-balance and pelagic ecosystem derivations of [@platt1977], [@silvert_platt1978], [@benoit_rochet2004], [@andersen_beyer2006], and [@cuesta2018].

In a landmark synthesis, [@hartvig2011] and [@scott2014] resolved species identity within the continuous size spectrum to create the **`mizer`** framework (*Multi-species size-spectrum ecological modelling in R*; applied to North Sea reference points by [@blanchard2014] and reviewed comprehensively by [@andersen2016], [@andersen2019], [@andersen2020], and [@pennino2025]).

Unlike age-structured models where somatic growth and natural mortality are external inputs, in `mizer` **both individual somatic growth $g_i(w,t)$ and predation mortality $\mu_{p,i}(w,t)$ emerge dynamically from individual-level foraging and strict conservation of mass/energy**.

---

### Complete Mathematical Formulation of Standard `mizer`

Let $N_i(w, t)$ ($\text{numbers} \cdot \text{g}^{-1} \cdot \text{m}^{-3}$ or total numbers per weight interval) denote the number density of species $i \in \{1, \dots, S\}$ at body mass $w \in [w_0, w_{\max, i}]$ and time $t$, and let $N_R(w, t)$ denote the background plankton/benthic resource spectrum for $w \in [w_{\text{cut}}, w_{\text{max}, R}]$.

#### 1. Governing McKendrick–von Foerster Conservation PDE
The time evolution of $N_i(w, t)$ is governed by the first-order hyperbolic continuity (advection-reaction) equation [@hartvig2011, Eq. 1; @scott2014, Eq. 3.1]:

$$\frac{\partial N_i(w, t)}{\partial t} + \frac{\partial \left(g_i(w, t) N_i(w, t)\right)}{\partial w} = -\mu_i(w, t) N_i(w, t)$$

coupled to the **boundary condition (recruitment flux)** at egg weight $w_0$:

$$g_i(w_0, t) N_i(w_0, t) = R_{\text{dd}, i}(t)$$

#### 2. Volumetric Search Rate, Predation Kernel, and Encounter Rate ($E_i(w,t)$)
A predator of species $i$ and mass $w$ searches a water volume at allometric rate $\gamma_i(w) = \gamma_i w^q$ ($q \approx 0.8$). Its size preference for prey of mass $w_p$ is described by a **predation kernel** $\phi_i(w / w_p)$, most commonly the Andersen–Ursin lognormal kernel [@scott2014, Eq. 3.4]:

$$\phi_i\!\left(\frac{w}{w_p}\right) = \exp\!\left( -\frac{\left(\ln\!\left(\frac{w}{\beta_i w_p}\right)\right)^2}{2\sigma_i^2} \right)$$

(or a power-law / box / double-logistic kernel). Given the species interaction (spatial overlap / vulnerability) matrix $\theta_{ij} \in [0, 1]$, the **total biomass of suitable prey encountered per unit time** $E_i(w, t)$ by predator $(i, w)$ is:

$$E_i(w, t) = \gamma_i(w) \int \phi_i\!\left(\frac{w}{w_p}\right) \left[ \theta_{iR} N_R(w_p, t) + \sum_{j=1}^{S} \theta_{ij} N_j(w_p, t) \right] w_p \, dw_p$$

#### 3. Feeding Level ($f_i(w,t)$) and Holling Type II Satiation
Given physiological maximum intake rate $I_{\max, i}(w) = h_i w^n$ ($n \approx 0.67\text{--}0.75$), the dimensionless **feeding level** $f_i(w, t) \in [0, 1)$ follows a Holling Type II functional response:

$$f_i(w, t) = \frac{E_i(w, t)}{E_i(w, t) + h_i w^n}$$

Thus actual biomass consumption per predator per unit time is $C_i(w, t) = f_i(w, t) h_i w^n$.

#### 4. Bioenergetic Allocation, Somatic Growth ($g_i(w,t)$), and Reproduction ($R_{p,i}(t)$)
Assimilated energy is a fraction $\alpha_i \in (0, 1)$ of ingested food ($\alpha_i \approx 0.6$), from which standard metabolic maintenance $k_{s,i} w^p$ ($p \approx 0.7\text{--}1.0$) is subtracted to give net available energy $E_{\text{avail}, i}(w, t)$:

$$E_{\text{avail}, i}(w, t) = \alpha_i \, f_i(w, t) \, h_i w^n - k_{s,i} w^p$$

A size-dependent maturation ogive $\psi_i(w)$ (transitioning from $0$ to $1$ around maturation weight $w_{\text{mat}, i}$ and scaling with $\left(w / w_{\max, i}\right)^{1-n}$ or $1 - (w/w_{\max,i})^{p-n}$ so growth ceases at $w_{\max, i}$) partitions net positive energy between **somatic growth** $g_i(w, t)$ and **gonadal reproduction** $E_{r,i}(w, t)$:

$$g_i(w, t) = \max\!\Big(0,\; E_{\text{avail}, i}(w, t)\Big) \Big(1 - \psi_i(w)\Big)$$

$$\psi_i(w) = \left[1 + \left(\frac{w}{w_{\text{mat}, i}}\right)^{-U}\right]^{-1} \left(\frac{w}{w_{\max, i}}\right)^{1 - n}$$

The **density-independent physiological egg production rate** $R_{\text{di}, i}(t)$ ($\text{eggs} \cdot \text{yr}^{-1}$) is obtained by integrating reproductive energy over all sizes, multiplying by reproductive efficiency $\epsilon_{\text{rep}, i}$ and female fraction $1/2$, and dividing by egg weight $w_0$:

$$R_{\text{di}, i}(t) = \frac{\epsilon_{\text{rep}, i}}{2 w_0} \int_{w_0}^{w_{\max, i}} N_i(w, t) \max\!\Big(0,\; E_{\text{avail}, i}(w, t)\Big) \psi_i(w) \, dw$$

In standard `mizer`, external density dependence at the larval stage (e.g., limited larval habitat or unmodeled larval mortality) is represented by passing $R_{\text{di}, i}(t)$ through a **Beverton–Holt stock-recruitment function** with maximum recruitment cap $R_{\max, i}$:

$$R_{\text{dd}, i}(t) = \frac{R_{\max, i} \, R_{\text{di}, i}(t)}{R_{\text{di}, i}(t) + R_{\max, i}}$$

#### 5. Total Mortality ($\mu_i(w,t)$): Background, Predation, Starvation, and Fishing
Total instantaneous mortality $\mu_i(w, t)$ is the sum of four mechanistic components:

$$\mu_i(w, t) = \mu_{0, i}(w) + \mu_{p, i}(w, t) + \mu_{s, i}(w, t) + F_i(w, t)$$

1. **Background / Senescence Mortality**: $\mu_{0, i}(w) = \mu_{0, i} w_{\max, i}^{n - 1}$ (allometrically scaled with $w_{\max,i}$).
2. **Endogenous Predation Mortality ($\mu_{p, i}(w_p, t) \equiv M2$)**: By strict mass balance, every gram of prey $(j, w_p)$ eaten by predators $(i, w)$ must appear as mortality on $(j, w_p)$ [@scott2014, Eq. 3.11]:

   $$\mu_{p, j}(w_p, t) = \sum_{i=1}^{S} \theta_{ij} \int_{w_0}^{w_{\max, i}} \phi_i\!\left(\frac{w}{w_p}\right) \Big(1 - f_i(w, t)\Big) \gamma_i(w) \, N_i(w, t) \, dw$$

   *(Notice the identity $\frac{f_i(w) h_i w^n}{E_i(w)} = 1 - f_i(w)$: when a predator is satiated [$f_i(w) \to 1$], its per-unit-prey clearance rate $(1 - f_i)\gamma_i(w)$ drops toward zero!)*
3. **Starvation Mortality ($\mu_{s, i}(w, t)$)**: Whenever available assimilated energy falls below metabolic maintenance ($E_{\text{avail}, i}(w, t) < 0$), individuals suffer starvation mortality proportional to the energy deficit divided by lipid energy stores $\xi_{\text{store}} w$ [@scott2014, Eq. 3.12]:

   $$\mu_{s, i}(w, t) = \begin{cases} \dfrac{-E_{\text{avail}, i}(w, t)}{\xi_{\text{store}} w} & \text{if } E_{\text{avail}, i}(w, t) < 0 \\[8pt] 0 & \text{otherwise} \end{cases}$$

4. **Fishing Mortality ($F_i(w, t)$)**: Summed over fishing gears/fleets $f$ with size selectivity $S_{f,i}(w)$, catchability $Q_{f,i}$, and fishing effort $E_{f}(t)$:

   $$F_i(w, t) = \sum_{f} S_{f, i}(w) \, Q_{f, i} \, E_{f}(t)$$

   generating total fishery yield in weight $Y_{f,i}(t) = \int_{w_0}^{w_{\max,i}} F_{f,i}(w, t) \, w \, N_i(w, t) \, dw$.

#### 6. Resource Spectrum Dynamics ($N_R(w,t)$)
The background plankton/benthos spectrum $N_R(w, t)$ evolves via semi-chemostat (or logistic) renewal at rate $r_R(w) = r_0 w^{n-1}$ toward carrying capacity $\kappa_R(w) = \kappa w^{-\lambda}$ ($\lambda \approx 2 + q - n \approx 2.05$):

$$\frac{\partial N_R(w, t)}{\partial t} = r_R(w) \Big(\kappa_R(w) - N_R(w, t)\Big) - \mu_{p, R}(w, t) N_R(w, t)$$

#### 7. Semi-Implicit Upwind Finite-Difference Discretization in `mizer`
On a logarithmically spaced weight grid $w_j = w_0 \cdot 10^{(j-1)\Delta x}$ with bin widths $\Delta w_j = w_{j+1} - w_j$, `mizer` solves the McKendrick–von Foerster PDE over time step $\Delta t$ by evaluating encounter $E_i$, feeding level $f_i$, growth $g_{i,j}^t$, and mortality $\mu_{i,j}^t$ explicitly at time $t$ (using Fast Fourier Transform [FFT] convolution on the log-weight grid for the integrals $E_i$ and $\mu_{p,i}$ when kernels are shift-invariant in $\ln(w/w_p)$), and then advancing $N_{i,j}^{t+\Delta t}$ via a **first-order upwind, backward-Euler (semi-implicit) scheme** [@scott2014, Appendix A]:

$$\frac{N_{i,j}^{t+\Delta t} - N_{i,j}^t}{\Delta t} + \frac{g_{i,j}^t N_{i,j}^{t+\Delta t} - g_{i,j-1}^t N_{i,j-1}^{t+\Delta t}}{\Delta w_j} = -\mu_{i,j}^t N_{i,j}^{t+\Delta t}$$

Solving for $N_{i,j}^{t+\Delta t}$ gives a lower-bidiagonal linear system solved in $O(n_w)$ forward substitution:

$$N_{i,j}^{t+\Delta t} = \frac{N_{i,j}^t + \frac{\Delta t}{\Delta w_j} g_{i,j-1}^t N_{i,j-1}^{t+\Delta t}}{1 + \Delta t \left(\frac{g_{i,j}^t}{\Delta w_j} + \mu_{i,j}^t\right)}$$

---

## The Jump-Growth Fokker–Planck PDE (Datta, Delius & Law 2010, 2011)

### Why the First-Order McKendrick–von Foerster PDE Fails for Dynamic Cohort Fitting
In the standard McKendrick–von Foerster PDE $\partial_t N + \partial_w (g N) = -\mu N$, somatic growth is treated as a purely deterministic continuous drift $dw/dt = g(w, t)$. Consequently:
1. **No Intrinsic Cohort Broadening**: All individuals born at the same instant $t_0$ at egg weight $w_0$ travel along the exact same characteristic trajectory $w(t)$ with zero variance in size-at-age! If recruitment is pulsed seasonally (or fluctuates interannually as in real stock assessments), a narrow delta-function pulse at $w_0$ remains an unrealistically sharp spike (broadened only by **numerical discretization error** $\propto \Delta w_j$ of the upwind scheme!).
2. **Grid-Dependent Numerical Diffusion**: In `mizer`'s first-order upwind scheme, changing the number of weight bins $n_w$ or time step $\Delta t$ alters the numerical diffusivity $D_{\text{num}} \approx \frac{1}{2} g \Delta w$, which artificially changes how recruitment pulses smear across length compositions!
3. **Deterministic Instability of the Community Size Spectrum**: Without growth dispersion or an external stock-recruitment clamp, the continuous McKendrick–von Foerster PDE exhibits traveling-wave instabilities (cohort cycles) for realistic predator-prey mass ratios ($\beta \sim 100$; [@datta2011]).

### Derivation of the Jump-Growth Equation (Datta, Delius & Law 2010)
[@datta2010] resolved this fundamental limitation from first principles by observing that fish do not grow by absorbing a continuous fluid: they grow by **discrete stochastic jumps** $\Delta w = \alpha w_p$ each time they capture and assimilate a discrete prey organism of body mass $w_p$.

Starting from the exact Markov **Master Equation** for discrete prey-capture events in continuous time [@datta2010, Eq. 5]:
$$\frac{\partial N_i(w,t)}{\partial t} = \int dw_p \Big[ \text{Transition in from } (w - \alpha_i w_p) - \text{Transition out from } w \Big] - \mu_i(w,t) N_i(w,t)$$
and performing a second-order **Kramers–Moyal Taylor expansion** in the prey mass increment $\alpha_i w_p \ll w$, [@datta2010] derived the **Jump-Growth Fokker–Planck (Advection–Diffusion–Reaction) PDE**:

$$\frac{\partial N_i(w, t)}{\partial t} + \frac{\partial \left(g_i(w, t) N_i(w, t)\right)}{\partial w} = \frac{1}{2} \frac{\partial^2 \left(d_i(w, t) N_i(w, t)\right)}{\partial w^2} - \mu_i(w, t) N_i(w, t)$$

where:
- The **first-order drift (mean somatic growth rate)** $g_i(w, t)$ is the first moment of assimilated prey mass:
  $$g_i(w, t) = \alpha_i \Big(1 - \psi_i(w)\Big) \Big(1 - f_i(w, t)\Big) \gamma_i(w) \int \phi_i\!\left(\frac{w}{w_p}\right) N_c(w_p, t) \, \mathbf{w_p} \, dw_p - k_{s,i} w^p \Big(1 - \psi_i(w)\Big)$$
- The **second-order macro-diffusion (growth variance) coefficient** $d_i(w, t)$ is the **second moment** of assimilated prey mass [@datta2010, Eq. 9]:
  $$d_i(w, t) = \alpha_i^2 \Big(1 - \psi_i(w)\Big)^2 \Big(1 - f_i(w, t)\Big) \gamma_i(w) \int \phi_i\!\left(\frac{w}{w_p}\right) N_c(w_p, t) \, \mathbf{w_p^2} \, dw_p$$

> [!IMPORTANT]
> **Mathematical & Practical Significance of the Jump-Growth Diffusion Term $d_i(w,t)$**:
> 1. **Zero Free Parameters**: In the pure mechanistic derivation, $d_i(w, t)$ requires **no additional free parameters** (or a single scaling multiplier if accounting for behavioral bout-feeding / schooling variation, analogous to Gadget's Beta-Binomial $\beta_s$).
> 2. **Physical Cohort Dispersion**: A recruiting year-class pulse at $w_0$ naturally broadens into a realistic bell-shaped size-at-age distribution as it advects to larger sizes, enabling direct likelihood fitting to time-series length-frequency distributions!
> 3. **Grid-Independent Convergence**: When solved on a sufficiently fine grid (or via Crank–Nicolson / implicit central differences for the $\partial_w^2(d_i N_i)$ term), the physical diffusivity $d_i(w, t)$ dominates numerical grid diffusivity, making model predictions invariant to grid resolution.

---

## Population Regulation Without Ad-Hoc Stock–Recruitment Relationships (Canales, Delius & Law 2020)

In standard fisheries models (including standard `mizer`, `SMS`, `Gadget`, and `LeMaRns`), population abundance is stabilized externally by imposing an ad-hoc **Beverton–Holt or Ricker stock-recruitment relationship (SRR)** at the boundary $w_0$, with one or two free species-specific parameters ($R_{\max, i}$ or steepness $h_i$) that have little direct physiological basis and strongly dictate MSY predictions.

[@canales2020] demonstrated that in a size-spectrum model with linear density-independent reproduction ($R_{\text{dd}, i} = R_{\text{di}, i}$, i.e., $R_{\max, i} \to \infty$), **three endogenous biological mechanisms** naturally regulate population abundance and coexistence without any ad-hoc SRR:
1. **Density-Dependent Somatic Growth via Prey Competition**: When a species' abundance increases, it depletes its prey size range, lowering its feeding level $f_i(w, t)$ and slowing juvenile growth $g_i(w, t)$. Because juveniles spend longer in small vulnerable size classes ($w < w_{\text{mat}, i}$) while experiencing predation mortality $\mu_{p,i}(w, t)$, cumulative survival to maturity $\exp\!\left(-\int_{w_0}^{w_{\text{mat}}} \frac{\mu_i(w)}{g_i(w)} dw\right)$ drops sharply!
2. **Early-Life Predation and Cannibalism**: Extending the predation kernel to planktonic/larval stages and accounting for intra- and inter-specific predation on larvae creates a self-regulating mortality bottleneck.
3. **Reproductive Allometry ($m_i > n_i$)**: When reproductive investment scales steeper with body mass than metabolic intake (matching empirical BOFFFF—Big Old Fat Fecund Female—observations), truncating the adult size structure via fishing directly reduces per-biomass egg output without needing an SRR curve.

---

## Mass-Balance, Foraging Arena, and State-Space Food-Web Models (`Ecopath`, `Ecosim`, `Rpath`, `EcoState`)

While size-spectrum models resolve continuous body mass, the most widely applied whole-ecosystem modeling framework in marine science is **Ecopath with Ecosim (EwE)**—originating with [@polovina1984] and [@christensen1992], extended dynamically and spatially via Foraging Arena Theory in `Ecosim` and `Ecospace` by [@walters1997], [@walters1999], [@christensen2004], [@christensen2005guide], and [@christensen2024ewe], critically evaluated by [@plaganyi2004], and automated via stepwise fitting by [@scott2016]—along with its reproducible R implementation **`Rpath`** [@lucey2020] and its recent TMB state-space stock assessment reformulation **`EcoState`** [@thorson2024].

### `Ecopath`: Static Mass Balance at Reference Equilibrium
`Ecopath` partitions the ecosystem into $i \in \{1, \dots, n\}$ functional groups (which may be single species, multi-species guilds, or juvenile/adult life stanzas). At a reference equilibrium period, production of group $i$ equals the sum of all mortality and export fluxes [@christensen2004, Eq. 1; @lucey2020, Eq. 1]:

$$B_i \left(\frac{P}{B}\right)_i \cdot EE_i = Y_i + E_i + BA_i + \sum_{j=1}^{n} B_j \left(\frac{Q}{B}\right)_j DC_{ij}$$

where:
- $B_i$ is biomass, $(P/B)_i = Z_i$ is the production-to-biomass ratio (total mortality rate under equilibrium), and $(Q/B)_j$ is the consumption-to-biomass ratio of predator $j$,
- $EE_i \in [0, 1]$ is **Ecotrophic Efficiency** (the fraction of production explained within the modeled system, so $M0_i = (P/B)_i (1 - EE_i)$ is residual unexplained mortality),
- $Y_i$ is fishery catch, $E_i$ is net migration, $BA_i$ is biomass accumulation, and
- $DC_{ij}$ is the **diet composition matrix** (fraction by weight of predator $j$'s diet consisting of prey $i$, with $\sum_i DC_{ij} = 1$).
For each group $i$, three of the four core parameters $\{B_i, (P/B)_i, (Q/B)_i, EE_i\}$ are supplied along with $DC_{ij}$ and $Y_i$, and the linear system is solved for the missing parameter (typically $EE_i$ or $B_i$).

---

### `Ecosim` and Walters' Foraging Arena Theory
When `Ecopath` is projected dynamically through time as a system of coupled ODEs (`Ecosim`; [@walters1997]):

$$\frac{dB_i}{dt} = g_i \sum_{j} Q_{ji}(t) - \sum_{k} Q_{ik}(t) + I_i - \left(M0_i + F_i(t) + e_i\right) B_i(t)$$

where $g_i = (P/Q)_i$ is net growth efficiency, a critical question arises: how should the consumption flux $Q_{ij}(t)$ of prey $i$ by predator $j$ depend on $B_i(t)$ and $B_j(t)$?

[@walters1997] showed that if one uses standard mass-action Lotka–Volterra ($Q_{ij} = a_{ij} B_i B_j$) or standard Holling Type II in a complex 30-species food web, **the ecosystem rapidly collapses via Gause competitive exclusion or violent predator-prey limit cycles**.

To resolve this, [@walters1997] formulated **Foraging Arena Theory**: at any given moment, total prey biomass $B_i$ is partitioned with respect to predator $j$ into an **invulnerable (refuge/sheltered/un-overlapping) state** $B_i - V_{ij}$ and a small **vulnerable foraging arena state** $V_{ij}$, with behavioural/spatial exchange rate $v_{ij}$ between them:

$$\frac{dV_{ij}}{dt} = v_{ij}\left(B_i - V_{ij}\right) - v'_{ij} V_{ij} - a_{ij} V_{ij} B_j$$

Assuming the behavioral exchange $V_{ij}$ equilibrates rapidly relative to annual population changes ($dV_{ij}/dt \approx 0$, with $v'_{ij} = v_{ij}$), solving for $V_{ij}$ yields the **Foraging Arena Functional Response** [@walters1997, Eq. 4; @christensen2004, Eq. 5; @thorson2024, Eq. 4]:

$$Q_{ij}(t) = a_{ij} V_{ij}(t) B_j(t) = \frac{a_{ij} \, v_{ij} \, B_i(t) \, B_j(t)}{2 v_{ij} + a_{ij} B_j(t)}$$

Calibrating to the Ecopath baseline $Q_{ij,0} = B_{j,0} (Q/B)_{j,0} DC_{ij,0}$ at reference biomasses $(B_{i,0}, B_{j,0})$ eliminates $a_{ij}$ in favor of a single dimensionless **vulnerability multiplier** $X_{ij} = \frac{v_{ij} B_{i,0}}{Q_{ij,0}} \in (1, \infty)$ (default $X_{ij} = 2$):

$$Q_{ij}(t) = Q_{ij,0} \, \frac{X_{ij}}{X_{ij} - 1 + \frac{B_j(t)}{B_{j,0}}} \, \frac{B_i(t)}{B_{i,0}}$$

- **Low Vulnerability ($X_{ij} \to 1$, "Bottom-Up / Donor Control")**: As $X_{ij} \to 1$, $Q_{ij}(t) \to Q_{ij,0} \frac{B_i(t)}{B_{i,0}}$ is **independent of predator biomass $B_j(t)$**! Doubling predator abundance $B_j$ simply means twice as many predators compete for the same trickle of prey entering the foraging arena, halving per-capita predator intake and preventing top-down prey extermination.
- **High Vulnerability ($X_{ij} \to \infty$, "Top-Down / Lotka–Volterra Control")**: As $X_{ij} \to \infty$, $Q_{ij}(t) \to Q_{ij,0} \frac{B_i(t)}{B_{i,0}} \frac{B_j(t)}{B_{j,0}}$ recovers mass-action top-down predation.

*(For critique of EwE's delay-difference multi-stanza assumptions and vulnerability overfitting when many $X_{ij}$ are tuned to short time series, see [@plaganyi2004] and the stepwise AIC fitting protocol of [@scott2016]).*

---

### `EcoState` (Thorson et al. 2024): State-Space Mass-Balance Assessment in TMB
Historically, EwE was fitted by shooting (deterministic trajectory matching) without process errors or formal state-space marginal likelihoods. In 2024, James Thorson, Kirstin Holsman, Gaichas, Aydin, and Punt introduced **`EcoState`** [@thorson2024], which embeds the Ecopath + Ecosim ODEs directly inside a **Template Model Builder (TMB) state-space stock assessment model**:

1. **Continuous-Time ODE with Annual Process Errors $\epsilon_{i,t}$**:
   Within each year $t \in \{t_0, \dots, T\}$, biomass $B_i(s)$ evolves over $s \in [t, t+1]$ via the Ecosim Foraging Arena ODE augmented with an estimated annual unexplained process-error rate $\epsilon_{i,t}$ [@thorson2024, Eq. 5]:
   $$\frac{dB_i(s)}{ds} = g_i \sum_j Q_{ji}(s) - \sum_k Q_{ik}(s) - \Big(M0_i + F_{i,t} - \epsilon_{i,t}\Big) B_i(s)$$
   integrated numerically in TMB using a differentiable Runge–Kutta / Adams–Bashforth ODE solver (`odeint`), where $\epsilon_{i,t} \sim \mathcal{N}(0, \sigma^2_{\text{proc}, i})$ are random effects marginalized via the **Laplace approximation**.
2. **Simultaneous Joint Likelihood Across All Data Streams**:
   `EcoState` estimates Ecopath baseline parameters ($\ln B_{i,0}, EE_i, DC_{ij,0}$), Ecosim vulnerabilities $\ln(X_{ij}-1)$, fishing mortalities $\ln F_{i,t}$, and process error variances $\sigma_{\text{proc},i}$ simultaneously by fitting to:
   - Survey biomass indices: $\ln I_{i,t} \sim \mathcal{N}\!\left(\ln(q_i \hat{B}_{i,t}),\, \sigma^2_{I,i,t}\right)$,
   - Fishery catches: $\ln C^{\text{obs}}_{i,t} \sim \mathcal{N}\!\left(\ln \hat{C}_{i,t},\, \sigma^2_{C,i}\right)$,
   - Time-varying stomach diet proportions $\mathbf{W}^{\text{obs}}_{j,t}$: modeled via a **Dirichlet likelihood** on predicted Ecosim consumption shares $\hat{p}_{ij,t} = \frac{\hat{Q}_{ij,t}}{\sum_k \hat{Q}_{kj,t}}$!
