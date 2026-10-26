# Age-Structured Multispecies Stock Assessment Models (`MSVPA`, `SMS`/`smsR`, `MSCAA`/`Hydra`, `CEATTLE`) {#sec-age-structured}

## Overview of the Age-Structured Multispecies Family

Age-structured multispecies models are direct extensions of single-species Virtual Population Analysis (VPA) and Statistical Catch-at-Age (SCAA) models. They partition total instantaneous mortality $Z_{s,a,y,q}$ of species $s$, age class $a$, year $y$, and season/quarter $q$ into three additive components (first formalized by [@helgason_gislason1979] and [@pope1979]):

$$Z_{s,a,y,q} = F_{s,a,y,q} + M1_{s,a,q} + M2_{s,a,y,q}$$ {#eq-total-mortality-decomposition}

where:
- $F_{s,a,y,q}$ is the instantaneous fishing mortality rate,
- $M1_{s,a,q}$ is the **residual (background) natural mortality rate** due to disease, senescence, abiotic stress, and predation by unmodeled predators,
- $M2_{s,a,y,q}$ is the **endogenous predation mortality rate** inflicted by the modeled predators in the ecosystem.

Because these models track discrete age cohorts $(s, a)$ rather than continuous body size, individual body weight $w_{s,a,y,q}$ and length $l_{s,a,y,q}$ are typically treated as known empirical inputs (measured from port and survey sampling), or mapped via pre-computed age-length keys ($\text{ALK}_{s,a,l,y,q}$).

Four major model lineages define this family:
- **Foundations and `MSVPA`**: [@ursin1967], [@ursin1973], and [@andersen_ursin1977] established the size-based bioenergetic and prey-suitability equations, evaluated for North Pacific stocks by [@livingston1983] and simplified into operational `MSVPA` by [@helgason_gislason1979], [@pope1979], and [@sparre1980].
- **Stochastic Multispecies Model (`SMS` / `smsR`)**: [@lewy2004] formulated the seasonal age-length statistical multispecies model and Dirichlet stomach-content likelihood used in ICES key runs [@ices2021; @ices2023].
- **North American Statistical Catch-at-Age (`MSASA` / `MSCAA` / `Hydra`)**: [@jurado_molina2005], [@kinzey2009], [@vankirk2010], [@curti2013], and [@trijoulet2019] developed and tested forward-projecting multispecies catch-at-age models across alternative predator functional responses and diet observation error regimes.
- **Climate-Enhanced Multispecies Assessment (`CEATTLE` / `Rceattle`)**: [@holsman_aydin2015], [@holsman2016], and [@adams2022] coupled temperature-dependent von Bertalanffy growth [@essington2001] and Wisconsin bioenergetics [@kitchell1977; @thornton_lessem1978] with age-structured stock assessment and multispecies biological reference points.

---

## Multispecies Virtual Population Analysis (`MSVPA`) and `MSFOR`

### Historical Development
Developed simultaneously by [@helgason_gislason1979] and [@pope1979] following the [@andersen_ursin1977] monograph, **MSVPA** was the standard operational multispecies model of the ICES Multispecies Working Group (MSWG) throughout the 1980s and 1990s, motivating the massive 1981 and 1991 ICES "Years of the Stomach" sampling programs in the North Sea [@stefansson1998; @lewy2004]. **MSFOR** is its forward-projection counterpart used to evaluate management scenarios given recruitment assumptions.

### Mathematical Formulation of MSVPA
Within each quarter $q \in \{1, 2, 3, 4\}$ of year $y$, cohort abundance $N_{i,a,y,q}$ of prey species $i$, age $a$ decays exponentially:

$$N_{i,a,y,q+1} = N_{i,a,y,q} \exp\!\left(-Z_{i,a,y,q}\right), \qquad \bar{N}_{i,a,y,q} = N_{i,a,y,q} \frac{1 - \exp(-Z_{i,a,y,q})}{Z_{i,a,y,q}}$$ {#eq-msvpa-cohort-decay}

Observed catch-at-age in numbers $C_{i,a,y,q}$ is assumed exact (zero observation error), obeying the Baranov catch equation:

$$C_{i,a,y,q} = \frac{F_{i,a,y,q}}{Z_{i,a,y,q}} N_{i,a,y,q} \left(1 - \exp(-Z_{i,a,y,q})\right) = F_{i,a,y,q} \bar{N}_{i,a,y,q}$$ {#eq-msvpa-baranov-catch}

The predation mortality $M2_{i,a,y,q}$ on prey $(i, a)$ by all predator species $j$ and predator ages $b$ is:

$$M2_{i,a,y,q} = \sum_{j} \sum_{b} \frac{\bar{N}_{j,b,y,q} \, R_{j,b,q} \, S_{i,a,j,b,q}}{\Phi_{j,b,y,q}}$$ {#eq-msvpa-predation-mortality}

where:
- $R_{j,b,q}$ is the prescribed per-capita quarterly food ration (in mass) consumed by predator $(j, b)$ in quarter $q$,
- $S_{i,a,j,b,q} \in [0, 1]$ is the **suitability coefficient** of prey $(i, a)$ for predator $(j, b)$ in quarter $q$,
- $\Phi_{j,b,y,q}$ is the **total available (suitable) food biomass** for predator $(j, b)$:

$$\Phi_{j,b,y,q} = \sum_{k} \sum_{c} S_{k,c,j,b,q} \, \bar{N}_{k,c,y,q} \, w_{k,c,y,q} + S_{\text{other}, j, b, q} \, B_{\text{other}, q}$$ {#eq-msvpa-suitable-biomass}

### Iterative Fixed-Point Solution Algorithm
In classic MSVPA, suitability $S_{i,a,j,b,q}$ is not parameterized as a smooth function of predator-prey size ratio; instead, it is estimated directly from observed average stomach content weight proportions $\bar{W}_{i,a,j,b,y,q}$ via fixed-point iteration:
1. **Initialize**: Set $M2^{(0)}_{i,a,y,q} = 0$ (or a prior guess) and run standard backward VPA on each species to obtain initial cohort abundances $\bar{N}^{(0)}_{i,a,y,q}$ and fishing mortalities $F^{(0)}_{i,a,y,q}$.
2. **Update Suitabilities**: For years $y \in \mathcal{Y}_{\text{stom}}$ with stomach content observations $\bar{W}_{i,a,j,b,y,q}$, compute empirical suitabilities:
   $$\tilde{S}_{i,a,j,b,y,q} = \frac{\bar{W}_{i,a,j,b,y,q}}{\bar{N}_{i,a,y,q} \, w_{i,a,y,q}}$$ {#eq-msvpa-suitability-update}
   Normalize so $\sum_{i,a} \tilde{S}_{i,a,j,b,y,q} = 1$, and average across stomach-sampling years $y \in \mathcal{Y}_{\text{stom}}$ (weighted by predator stomach sample size) to obtain time-invariant quarter-specific suitabilities $S^{(m)}_{i,a,j,b,q}$.
3. **Update $M2$**: Substitute $S^{(m)}_{i,a,j,b,q}$ and $\bar{N}^{(m-1)}$ into the $M2$ equation to compute $M2^{(m)}_{i,a,y,q}$.
4. **Update VPA**: Re-run backward VPA for all species using $Z^{(m)}_{i,a,y,q} = F^{(m)}_{i,a,y,q} + M1_{i,a,q} + M2^{(m)}_{i,a,y,q}$ to obtain $\bar{N}^{(m)}_{i,a,y,q}$.
5. **Iterate**: Repeat Steps 2–4 until $\max |M2^{(m)} - M2^{(m-1)}| < \epsilon$.

### Limitations of MSVPA
As analyzed by [@stefansson1998] and [@lewy2004]:
- **No Observation Error**: MSVPA treats commercial catch-at-age and stomach content proportions as exact deterministic quantities. Sampling noise in rare stomach items propagates directly into suitability $S_{i,a,j,b,q}$.
- **Over-Parameterization of Suitability**: Estimating a free suitability parameter for every $(i, a, j, b, q)$ combination requires thousands of parameters and leaves zero degrees of freedom for prey age classes not observed in a particular predator age's stomachs.
- **No Uncertainty Estimates**: Because there is no statistical likelihood function, MSVPA cannot produce confidence intervals on $M2$, SSB, or recruitment, nor can it perform formal hypothesis testing (e.g., AIC/LRT between functional responses).

---

## The Stochastic Multispecies Model (`SMS` and `smsR`)

To overcome the deterministic limitations of MSVPA, Peter Lewy and Morten Vinther at DTU Aqua developed **SMS** (Stochastic Multispecies Model; [@lewy2004]), implemented in **ADMB** (Auto Differentiation Model Builder) and recently ported to **TMB/R** as `smsR` (Nis Sand Jacobsen et al.; [@ices2021; @ices2023]).

SMS is the flagship operational multispecies stock assessment model in ICES, providing the official time-varying natural mortality matrices ($M = M1 + M2$) used in single-species stock assessments (e.g., `SAM` and `SS3`) for the North Sea and Baltic Sea [@ices2023].

### State-Space Indexing and Population Dynamics
SMS is a forward-projecting, seasonal (typically quarterly, $q \in \{1, 2, 3, 4\}$), age-structured statistical catch-at-age model. In addition to tracked dynamic VPA species (which can act as both predators and prey), SMS includes **"other predators"** (e.g., seabirds, grey seals, harbour porpoises, horse mackerel, hake) whose abundance time series $N_{A,a,y,q}$ are entered as known external covariates, and an **"other food"** biomass pool $B_{\text{OF}}$ [@lewy2004, p. 7].

Following the notation of [@lewy2004]:
- Let $(A, a)$ denote a predator species $A$ of age (or size class) $a$, with mean weight $w_{A,a,y,q}$.
- Let $(B, b)$ denote a prey species $B$ of age (or size class) $b$, with mean weight $w_{B,b,y,q}$.

Stock survival across seasons $q$ within year $y$ follows:
$$N_{s,a,y,q+1} = N_{s,a,y,q} \exp\!\left(-Z_{s,a,y,q}\right), \qquad q = 1, 2, 3$$ {#eq-sms-cohort-quarter-step}
At the transition from the final quarter $q=4$ of year $y$ to $q=1$ of year $y+1$:
$$N_{s,a+1,y+1,1} = N_{s,a,y,4} \exp\!\left(-Z_{s,a,y,4}\right), \qquad a < a_{\max}-1$$ {#eq-sms-cohort-annual-step}
and for the plus-group $a_{\max}$:
$$N_{s,a_{\max},y+1,1} = N_{s,a_{\max}-1,y,4} \exp\!\left(-Z_{s,a_{\max}-1,y,4}\right) + N_{s,a_{\max},y,4} \exp\!\left(-Z_{s,a_{\max},y,4}\right)$$ {#eq-sms-plusgroup-step}

Recruitment $N_{s,a_{\text{rec}},y,q_{\text{rec}}}$ is parameterized either as free log-parameters $\ln R_{s,y}$ or via a stochastic stock-recruitment relationship (Ricker, Beverton-Holt, or hockey-stick) with lognormal process error:
$$N_{s,a_{\text{rec}},y,q_{\text{rec}}} = f_{\text{SRR}}\!\left(\text{SSB}_{s,y-a_{\text{rec}}}\right) \exp\!\left(\epsilon_{R,s,y}\right), \qquad \epsilon_{R,s,y} \sim \mathcal{N}\!\left(0, \sigma^2_{R,s}\right)$$ {#eq-sms-stock-recruitment}

### Separable Fishing Mortality Model
To prevent overfitting catch-at-age noise, SMS models instantaneous fishing mortality $F_{s,a,y,q}$ using a partially separable structure across age $a$, year $y$, and season $q$ within multi-year regime blocks $Y_j$ [@lewy2004, Eq. 2]:

$$F_{s,a,y,q} = v_{s,a,Y_j} \, F^{(y)}_{s,y} \, F^{(q)}_{s,a,Y_j,q}$$ {#eq-sms-fishing-mortality-separable}

where:
- $v_{s,a,Y_j}$ is the age effect (fishery selectivity at age $a$ in year-block $Y_j$),
- $F^{(y)}_{s,y}$ is the annual fishing effort/mortality multiplier in year $y$,
- $F^{(q)}_{s,a,Y_j,q}$ is the seasonal allocation fraction of fishing mortality in quarter $q$ (constrained so the mean across seasons equals $0.25$).

Expected commercial catch-at-age $\hat{C}_{s,a,y,q}$ is given by the Baranov equation:
$$\hat{C}_{s,a,y,q} = \frac{F_{s,a,y,q}}{Z_{s,a,y,q}} N_{s,a,y,q} \left(1 - \exp(-Z_{s,a,y,q})\right) = F_{s,a,y,q} \bar{N}_{s,a,y,q}$$ {#eq-sms-expected-catch}

### Parametric Suitability and the Andersen–Ursin Size-Selection Submodel
Instead of MSVPA's non-parametric empirical suitabilities, SMS parameterizes the suitability $S(B, b, A, a, q)$ of prey $(B, b)$ for predator $(A, a)$ in season $q$ as a multiplicative product of three biological components [@lewy2004, Eq. 5]:

$$S(B, b, A, a, y, q) = \text{vul}(B, A) \times \text{size}(B, b, A, a, y, q) \times \text{overlap}(B, A, q)$$ {#eq-sms-suitability-components}

1. **Species Vulnerability ($\text{vul}(B, A)$)**: Time-invariant parameter measuring the behavioral/morphological vulnerability of prey species $B$ to predator species $A$ (with $\text{vul}(\text{Other Food}, A) \equiv 1$ as the reference baseline).
2. **Spatial/Seasonal Overlap ($\text{overlap}(B, A, q)$)**: Parameter representing the spatial co-occurrence of prey $B$ and predator $A$ in season $q$ (fixed to $1$ in seasons with full overlap, estimated in seasons with migratory separation).
3. **Asymmetric Lognormal Size Preference ($\text{size}(B, b, A, a, y, q)$)**: Following [@andersen_ursin1977], SMS models size preference as a Gaussian (or asymmetric split-Gaussian) function of the log predator-to-prey weight ratio $\xi = \ln\!\left(\frac{w_{A,a,y,q}}{w_{B,b,y,q}}\right)$, restricted to the physiologically feasible range $[\xi_{\min, A}, \xi_{\max, A}]$ [@lewy2004, Eq. 6]:

$$\text{size}(B, b, A, a, y, q) = \begin{cases} \exp\!\left(-\dfrac{\left(\ln\frac{w_{A,a}}{w_{B,b}} - \eta_A\right)^2}{2\sigma_{1,A}^2}\right) & \text{if } \xi_{\min,A} \le \ln\frac{w_{A,a}}{w_{B,b}} \le \eta_A \\[12pt] \exp\!\left(-\dfrac{\left(\ln\frac{w_{A,a}}{w_{B,b}} - \eta_A\right)^2}{2\sigma_{2,A}^2}\right) & \text{if } \eta_A < \ln\frac{w_{A,a}}{w_{B,b}} \le \xi_{\max,A} \\[10pt] 0 & \text{otherwise} \end{cases}$$ {#eq-sms-size-preference-kernel}

When SMS is run in **length-class stomach mode** (because ICES stomach databases record predator and prey *lengths* $l$, not ages $a$), SMS converts age abundances $N_{s,a,y,q}$ into length classes $l$ using seasonal empirical Age-Length Keys ($\text{ALK}_{s,a,l,y,q}$) before evaluating $\text{size}(B, l_B, A, l_A)$ using constant length-weight relationships $w_s(l) = \alpha_s l^{\beta_s}$, and then maps $M2$ back from length classes to age classes [@lewy2004, pp. 8--9].

### Gastric-Evacuation-Based Consumption Rate ($R_{A,a,y,q}$)
Total quarterly food intake (ration) $R_{A,a,y,q}$ of predator $(A,a)$ is input as an empirical covariate derived from average observed stomach fullness $\bar{S}_{A,a,q}$ and temperature-dependent gastric evacuation models [@andersen_beyer2005; @lewy2004, Eqs. 7--8]:

$$R_{A,a,y,q} = \alpha_{R,A} \, w_{A,a,y,q}^{\beta_{R,A}} \quad \text{or} \quad R_{A,a,y,q} = \rho_A \, \bar{S}_{A,a,q}^{\gamma_A} \, \exp(\delta_A T_{y,q}) \, L_{A,a}^{\lambda_A} \, D_q$$ {#eq-sms-predator-ration}

where $T_{y,q}$ is ambient water temperature and $D_q$ is the number of hours in quarter $q$.

### Predation Mortality ($M2$) and Exact Integration Over the Quarter
Let $\bar{B}_{B,b,y,q} = \bar{N}_{B,b,y,q} w_{B,b,y,q}$ be the mean quarterly biomass of prey $(B, b)$, and let total available food for predator $(A, a)$ be:

$$\Phi(A, a, y, q) = \sum_{B} \sum_{b} \bar{B}_{B,b,y,q} \, S(B, b, A, a, y, q) + B_{\text{OF}, A} \, S(\text{OF}, A, a, y, q)$$ {#eq-sms-total-food-biomass}

In differential form at instant $t \in [q, q+1]$, the instantaneous predation rate is $M2_{B,b}(t) = \sum_{A,a} \frac{N_{A,a}(t) r_{A,a} S(B,b,A,a)}{\Phi(A,a,t)}$. Assuming $\Phi(A,a,t)$ is approximately constant within a single quarter (since $B_{\text{OF}}$ is large and individual prey decline is offset across cohorts), integrating over the quarter yields the exact discrete-time SMS $M2$ equation [@lewy2004, Eq. 11]:

$$M2_{B,b,y,q} = \sum_{A} \sum_{a} \frac{\bar{N}_{A,a,y,q} \, R_{A,a,y,q} \, S(B, b, A, a, y, q)}{\Phi(A, a, y, q)}$$ {#eq-sms-predation-mortality-m2}

Since $\bar{N}_{B,b,y,q} = N_{B,b,y,q}\frac{1 - e^{-Z_{B,b,y,q}}}{Z_{B,b,y,q}}$ depends on $M2_{B,b,y,q}$ (and predator survival also depends on $M2$ for intermediate predators like whiting or young cod), SMS performs $2\text{--}4$ internal Newton/fixed-point iterations per quarter evaluation within the ADMB/TMB objective function, or uses the explicit beginning-of-quarter biomass approximation when quarterly time steps are short.

### Expected Stomach Content Proportions
Under the SMS equations, the expected proportion by weight $\hat{W}(B, b, A, a, y, q)$ of prey $(B, b)$ in the stomach of predator $(A, a)$ in year $y$, quarter $q$ is [@lewy2004, Eq. 12]:

$$\hat{W}(B, b, A, a, y, q) = \frac{\bar{B}_{B,b,y,q} \, S(B, b, A, a, y, q)}{\Phi(A, a, y, q)}$$ {#eq-sms-stomach-diet-proportion}

Note that $\sum_{B,b} \hat{W}(B, b, A, a, y, q) + \hat{W}(\text{OF}, A, a, y, q) = 1$.

### Statistical Likelihood Function of SMS
SMS estimates all parameters $\boldsymbol{\theta} = (\ln N_{\text{init}}, \ln R_{s,y}, v_{s,a}, F^{(y)}_{s,y}, F^{(q)}_{s,q}, q_{\text{survey}}, \text{vul}(B,A), \eta_A, \sigma_A, \text{overlap}, \dots)$ by minimizing the negative joint log-likelihood across four data types [@lewy2004, pp. 10--16]:

$$\mathcal{L}_{\text{total}}(\boldsymbol{\theta}) = \mathcal{L}_{\text{catch}}(\boldsymbol{\theta}) + \mathcal{L}_{\text{CPUE}}(\boldsymbol{\theta}) + \mathcal{L}_{\text{stomach}}(\boldsymbol{\theta}) + \mathcal{L}_{\text{SRR}}(\boldsymbol{\theta})$$ {#eq-sms-objective-function}

#### 1. Commercial Catch-at-Age Likelihood (Lognormal)
Observed catch-at-age $C_{s,a,y,q}$ is assumed lognormally distributed around expected catch $\hat{C}_{s,a,y,q}$:
$$\ln C_{s,a,y,q} \sim \mathcal{N}\!\left(\ln \hat{C}_{s,a,y,q},\; \sigma^2_{\text{catch}, s, a, q}\right)$$ {#eq-sms-catch-observation-model}
contributing to the negative log-likelihood:
$$-\ln L_{\text{catch}} = \frac{1}{2} \sum_{s,a,y,q} \left[ \ln\!\left(2\pi \sigma^2_{\text{catch},s,a,q}\right) + \frac{\left(\ln C_{s,a,y,q} - \ln \hat{C}_{s,a,y,q}\right)^2}{\sigma^2_{\text{catch},s,a,q}} \right]$$ {#eq-sms-catch-loglik}

#### 2. Survey Abundance Index (CPUE) Likelihood (Lognormal)
For survey fleet $f$ operating from season fraction $\alpha_{f}$ to $\beta_{f}$ in quarter $q$, the expected survey index $\hat{I}_{f,s,a,y}$ is:
$$\hat{I}_{f,s,a,y} = Q_{f,s,a} \, N_{s,a,y,q} \frac{\exp(-\alpha_f Z_{s,a,y,q}) - \exp(-\beta_f Z_{s,a,y,q})}{(\beta_f - \alpha_f) Z_{s,a,y,q}}$$ {#eq-sms-survey-cpue-index}
with lognormal negative log-likelihood:
$$-\ln L_{\text{CPUE}} = \frac{1}{2} \sum_{f,s,a,y} \left[ \ln\!\left(2\pi \sigma^2_{\text{CPUE},f,s,a}\right) + \frac{\left(\ln I_{f,s,a,y} - \ln \hat{I}_{f,s,a,y}\right)^2}{\sigma^2_{\text{CPUE},f,s,a}} \right]$$ {#eq-sms-cpue-loglik}

#### 3. Stomach Content Composition Likelihood (Dirichlet vs. Multinomial)
Let $U = (A, a, y, q)$ index a specific predator-age-year-quarter sampling stratum, and let $k \in \{1, \dots, K_U\}$ index the $K_U$ prey entities $(B, b)$ (plus "other food") observed in that predator's diet, with observed weight proportions $\mathbf{W}_U = (W_{1,U}, \dots, W_{K_U,U})$ satisfying $\sum_{k=1}^{K_U} W_{k,U} = 1$, and model expectations $\hat{p}_{k,U} = \hat{W}(k, U)$.

[@lewy2004, pp. 11--14] derived two alternative likelihood formulations for continuous stomach weight proportions:

- **Option A: Multinomial Approximation on Effective Prey Counts**
  Defining an effective number of independent prey observations $N_{\text{stom}, U}$ (e.g., number of stomachs or hauls sampled), $N_{\text{stom}, U} \mathbf{W}_U$ is treated as a multinomial vector with log-likelihood kernel:
  $$\ln L_{\text{stom, Mult}} = \sum_{U} N_{\text{stom}, U} \sum_{k=1}^{K_U} W_{k,U} \ln \hat{p}_{k,U}$$ {#eq-sms-multinomial-stomach-loglik}

- **Option B: Exact Dirichlet Distribution (Default in SMS)**
  Because $W_{k,U} \in (0, 1)$ is a continuous biomass proportion rather than a discrete integer count, [@lewy2004, Eqs. 19--24] modeled $\mathbf{W}_U$ using a **Dirichlet distribution** with concentration parameters $\alpha_{k,U} = v_U \hat{p}_{k,U}$ and total precision sum $v_U = \sum_{k=1}^{K_U} \alpha_{k,U}$:

  $$f_{\text{Dir}}(\mathbf{W}_U \mid v_U \hat{\mathbf{p}}_U) = \frac{\Gamma(v_U)}{\prod_{k=1}^{K_U} \Gamma(v_U \hat{p}_{k,U})} \prod_{k=1}^{K_U} W_{k,U}^{v_U \hat{p}_{k,U} - 1}$$ {#eq-sms-dirichlet-density}

  Under the Dirichlet distribution, the expectation and variance of each prey share $W_{k,U}$ are:
  $$\mathbb{E}[W_{k,U}] = \hat{p}_{k,U}, \qquad \text{Var}(W_{k,U}) = \frac{\hat{p}_{k,U}(1 - \hat{p}_{k,U})}{v_U + 1}$$ {#eq-sms-dirichlet-moments}
  Defining the variance-scaling factor $V_U = \frac{1}{v_U + 1} \iff v_U = \frac{1}{V_U} - 1$, SMS relates $V_U$ to the predator species' overall diet dispersion parameter $\theta_{\text{stom}, A}$ and the number of stomachs/hauls sampled $N_{\text{haul}, U}$ via:
  $$V_U = \frac{\theta_{\text{stom}, A}}{N_{\text{haul}, U}} \implies v_U = \frac{N_{\text{haul}, U}}{\theta_{\text{stom}, A}} - 1$$ {#eq-sms-dirichlet-precision}
  This allows $\theta_{\text{stom}, A}$ to be estimated directly inside the maximum likelihood optimization via the exact Gamma-function log-likelihood! (Wait: to handle observed zeros $W_{k,U} = 0$, SMS pools unobserved interior age classes with adjacent age classes or applies a small threshold; see @sec-stomach-data for full details).

---

## Statistical Multispecies Catch-at-Age (`MSASA`, `MSCAA`, and `Hydra`)

In North America (particularly NOAA AFSC for the Aleutian Islands, Gulf of Alaska, and Bering Sea, and NOAA NEFSC for Georges Bank and the Gulf of Maine), a family of **Multispecies Statistical Catch-at-Age (`MSASA` / `MSCAA`)** models was developed in ADMB and TMB by [@kinzey2009], [@vankirk2010], [@curti2013], Gaichas et al. (`Hydra` MSE operating model), and evaluated in simulation studies by [@trijoulet2019].

### Kinzey & Punt (2009): Aleutian Islands Age–Length Hybrid `MSCAA` and the Seven Predator Functional Responses
[@kinzey2009] developed a foundational integrated multispecies statistical catch-at-age framework in ADMB and applied it to the Aleutian Islands (1960–2003) for walleye pollock (*Theragra chalcogramma*, ages 0–15+, $1\text{--}80\text{ cm}$ in $5\text{ cm}$ bins), Atka mackerel (*Pleurogrammus monopterygius*, ages 0–15+, $1\text{--}65\text{ cm}$), and Pacific cod (*Gadus macrocephalus*, ages 0–12+, $1\text{--}120\text{ cm}$, harvested by 3 fleets: trawl, pot, and longline).

1. **Predation Mortality ($P_{k,a,y}$), Gamma Length-Ratio Selectivity ($\gamma^{r,u}_{k,a}$), and Age–Length Transition Mapping ($\chi^{r,u,l}$)**:
   Total mortality on prey species $k$, age $a$ in year $y$ is $Z_{k,a,y} = \sum_f F^f_{k,a,y} + P_{k,a,y} + M_k$, where predation mortality $P_{k,a,y}$ inflicted by all predator species $r$ and predator ages $u$ is [@kinzey2009, Eqs. 6--7, p. 77]:
   $$P_{k,a,y} = \sum_{r,u} V^{r,u}_{k,a,y}, \qquad V^{r,u}_{k,a,y} = N_{r,u,y} \, \phi^{r,u}_{k,a,y} \, \gamma^{r,u}_{k,a}$$ {#eq-jurado-predation-capacity}
   Instead of a symmetric lognormal weight-ratio kernel, [@kinzey2009, Eq. 18, p. 80] parameterize predator $(r, u)$ size selectivity $\gamma^{r,u}_{k,a} \in (0, 1]$ for prey $(k, a)$ as a **skewed Gamma function of the log predator-to-prey expected length ratio** $G^{r,u}_{k,a} = \ln(\bar{L}_{r,u} / \bar{L}_{k,a})$, normalized to equal $1$ at its mode $\tilde{G}^r = (\alpha^r - 1)\beta^r$:
   $$\gamma^{r,u}_{k,a} = \left(\frac{G^{r,u}_{k,a}}{\tilde{G}^r}\right)^{\alpha^r - 1} \exp\!\left[-\frac{G^{r,u}_{k,a} - \tilde{G}^r}{\beta^r}\right], \qquad \tilde{G}^r = (\alpha^r - 1)\beta^r$$ {#eq-jurado-gamma-size-preference}
   The total number of prey $(k, a)$ eaten by predator $(r, u)$ in year $y$ follows the Baranov predation equation:
   $$E^{r,u}_{k,a,y} = \frac{V^{r,u}_{k,a,y}}{Z_{k,a,y}} N_{k,a,y} \left(1 - e^{-Z_{k,a,y}}\right)$$ {#eq-jurado-baranov-consumption}
   Because stomach samples ($1,987$ pollock, $1,033$ mackerel, and $3,278$ cod stomachs) are recorded by **$5\text{ cm}$ predator length bin $l$** and **$5\text{ cm}$ prey length bin $m$**, Kinzey & Punt map age-specific predation into length space using prespecified empirical age-to-length transition matrices $\chi^{r,u,l}$ [@kinzey2009, Eqs. 9--14 & Table 3, pp. 77--84]:
   - **Prey numbers eaten by predator length class $(r, l)$**: $E^{r,l}_{k,a,y} = \sum_u E^{r,u}_{k,a,y} \chi^{r,u,l}$.
   - **Predicted prey length composition $\tilde{T}^{r,l}_{k,m}$ inside predator $(r, l)$ stomachs** (fitted in multinomial likelihood $L_6$):
     $$\tilde{T}^{r,l}_{k,m} = \frac{\sum_y n^{r,l}_y \sum_a \chi^{k,a,m} E^{r,l}_{k,a,y}}{\sum_y n^{r,l}_y \sum_{a'} E^{r,l}_{k,a',y}}$$ {#eq-jurado-length-stomach-proportion}
   - **Predicted prey species mass $Q^{r,u}_{k,y}$ and diet weight proportions $\tilde{Q}^{r,l}_{k,y}$** (fitted in multinomial likelihood $L_7$):
     $$Q^{r,u}_{k,y} = \begin{cases} \sum_{a=0}^{A_{L,k}} E^{r,u}_{k,a,y} w_{k,a} & \text{modeled prey } k \\[6pt] Q^r_0 \left(1 - e^{-\phi^{r,u}_o N_{r,u,y}}\right) & \text{"other" prey} \end{cases}, \qquad Q^{r,l}_{k,y} = \sum_u Q^{r,u}_{k,y} \chi^{r,u,l}, \qquad \tilde{Q}^{r,l}_{k,y} = \frac{Q^{r,l}_{k,y}}{\sum_{k'} Q^{r,l}_{k',y}}$$ {#eq-jurado-annual-consumption-biomass}
   - **Predicted vs. Bioenergetic Daily Ration ($\hat{\Omega}^{r,u}_y$, fitted in likelihood $L_5$)**:
     $$\hat{\Omega}^{r,u}_y = \frac{\sum_k Q^{r,u}_{k,y}}{365 \, N_{r,u,y} e^{-Z_{r,u,y}/2}}, \qquad L_5 = \sum_{r,u} \frac{\left(\Omega^{r,u} - \bar{\hat{\Omega}}^{r,u}\right)^2}{2\sigma_d^2}$$ {#eq-jurado-daily-ration-penalty}
     where $\Omega^{r,u} = \frac{1}{365}\frac{K^r}{A^r}(w^{r,u})^{d^r}(w^{r,\infty})^{d^r - 1}$ is the [@essington2001] bioenergetic daily ration!

2. **Unified Equilibrium-Normalized Parameterization of Seven Functional Responses (Table 1, p. 79)**:
   To compare prey-dependent and predator-dependent functional responses on an equal footing, [@kinzey2009, Eqs. 15--17 & Table 1] define the predator age-effect $\theta^{r,u} = 1 + \frac{\omega^r \tilde{\omega}^r}{u + \tilde{\omega}^r}$ and normalize available prey abundance $\Phi^{r,u}_y$ and competing predator abundance $\Psi^{k,a}_y$ **relative to their unfished equilibrium values in year $y_P = 1960$**:
   $$\Phi^{r,u}_y = \frac{\sum_{k,a} N_{k,a,y} \gamma^{r,u}_{k,a}}{\sum_{k,a} N_{k,a,y_P} \gamma^{r,u}_{k,a}}, \qquad \Psi^{k,a}_y = \frac{\sum_{r,u} N_{r,u,y} \gamma^{r,u}_{k,a}}{\sum_{r,u} N_{r,u,y_P} \gamma^{r,u}_{k,a}}$$ {#eq-jurado-available-food-index}
   Because $\Phi^{r,u}_{y_P} = \Psi^{k,a}_{y_P} = 1$ at unfished equilibrium, all seven functional responses for $\phi^{r,u}_{k,a,y}$ share the exact same unfished baseline $\phi^{r,u}_{k,a,y_P} = v^r_k \theta^{r,u}$ (see @tbl-mscaa-functional-responses):

   | Type | Functional Response Name | Mathematical Formula for $\phi^{r,u}_{k,a,y}$ | Params | Aleutian $\text{AIC}$ ($\Delta\text{AIC}$) |
   | :--- | :--- | :--- | :---: | :---: |
   | **I** | **Linear (Lotka–Volterra)** | $\theta^{r,u} v^r_k$ | $1$ | $7411$ ($+716$) |
   | **II** | **Asymptotic (Holling Type II)** | $\dfrac{v^r_k \theta^{r,u} \left[1 + \tilde{v}^r_k\right]}{1 + \tilde{v}^r_k \Phi^{r,u}_y}$ | $2$ | $7080$ ($+384$) |
   | **III** | **Sigmoid (Holling Type III)** | $\dfrac{v^r_k \theta^{r,u} \left(1 + \tilde{v}^r_k\right) \left(\Phi^{r,u}_y\right)^{\tilde{\tilde{v}}^r_k - 1}}{1 + \tilde{v}^r_k \left(\Phi^{r,u}_y\right)^{\tilde{\tilde{v}}^r_k}}$ | $3$ | $7306$ ($+611$) |
   | **IV** | **Predator Interference (Beddington–DeAngelis)** | $\dfrac{v^r_k \theta^{r,u} \left[1 + \tilde{v}^r_k\right]}{1 + \tilde{v}^r_k \Phi^{r,u}_y + \tilde{\tilde{v}}^r_k \left(\Psi^{k,a}_y - 1\right)}$ | $3$ | $7098$ ($+403$) |
   | **V** | **Predator Pre-emption (Multiplicative Interference)** | $\dfrac{v^r_k \theta^{r,u} \left[1 + \tilde{v}^r_k\right]}{\left(1 + \tilde{v}^r_k \Phi^{r,u}_y\right)\left[1 + \tilde{\tilde{v}}^r_k \left(\Psi^{k,a}_y - 1\right)\right]}$ | $3$ | **$6696$ ($0$)** |
   | **VI** | **Hassell–Varley** | $\dfrac{v^r_k \theta^{r,u} \left[1 + \tilde{v}^r_k\right]}{\tilde{v}^r_k \Phi^{r,u}_y + \left(\Phi^{r,u}_y\right)^{\tilde{\tilde{v}}^r_k}}$ | $3$ | $7169$ ($+473$) |
   | **VII** | **Ecosim (Foraging Arena Vulnerability)** | $\dfrac{v^r_k \theta^{r,u}}{1 + \tilde{v}^r_k \left(\Psi^{k,a}_y - 1\right)}$ | $2$ | $7046$ ($+351$) |

   : Comparison and Empirical AIC of Seven Functional Responses in Aleutian Islands Walleye Pollock and Atka Mackerel [@kinzey2009] {#tbl-mscaa-functional-responses}

3. **Key Empirical Findings [@kinzey2009, pp. 89--101]**:
   - **Dominance of Predator-Dependent Functional Responses (Types V & VII)**: The two best-fitting functional responses by AIC were **Type V (Predator Pre-emption, $\Delta\text{AIC} = 0$)** and **Type VII (Ecosim Foraging Arena, $\Delta\text{AIC} = 351$)**, both of which outperformed standard prey-only Holling Type II ($\Delta\text{AIC} = 384$) and Type III ($\Delta\text{AIC} = 611$). In Type V, predator interference $\Psi^{k,a}_y$ multiplicatively dampens per-capita attack rates at low and intermediate prey densities, preventing unrealistic over-depletion when predator stocks are large.
   - **Divergence Between Single-Species and Multispecies Trajectories**: The multispecies model estimated that Aleutian pollock spawning biomass declined **over $3\times$ more steeply** from 1964 to 2003 than estimated by the single-species model (with non-overlapping 95% confidence bands in the 1960s–1980s, Fig. 9), and estimated age-0 pollock predation mortality $P_{1,0,y} \approx 5.0\text{ yr}^{-1}$ (only $\sim 0.5\%$ survival from age 0 to age 1).

---

### Van Kirk et al. (2010): Gulf of Alaska `MSASA` and Simultaneity vs. Identifiability
[@vankirk2010] constructed an ADMB-based **Multispecies Age-Structured Assessment (`MSASA`)** model for three strongly interacting species in the Gulf of Alaska (1981–2001): walleye pollock (*Theragra chalcogramma*, ages 1–10+), Pacific cod (*Gadus macrocephalus*, ages 1–12+), and arrowtooth flounder (*Atheresthes stomias*, ages 1–15+), fitted to commercial catch, NMFS triennial bottom-trawl survey data, and $2,704$ predator stomachs from the AFSC REEM database.

1. **Predation Mortality ($P_{i,a,t}$) and Beginning-of-Year Approximation**:
   Total mortality is $Z_{i,a,t} = M^r_{i,a,t} + F_{i,a,t} + P_{i,a,t}$. Because true instantaneous predation mortality $P^{\text{inst}}_{i,a,t} = \frac{1}{\bar{B}_{i,a,t}}\sum_j \sum_b I_{j,b} \bar{N}_{j,b,t} \frac{\varphi_{i,a,j,b,t}}{\varphi_{j,b,t}}$ is implicitly recursive in $\bar{N}_{i,a,t}$ and $\bar{B}_{i,a,t}$, `MSASA` uses the explicit beginning-of-year abundance approximation [@vankirk2010, Eqs. 1--6]:
   $$P_{i,a,t} = \frac{1}{B_{i,a,t}} \sum_{j} \sum_{b} I_{j,b} \, N_{j,b,t} \, \frac{\varphi_{i,a,j,b,t}}{\varphi_{j,b,t}}$$ {#eq-curti-predation-mortality}
   where $I_{j,b}$ ($\text{kg} \cdot \text{yr}^{-1}$) is the annual predator ration derived from [@essington2001] / [@aydin2005] von Bertalanffy bioenergetics, and available prey biomass is:
   $$\varphi_{i,a,j,b,t} = \nu_{i,a,j,b} \, B_{i,a,t}, \qquad \varphi_{j,b,t} = B_{\text{oth}} + \sum_{i} \sum_{a} \nu_{i,a,j,b,t} \, B_{i,a,t}$$ {#eq-curti-available-biomass}
   with constant non-modeled prey biomass $B_{\text{oth}} = e^{15}$ ($\nu_{\text{oth}} \equiv 1$) and Andersen–Ursin suitability:
   $$\nu_{i,a,j,b} = \rho_{i,j} \, g_{i,a,j,b}, \qquad g_{i,a,j,b} = \exp\!\left[-\frac{1}{2\sigma_j^2}\left(\ln\frac{w_{j,b}}{w_{i,a}} - \eta_j\right)^2\right]$$ {#eq-curti-suitability-weight-ratio}

2. **Joint Objective Function and $\sqrt{n}$-Weighted Stomach Least Squares**:
   Unlike `SMS` (which uses a Dirichlet stomach likelihood), `MSASA` combines lognormal total catch and survey biomass, multinomial catch-at-age and survey-at-age/length (with fixed effective sample size $w = 100$), and a **weighted least-squares (normal) objective function on stomach content weight proportions** [@vankirk2010, Table 4, p. 1139]:
   $$\mathcal{L}_{\text{stom}} = \sum_{p=1}^{3} w_{\text{stom}, p} \sum_{i,a,j,b,t \in p} \left[ \left(\frac{\varphi_{i,a,j,b,t}}{\varphi_{j,b,t}} + 0.0001\right) - \left(\frac{\hat{\varphi}_{i,a,j,b,t}}{\hat{\varphi}_{j,b,t}} + 0.0001\right) \right]^2$$ {#eq-curti-stomach-least-squares}
   where stomach proportions are pooled into three 7-year blocks ($p \in \{1981\text{--}1987,\, 1988\text{--}1994,\, 1995\text{--}2001\}$) with weights set to the square root of predator stomach sample sizes ($w_{\text{stom}} = \{11, 26, 44\}$).

3. **Key Lessons on Parameter Estimability [@vankirk2010, pp. 1140--1145]**:
   - **Temporal Pooling of Stomach Data is Essential for Convergence**: Attempting to fit unpooled annual stomach proportions failed to converge in ADMB, and pooling all 21 years into a single time-invariant average (`M1`) produced a singular Hessian unless stomach weights were inflated so high that catch and survey fits degraded. Dividing diets into 7-year regime blocks (`M2`–`M4`) with time-varying species preference $\rho_{i,j}$ (`M3`, $\text{AIC} = 2121.68$) and size preference $\eta_j, \sigma_j$ (`M4`, $\text{AIC} = 2050.65$) achieved positive-definite Hessians and captured the tripling of arrowtooth flounder predation on age-1 pollock.
   - **Confounding of Survey Catchability ($Q_i$) and Residual Mortality ($M^r_i$) with $P_{i,a,t}$**: When Van Kirk et al. attempted to estimate survey catchability $Q_i$ or residual natural mortality $M^r_i$ freely alongside predation mortality $P_{i,a,t}$, ADMB pushed prey abundance upward, drove $Q_i \to 0$, or counter-intuitively *increased* $M^r_i$ to absorb the extra older cohorts generated by high juvenile recruitment. Thus, $Q_i$ and $M^r_i$ had to be fixed externally.

---

### Curti et al. (2013): Georges Bank `MSCAA`, Ecosystem Biomass Closure, and Monte Carlo Error Propagation
[@curti2013] developed an ADMB `MSCAA` model for Georges Bank (1978–2007) covering Atlantic cod (*Gadus morhua*), silver hake (*Merluccius bilinearis*), and Atlantic herring (*Clupea harengus*), introducing three important methodological refinements:

1. **Sparre (1980) Stationary Total Ecosystem Biomass Closure**:
   Instead of assuming a constant biomass of "other food" $B_{\text{other}}$, [@curti2013, Eqs. 6--10] followed [@sparre1980] and [@tsou_collie2001] in assuming a stationary **total ecosystem biomass** $\text{EB} = 15 \times 10^6\text{ metric tons}$ (and tested $\text{EB} = 7.34 \times 10^6\text{ t}$ from [@rochet2011] in a sensitivity analysis), so that non-modeled "other food" varies dynamically each year as the complement of modeled fish biomass:
   $$B_{\text{other}, t} = \text{EB} - \sum_{i} \sum_{a} B_{i,a,t}$$ {#eq-curti-sparre-other-food}
   Suitabilities $v_{i,a,j,b,t} = \rho_{i,j} g_{i,a,j,b,t}$ (with $v_{\text{other}} = \rho_{\text{other}} g_{\text{other}} = 1$) are normalized across all prey categories for each predator $(j, b)$:
   $$\tilde{v}_{i,a,j,b,t} = \frac{v_{i,a,j,b,t}}{\sum_{k}\sum_{c} v_{k,c,j,b,t} + v_{\text{other}}}, \qquad \phi_{i,a,j,b,t} = \tilde{v}_{i,a,j,b,t} B_{i,a,t}, \qquad \phi_{j,b,t} = \tilde{v}_{\text{other}} B_{\text{other}, t} + \sum_{k,c} \phi_{k,c,j,b,t}$$ {#eq-curti-normalized-vulnerability}
   and predation mortality is driven by age-specific consumption-to-biomass ratios $(CB)_{j,b}$:
   $$M2_{i,a,t} \approx \frac{1}{N_{i,a,t} W_{i,a,t}} \sum_{j} \sum_{b} (CB)_{j,b} \, B_{j,b,t} \, \frac{\phi_{i,a,j,b,t}}{\phi_{j,b,t}}$$ {#eq-curti-predation-mortality-sparre}

::: {.callout-important title="Ecosystem Biomass Scaling Invariance"}
[@curti2013, p. 478] demonstrated that halving assumed total ecosystem biomass $\text{EB}$ from $15\text{ Mt}$ to $7.34\text{ Mt}$ changed estimated recruitment, abundance, and $M2_{i,a,t}$ by **less than $3.5\%$** ($0.53\%$, $0.59\%$, and $3.49\%$ maximum differences, respectively!), because the estimated vulnerability parameters $\rho_{i,j}$ (defined relative to $\rho_{\text{other}} = 1$) automatically scaled down by $51\text{--}53\%$ to preserve the observed diet proportions $\phi_{i,a,j,b,t}/\phi_{j,b,t}$.
:::


2. **Empirical Pre-Estimation of Size Preference $(\eta_{i,j}, \sigma_{i,j})$ vs. Internal Estimation of $\ln(\rho_{i,j})$**:
   Because NEFSC Food Web Dynamics stomach records routinely measure predator length and prey species mass, but only a subset ($1,900$ records) measured individual prey lengths, Curti et al. found that attempting to estimate $(\eta_{i,j}, \sigma_{i,j})$ inside the integrated model caused severe parameter confounding. They therefore computed $\eta_{i,j}$ and $\sigma_{i,j}$ externally from the $1,900$ paired predator-prey weight records (Table 5, p. 477) and estimated only the log species-preference coefficients $\ln(\rho_{i,j})$ inside ADMB against 5-year-binned multinomial food-habits proportions weighted by McAllister–Ianelli effective sample sizes [@mcallister_ianelli1997].

3. **Regularizing Penalty Functions [@curti2013, Table 3, p. 474]**:
   To prevent biologically unrealistic boundary solutions during optimization and Monte Carlo runs, Curti et al. added three smooth quadratic penalties to the negative log-likelihood:
   - **Initial Cohort Decay Penalty ($Yp_i$)**: Penalizes deviations of first-year abundances $N_{i,a,1}$ from exponential decay $\text{Yr1}_{i,a,1} = N_{i,a-1,1} e^{-\bar{Z}_{i,a-1}}$:
     $$Yp_i = \sum_{a > 1} \left(N_{i,a,1} - \text{Yr1}_{i,a,1}\right)^2$$ {#eq-curti-penalty-first-year}
   - **Log-Recruitment CV Threshold Penalty ($Rp_i$)**: Prevents the coefficient of variation of log-recruitment from exceeding species thresholds $\text{Rth}_i \in \{1.0, 0.7, 0.4\}$:
     $$Rp_i = 0.01 \left[\text{CV}(\ln N_{i,1,t}) - \text{Rth}_i\right]^2 \cdot \mathbb{I}\!\left(\text{CV}(\ln N_{i,1,t}) > \text{Rth}_i\right)$$ {#eq-curti-penalty-recruitment-cv}
   - **Minimum Age-Specific Biomass Threshold Penalty ($Bp_i$)**: Prevents any cohort biomass $B_{i,a,t}$ from collapsing below $\text{Bth}_i = 10\text{ kg}$ during trial steps:
     $$Bp_i = 0.01 \sum_{a,t} \left(B_{i,a,t} - \text{Bth}_i\right)^2 \cdot \mathbb{I}\!\left(B_{i,a,t} < \text{Bth}_i\right)$$ {#eq-curti-penalty-biomass-threshold}

4. **Monte Carlo Error-Propagation Findings [@curti2013, pp. 474--482, Figs. 4--7]**:
   Running $400\text{--}500$ Monte Carlo replicates across progressively higher observation error levels in commercial catch ($\sigma_{\text{TC}} \in \{0.1, 0.2, 0.4, 0.6, 0.8\}$, $N_{\text{eff,CP}} \in \{1000, 500, 75, 40, 20\}$), survey indices ($\sigma_{\text{TS}} \in \{0.2, 0.4, 0.6, 0.8\}$, $N_{\text{eff,SP}} \in \{1000, 200, 100, 60\}$), and predator diets ($N_{\text{eff,FH}} \in \{2500, 500, 100, 50, 10\}$) revealed a clear separation of information sources:
   - **Fishing mortality $F_{i,t}$** is most sensitive to measurement error in **commercial catches**.
   - **Recruitment $N_{i,1,t}$ and predation mortality $M2_{i,a,t}$** are most sensitive to measurement error in **trawl survey indices** (breaking down when survey lognormal $\sigma_{\text{TS}} \ge 0.6$).
   - **Robustness to Noisy Diet Data**: Across all parameters, increasing noise in predator stomach contents (even down to an effective sample size of $N_{\text{eff,FH}} = 10$ per 5-year block) had the **smallest impact** on recruitment, $F$, and $M2$ precision, because the temporal trajectory of $M2_{i,a,t}$ is primarily driven by predator and prey abundance trajectories anchored by the catch and survey datasets!

---

### Alternative Functional Responses in `MSCAA` and `Hydra`
Whereas SMS uses the fixed-ration Andersen–Ursin equation (which corresponds to a Holling Type II functional response in the limit where total consumption rate $R_{A,a}$ is held constant by adjusting "other food"), [@kinzey2009] and Gaichas et al. (`Hydra`) formulated `MSCAA` to allow explicit comparison of **alternative functional response forms** for the per-predator consumption rate $I_{i,a,j,b,t}$ of prey $(i, a)$ by predator $(j, b)$:

Let $\tilde{B}_{i,a,j,b,t} = S_{i,a,j,b} \, \bar{N}_{i,a,t} \, w_{i,a,t}$ be the suitable prey biomass of $(i, a)$ for predator $(j, b)$, and let $\Phi_{j,b,t} = \sum_{k,c} \tilde{B}_{k,c,j,b,t} + S_{\text{OF},j,b} B_{\text{OF}}$ be total suitable food.

::: {.panel-tabset}
## Type I (Linear)
Linear functional response with constant search / clearance rate $\alpha_{j,b}$:
$$M2_{i,a,t} = \sum_{j,b} \alpha_{j,b} \, S_{i,a,j,b} \, \bar{N}_{j,b,t}$$ {#eq-mscaa-type1-predation}

## Type II (Hyperbolic / SMS)
Hyperbolic satiation with half-saturation constant $H_{j,b}$:
$$M2_{i,a,t} = \sum_{j,b} \frac{R_{j,b,t} \, \bar{N}_{j,b,t} \, S_{i,a,j,b}}{\Phi_{j,b,t} + H_{j,b}}$$ {#eq-mscaa-type2-predation}

*(Note: when $H_{j,b} = 0$ and $S_{\text{OF},j,b} B_{\text{OF}} > 0$, this reduces identically to the SMS equation).*

## Type III (Sigmoidal Switching)
Sigmoidal prey-switching with exponent $\theta_{j} > 1$:
$$M2_{i,a,t} = \sum_{j,b} \frac{R_{j,b,t} \, \bar{N}_{j,b,t} \, S_{i,a,j,b} \, \left(\bar{N}_{i,a,t} w_{i,a,t}\right)^{\theta_j - 1}}{\sum_{k,c} S_{k,c,j,b} \left(\bar{N}_{k,c,t} w_{k,c,t}\right)^{\theta_j} + \left(S_{\text{OF},j,b} B_{\text{OF}}\right)^{\theta_j}}$$ {#eq-mscaa-type3-predation}

When $\theta_j = 2$, rare prey experience disproportionately low predation mortality ($M2_{i,a,t} \to 0$ as $B_{i,a,t} \to 0$), providing an automatic stabilizing refuge at low stock size.

## Predator Interference
Beddington–DeAngelis / Crowley–Martin pre-emption with predator interference parameter $\gamma_{j}$:
$$M2_{i,a,t} = \sum_{j,b} \frac{R_{j,b,t} \, \bar{N}_{j,b,t} \, S_{i,a,j,b}}{\Phi_{j,b,t} + H_{j,b}\left(1 + \gamma_j \bar{N}_{j,b,t}\right)}$$ {#eq-mscaa-predator-interference-predation}

This dampens top-down control when predator abundance spikes, mimicking the behavioral foraging-arena limitation of Ecosim [@walters1997].
:::

::: {.callout-tip title="Simulation-Estimation Takeaways: Trijoulet et al. (2019)"}
Using a self-test and cross-test simulation framework in the TMB-based multispecies model `ms-keyrun` / `MSCAA`, [@trijoulet2019] demonstrated that:

- Ignoring trophic interactions (fitting single-species SCAA with constant $M$ to data generated from a multispecies food web) produces severe, systematic biases in SSB ($>30\%$), recruitment ($>50\%$), and $F_{\text{MSY}}$ reference points for both forage fish and cannibalistic predators.
- Even when diet composition data are noisy (effective sample size $N_{\text{stom}} = 20\text{--}50$), `MSCAA` models accurately recover the temporal trajectory of $M2_{i,a,t}$ and eliminate the retrospective bias seen in constant-$M$ single-species assessments.
:::

---

## Climate-Enhanced Multispecies Stock Assessment (`CEATTLE` / `Rceattle`)

Developed by Kirstin Holsman, Jim Ianelli, Kerim Aydin, André Punt, Grant Adams, and colleagues at NOAA AFSC and the University of Washington [@holsman2016; @adams2022; @ices2021; @ices2023], **CEATTLE** (*Climate-Enhanced, Age-based model with Temperature-specific Trophic Linkages and Energetics*) and its open-source TMB/R package **`Rceattle`** represent the state-of-the-art in integrating **temperature-dependent growth, Wisconsin bioenergetics, and multispecies harvest control rules** into an age-structured statistical stock assessment.

Originally built in ADMB for the Eastern Bering Sea complex—walleye pollock (*Gadus chalcogrammus*, ages 1–12), Pacific cod (*Gadus macrocephalus*, ages 1–12), and arrowtooth flounder (*Atheresthes stomias*, ages 1–21, sex-structured)—and subsequently ported to TMB (`Rceattle`) and expanded to the Gulf of Alaska (adding Pacific halibut *Hippoglossus stenolepis*; [@adams2022]), California Current, and Baltic Sea, `CEATTLE` is used operationally by the North Pacific Fishery Management Council (NPFMC).

### Temperature-Dependent Generalized von Bertalanffy Growth (Holsman et al. 2016)
Unlike standard age-structured models that treat weight-at-age $W_{ij,y}$ solely as an empirical look-up table, [@holsman2016, Eqs. 1--2 & Table 1 Eqs. T1.10a--c, pp. 12--14, 50] derived a **temperature-dependent generalized von Bertalanffy growth function (VBGF)** based on [@essington2001] and [@holsman_aydin2015], fitted to $21,388$ pollock, $14,362$ Pacific cod, and $772$ arrowtooth flounder otolith weight-at-age records:

$$W_{ij,y} = W_{\infty, iy} \left(1 - \exp\!\left(-K_i \left(1 - d_{i,y}\right)\left(j - t_{0,i}\right)\right)\right)^{\frac{1}{1 - d_{i,y}}} e^{\varepsilon}, \qquad \varepsilon \sim \mathcal{N}\!\left(0, \sigma_{d,i}^2\right)$$ {#eq-ceattle-von-bertalanffy-weight}

where:
- $W_{\infty, iy} = \left(\frac{H_i}{K_i}\right)^{\frac{1}{1 - d_{i,y}}}$ is the time-varying asymptotic body mass of species $i$ in year $y$, determined by the bioenergetic assimilation constant $H_i$ ($\text{kg}\cdot\text{kg}^{-d}\cdot\text{yr}^{-1}$; e.g., $16.34$ for pollock, $9.30$ for cod) and catabolic energy loss constant $K_i$ ($\text{kg}\cdot\text{kg}^{-1}\cdot\text{yr}^{-1}$; $0.22$ for pollock, $0.45$ for cod; Table 3, p. 54),
- $d_{i,y}$ is the **realized allometric slope of consumption**, which varies annually as an exponential-linear function of bottom temperature $T_y$ ($^\circ\text{C}$):
   $$d_{i,y} = \exp\!\left(\alpha_{d,i,y} + \alpha 0_{d,i} + \beta_{d,i} T_y\right)$$ {#eq-ceattle-metabolic-exponent}
   allowing both the growth rate and asymptotic weight $W_{\infty, iy}$ to respond dynamically to historical and IPCC-projected bottom temperatures!
- For sexually dimorphic species (arrowtooth flounder, where males have higher natural mortality $M_{\text{male}} = 0.35$ vs. $M_{\text{fem}} = 0.20$ and smaller asymptotic size), sex-specific weights $W^{\text{fem}}_{ij,y}, W^{\text{male}}_{ij,y}$ and residual mortalities are combined using the surviving female cohort proportion [@holsman2016, Eqs. T1.18--T1.21, p. 51]:
   $$\omega_{ij} = \frac{e^{-j M_{\text{fem}}}}{e^{-j M_{\text{fem}}} + e^{-j M_{\text{male}}}}, \qquad W_{ij,y} = \omega_{ij} W^{\text{fem}}_{ij,y} + (1 - \omega_{ij}) W^{\text{male}}_{ij,y}$$ {#eq-ceattle-sex-ratio-weight}

### Wisconsin Bioenergetics Consumption and Thermal Scaling ($f(T_y)_p$)
In `CEATTLE`, the individual specific annual ration $\delta_{pa,y}$ ($\text{kg prey} \cdot \text{kg predator}^{-1} \cdot \text{yr}^{-1}$) of predator species $p$, age $a$ in year $y$ is modeled via the Wisconsin Bioenergetics allometric equation [@holsman2016, Table 2, Eqs. T2.4--T2.5, p. 52; @adams2022, Eqs. 10--11]:

$$\delta_{pa,y} = \hat{\varphi}_p \, \alpha_\delta \, W_{pa,y}^{(1 + \beta_\delta)} \, f(T_y)_p$$ {#eq-ceattle-annual-ration}

where:
- $\hat{\varphi}_p$ is the annual relative foraging rate ($\text{d}\cdot\text{yr}^{-1}$ scalar),
- $\alpha_\delta$ and $\beta_\delta < 0$ are the intercept and allometric slope of maximum specific consumption ($\beta_\delta = -0.460$ for pollock, $-0.122$ for cod, $-0.245$ for arrowtooth flounder; Table 3, p. 53),
- $f(T_y)_p$ is the [@kitchell1977] **dome-shaped temperature scaling algorithm**:
   $$f(T_y)_p = V^X \exp\!\left(X (1 - V)\right)$$ {#eq-ceattle-temperature-effect-thornton}
   with:
   $$V = \frac{T^{\text{cm}}_p - T_y}{T^{\text{cm}}_p - T^{\text{co}}_p}, \qquad X = \frac{Z^2 \left(1 + \left(1 + 40/Y\right)^{0.5}\right)^2}{400}$$ {#eq-ceattle-thornton-temperature-shape}
   $$Z = \ln(Q^c_p)\left(T^{\text{cm}}_p - T^{\text{co}}_p\right), \qquad Y = \ln(Q^c_p)\left(T^{\text{cm}}_p - T^{\text{co}}_p + 2\right)$$ {#eq-ceattle-thornton-temperature-constants}
   where $T^{\text{co}}_p$ is the physiological optimum consumption temperature ($10.0^\circ\text{C}$ for pollock, $13.7^\circ\text{C}$ for cod, $19.6^\circ\text{C}$ for arrowtooth), $T^{\text{cm}}_p$ is the maximum lethal physiological temperature ($15.0^\circ\text{C}$, $21.0^\circ\text{C}$, $34.13^\circ\text{C}$), and $Q^c_p$ is the $Q_{10}$ thermal slope ($2.60, 2.41, 2.18$). (`Rceattle` also supports the [@thornton_lessem1978] two-logistic thermal curve and energy-density caloric scaling $\mathcal{E}_{\text{pred}} / \sum \hat{W} \mathcal{E}_{\text{prey}}$).

### Predation Mortality ($M2_{ij,y}$) and Analytical/Iterative Suitability ($\bar{S}_{paij}$)
Following [@jurado_molina2005], `CEATTLE` expresses predation mortality $M2_{ij,y}$ on prey species $i$, age $j$ in year $y$ as a Holling Type II functional response where suitability $\bar{S}_{paij} \in [0, 1]$ of modeled prey $(i, j)$ and suitability $1 - \sum_{ij} \bar{S}_{paij}$ of "other prey" $B^{\text{other}}_p$ sum to unity [@holsman2016, Table 2, Eqs. T2.1--T2.3, p. 52]:

$$M2_{ij,y} = \sum_{p} \sum_{a} \left( \frac{N_{pa,y} \, \delta_{pa,y} \, \bar{S}_{paij}}{\sum_{i,j} \left(\bar{S}_{paij} B_{ij,y}\right) + B^{\text{other}}_p \left(1 - \sum_{i,j} \bar{S}_{paij}\right)} \right)$$ {#eq-ceattle-m2-predation}

Instead of adding thousands of noisy diet cells to the likelihood, [@holsman2016] and [@adams2022] calculate $\bar{S}_{paij}$ analytically from observed gravimetric stomach proportions $U_{paij,y}$ (converted from length-based stomach samples to age-based proportions using annual empirical age-length transition matrices $\text{ALK}_y$) and model-estimated biomasses $B_{ij,y}$:

$$\bar{S}_{paij} = \frac{1}{n_y} \sum_{y} \left( \frac{\dfrac{\bar{U}_{paij}}{B_{ij,y}}}{\sum_{i,j} \left(\dfrac{\bar{U}_{paij}}{B_{ij,y}}\right) + \dfrac{1 - \sum_{i,j} \bar{U}_{paij}}{B^{\text{other}}_p}} \right)$$ {#eq-ceattle-mean-suitability}

In `Rceattle` [@adams2022, Section 2.4, p. 4], the single-species assessment models are first optimized (tuning McAllister–Ianelli effective sample sizes for age/length compositions; [@mcallister_ianelli1997]), and then the multispecies model is iteratively re-optimized, updating $\bar{S}_{paij}$ at each outer step until convergence (which occurs in **3 iterations**).

### Multispecies Biological Reference Points (MBRPs) and Harvest Control Rules (Holsman et al. 2016)
A major contribution of [@holsman2016, Section 2.3 & Table 5, pp. 15--20, 57] is the rigorous mathematical formalization of **12 Multispecies Biological Reference Point (MBRP) harvest scenarios** under NPFMC Tier 3 control rules (where target Acceptable Biological Catch proxy biomass is $\overline{\text{SSB}}_{F,i} = 0.40 \, \overline{\text{SSB}}_{0,i}$ and the minimum MSY biomass limit is $B_{\text{MSY}, i} = 0.35 \, \overline{\text{SSB}}_{0,i}$):

1. **Non-Uniqueness of "Unfished Biomass" ($\overline{\text{SSB}}_{0,i}$) and Fishing-Order Dependence (Scenarios 1.1–1.4)**:
   - Under **Scenario 1.1**, $\overline{\text{SSB}}_{0,i}$ is defined by setting $F_i = 0$ for **all species simultaneously**:
     $$\overline{\text{SSB}}_{F,i} = 0.40 \, \overline{\text{SSB}}_{0,i} \Big|_{F_1 = F_2 = F_3 = 0}$$ {#eq-ceattle-multispecies-biological-reference-point}
   - Under **Scenarios 1.3 & 1.4**, $\overline{\text{SSB}}'_{0,1}$ for prey (pollock) is defined **conditionally on predators (cod and arrowtooth) being fished at their target $F^*_{\text{ABC}, 2}, F^*_{\text{ABC}, 3}$**:
     $$\overline{\text{SSB}}'_{0,1} = \left(\overline{\text{SSB}}_{0,1} \;\middle|\; F_1 = 0,\, \overline{\text{SSB}}_{F,2},\, \overline{\text{SSB}}_{F,3}\right)$$ {#eq-ceattle-dynamic-ssb0-conditional}
   - *Mathematical Insight*: When predators are unfished ($F_2 = F_3 = 0$), predator biomass expands and inflicts massive $M2_{1,a}$ on juvenile pollock (plus intense pollock cannibalism, which triggers **over-compensatory damped oscillations** in projected SSB; Fig. 10, p. 60). Consequently, simultaneous unfished pollock biomass $\overline{\text{SSB}}_{0,1}$ can actually be **lower** than harvested pollock biomass when predators are fished!
2. **Spawning Biomass-Per-Recruit Approach (Scenario 2.1, Eqs. T5.2.1a–b)**:
   Solves for $F_{x,i}$ such that $\text{SPR}_{F,i} = 0.40 \left(\text{SPR}_{F,i} \mid F_{x,i} = 0\right)$, where:
   $$\text{SPR}_{F,i} = \frac{1}{R_{i,1}} \sum_{j=1}^{A_i} \rho_{ij} \, \bar{W}_{ij} \, N_{ij} \exp\!\left(-\left(F_{x,i} s^f_{ij} + M1_{ij} + \overline{M2}_{ij}\right)\right)$$ {#eq-ceattle-multispecies-spr}
3. **System-Wide Aggregate Biomass Proxies (Scenarios 3.1–3.4, Eqs. T5.3.1–T5.3.4)**:
   Finds a system-wide multiplier $F_{\text{ALL}}$ such that total community female spawning biomass satisfies $\sum_i \overline{\text{SSB}}_{F,i} = 0.40 \sum_i \overline{\text{SSB}}_{0,i}$, scaling species-specific fishing mortality either by total natural mortality $F_{x,i} = F_{\text{ALL}} (\bar{M}_i / \sum_k \bar{M}_k)$ or by unfished biomass $F_{x,i} = F_{\text{ALL}} (\overline{\text{SSB}}_{0,i} / \sum_k \overline{\text{SSB}}_{0,k})$, with (Scenarios 3.2, 3.4) or without (Scenarios 3.1, 3.3) the conservation constraint $\text{SSB}_{F,i,y} \ge 0.35 \, \text{SSB}_{0,i}$.
4. **Simplex Optimization of Individual vs. Aggregate MSY (Scenarios 4.1–4.3, Eqs. T5.4.1–T5.4.3)**:
   Uses 8 iterations of Nelder–Mead simplex optimization to find $\mathbf{F}_x = (F_{x,1}, F_{x,2}, F_{x,3})$ maximizing either individual species yield $\max(\overline{\text{ABC}}_{x,i} \mid F_{x,i})$ (Scenario 4.1) or total combined community yield $\max\left(\sum_i \overline{\text{ABC}}_{x,i} \mid \mathbf{F}_x\right)$ without (Scenario 4.2) or with (Scenario 4.3) the minimum biomass threshold $\text{SSB}_{F,i,y} \ge 0.35 \, \text{SSB}_{0,i}$:
   - *Predator Eradication Paradox*: Unconstrained maximization of aggregate yield (Scenario 4.2) achieved $3.2\text{ million t}$ by **fishing Pacific cod and arrowtooth flounder to extinction** in order to eliminate predation on pollock (since annual predation on juvenile pollock is $2\text{--}5\times$ the annual fishery harvest!). Imposing the biodiversity constraint $\text{SSB}_{F,i} \ge 0.35 \, \text{SSB}_{0,i}$ (Scenario 4.3) prevented predator collapse and demonstrated that **MBRPs are not inherently more or less conservative than single-species BRPs**; rather, their conservatism depends on whether single-species or aggregate targets are used and whether minimum biomass thresholds are enforced.

---

### Gulf of Alaska `Rceattle` Ensemble and Hybrid "Predators-as-Fleets" Scaling (Adams et al. 2022)
In adapting `Rceattle` (implemented in TMB with recruitment and initial cohort deviations treated as Laplace-integrated random effects) to the Gulf of Alaska for pollock, Pacific cod, arrowtooth flounder, and Pacific halibut (*Hippoglossus stenolepis*), [@adams2022] addressed a universal practical problem in multispecies assessment: **how to include an important apex predator (Pacific halibut) whose stock assessment is conducted externally on a coastwide spatial scale** (by the International Pacific Halibut Commission, IPHC, using a 4-model Stock Synthesis ensemble).

[@adams2022] evaluated 16 model configurations across 3 time-series durations ($1977\text{--}2018$, $1993\text{--}2018$, $1996\text{--}2018$) grouped into 4 ensemble averages:
1. **Ensemble Avg 1 (Single-Species)** vs. **Ensemble Avg 2 (3-Species MS without Halibut)**: Including predation by arrowtooth flounder (the dominant predator in the GOA), cod, and pollock increased average estimated pollock total biomass by **$2.03\text{--}2.21\times$** and pollock SSB by **$1.70\text{--}1.84\times$** relative to single-species models.
2. **Ensemble Avg 3 (Fixed External Predator Abundance)**: Inputs sex-specific mid-year halibut numbers-at-age directly from each of the 4 IPHC Stock Synthesis ensemble models (scaled by the proportion of age-5+ halibut survey biomass in IPHC Regulatory Areas 3A & 3B) without fitting halibut data in the `Rceattle` likelihood.
3. **Ensemble Avg 4 (Relative Predator Survey Index with Estimated Scaling $\tilde{q}_i$)**: Inputs relative halibut abundance-at-age from the IPHC setline survey multiplied by a freely estimated scaling parameter $\tilde{q}_{\text{halibut}}$.
   - *Critical Identifiability Warning* [@adams2022, Sections 3.1 & 4.2, pp. 6, 11]: When $\tilde{q}_{\text{halibut}}$ was freely estimated inside the multispecies likelihood (Models 11 & 16), the optimizer drove $\tilde{q}_{\text{halibut}} \to 0$ ($0.0021\text{--}0.0025$ with large standard errors), **shrinking halibut predation mortality to zero** and reproducing the 3-species no-halibut model! Even though Models 11 & 16 achieved a slightly lower AIC on the pollock/cod/arrowtooth survey and catch data by removing the external constraint of halibut predation, doing so was biologically spurious—a textbook demonstration of why predator abundance indices in "predators-as-fleets" extensions cannot have unanchored scaling parameters $\tilde{q}_i$ without informative priors on predator scale.
