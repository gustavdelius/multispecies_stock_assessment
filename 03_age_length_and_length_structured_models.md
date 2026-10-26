# Age-Length and Length-Structured Assessment Models (`Gadget`/`gadget3`, `SS3`, `MULTIFAN-CL`, `LeMans`/`LeMaRns`, `CSMS`/`CSM`) {#sec-length-structured}

## Why Track Length in Addition to (or Instead of) Age?

In marine ecosystems, nearly every ecological and fishery interaction—gape-limited predation, trawl mesh selection, gillnet dome-shaped selectivity, metabolic rate, and fecundity—is a direct physical function of **body size (length $L$ or weight $W$)** rather than chronological age $a$ [@stefansson1998; @punt2016]. Purely age-structured models (`MSVPA`, `SMS`, `CEATTLE`) suffer from three structural limitations when size-dependent processes are strong:
1. **Rosa Lee's Phenomenon (Size-Selective Mortality Within an Age Class)**: When fishing or predation selectively removes the fastest-growing (or slowest-growing) individuals within an age cohort, the surviving individuals of age $a$ have a distorted length and weight distribution that cannot be represented by a fixed empirical weight-at-age $w_{s,a,t}$ [@methot2013; @begley2004overview].
2. **Consumption-Driven Growth Plasticity**: When forage fish collapse, predators grow more slowly in length and lose body condition (weight-at-length). A joint age-length-weight model can track both stunted length growth and starvation-induced weight loss dynamically [@stefansson1998].
3. **Direct Observation of Length Frequencies**: Measuring fish lengths at sea or port is fast and cheap, whereas reading otolith ages is expensive, prone to reader bias, and unavailable for crustaceans (shrimp, Nephrops, crabs), cephalopods, and many tropical or data-limited stocks [@fournier1998; @punt2016].

---

## `Gadget` (`Gadget2` and `gadget3`): Statistical Age-Length-Structured Multispecies Framework

### Lineage and Architecture
**Gadget** (*Globally Applicable Area Disaggregated General Ecosystem Toolbox*) evolved from two Nordic multispecies research models of the late 1980s and 1990s [@stefansson1998; @begley2004overview], synthesized under the EU `dst2` project [@stefansson2003; @stefansson2004dst2; @begley2004guide; @begley2004overview]:
- **MULTSPEC** [@bogstad1997, IMR Norway]: An area-, age-, and length-structured model of Barents Sea cod, capelin, herring, harp seals, and minke whales.
- **BORMICON** (*Boreal Migration and Consumption model*; [@stefansson_palsson1997], MRI Iceland): A statistical multi-area age-length assessment model of Icelandic cod, capelin, and northern shrimp (*Pandalus borealis*) [@taylor2004shrimp; @taylor2004cod; @taylor2007].

Originally written in C++ (`Gadget2`, [@begley2004guide]) and optimized via Hooke & Jeeves pattern search, Simulated Annealing, and BFGS, with applications ranging from Icelandic cod–capelin–shrimp interactions [@taylor2004shrimp; @taylor2004cod; @taylor2007] and Mediterranean hake [@bartolino2011] to Eastern Baltic cod ontogenetic diet shifts [@kulatska2019], the framework has recently been re-engineered as **`gadget3` (`g3`)** [@lentz_gadget3; @ices2021], which compiles modular R model definitions (`g3a_grow_*`, `g3a_predate_*`, `g3l_distribution_*`) directly into **TMB C++ templates with exact automatic differentiation**.

### State-Space Tensor ($N_{s,a,l,r,t}$ and $W_{s,a,l,r,t}$)
Gadget tracks sub-stocks $s$ (which may represent species, sex, or immature vs. mature life stages, e.g., immature hake vs. mature hake in [@bartolino2011]) across discrete time steps $t$ (typically monthly or quarterly, $\Delta t = 1/12$ or $1/4\text{ yr}$), spatial areas $r$, discrete ages $a \in \{a_{\min}, \dots, a_{\max}\}$, and discrete length bins $l \in \{1, \dots, n_L\}$ with bin boundaries $[L_{\min, l}, L_{\max, l})$ and midpoint $L_l$.

Crucially, in every $(s, a, l, r, t)$ cell, Gadget tracks **two coupled state variables** [@begley2004overview, p. 4; @lentz_gadget3]:
1. $N_{s,a,l,r,t}$: The **number of individuals** in cell $(s, a, l, r)$ at time $t$.
2. $W_{s,a,l,r,t}$: The **mean body weight** of individuals in cell $(s, a, l, r)$ at time $t$ (which can either follow a fixed length-weight curve $W_l = \alpha_s L_l^{\beta_s}$ or vary dynamically in response to food intake and growth).

Within each time step $t$, biological and fishery processes are executed sequentially in a strict order:
$$\text{1. Migration } (\mathbf{M}) \longrightarrow \text{2. Predation \& Fishing } (\mathbf{P}) \longrightarrow \text{3. Somatic Growth \& Maturation } (\mathbf{G}) \longrightarrow \text{4. Ageing } (\mathbf{A}) \longrightarrow \text{5. Spawning \& Recruitment } (\mathbf{R})$$ {#eq-gadget-substep-sequence}

---

### Somatic Growth and the Beta-Binomial Length-Transition Distribution
In a discrete length-bin model, somatic growth during time step $\Delta t$ must specify:
1. The **mean length increment** $\Delta L_l$ and **mean weight increment** $\Delta W_l$ of fish currently in length bin $l$, and
2. The **probability distribution** $p(\Delta j \mid l)$ that a fish in length bin $l$ jumps up by $\Delta j \in \{0, 1, \dots, \Delta j_{\max}\}$ length bins during $\Delta t$.

#### A. Mean Length & Weight Increments (`g3a_grow_*`)
Gadget provides two classes of growth equations [@begley2004overview, pp. 8--10; @lentz_gadget3]:

1. **Consumption-Independent von Bertalanffy Growth (`g3a_grow_lengthvbsimple` / `g3a_grow_weightsimple`)**:
   $$\Delta L_l = \left(L_{\infty, s} - L_l\right)\left(1 - \exp(-\kappa_s \Delta t)\right)$$ {#eq-gadget-growth-length-vb}
   $$\Delta W_l = \alpha_s \left(L_l + \Delta L_l\right)^{\beta_s} - \alpha_s L_l^{\beta_s}$$ {#eq-gadget-growth-weight-allometric}

2. **Consumption-Dependent (Two-Way Bioenergetic) Growth (`g3a_grow_length_multspec` and `g3a_grow_weight_weightjones`)**:
   Let $\psi_{s,l,r,t} \in [0, 1]$ be the predator's realized **feeding level** in length bin $l$ during time step $t$ (ratio of actual biomass consumed to maximum possible consumption $M_{s,l,t}$; see @sec-gadget-predation), and let $T_{r,t}$ be ambient temperature.
   - **MULTSPEC Length & Weight Increment (`g3a_grow_length_multspec` / `g3a_grow_weight_multspec`)**:
     $$\Delta L_l = \Delta t \cdot p_{0,s} \, L_l^{p_{1,s}} \, \psi_{s,l,r,t} \left(p_{2,s} T_{r,t} + p_{3,s}\right)$$ {#eq-gadget-growth-length-jones}
     $$\Delta W_l = \Delta t \cdot q_{0,s} \, W_{s,a,l}^{q_{1,s}} \left(\psi_{s,l,r,t} - q_{2,s}\right)\left(q_{3,s} T_{r,t} + q_{4,s}\right)$$ {#eq-gadget-growth-weight-jones}
     where $q_{2,s}$ is the maintenance feeding level below which fish **lose weight** ($\Delta W_l < 0$, while $\Delta L_l \ge 0$ since fish skeletons do not shrink!).
   - **Jones (1976) Bioenergetic Weight Increment (`g3a_grow_weight_weightjones`, [@jones1976])**:
     $$\Delta W_l = \frac{\Delta t \left( \frac{C_{s,l,t}}{\Delta t} - r_{0,s} W_{s,a,l}^{r_{1,s}} \exp(r_{2,s} T_{r,t}) \right)}{r_{3,s}}$$ {#eq-gadget-growth-dynamic-energy}
     coupled to length growth via a condition-factor feedback function $f(\delta W_l)$ that suppresses length growth $\Delta L_l$ whenever actual body weight $W_{s,a,l}$ falls below reference weight $W_{\text{ref}, l} = \alpha_s L_l^{\beta_s}$ [@begley2004overview, p. 9; @lentz_gadget3].

#### B. The Beta-Binomial Length-Bin Jump Distribution (`g3a_grow_impl_bbinom`)
Once the mean length increment $\Delta L_l$ is computed (and normalized by the bin width $\delta L$ to give mean bin growth $\mu_l = \Delta L_l / \delta L$), Gadget distributes growing individuals across integer bin jumps $\Delta j \in \{0, 1, \dots, n_{\max}\}$ using a **Beta-Binomial distribution** [@stefansson1998; @taylor2004cod, Eq. 5; @lentz_gadget3]:

Given:
- $n_{\max}$: Maximum allowed number of length-bin jumps in a single time step $\Delta t$,
- $\mu_l = \frac{\Delta L_l}{\delta L} \in [0, n_{\max}]$: Expected number of bin jumps,
- $\beta_s > 0$: Single free dispersion parameter controlling variance in individual growth rates,

define $\alpha_l$ via the mean constraint $\mu_l = \frac{n_{\max} \alpha_l}{\alpha_l + \beta_s}$:
$$\alpha_l = \frac{\beta_s \, \mu_l}{n_{\max} - \mu_l}$$ {#eq-gadget-beta-binomial-alpha}
Then the probability $G_{l, l + \Delta j} = \Pr(\text{Jump} = \Delta j \mid l)$ of growing by $\Delta j \in \{0, 1, \dots, n_{\max}\}$ bins in $\Delta t$ is:

$$G_{l, l + \Delta j} = \binom{n_{\max}}{\Delta j} \frac{\text{B}(\Delta j + \alpha_l,\; n_{\max} - \Delta j + \beta_s)}{\text{B}(\alpha_l, \beta_s)} = \frac{\Gamma(n_{\max} + 1)}{\Gamma(\Delta j + 1)\Gamma(n_{\max} - \Delta j + 1)} \frac{\Gamma(\alpha_l + \beta_s)}{\Gamma(n_{\max} + \alpha_l + \beta_s)} \frac{\Gamma(\Delta j + \alpha_l)\Gamma(n_{\max} - \Delta j + \beta_s)}{\Gamma(\alpha_l)\Gamma(\beta_s)}$$ {#eq-gadget-beta-binomial-transition}

In `gadget3` [@lentz_gadget3], this is evaluated efficiently via an exact forward recurrence relation in $\Delta j$ so that no `lgamma` calls are needed inside the inner loop, and both cohort numbers $N_{s,a,l'}$ and cell mean weights $W_{s,a,l'}$ are updated cleanly via matrix multiplication:
$$N_{s,a,l',t+\Delta t} = \sum_{\Delta j = 0}^{n_{\max}} G_{l' - \Delta j, l'} \, N_{s,a,l'-\Delta j, t}, \qquad (NW)_{s,a,l',t+\Delta t} = \sum_{\Delta j = 0}^{n_{\max}} G_{l' - \Delta j, l'} \, N_{s,a,l'-\Delta j, t} \left(W_{s,a,l'-\Delta j, t} + \Delta W_{l'-\Delta j, \Delta j}\right)$$ {#eq-gadget-length-bin-update}

> [!TIP]
> **Why the Beta-Binomial Growth Implementation is a Major Lesson for `mizer`**:
> With **only one extra parameter ($\beta_s$)** beyond the deterministic growth rate $\Delta L_l$, `g3a_grow_impl_bbinom` generates realistic, strictly non-negative, bounded discrete-bin dispersion around the mean growth trajectory, completely eliminating the numerical vs. biological diffusion problem of deterministic first-order PDE schemes while remaining $C^\infty$-differentiable in TMB!

---

### Predation, Functional Response, and Overconsumption Protection (`g3a_predate`) {#sec-gadget-predation}
In Gadget, both **commercial/survey fishing fleets** and **biological predators** are unified under a single mathematical abstraction (`g3a_predate` in [@lentz_gadget3]).

For a biological predator $p$ in length group $L$ preying on prey stock $i$ in length group $l$ [@begley2004overview, pp. 13--15; @taylor2004shrimp, Eqs. 7--11; @lentz_gadget3]:

1. **Predator-Prey Length Suitability ($S_{p,i}(L, l)$)**:
   Gadget offers several parametric suitability functions of predator length $L$ and prey length $l$, most commonly:
   - **Andersen & Ursin Log-Ratio Suitability (`g3_suitability_andersen`)**:
     $$S_{p,i}(L, l) = p_0 + p_2 \exp\!\left(-\frac{\left(\ln\frac{L}{l} - p_1\right)^2}{p_3 \text{ (or } p_4\text{)}}\right)$$ {#eq-gadget-suitability-lognormal}
   - **Richards / Asymmetric Logistic Suitability (`g3_suitability_richards`)**:
     $$S_{p,i}(L, l) = \left(\frac{p_3}{1 + \exp(-p_0 - p_1 l - p_2 L)}\right)^{1/p_4}$$ {#eq-gadget-suitability-sigmoidal}

2. **Maximum Possible Consumption ($M_{p,L,t}$, `g3a_predate_maxconsumption`)**:
   The maximum biomass a single predator of length $L$ can consume in time step $\Delta t$ at temperature $T_t$ is:
   $$M_{p,L,t} = \Delta t \cdot m_{0,p} \, L^{m_{1,p}} \exp\!\left(m_{2,p} T_t - m_{3,p} T_t^3\right)$$ {#eq-gadget-max-consumption}

3. **Desired Consumption and Holling Type II/III Functional Response (`g3a_predate_catchability_predator`)**:
   Let $B_{i,l,t} = N_{i,l,t} W_{i,l,t}$ be the available prey biomass of $(i, l)$, and let the suitability-weighted prey energy/biomass be:
   $$E_{p,L,i,l,t} = \left(S_{p,i}(L, l) \, \mathcal{E}_i \, N_{i,l,t} \, W_{i,l,t}\right)^{d_{i}}$$ {#eq-gadget-suitable-biomass}
   where $d_i = 1$ yields a **Holling Type II** functional response and $d_i > 1$ yields a **Holling Type III (prey-switching)** functional response.
   The predator's **feeding level** $\psi_{p,L,t} \in [0, 1]$ is given by the half-saturation parameter $H_p$:
   $$\psi_{p,L,t} = \frac{\sum_{i,l} E_{p,L,i,l,t}}{\sum_{i,l} E_{p,L,i,l,t} + H_p \, \Delta t + E_{\text{other}, p}}$$ {#eq-gadget-feeding-level}
   and the total desired consumption $C^{\text{des}}_{p,L,i,l,t}$ of prey $(i, l)$ by all $N_{p,L,t}$ predators in length group $L$ is:
   $$C^{\text{des}}_{p,L,i,l,t} = N_{p,L,t} \, M_{p,L,t} \, \psi_{p,L,t} \, \frac{E_{p,L,i,l,t}}{\sum_{k,m} E_{p,L,k,m,t} + E_{\text{other}, p}}$$ {#eq-gadget-desired-consumption}

4. **Differentiable Overconsumption Cap and Understocking Penalty**:
   During optimization, parameter trial values may cause total desired consumption across all predators and fleets $D_{i,l,t} = \sum_{p,L} C^{\text{des}}_{p,L,i,l,t} + \sum_{f} C^{\text{des}}_{f,i,l,t}$ to exceed the available prey biomass $B_{i,l,t} = N_{i,l,t} W_{i,l,t}$.
   In `gadget3` [@lentz_gadget3], total consumption is smoothly capped at $95\%$ of available biomass ($0.95 B_{i,l,t}$) using a differentiable logistic minimum function (`dif_pmin`), scaling each predator's actual consumption down proportionally:
   $$C^{\text{act}}_{p,L,i,l,t} = C^{\text{des}}_{p,L,i,l,t} \cdot \text{dif\_pmin}\!\left(1,\; \frac{0.95 \, B_{i,l,t}}{D_{i,l,t}}\right)$$ {#eq-gadget-actual-consumption}
   Simultaneously, any "unmet demand" (understocking $U_{i,l,t} = \max(0, D_{i,l,t} - 0.95 B_{i,l,t})$) adds a smooth quadratic penalty $\ell_{\text{understock}} = w_U \sum_{i,l,t} U_{i,l,t}^p$ to the objective function, guiding the gradient optimizer away from biologically infeasible regions where prey stocks are driven negative.

---

### Composite Likelihoods and Two-Stage Iterative Reweighting in `Gadget`
Gadget fits simultaneously to arbitrary combinations of raw, unaggregated observations [@begley2004overview, pp. 18--23; @lentz_gadget3]:
1. **Survey Abundance Indices (`g3l_distribution_surveyindices_log` / `linear`)**:
   $$\ell_{\text{index}} = \sum_{t} \left(\ln I_{t} - \left(\alpha + \beta \ln \hat{I}_{t}\right)\right)^2$$ {#eq-gadget-survey-loglik}
   where $\alpha = \ln q$ is log catchability and $\beta$ is either fixed at $1$ (linear proportionality) or estimated (hyperstability/hyperdepletion).
2. **Length-Frequency & Age-Length Compositions (`g3l_distribution_sumofsquares`, `multinomial`, `mvlogistic`)**:
   - *Sum of Squares of Proportions*: $\ell_{\text{SS}} = \sum_{t,a,l} \left(\pi_{t,a,l} - \hat{\pi}_{t,a,l}\right)^2$, where $\pi_{t,a,l} = N^{\text{obs}}_{t,a,l} / \sum_{a',l'} N^{\text{obs}}_{t,a',l'}$.
   - *Multinomial*: $\ell_{\text{mult}} = 2 \sum_{t} \left[\ln\Gamma(N_t + 1) - \sum_{a,l} \ln\Gamma(N^{\text{obs}}_{t,a,l} + 1) + \sum_{a,l} N^{\text{obs}}_{t,a,l} \ln \hat{\pi}_{t,a,l}\right]$.
3. **Stomach Content Proportions (`StomachContent` / `g3l_distribution_sumofsquares`)**:
   Compares observed prey species/length weight ratios $W^{\text{obs}}_{p,L,i,l,t}$ in predator $(p, L)$ stomachs against model-predicted consumption fractions:
   $$\hat{W}_{p,L,i,l,t} = \frac{C^{\text{act}}_{p,L,i,l,t}}{\sum_{k,m} C^{\text{act}}_{p,L,k,m,t}}$$ {#eq-gadget-expected-diet-proportion}

#### Stefánsson's (1998) Two-Stage Iterative Likelihood Reweighting
Because a typical Gadget model combines $10\text{--}30$ disparate likelihood components $\ell_c$ ($c = 1, \dots, C$) with different units and sample sizes, [@stefansson1998], [@stefansson2003], and [@taylor2004cod] developed an objective two-stage reweighting procedure (`gadget.iterative`):
1. **Stage 1 (Independent Component Minima)**: For each likelihood group $g$, run an optimization assigning a very high weight to group $g$ to determine the minimum achievable residual sum of squares $\ell_g^*$ (representing the intrinsic observation variance $\hat{\sigma}_g^2 = \ell_g^* / \text{df}_g$ of dataset $g$ without model-conflict bias).
2. **Stage 2 (Inverse-Variance Weighting)**: Set the weight $w_g$ of each likelihood component equal to its inverse estimated variance:
   $$w_g = \frac{\text{df}_g}{\ell_g^*}$$ {#eq-gadget-iterative-reweighting}
   so that in the final joint optimization $\mathcal{L}_{\text{tot}} = \sum_{g} w_g \ell_g(\boldsymbol{\theta})$, every well-specified likelihood component has an expected weighted score of $\approx \text{df}_g$ (chi-square degrees of freedom).

---

### Case Study: Baltic Sea Multispecies `Gadget` Model and Ontogenetic Diet Shifts (Kulatska et al. 2019)
[@kulatska2019] constructed a quarterly (1974–2013) age-length multispecies `Gadget` model (using `Rgadget` and the Marine Fisheries Database `MFDB`) for the Baltic Sea (ICES Subdivisions 25–32), providing a template for how benthic and pelagic prey functional groups are parameterized in a length-structured assessment:

1. **Dynamic Fish Stocks vs. Constant-Biomass Invertebrate Pools**:
   Eastern Baltic cod (*Gadus morhua*, split into juvenile substock ages 0–1, $4\text{--}35\text{ cm}$ and adult substock ages 2–12+, $10\text{--}120\text{ cm}$ at $1\text{ cm}$ resolution) is modeled as the active predator on two dynamic pelagic clupeid stocks—Central Baltic herring (*Clupea harengus*, ages 0–10+, $3.5\text{--}32\text{ cm}$) and Baltic sprat (*Sprattus sprattus*, ages 0–8+, $3\text{--}17\text{ cm}$)—plus three constant-biomass benthic/invertebrate prey pools: the benthic isopod *Saduria entomon* ($0.5\text{--}9\text{ cm}$), pelagic/hyperbenthic mysids (*Mysis mixta*, $0.5\text{--}3\text{ cm}$), and "other food".

2. **Empirical Identification of Ontogenetic Predator Size Classes via Bray–Curtis Clustering**:
   Rather than choosing arbitrary predator length bins for the stomach likelihood (`StomachContent`), [@kulatska2019, Section 2.3.1, p. 341] applied hierarchical cluster analysis with **Bray–Curtis dissimilarity** ($\text{BC} \ge 0.2$ threshold) to $1974\text{--}2013$ EU/ICES Baltic cod stomach records (Quarters 1 and 4, using only digestion stages 1–2 and excluding stage 0 intact prey to eliminate trawl-net feeding!), identifying **five ontogenetic cod diet regimes**:
   $$L_{\text{cod}} \in \big\{[16, 25),\; [25, 35),\; [35, 60),\; [60, 80),\; [80, 112]\text{ cm}\big\}$$ {#eq-kulatska-cod-length-bins}
   reflecting the ontogenetic transition from mysids and *Saduria* ($<25\text{ cm}$) to sprat ($25\text{--}60\text{ cm}$) and large herring ($>60\text{ cm}$).

3. **Prey-Specific Suitability Functions ($S(L, l)$) and Evacuation Ration ($M_L$)**:
   - **Size-Structured Prey (Herring, Sprat, and *Saduria entomon*)**: Modeled via the symmetric Andersen–Ursin logarithmic length-ratio suitability curve [@kulatska2019, Eq. 4, p. 340]:
     $$S_{\text{prey}}(L, l) = p_0 + p_2 \exp\!\left(-\frac{\left(\ln\frac{L}{l} - p_1\right)^2}{p_3}\right)$$ {#eq-kulatska-fish-suitability}
     yielding optimal prey lengths $l^* = L e^{-p_1}$ equal to **$19\%$ of cod length for sprat** and **$28\%$ of cod length for herring**.
   - **Mysids (*Mysis mixta*)**: Because mysids span a narrow size range ($0.5\text{--}3\text{ cm}$) eaten primarily by juvenile cod, suitability is parameterized as a decreasing exponential function of predator length $L$ alone:
     $$S_{\text{mysid}}(L) = \frac{\gamma}{1 + \exp(-\alpha - \beta L)}$$ {#eq-kulatska-mysid-suitability}
   - **Maximum Consumption Rate ($M_L$)**: Estimated externally from observed stomach content weights $S_k$ (including empty stomachs) via the [@jones1978] evacuation rate model, yielding the quarterly allometric ration:
     $$M_L = \Delta t_{\text{days}} \cdot m_0 \, L^{m_1}, \qquad m_0 = 3.167 \times 10^{-4}\text{ kg}\cdot\text{d}^{-1}\cdot\text{cm}^{-m_1},\quad m_1 = 2.66$$ {#eq-kulatska-max-consumption}

4. **Diagnosing Regime Shifts and Spatial Mismatch [@kulatska2019, pp. 343--347]**:
   Optimized across 21 likelihood components using Simulated Annealing $\to$ Hooke & Jeeves $\to$ BFGS with [@stefansson2003] / [@taylor2007] iterative inverse-variance reweighting, the model accurately captured both ontogenetic diet shifts and interdecadal clupeid switches during $1974\text{--}1988$ (Bray–Curtis dissimilarity $0.08\text{--}0.28$), but exhibited higher diet residuals in $2007\text{--}2013$ ($\text{BC} = 0.24\text{--}0.52$). Kulatska et al. showed that this residual pattern was diagnostic of two unmodeled physical processes in the recent Baltic regime: (a) expanding bottom hypoxia decimating benthic *Saduria entomon* biomass, and (b) **spatial predator–prey mismatch** (Eastern Baltic cod contracting southwestward into ICES SD 25 while sprat and herring shifted northeastward into SD 28–32).

---

## Single-Species Length-Based & Integrated Assessments as a Methodological Bridge

Building on classical cohort and yield theory [@beverton_holt1956; @beverton_holt1957; @cadima2003; @lart2022], single-species integrated assessment models (`MULTIFAN-CL` [@fournier1998], `Stock Synthesis` `SS3` [@methot2013], Reviewed by [@maunder2013]), size-transition matrix assessments [@punt2016], and non-equilibrium length-based mixed-effects models (`LIME` [@rudd2018]) solved several fundamental statistical problems in fitting length-composition data that are directly applicable to size-spectrum models like `mizer`.

### `MULTIFAN-CL` and Fournier's Robustified Length-Frequency Likelihood
Developed by David Fournier, John Hampton, and John Sibert [@fournier1998] for Western and Central Pacific tuna assessments, **MULTIFAN-CL** was the first large-scale ADMB integrated model to estimate age-structured dynamics directly from length-frequency data.

1. **Length-at-Age Mixture Model**:
   At any time $t$, the lengths of fish in age class $a \in \{1, \dots, A\}$ are assumed normally distributed with mean $\mu_{a,t}$ (following a von Bertalanffy curve with seasonal oscillations and cohort-specific growth deviations) and standard deviation $\sigma_a$ (scaling linearly with mean length):
   $$\sigma_a = \sigma_1 \exp\!\left(\lambda_{\sigma} \left(\frac{\mu_a - \mu_1}{\mu_A - \mu_1}\right)\right)$$ {#eq-mfcl-length-sd}
   Crucially, fishery selectivity $s_{f,a}$ is parameterized as a smooth cubic-spline function of **mean length** $\mu_a$ rather than age (`length-dependent selectivity`), reducing the number of selectivity parameters and linking cohorts across growth variations.

2. **Fournier's Robustified Normal Likelihood for Length Compositions**:
   Standard multinomial likelihoods $\sum_l \tilde{P}_l \ln P_l$ are notoriously sensitive to **outliers in the tails** (a single large fish caught in a bin where $\hat{P}_l = 10^{-6}$ creates a massive log-likelihood penalty) and **overdispersion from schooling** (fish of similar sizes school together, so haul length samples have much smaller effective sample size $\tau_t$ than the number of fish measured).
   [@fournier1998, Eq. 3.6] derived the **robustified normal likelihood** for observed length proportions $\tilde{p}_{f,t,l}$ vs. predicted proportions $\hat{p}_{f,t,l}$ across $L$ length bins:

   $$\ln L_{\text{MFCL}} = -\frac{1}{2} \sum_{f,t} \sum_{l=1}^{L} \ln\!\left(2\pi \left(\xi_{f,t,l} + \frac{0.1}{L}\right)\right) + \sum_{f,t} \sum_{l=1}^{L} \ln\!\left[ \exp\!\left(-\frac{\left(\tilde{p}_{f,t,l} - \hat{p}_{f,t,l}\right)^2}{2\left(\xi_{f,t,l} + \frac{0.1}{L}\right) / \tau_{f,t}}\right) + 0.01 \right]$$ {#eq-mfcl-length-loglik}

   where $\xi_{f,t,l} = (1 - \tilde{p}_{f,t,l})\tilde{p}_{f,t,l}$, $\tau_{f,t} = \min(N_{\text{fish}, f, t}, 1000)$ is the effective sample size, $+ \frac{0.1}{L}$ prevents the variance from vanishing in bins where $\tilde{p}_{f,t,l} \approx 0$ so empty tail bins do not dominate the variance weight, and $+ 0.01$ inside the logarithm **caps the maximum penalty of any single outlier bin at $\ln(0.01) \approx -4.605$** (equivalent to a $3\sigma$ normal residual!).

### `Stock Synthesis` (`SS3`): Age-Length Kernels and Growth Morphs/Platoons
In **Stock Synthesis** [@methot2013, pp. 86--90], population numbers $N_{t,g,a}$ are tracked by age $a$ and growth morph/platoon $g$, and mapped at every time step to discrete population length bins $l \in \{1, \dots, n_L\}$ via an **age-length transition matrix** $\varphi_{l, g, a}$:

$$\varphi_{l, g, a} = \Phi\!\left(\frac{L_{\max, l} - \tilde{L}_{g,a}}{\sigma_{g,a}}\right) - \Phi\!\left(\frac{L_{\min, l} - \tilde{L}_{g,a}}{\sigma_{g,a}}\right)$$ {#eq-mfcl-length-bin-proportion}

where $\Phi(\cdot)$ is the standard normal CDF, $\tilde{L}_{g,a}$ is the mean length at exact mid-season age for morph $g$, and $\sigma_{g,a}$ is its standard deviation.
- **Combined Age- and Length-Specific Selectivity**: For fleet $f$ with age selectivity $S^{(A)}_{f,a}$ and length selectivity $S^{(L)}_{f,l}$, the effective survival and catch-at-age are integrated over the length distribution within each age class [@methot2013, Eq. A.1.24]:
  $$F_{t,f,g,a} = E_{t,f} \, S^{(A)}_{f,a} \sum_{l=1}^{n_L} \varphi_{l,g,a} \, S^{(L)}_{f,l}$$ {#eq-mfcl-fishing-mortality}
- **Capturing Rosa Lee's Phenomenon via Platoons**: To prevent the assumption that the normal length distribution $\varphi_{l,g,a}$ resets every year regardless of size-selective fishing on the right tail of the cohort, SS3 splits each recruiting cohort into $3$ or $5$ **growth platoons** $g$ (slow, medium-slow, average, medium-fast, fast growers with initial Gaussian weights). Because fast-growing platoons enter the length-selective fishery earlier, they are depleted faster ($F_{t,f,\text{fast},a} > F_{t,f,\text{slow},a}$), naturally shifting the surviving cohort's mean length downward!

### Pure Size-Structured Assessments and Size-Transition Matrices ($\mathbf{X}_{l,l'}$)
In their comprehensive review of size-structured stock assessments (used worldwide for lobsters, crabs, abalone, scallops, prawns, and tropical finfish), [@punt2016] formalized the general discrete-time, discrete-size population equation:

$$\mathbf{N}_{t+1} = \mathbf{X}_t \, \mathbf{S}_t \, \mathbf{N}_t + \mathbf{R}_{t+1}$$ {#eq-punt-size-dynamics}

where:
- $\mathbf{N}_t = (N_{1,t}, N_{2,t}, \dots, N_{n_L,t})^\top$ is the vector of numbers across $n_L$ size bins,
- $\mathbf{S}_t = \text{diag}\!\left(\exp(-Z_{1,t}), \dots, \exp(-Z_{n_L,t})\right)$ is the diagonal survival matrix during time step $t$,
- $\mathbf{X}_t$ is the lower-triangular (since animals do not shrink in length) **$n_L \times n_L$ Size-Transition Matrix**, where $X_{l', l}$ is the probability that an individual in size bin $l$ grows into size bin $l' \ge l$ at the end of the time step ($\sum_{l' \ge l} X_{l', l} = 1$).

#### Parameterizing the Size-Transition Matrix $\mathbf{X}_{l', l}$ [@punt2016, Eqs. 2--5]
Instead of estimating $n_L(n_L-1)/2$ free transition probabilities, size-structured assessments model the continuous growth increment $\Delta L = L_{t+1} - L_l \ge 0$ of an individual at midpoint $L_l$ using a **Gamma distribution** (or truncated Normal / Log-normal distribution) with mean increment $\bar{\Delta}_l = (L_\infty - L_l)(1 - e^{-k \Delta t})$ and variance $\sigma^2_{\Delta, l} = \beta_{\Delta} \bar{\Delta}_l$ (so shape $\alpha_l = \bar{\Delta}_l / \beta_{\Delta}$ and scale $\beta_{\Delta}$):

$$X_{l', l} = \int_{L_{\min, l'} - L_l}^{L_{\max, l'} - L_l} \frac{\Delta^{\alpha_l - 1} \exp(-\Delta / \beta_{\Delta})}{\Gamma(\alpha_l) \, \beta_{\Delta}^{\alpha_l}} \, d\Delta, \qquad l' \ge l$$ {#eq-punt-gamma-transition}

> [!IMPORTANT]
> **Critical Insight from Punt et al. (2016) for Fitting `mizer` to Length Data**:
> [@punt2016, pp. 50--53] prove via simulation that in any size-structured assessment model, **asymptotic size $L_\infty$ (or growth rate $g(w)$), natural mortality $M(w)$, and dome-shaped vs. asymptotic gear selectivity $S_f(w)$ are strongly confounded** when fitted to commercial length compositions alone:
> - Absence of large fish in the catch can be explained equally well by: (a) slower growth / lower $L_\infty$, (b) higher predation/natural mortality $M$, (c) higher fishing mortality $F$, or (d) dome-shaped selectivity $S_f(w)$.
> - **How to break the confounding**: (1) Use scientific trawl survey length distributions with known asymptotic selectivity, (2) include **mean weight-at-age / conditional age-at-length** data where available, or (3) include **stomach content / predator-prey size data** which independently anchors $M2(w)$ and growth!

---

### `LIME`: Length-Based Integrated Mixed Effects for Non-Equilibrium Length Data (Rudd & Thorson 2018)
Traditional length-only assessment methods for data-limited stocks—such as Length-Based Spawning Potential Ratio (`LB-SPR`, [@hordyk2015]) and mean-length mortality estimators [@beverton_holt1956; @gedamke_hoenig2006]—assume **strict steady-state equilibrium** (constant recruitment $R_t = R_0$ and constant fishing mortality $F_t = F$). In reality, recruitment pulses and time-varying $F_t$ distort single-year length compositions, causing severe bias in equilibrium length estimators.

To relax these equilibrium assumptions while requiring as little as **1 year of length-composition data** plus known basic life-history parameters ($L_\infty, k, t_0, M, L_{m}^{50}, \alpha, \beta, \text{CV}_L$), Merrill Rudd and James Thorson [@rudd2018] developed **LIME** (*Length-based Integrated Mixed Effects*), implemented in **R and TMB**:

1. **Age-Length Transition Matrix ($\varphi_{a,l}$) and Length-Based Selectivity ($S_l$)**:
   Given von Bertalanffy mean length-at-age $L_a = L_\infty\left(1 - e^{-k(a - t_0)}\right)$ and length standard deviation $\sigma_{L,a} = L_a \cdot \text{CV}_L$ (typically $\text{CV}_L = 0.10$), the probability $\varphi_{a,l}$ of an individual of age $a$ falling into length bin $l$ with bounds $[l_{\text{low}}, l_{\text{high}})$ is [@rudd2018, Table 1, Eq. 5, p. 29]:
   $$\varphi_{a,l} = \begin{cases} \Phi\!\left(\dfrac{l_{\text{high}} - L_a}{\sigma_{L,a}}\right) & l = 1 \\[10pt] \Phi\!\left(\dfrac{l_{\text{high}} - L_a}{\sigma_{L,a}}\right) - \Phi\!\left(\dfrac{l_{\text{low}} - L_a}{\sigma_{L,a}}\right) & 1 < l < n_l \\[10pt] 1 - \Phi\!\left(\dfrac{l_{\text{low}} - L_a}{\sigma_{L,a}}\right) & l = n_l \end{cases}$$ {#eq-punt-length-at-age-distribution}
   Gear selectivity is parameterized as a two-parameter logistic function of length $x_l$ (with parameters $L_s^{50}$ and $L_s^{95}$) and mapped to age $a$ by integrating across the age-length matrix:
   $$S_l = \frac{1}{1 + \exp\!\left(-\dfrac{\ln(19)\left(x_l - L_s^{50}\right)}{L_s^{95} - L_s^{50}}\right)}, \qquad S_a = \sum_{l=1}^{n_l} \varphi_{a,l} \, S_l$$ {#eq-punt-logistic-selectivity}

2. **Mixed-Effects State-Space Dynamics (Random-Walk $F_t$ and Laplace-Integrated Recruitment $R_t$)**:
   - **Autocorrelated Recruitment Random Effects**: Annual recruitment $R_t$ is modeled as a lognormal random effect (with optional $\text{AR}(1)$ autocorrelation $\rho_R$ and TMB epsilon bias-correction; [@rudd2018, Table 2, Eqs. 7, 17--19, p. 30]):
     $$R_t \sim \text{Lognormal}\!\left(\ln R_0 - \frac{\sigma_R^2}{2},\; \sigma_R^2\right)$$ {#eq-punt-recruitment-prior}
     where recruitment standard deviation $\sigma_R$ is estimated with a hierarchical lognormal prior centered on the meta-analytic fish average ($\sigma_R \sim \text{Lognormal}(\ln(0.737),\, 0.353^2)$; [@thorson2014]).
   - **Random-Walk Penalty on Fishing Mortality $F_t$**: Annual fishing mortality $F_t$ is estimated as a fixed effect constrained by a first-difference random-walk penalty in log-space ($\sigma_F = 0.20$):
     $$\ln F_t \sim \mathcal{N}\!\left(\ln F_{t-1},\; \sigma_F^2\right)$$ {#eq-punt-fishing-mortality-walk}
   - **Predicted Catch Length Proportions ($\hat{p}_{t,l}$)**: Predicted catch-at-age $\hat{C}_{t,a} = \frac{F_t S_a}{M + F_t S_a} N_{t,a}\left(1 - e^{-M - F_t S_a}\right)$ is projected into length bins via $\hat{C}_{t,l} = \sum_{a} \hat{C}_{t,a} \varphi_{a,l}$, giving predicted length proportions $\hat{p}_{t,l} = \hat{C}_{t,l} / \sum_{l'} \hat{C}_{t,l'}$.

3. **Self-Weighting Dirichlet-Multinomial Length Likelihood and Derived Reference Points**:
   LIME fits observed catch length counts $\tilde{\pi}_{t,l}$ using the **Dirichlet-Multinomial (DM) likelihood** of [@thorson2017] (see @sec-dm-likelihood), estimating the DM linear overdispersion parameter $\theta$ (`log_theta`) simultaneously with $F_t, L_s^{50}, L_s^{95}, \sigma_R$, and integrating over the latent recruitment vector $\mathbf{r} = \{R_t\}$ via TMB's Laplace approximation [@rudd2018, Table 2, Eqs. 16, 22].
   From the estimated parameters, LIME computes annual Spawning Potential Ratio:
   $$\text{SPR}_t = \frac{\text{SBPR}_t}{\text{SBPR}_0} = \frac{\sum_{a} \text{Mat}_a \, W_a \, N_{F_t, a}}{\sum_{a} \text{Mat}_a \, W_a \, N_{F=0, a}}$$ {#eq-punt-spr}
   Notice that because $\varphi_{a,l}$ allows cohorts of different ages to overlap in length, LIME can separate a recruitment pulse $\epsilon_{R,t}$ (which creates a mode in juvenile length bins that propagates across years at rate $g(L)$) from an increase in fishing mortality $F_t$ (which truncates the right tail of the length distribution across all adult cohorts simultaneously)—a key mechanism for dynamic length-based calibration in `mizer`!

---

## Discrete Length-Structured Multispecies Models (`LeMans`, `LeMaRns`, `CSMS`, `CSM`)

### `LeMans` and `LeMaRns` (Hall et al. 2006; Rochet et al. 2011; Thorpe et al. 2015; Spence et al. 2020, 2021)
Originally created by Stephen Hall, Collie, Jennings, Bravington, and Pope [@hall2006] for Georges Bank, extended to evaluate selective fishing across size-based models [@rochet2011], and adapted to the North Sea by [@thorpe2015] and Michael Spence et al. as the R/C++ package **`LeMaRns`** [@spence2020; @spence2021advice; @spence2021uncertainty], **LeMans** (*Length-based Multispecies analysis by numerical simulation*) is a discrete-time, discrete-length-bin multispecies model with fixed von Bertalanffy growth and one-way predation mortality feedback.

#### Mathematical Equations of `LeMaRns` [@spence2020, pp. 3--7]
Let $N_{i,j,t}$ be the number of individuals of species $i \in \{1, \dots, n_s\}$ in length class $j \in \{1, \dots, n_l\}$ with bounds $[L_j, L_{j+1})$ and midpoint $\bar{L}_j = (L_j + L_{j+1})/2$, body mass $W_{i,j} = \alpha_i \bar{L}_j^{\beta_i}$, and maturity ogive $m_{i,j} = \left(1 + \exp(-\kappa_{m,i}(\bar{L}_j - L_{\text{mat},i}))\right)^{-1}$.

1. **Reproduction & Recruitment (Step 1)**:
   Spawning stock biomass is $\text{SSB}_{i,t} = \sum_{j=1}^{n_l} m_{i,j} N_{i,j,t} W_{i,j}$. At the first time step of each year, recruits $R_{i,t}$ enter the smallest length bin $j = 1$ via a Hockey-Stick, Ricker, Beverton-Holt, or linear stock-recruitment relationship (with optional lognormal process error in [@spence2021advice]):
   $$N^{(1)}_{i,1,t} = N_{i,1,t} + R_{i,t}, \qquad N^{(1)}_{i,j,t} = N_{i,j,t} \quad (j > 1)$$ {#eq-lemarns-recruitment-step}

2. **Mortality: $M1$, Fishing $F$, and Predation $M2$ (Step 2)**:
   Total mortality in time step $\Delta t$ is $Z_{i,j,t} = \Delta t \left(M1_{i,j} + F_{i,j,t} + M2_{i,j,t}\right)$, so surviving numbers are:
   $$N^{(2)}_{i,j,t} = N^{(1)}_{i,j,t} \exp(-Z_{i,j,t})$$ {#eq-lemarns-survival-step}
   where **predation mortality $M2_{i,j,t}$** on prey $(i, j)$ by all predators $(m, k)$ is governed by a Holling Type II (or Type III) functional response:
   $$M2_{i,j,t} = \sum_{m=1}^{n_s} \sum_{k=1}^{n_l} I_{m,k} \, N^{(1)}_{m,k,t} \, \frac{\theta_{m,k,i,j}}{\sum_{a=1}^{n_s} \sum_{b=1}^{n_l} \theta_{m,k,a,b} \, N^{(1)}_{a,b,t} W_{a,b} + W_{\text{other}}}$$ {#eq-lemarns-m2-predation}
   Here:
   - $I_{m,k} = W_{m,k} \left(L_{\infty, m} / \bar{L}_k\right)^{0.4} \cdot \text{ration\_scalar}_m$ is the maximum per-capita consumption rate of predator $(m, k)$,
   - $\theta_{m,k,i,j} = \tau_{m,i} \cdot \psi_{m,k,i,j} \cdot o_{m,i}$ is the suitability product of taxonomic preference $\tau_{m,i} \in \{0, 1\}$, spatial overlap $o_{m,i}$, and Ursin lognormal weight-ratio preference:
     $$\psi_{m,k,i,j} = \exp\!\left(-\frac{\left(\log_{10}\!\frac{W_{m,k}}{W_{i,j}} - \eta_m\right)^2}{2\sigma_m^2}\right)$$ {#eq-lemarns-diet-preference}

3. **Deterministic Fractional Somatic Growth (Step 3)**:
   Under deterministic von Bertalanffy growth, the time required for an individual of species $i$ to grow from the lower boundary $L_j$ to the upper boundary $L_{j+1}$ of length bin $j$ is:
   $$\delta t_{i,j} = \frac{1}{k_i} \ln\!\left(\frac{L_{\infty, i} - L_j}{L_{\infty, i} - L_{j+1}}\right)$$ {#eq-lemarns-bin-growth-time}
   Therefore, in a discrete time step $\Delta t$ (chosen small enough that $\Delta t \le \min_{i,j} \delta t_{i,j}$, e.g., $\Delta t = 0.1\text{ yr}$), a fixed fraction $\phi_{i,j}$ of individuals in length bin $j$ advance to length bin $j+1$:
   $$\phi_{i,j} = \frac{\Delta t}{\delta t_{i,j}} = \frac{\Delta t \cdot k_i}{\ln\!\left(\frac{L_{\infty, i} - L_j}{L_{\infty, i} - L_{j+1}}\right)}$$ {#eq-lemarns-transition-proportion}
   giving the state update:
   $$N_{i,j,t+1} = (1 - \phi_{i,j}) N^{(2)}_{i,j,t} + \phi_{i,j-1} N^{(2)}_{i,j-1,t}$$ {#eq-lemarns-population-update}
   *(Note: Mathematically, this update equation is identical to the **explicit first-order upwind finite-difference scheme** of the McKendrick–von Foerster PDE in `mizer`! Consequently, it exhibits the exact same numerical Courant-Friedrichs-Lewy [CFL] artificial cohort diffusion $\approx \frac{1}{2} g \delta L (1 - \phi)$ analyzed by [@datta2010]).*

---

### `CSMS` and `CSM` (Pope et al. 2006, 2019): Moment-Based Reduction of Length Spectra
John Pope and colleagues [@pope2006; @gislason2010; @pope2019] developed two remarkably parsimonious length-based multispecies models to explore multispecies MSY and community structure:
1. **CSMS (Charmingly Simple Multispecies Size-based model; [@pope2006])**:
   Parameterizes an entire 13-species North Sea fish community using only **15 life-history scaling parameters** (linking von Bertalanffy $k_i = c_k L_{\infty,i}^{-0.31}$, maturity $L_{\text{mat},i} = c_m L_{\infty,i}$, Ricker recruitment parameters, and Andersen-Ursin predation solely to each species' asymptotic length $L_{\infty,i}$), reproducing the community size spectrum slope and steeper slopes under fishing.
2. **CSM (Charmingly Simple Model; [@pope2019])**:
   Instead of tracking 50 discrete length bins per cohort, CSM notes that the log-length distribution of a cohort is approximately Gaussian and therefore **tracks only the 0th through 4th moments** $\left(M_0 = N_{s,a},\; M_1 = \mathbb{E}[\ln L],\; M_2 = \text{Var}(\ln L),\; \dots\right)$ of each cohort's size distribution! Because the product of a Gaussian cohort log-length density and a Gaussian Ursin log-size preference kernel has an **exact closed-form Gaussian integral**, CSM computes size-based predation mortality $M2_{s,a}$ analytically without any discrete length bins, achieving a $50\times$ speedup over bin-based models.
