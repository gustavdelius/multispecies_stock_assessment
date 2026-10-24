# Stomach Content Data, Gastric Evacuation Physics, and Diet Likelihood Formulations {#sec-stomach-data}

## Why Stomach Content Data Cannot Be Treated as Direct Diet Proportions

A common pitfall when calibrating ecological models is to treat the wet weight $S_{i,k}$ of prey $i$ observed inside stomach $k$ as a direct measure of the predator's consumption rate $C_i$. In reality, what is observed inside a sampled fish stomach is the **standing stock of partially digested prey**, which is related to the **ingestion flux (consumption rate)** only through the physics of **gastric evacuation** and **spatial sampling design** [@rindorf2004; @andersen_beyer2005; @buckland2017; @jacobsen2023]:

$$\text{Ingestion Rate } C_i(t) \xrightarrow{\text{Digestion / Evacuation Rate } \frac{dS_i}{dt}} \text{True Stomach Content } S_i(t) \xrightarrow{\text{Trawl Survey \& Lab Protocol}} \text{Observed Stomach Data}$$

---

## Biophysical Modelling of Gastric Evacuation (Andersen & Beyer 2005)

### The Square-Root Cylinder Evacuation Equation
Through extensive laboratory feeding experiments on gadoids and other predatory fishes, @andersen_beyer2005 established the mechanistic **cylinder model of gastric evacuation**.

Because digestive enzymes and hydrochloric acid are secreted by the stomach mucosa and act only on the **outer surface area** of the food bolus in the stomach (modeled as a cylinder of fixed length $l_s \propto L_{\text{pred}}$ and radius $r$, where volume/mass is $S \propto l_s r^2$ and surface area is $A_s \propto l_s r \propto S^{0.5}$), the instantaneous gastric evacuation rate (GER) of total stomach mass $S(t)$ (in $\text{g} \cdot \text{h}^{-1}$) follows a power law with exponent $\alpha = 0.5$ [@andersen_beyer2005, Eqs. 1--3]:

$$\frac{dS(t)}{dt} = -\rho \, S(t)^{\alpha}, \qquad \alpha \approx 0.5$$

where the rate parameter $\rho$ scales with predator total length $L$ (cm), ambient water temperature $T$ ($^\circ\text{C}$), and meal energy density $E$ ($\text{kJ} \cdot \text{g}^{-1}$ wet weight):

$$\rho(L, T, E) = \rho_0 \, L^{\lambda} \, \exp(\delta T) \, E^{-\mu}$$

Typical empirical values for North Sea and Baltic piscivores (cod, whiting, saithe; @andersen_beyer2005 [Table 1]) are:
- Surface-to-volume exponent: $\alpha = 0.50$,
- Predator length scaling: $\lambda \approx 1.30\text{--}1.44$ (since $\lambda + 3\alpha \approx 2.8\text{--}2.95$, matching the allometric scaling of maximum consumption $W^{0.75\text{--}0.8}$),
- Temperature coefficient: $\delta \approx 0.08\text{--}0.10\text{ }^\circ\text{C}^{-1}$,
- Prey energy-density resistance exponent: $\mu \approx 0.85\text{--}0.90$ (lipid-rich, high-energy prey such as herring or sprat [$E \approx 6\text{--}9\text{ kJ}\cdot\text{g}^{-1}$] are evacuated significantly more slowly in wet weight per hour than lean fish or crustaceans [$E \approx 3.5\text{--}4.5\text{ kJ}\cdot\text{g}^{-1}$], so **rate of energy delivery** $E \cdot |dS/dt| \propto E^{1-\mu} \approx E^{0.12}$ is nearly invariant across prey types!).

### Evacuation of Individual Prey Items in Mixed Meals
When a predator stomach contains $K$ distinct prey items/categories with individual wet masses $S_1, S_2, \dots, S_K$ and total stomach mass $S_{\text{tot}} = \sum_{k=1}^K S_k$, the evacuation rate of prey category $i$ is proportional to its mass share $S_i / S_{\text{tot}}$ in the cylinder multiplied by the total evacuation rate [@andersen_beyer2005, Eq. 6; @rindorf2004, Eq. 2]:

$$\text{GER}_i = -\frac{dS_i}{dt} = \rho_i \left(\frac{S_i}{S_{\text{tot}}}\right) S_{\text{tot}}^{\alpha} = \rho_i \, S_i \, S_{\text{tot}}^{\alpha - 1}$$

where $\rho_i = \rho_0 L^\lambda e^{\delta T} E_i^{-\mu}$ accounts for the specific energy density and cuticle/exoskeleton resistance of prey $i$.

---

## The Nonlinear Stomach-Pooling Bias (Rindorf & Lewy 2004)

In over a century of stomach sampling (including the 1981 and 1991 ICES North Sea Stomach Sampling Projects), stomachs of predators of the same species and length class caught in the same trawl haul were **physically pooled in a bucket aboard ship** to save laboratory processing time, and only the **mean stomach content weight** $\bar{S}_i = \frac{1}{n} \sum_{k=1}^n S_{i,k}$ across the $n$ sampled predators (including empty stomachs) was recorded [@rindorf2004].

### Proof of Jensen's Inequality Bias
Over a stationary feeding period of duration $D$, the true mean population consumption rate $\bar{C}_i$ per predator is the expectation of the individual evacuation rate across all individual stomachs $k \in \{1, \dots, n\}$:

$$\bar{C}_{\text{true}, i} = D \cdot \mathbb{E}\!\left[ \rho_i \, S_{i,k} \, S_{\text{tot}, k}^{\alpha - 1} \right] = \frac{D}{n} \sum_{k=1}^{n} \rho_i \, S_{i,k} \, S_{\text{tot}, k}^{\alpha - 1}$$

However, when stomachs are pooled prior to calculation, only the sample means $\bar{S}_i = \frac{1}{n}\sum_{k=1}^n S_{i,k}$ and $\bar{S}_{\text{tot}} = \frac{1}{n}\sum_{k=1}^n S_{\text{tot},k}$ are known, and the naive pooled estimator is:

$$\hat{C}_{\text{pooled}, i} = D \, \rho_i \, \bar{S}_i \, \bar{S}_{\text{tot}}^{\alpha - 1}$$

For total consumption ($S_i = S_{\text{tot}}$), $\bar{C}_{\text{true}} = D \rho \, \mathbb{E}[S_{\text{tot}}^\alpha]$ whereas $\hat{C}_{\text{pooled}} = D \rho \, (\mathbb{E}[S_{\text{tot}}])^\alpha$. Because $f(S) = S^\alpha$ is a **strictly concave function** for $0 < \alpha < 1$ ($\alpha \approx 0.5$), **Jensen's Inequality** guarantees that:

$$\mathbb{E}\!\left[S_{\text{tot}}^{\alpha}\right] < \Big(\mathbb{E}[S_{\text{tot}}]\Big)^{\alpha} \implies \hat{C}_{\text{pooled}} > \bar{C}_{\text{true}}$$

@rindorf2004 demonstrated analytically and via bootstrap simulations on individual cod stomachs that:
1. **Total Consumption Overestimation**: Because individual fish stomachs are highly skewed (many empty or near-empty stomachs, and a few stomachs with recently ingested large fish), applying $\bar{S}_{\text{tot}}^{0.5}$ to pooled haul averages **overestimates total food consumption by $25\%$ to $>100\%$**!
2. **Prey-Composition Bias Toward Large-Meal Piscivory**: When a predator occasionally eats a large fish prey item $A$ (high $S_{\text{tot}, k}$, so $S_{\text{tot}, k}^{\alpha - 1} = S_{\text{tot}, k}^{-0.5}$ is small) and frequently eats small benthic invertebrate meals $B$ (low $S_{\text{tot}, k}$, so $S_{\text{tot}, k}^{-0.5}$ is large), the naive pooled formula $\hat{C}_{\text{pooled}, i} \propto \bar{S}_i \bar{S}_{\text{tot}}^{-0.5}$ applies the **same** average total stomach weight $\bar{S}_{\text{tot}}^{-0.5}$ to both prey types. This **overestimates the diet share of large fish prey** and **underestimates the diet share of small invertebrate prey** [@rindorf2004, pp. 2491--2494]!

### Second-Order Taylor Bias Correction (Rindorf & Lewy 2004)
When individual stomachs cannot be analyzed separately, @rindorf2004 [Eqs. 11--14] derived a second-order Taylor series correction for total consumption among non-empty stomachs (where $p_+$ is the proportion of non-empty stomachs in the sample and $\text{CV}_+$ is the coefficient of variation of stomach content weight among feeding fish):

$$\bar{C}_{\text{corrected}} \approx D \, \rho \, p_{+}^{1 - \alpha} \, \bar{S}_{\text{tot}}^{\alpha} \left[ 1 + \frac{\alpha(\alpha - 1)}{2} \text{CV}_{+}^2 \right]$$

For $\alpha = 0.5$, $\frac{\alpha(\alpha - 1)}{2} = -\frac{1}{8} = -0.125$, reducing the pooling bias by more than an order of magnitude when $p_+$ (fraction of non-empty stomachs) is recorded.

---

## Prey Condition, Differential Digestion, and Dietary Metrics (Hyslop 1980; Buckland et al. 2017)

In empirical stomach content analysis, prey are quantified using one or more of four standard metrics [@hyslop1980; @buckland2017, Table 1]:
1. **Gravimetric ($\%W$)**: Wet or dry mass of prey category $i$ as a percentage of total stomach content mass.
2. **Volumetric ($\%V$)**: Displaced volume, grid area, or visual points of prey $i$ as a percentage of total volume.
3. **Numerical ($\%N$)**: Count of prey individuals of category $i$ as a percentage of total counted prey items.
4. **Frequency of Occurrence ($\%F_i$)**: Percentage of sampled non-empty stomachs in which prey category $i$ is present (binary presence/absence per stomach).

### Biases Revealed by Prey Condition Grading (Buckland et al. 2017)
By grading stomach contents across **7 digestion states** (from State 1 = intact freshly ingested prey to State 6 = isolated hard parts/otoliths/chaetae and State 7 = unidentifiable organic mush and mucus; @buckland2017 [Table A1]), Buckland et al. demonstrated that:
- **Unidentifiable Matter & Mucus Adhesion**: Unidentifiable organic tissue and stomach mucus were present in $33\text{--}77\%$ of stomachs and physically adhered to small prey in $>61\%$ of stomachs, making direct separation for $\%W$ and $\%V$ imprecise.
- **Differential Retention of Hard Parts**: Soft-bodied prey (polychaetes, small teleost larvae, nematodes) digest rapidly into unidentifiable mush, whereas crustacean carapaces, bivalve shells, and large teleost otoliths persist in stomachs for days. Raw $\%W$ of partially digested remains can either severely underestimate teleost importance (if only a tiny otolith remains) or overestimate shelled prey (if indigestible shells are weighed).
- **Prey Weight Reconstruction vs. Digested Stomach Weight**:
   - For size-preference estimation (e.g., fitting the `mizer` or `SMS` predation kernel $\phi(w_{\text{pred}} / w_{\text{prey}})$), one requires the **fresh ingested body weight $w_{\text{prey}}$** of the prey at the moment of capture, which must be reconstructed from intact prey lengths or otolith/carapace allometric regressions ($w_{\text{prey}} = a L_{\text{otolith}}^b$).
   - Conversely, for evacuation-rate calculation (`SMS` / `FishStomachs`), one requires the **actual digested wet mass $S_i$ currently in the stomach** to plug into $\rho_i S_i S_{\text{tot}}^{\alpha - 1}$! Recording *both* digested mass $S_i$ and reconstructed ingested size/mass $w_{\text{prey}}$ (as in the ICES and DAPSTOM databases) is therefore essential.
- **Robustness of Hierarchical $\%F$**: When comparing across studies with variable digestion states, `%F` (reported in a nested taxonomic hierarchy from Class $\to$ Order $\to$ Family $\to$ Genus $\to$ Species; @buckland2017 [Table B1]) is invariant to mucus/evaporation weight distortion and partial fragmentation.

---

## Stomach Content Databases and Spatiotemporal Raising Pipelines (Jacobsen et al. 2023)

### Major Open-Access Stomach Content Databases
Two primary databases underpin North Atlantic and European multispecies stock assessments [@ices2010stomach; @jacobsen2023], complemented by global predator–prey body-size compilations [@barnes2008; @barnes2010]:
1. **ICES Stomach Content Database** (`stomachdata.ices.dk`; @ices2010stomach):
   Contains $>800,000$ prey records from the North Sea (1981 & 1991 "Years of the Stomach" plus IBTS collections), Baltic Sea (cod and flounder stomachs 1963–2021+), and Barents/Celtic Seas. Organized into four normalized relational tables linked by `tblUploadID` / `tblHaulID` / `tblPredatorID` [@jacobsen2023, pp. 13--17]:
   - `File_Information (FI)`: Upload metadata, country, institute, evacuation parameters.
   - `Haul_Information (HH)`: ICES survey, gear, ship, haul duration, latitude/longitude, ICES statistical rectangle, bottom depth, bottom temperature $T$, year, month, time of day.
   - `Predator_Information (PI)`: Predator species (`AphiaID`), predator length class $L$, individual or pooled sample flag, number of feeding stomachs $n_+$, number of regurgitated stomachs $n_{\text{reg}}$, number of empty stomachs $n_0$, mean predator weight.
   - `Prey_Information (PP)`: Prey species (`AphiaID`), digestion stage, prey length class $l$, prey count $N_{\text{prey}}$, total observed wet weight $S_i$ in stomach, and gravimetric method.
2. **CEFAS `DAPSTOM` Database** [@pinnegar2014; @jacobsen2023, pp. 10--12; @morgan2022]:
   Contains $>225,000$ stomach records from 1836 to 2016 across Northwest European shelf seas, providing individual predator-prey length and mass ratios (PPMR) across a 180-year historical baseline.

---

### Design-Based Stratified Raising (`FishStomachs` R Package)
Because research vessels sample a fixed number of stomachs per length class per haul (e.g., at most 5 cod per 5-cm length bin per haul) rather than sampling proportionally to true predator density in the sea, simple arithmetic averages across all sampled stomachs are **spatially biased** toward areas with low predator density.

The ICES standard raising pipeline—implemented in Morten Vinther's **`FishStomachs`** R package [@jacobsen2023, pp. 34--38]—raises individual/haul stomach observations to the population level in three hierarchical steps:
1. **Step 1 (Regurgitation Correction & Evacuation Conversion at Haul Level $h$)**:
   Within each haul $h$ and predator length class $L$, regurgitated stomachs ($n_{\text{reg}, h, L}$) are assumed to have had the same mean diet as feeding stomachs ($n_{+, h, L}$), so observed weights are scaled by $\frac{n_{+,h,L} + n_{\text{reg},h,L}}{n_{+,h,L} \left(n_{+,h,L} + n_{\text{reg},h,L} + n_{0,h,L}\right)}$ and converted from stomach mass $S_{i,h,L}$ to evacuated ration rate $R_{i,h,L}$ using the @andersen_beyer2005 cylinder equation and @rindorf2004 pooling bias correction.
2. **Step 2 (Haul-to-Rectangle Raising via √Stomachs or CPUE)**:
   Within each ICES statistical rectangle $r$, haul rations $R_{i,h,L}$ are weighted by $\sqrt{n_{\text{stom}, h, L}}$ (or haul predator CPUE) to obtain the rectangle mean ration $R_{i,r,L} = \frac{\sum_{h \in r} w_{h,L} R_{i,h,L}}{\sum_{h \in r} w_{h,L}}$.
3. **Step 3 (Rectangle-to-Stock Raising via Predator Survey Abundance $I_{r,L}$)**:
   Across all rectangles $r \in \mathcal{R}$ in the stock assessment area, rectangle rations are weighted by the **survey abundance index $I_{r,L}$ of predators of length $L$ in rectangle $r$**:
   $$\bar{R}_{\text{pop}, i, L} = \frac{\sum_{r \in \mathcal{R}} I_{r, L} \, R_{i, r, L}}{\sum_{r \in \mathcal{R}} I_{r, L}}, \qquad \hat{W}_{\text{pop}, i, L} = \frac{\bar{R}_{\text{pop}, i, L}}{\sum_{k} \bar{R}_{\text{pop}, k, L}}$$

---

### Model-Based Spatiotemporal Standardization (`sdmTMB`)
When stomach sampling is patchy in space and time (missing ICES rectangles in certain quarters/years), design-based raising suffers from missing-stratum bias. @jacobsen2023 [pp. 26--33] developed a **model-based spatiotemporal standardization workflow** using **`sdmTMB`** [@anderson2022]:
1. Define the **Feeding Ratio (FR)** of prey $i$ in individual predator $k$ of length $L_k$ at spatial coordinates $\mathbf{s}_k$ and year-quarter $t$ as:
   $$\text{FR}_{i,k} = \frac{S_{i,k}}{W_{\text{pred}, k}}$$
2. Fit a spatial/spatiotemporal Generalized Linear Mixed Model (GLMM) with a **Tweedie** distribution (power variance $V(\mu) = \phi \mu^p$ with $1 < p < 2$, which naturally models continuous positive ratios with a point mass at exact zero $\text{FR}_{i,k} = 0$!) and log link:
   $$\mathbb{E}[\text{FR}_{i,k}] = \mu_{i,k}, \qquad \ln \mu_{i,k} = \mathbf{X}_k \boldsymbol{\beta} + s(L_k) + \omega_i(\mathbf{s}_k) + \epsilon_i(\mathbf{s}_k, t)$$
   where $s(L_k)$ is a smooth spline of predator length, $\omega_i(\mathbf{s}) \sim \text{MVN}(\mathbf{0}, \boldsymbol{\Sigma}_\omega)$ is a time-invariant spatial Gaussian Markov Random Field (GMRF) approximated via SPDE on a finite-element mesh, and $\epsilon_i(\mathbf{s}, t)$ is a spatiotemporal random field.
3. Predict $\hat{\mu}_{i}(\mathbf{s}_g, L, t)$ across a fine spatial prediction grid $g \in \mathcal{G}$, weight by spatially predicted predator density $\hat{D}_{\text{pred}}(\mathbf{s}_g, L, t)$ (also estimated via `sdmTMB` on trawl survey catches), and normalize across prey types $i$ to obtain spatially unbiased population diet proportions and their standard errors!

---

## Statistical Observation Models and Likelihoods for Stomach Content Data

Once population-level diet proportions $\mathbf{W}_{U} = (W_{1,U}, \dots, W_{K,U})^\top$ (satisfying $\sum_{k=1}^K W_{k,U} = 1$) or individual predator-prey mass ratios (PPMR) are prepared for predator stratum $U = (A, L, y, q)$, how should they enter the stock assessment objective function?

### Summary Comparison of Compositional Diet Likelihoods

| Likelihood Distribution | Data Support | Log-Likelihood $\ln L(\mathbf{W}_U \mid \hat{\mathbf{p}}_U, \theta)$ | Variance Structure $\text{Var}(W_k)$ | Handling of Zeros ($W_k = 0$) | Used In |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1. Multinomial** | Proportions $\mathbf{W}_U \in [0,1]^K$ scaled by sample size $N_U$ | $N_U \sum_{k=1}^K W_{k,U} \ln \hat{p}_{k,U} + \text{const}$ | $\frac{\hat{p}_k(1-\hat{p}_k)}{N_U}$ (Fixed by $N_U$; no free dispersion parameter) | **Natural**: $W_{k,U} = 0$ simply contributes $0 \cdot \ln \hat{p}_{k,U} = 0$! | `MSCAA`, `CEATTLE`, `Gadget`, `SMS` (Option A) |
| **2. Dirichlet** | Strictly positive simplex $W_{k,U} \in (0, 1)$, $\sum_k W_{k,U} = 1$ | $\ln\Gamma(v_U) - \sum_{k=1}^K \ln\Gamma(v_U \hat{p}_{k,U}) + \sum_{k=1}^K (v_U \hat{p}_{k,U} - 1)\ln W_{k,U}$ | $\frac{\hat{p}_k(1-\hat{p}_k)}{v_U + 1}$ (Estimable precision $v_U = N_U/\theta - 1$) | **Undefined at $W_k = 0$** ($\ln 0 = -\infty$); requires pooling zero bins or tail compression | `SMS` (Default), `EcoState`, `MSCAA` |
| **3. Dirichlet-Multinomial (DM)** | Count vector $N_U \mathbf{W}_U$ with overdispersion $\beta = \theta N_U$ | $\ln\Gamma(N_U+1) + \ln\Gamma(\beta) - \ln\Gamma(N_U+\beta) + \sum_{k=1}^K \left[\ln\frac{\Gamma(N_U W_k + \beta \hat{p}_k)}{\Gamma(N_U W_k + 1)\Gamma(\beta \hat{p}_k)}\right]$ | $\frac{\hat{p}_k(1-\hat{p}_k)}{N_U} \left(\frac{N_U + \beta}{1 + \beta}\right)$ (Linear or saturating parameterization; @thorson2017) | **Natural**: Exact handling of $W_{k,U} = 0$ AND estimable overdispersion $\theta$! | `SS3`, `LIME`, `WHAM`, modern `Rceattle` |
| **4. Additive Logistic-Normal (Aitchison)** | Log-ratios $y_k = \ln(W_k / W_K) \in \mathbb{R}^{K-1}$ | $-\frac{1}{2}\ln|2\pi\boldsymbol{\Sigma}_U| - \frac{1}{2}(\mathbf{y}_U - \hat{\boldsymbol{\mu}}_U)^\top \boldsymbol{\Sigma}_U^{-1} (\mathbf{y}_U - \hat{\boldsymbol{\mu}}_U) - \sum_{k=1}^K \ln W_{k,U}$ | Arbitrary $(K-1)\times(K-1)$ covariance $\boldsymbol{\Sigma}_U$ (allows positive/negative correlations!) | Requires zero replacement or censored Tobit normal integration | `SAM` (compositions), `gadget3` (`mvlogistic`) |

### Mathematical Derivation of the Compound Dirichlet-Multinomial Likelihood (Thorson et al. 2017) {#sec-dm-likelihood}
Because compositional observations in both age/length sampling and stomach content analysis exhibit strong **overdispersion** relative to a pure multinomial draw of nominal sample size $n$, @thorson2017 derived a model-based compound **Dirichlet-Multinomial (DM)** likelihood that estimates effective sample size $n_{\text{eff}}$ simultaneously with all stock assessment parameters.

1. **Hierarchical Compound Integral Formulation**:
   Let $\mathbf{p} = (p_1, \dots, p_A)$ denote the assessment model's expected proportion vector across $A$ categories ($\sum_{a=1}^A p_a = 1$), and suppose that due to **process error** (unmodeled spatial/temporal variation in cohort availability or predator foraging) or **clustered sampling** (schooling of fish within trawl hauls), the true local proportion vector $\boldsymbol{\mu} = (\mu_1, \dots, \mu_A)$ encountered during sampling is drawn from a Dirichlet distribution with expectation $\mathbf{p}$ and concentration sum $\beta = \sum_{a=1}^A \beta p_a > 0$, followed by a multinomial sample of size $n$ yielding observed proportions $\tilde{\boldsymbol{\pi}} = (\tilde{\pi}_1, \dots, \tilde{\pi}_A)$ [@thorson2017, Eqs. 2--5, p. 86]:

   $$\Pr(\tilde{\boldsymbol{\pi}} \mid \mathbf{p}, n, \beta) = \int \text{Multinomial}(n\tilde{\boldsymbol{\pi}} \mid \boldsymbol{\mu}, n) \, \text{Dirichlet}(\boldsymbol{\mu} \mid \beta \mathbf{p}) \, d\boldsymbol{\mu}$$

   Evaluating this Dirichlet integral analytically yields the closed-form **Dirichlet-Multinomial log-likelihood**:

   $$\ln f_{\text{DM}}(\tilde{\boldsymbol{\pi}} \mid \mathbf{p}, n, \beta) = \ln\Gamma(n + 1) - \sum_{a=1}^{A} \ln\Gamma(n\tilde{\pi}_a + 1) + \ln\Gamma(\beta) - \ln\Gamma(n + \beta) + \sum_{a=1}^{A} \Big[\ln\Gamma(n\tilde{\pi}_a + \beta p_a) - \ln\Gamma(\beta p_a)\Big]$$

   Notice that when $\tilde{\pi}_a = 0$, the term inside the final sum is $\ln\Gamma(\beta p_a) - \ln\Gamma(\beta p_a) = 0$, so **exact observed zeros require no ad hoc constant addition or bin pooling**!

2. **Two Parameterizations of Effective Sample Size ($n_{\text{eff}}$)**:
   The marginal variance of an observed proportion $\tilde{\pi}_a$ under the DM distribution is $\text{Var}(\tilde{\pi}_a) = \frac{p_a(1 - p_a)}{n}\left(\frac{n + \beta}{1 + \beta}\right)$. Equating this to the equivalent multinomial variance $\frac{p_a(1 - p_a)}{n_{\text{eff}}}$ yields $n_{\text{eff}} = \frac{n + n\beta}{n + \beta}$. @thorson2017 [Eqs. 9--12, p. 86] introduced two distinct parameterizations of how $n_{\text{eff}}$ scales with input sample size $n_y$ across years:
   - **Parameterization #1 (Linear Scaling, $\beta_y = \theta n_y$, Default)**:
     $$n_{\text{eff}, y} = \frac{1 + \theta n_y}{1 + \theta} = \frac{1}{1 + \theta} + n_y \left(\frac{\theta}{1 + \theta}\right)$$
     As $\theta \to \infty$, $n_{\text{eff}, y} \to n_y$ (recovering the unweighted multinomial); for small $\theta$ and large $n_y$, $\frac{\theta}{1 + \theta} \in (0, 1)$ acts directly as a **multiplicative data-weighting factor** $\lambda$ (analogous to McAllister–Ianelli or Francis `TA1.8` reweighting), suitable when overdispersion arises from **clustered/schooling sampling**.
   - **Parameterization #2 (Saturating Scaling, $\beta_y = \beta$ constant)**:
     $$n_{\text{eff}, y} = \frac{n_y + n_y \beta}{n_y + \beta} \xrightarrow{n_y \gg \beta} 1 + \beta$$
     Here $n_{\text{eff}, y}$ increases linearly with $n_y$ at small sample sizes ($n_y \ll \beta$) but **asymptotes at an upper ceiling $1 + \beta$** as $n_y \to \infty$, which is theoretically exact when overdispersion arises from **process error / model misspecification** in selectivity or availability (because measuring infinitely many fish in year $y$ cannot reduce year-specific process error $\text{Var}(\mu_a) = \frac{p_a(1-p_a)}{1+\beta}$!).

> [!TIP]
> **Choosing Between Dirichlet, Dirichlet-Multinomial, and Logistic-Normal for Multispecies Assessments**:
> While the **Dirichlet** distribution (used in `SMS` and `EcoState`) works well for aggregated species-level diet shares where all modeled prey are present ($W_{k,U} > 0$), when diet data are disaggregated simultaneously by **prey species AND prey size bins** (as required for `Gadget` or dynamic `mizer`), many $(i, l)$ prey bins have observed $W_{k,U} = 0$. In that regime, the **Dirichlet-Multinomial** [@thorson2017] is strictly superior to the plain Dirichlet because $W_{k,U} = 0$ is in the natural support of the distribution while still estimating $\theta_{\text{stom}}$ via maximum likelihood. Note, however, that both Dirichlet and DM imply strictly negative correlations $\text{Corr}(\tilde{\pi}_a, \tilde{\pi}_{a'}) < 0$ between bins; if positive correlations between adjacent length bins are strong, the **Logistic-Normal** or **Tweedie spatial GLMM** (`sdmTMB`) should be used.

---

## Fitting Exponential Family & Gamma Distributions to Predator–Prey Mass Ratios (Morgan 2022)

For continuous size-spectrum models like `mizer`, stomach content data provide a second, uniquely powerful source of information beyond species diet shares: the **empirical distribution of individual Predator–Prey Mass Ratios (PPMR)**, $\text{PPMR} = w_{\text{pred}} / w_{\text{prey}}$.

In his University of York dissertation analyzing $>300,000$ stomach records from the combined CEFAS `DAPSTOM` and ICES databases, Owain Morgan [@morgan2022] derived an elegant analytical link between **realized stomach log-PPMR distributions** and the **intrinsic `mizer` feeding kernel $\varphi_{w'}(w)$** using the multiplicative closure properties of the exponential family.

### Linking Realized Diet Size Distribution to the Feeding Kernel via the Sheldon Spectrum
Let $w'$ denote the body mass of the predator and $w$ denote the body mass of an individual prey item in its stomach, so the logarithmic predator-prey mass ratio is:
$$x = \ln\!\left(\frac{w'}{w}\right) \iff w = w' e^{-x}$$

In `mizer`, the realized number density $I_{w'}(w)$ of prey of weight $w$ consumed by a predator of weight $w'$ is proportional to the product of the predator's intrinsic **feeding preference kernel** $\varphi_{w'}(w) = \phi(w'/w) = \phi(e^x)$ and the **ambient prey abundance spectrum** $N(w)$ in the sea [@morgan2022, Eq. 1]:

$$I_{w'}(w) \propto \varphi_{w'}(w) \, N(w)$$

By the Sheldon spectrum power law [@morgan2022, Eq. 2], background marine abundance density scales as $N(w) = n_0 w^{-\lambda}$ (with $\lambda \approx 2$ for number density per unit weight, or $w N(w) \propto w^{1-\lambda} = w^{-1}$ per logarithmic weight bin $dx = -dw/w$). Substituting $w = w' e^{-x}$ into $N(w)$ gives:

$$N(w' e^{-x}) = n_0 (w')^{-\lambda} \exp(\lambda x)$$

Therefore, in terms of the log-PPMR variable $x = \ln(w'/w)$, the **realized stomach log-PPMR probability density** $f_{\text{real}}(x)$ and the **intrinsic feeding preference kernel** $\phi(x)$ are related by an **exponential tilt factor** $e^{\tilde{\lambda} x}$ (where $\tilde{\lambda} = \lambda - 1 \approx 1$ for number counts in log-bins, or $\tilde{\lambda} = \lambda - 2 \approx 0$ for biomass shares in log-bins):

$$f_{\text{real}}(x) \propto \phi(x) \, \exp(\tilde{\lambda} x) \iff \phi(x) \propto f_{\text{real}}(x) \, \exp(-\tilde{\lambda} x)$$

### Exact Multiplicative Closure of Exponential Family & Gamma Distributions (Morgan 2022, Section 4.6)
@morgan2022 [pp. 7--9] showed that if the realized log-PPMR $X = \ln(w'/w) > 0$ follows any distribution in the natural **Exponential Family** with sufficient statistic vector $T(x)$ containing $x$:
$$f_{\text{real}}(x \mid \boldsymbol{\eta}) = h(x) \exp\!\left(\boldsymbol{\eta}^\top T(x) - A(\boldsymbol{\eta})\right)$$
then multiplying by the Sheldon abundance correction $\exp(-\tilde{\lambda} x)$ simply shifts the canonical parameter corresponding to $x$ by $-\tilde{\lambda}$, keeping the intrinsic feeding kernel $\phi(x)$ within the **exact same parametric distribution family**!

Specifically:
1. **Normal Distribution on $x = \ln(w'/w)$ (Lognormal on $w'/w$, as in standard `mizer`)**:
   If $f_{\text{real}}(x) = \mathcal{N}(x \mid \mu_{\text{real}}, \sigma^2)$, then multiplying by $e^{-\tilde{\lambda} x}$ yields another Normal distribution with the **same variance $\sigma^2$** and a shifted mean:
   $$\phi(x) = \mathcal{N}\!\left(x \;\middle|\; \mu_{\text{pref}} = \mu_{\text{real}} - \tilde{\lambda} \sigma^2,\; \sigma^2\right)$$
2. **Gamma Distribution on $x = \ln(w'/w)$ [@morgan2022, Eqs. 19--24]**:
   Because empirical log-PPMR distributions in `DAPSTOM` and ICES stomachs are positively supported ($x = \ln(w'/w) > 0$, since predators rarely eat prey heavier than themselves) and right-skewed (sharp gape-limitation cutoff at small $x \approx 2\text{--}3$, long tail toward small benthic/planktonic prey at large $x \approx 10\text{--}15$), @morgan2022 fitted a **Gamma distribution** to $x = \ln(w'/w)$:
   $$f^{\Gamma}_{\alpha, \beta}(x) = \frac{\beta^\alpha}{\Gamma(\alpha)} x^{\alpha - 1} e^{-\beta x}, \qquad x > 0$$
   When multiplied by $e^{-\tilde{\lambda} x}$, the Gamma density satisfies the exact algebraic identity [@morgan2022, Eq. 24]:
   $$e^{-\tilde{\lambda} x} f^{\Gamma}_{\alpha, \beta}(x) = \left(\frac{\beta}{\beta + \tilde{\lambda}}\right)^\alpha f^{\Gamma}_{\alpha,\; \beta + \tilde{\lambda}}(x)$$
   Thus, if the observed log-PPMR in stomachs is $\text{Gamma}(\hat{\alpha}, \hat{\beta})$ (with empirical estimates across North Sea cod, haddock, whiting, plaice, and sole clustering around $\hat{\alpha} \approx 6.8\text{--}12.3$, $\hat{\beta} \approx 1.15\text{--}1.85$; @morgan2022 [Table 1]), the underlying `mizer` feeding kernel is analytically a **$\text{Gamma}(\hat{\alpha},\; \hat{\beta} + \tilde{\lambda})$** distribution on $\ln(w'/w)$!
