# Mathematical & Statistical Literature Review of Multispecies Stock Assessment Models

> **A Comprehensive Reference Monograph on Multispecies Stock Assessment, Length-Structured Integrated Assessment, Stomach Content Likelihoods, and Dynamic Calibration for Size-Spectrum Models (`mizer`)**

---

## Navigation & Chapter Structure

This repository contains a six-chapter mathematical and statistical literature review rendered as a Quarto book (configured in [`_quarto.yml`](_quarto.yml)) with a comprehensive BibTeX bibliography ([`references.bib`](references.bib)) and reference section ([`references.md`](references.md)).

| Chapter | Title & Core Scope | Key Models & Methodologies Covered |
| :--- | :--- | :--- |
| **[Overview](index.md)** | **Overview & Thematic Guide** | Thematic guide to the literature across all eight methodological areas. |
| **[Chapter 1](01_taxonomy_and_foundations.md)** | **Taxonomy and Theoretical Foundations of Multispecies Stock Assessment Models** | Strategic vs. Tactical MICE "sweet spot"; State-space coordinate systems (Age, Age-Length, Length, Continuous Mass, Biomass Pools); One-Way ($M2$-only) vs. Two-Way (Bioenergetic) trophic feedback; Andersen & Ursin (1977) differential equations and Livingston (1983) nonlinear sensitivity analysis, and how `SMS` vs. `mizer` diverged; Jacobian community matrices ($\mathbf{J}$) and quadratic multispecies yield approximations. |
| **[Chapter 2](02_age_structured_multispecies_models.md)** | **Age-Structured Multispecies Stock Assessment Models** | **`MSVPA`/`MSFOR`** iterative cohort equations; **`SMS` / `smsR`** (Lewy & Vinther 2004) seasonal survival, separable $F$, Andersen–Ursin size preference, evacuation rations, exact $M2$, and Lognormal + Dirichlet likelihoods; **`MSASA` / `MSCAA` / `Hydra`** (Kinzey & Punt 2009; Van Kirk et al. 2010; Curti et al. 2013; Trijoulet et al. 2019; 7 predator functional responses including Pre-emption and Ecosim Foraging Arena); **`CEATTLE` / `Rceattle`** (Holsman et al. 2016; Adams et al. 2022; temperature-dependent generalized VBGF, Wisconsin bioenergetics, 12 MBRP harvest control rules, and hybrid predator scaling). |
| **[Chapter 3](03_age_length_and_length_structured_models.md)** | **Age-Length and Length-Structured Assessment Models** | **`Gadget2` & `gadget3`** (5D state tensor $N_{s,a,l,r,t}, W_{s,a,l,r,t}$, MULTSPEC/Jones consumption-coupled growth, **Beta-Binomial** discrete bin transitions `g3a_grow_impl_bbinom`, Type II/III predation with `dif_pmin` overconsumption cap, Stefánsson's two-stage iterative reweighting, and Kulatska et al. 2019 Baltic cod-herring-sprat-*Saduria* case study); **`MULTIFAN-CL`** (Fournier's robustified normal length likelihood); **`Stock Synthesis` (`SS3`)** (growth morphs/platoons & Rosa Lee's phenomenon); **Size-Transition Matrices $\mathbf{X}_{l,l'}$** (Punt et al. 2016); **`LIME`** (Rudd & Thorson 2018 length-based mixed effects); **`LeMans` / `LeMaRns`** and **`CSMS` / `CSM`** (moment-based length spectra). |
| **[Chapter 4](04_size_spectrum_and_food_web_models.md)** | **Continuous Size-Spectrum (`mizer`) and State-Space Food-Web Models** | **`mizer`** complete PDEs, bioenergetics, reproduction, starvation, and semi-implicit upwind discretization; **Jump-Growth Fokker–Planck PDE** (Datta, Delius & Law 2010) and second-order growth diffusion $d_i(w,t)$; **Regulation without SRR** (Canales, Delius & Law 2020); **`Ecopath`, `Ecosim`, `Rpath`**, Walters' **Foraging Arena Theory** ($v_{ij}$ exchange), and **`EcoState`** (Thorson et al. 2024 TMB state-space ODE model). |
| **[Chapter 5](05_stomach_content_data_and_likelihoods.md)** | **Stomach Content Data, Gastric Evacuation Physics, and Diet Likelihood Formulations** | Biophysical **cylinder model of gastric evacuation** ($\frac{dS}{dt} = -\rho L^\lambda e^{\delta T} E^{-\mu} S^{0.5}$, Andersen & Beyer 2005); **Jensen's inequality pooling bias** and Taylor correction (Rindorf & Lewy 2004); Prey condition & `%F` vs. `%W` (Buckland et al. 2017); **ICES & `DAPSTOM`** databases, `FishStomachs` CPUE raising, and **`sdmTMB`** Tweedie spatiotemporal diet standardization (Jacobsen et al. 2023); Multinomial, Dirichlet, **Dirichlet-Multinomial** (Thorson et al. 2017 linear & saturating parameterizations), Logistic-Normal, and **Exponential Family / Gamma log-PPMR likelihoods** (Morgan 2022). |
| **[Chapter 6](06_statistical_fitting_uncertainty_and_mizer_synthesis.md)** | **Statistical Estimation, Data Weighting, Uncertainty Quantification, and `mizer` Roadmap** | Master comparative table across 8 model families; **Francis (2011) `TA1.8`** vs. McAllister–Ianelli vs. Dirichlet-Multinomial (Thorson et al. 2017) data weighting; **TMB Laplace approximation** (Kristensen et al. 2016); GA + Parallel Tempering MCMC (Zhang et al. 2015) and History Matching + SMC + Multi-Model Ensembles (Spence et al. 2016, 2018, 2021); **Five-Pillar Architectural Blueprint for Dynamic Statistical Stock Assessment in `mizer`**. |
| **[References](references.md)** | **Bibliography** | Complete bibliography rendered from [`references.bib`](references.bib) with DOIs. |

---

## Building the Book

To render the HTML book into `docs/`:

```bash
quarto render
```
