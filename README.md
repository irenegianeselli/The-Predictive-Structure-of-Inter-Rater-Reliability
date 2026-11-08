# The Predictive Structure of Inter-Rater Reliability

**Reproducibility materials for the manuscript:**  
*The Predictive Structure of Inter-Rater Reliability: Implications of Agreement and Alignment for a Future Rater*

## Overview

This repository contains the R implementation and numerical outputs of a Monte Carlo study examining the predictive information provided by observed inter-rater agreement and estimated decision-threshold dispersion for an additional rater.

The study uses an exchangeable Gaussian signal detection theory (SDT) model and evaluates prospective criterion displacement across 86,400 Monte Carlo replications.

## Repository Contents

- **R code:** Complete implementation of the simulation design, threshold estimation, predictive models, and numerical analyses.
- **Simulation outputs:** Replication-level results across sample size, panel size, prevalence, and threshold dispersion.
- **Predictive analyses:** Held-out and prevalence-transport results for Cohen's κ, flexible κ specifications, SDT-calibrated κ, and the SCI-derived prospective scale.
- **Oracle benchmarks:** Predictive results using generating threshold parameters, including marginal and conditional oracle specifications.
- **Bootstrap results:** Paired RMSE comparisons and confidence intervals.
- **Calibration diagnostics:** Numerical inversion checks, boundary frequencies, quadrature comparisons, and sensitivity analyses.

## Simulation Design

The simulation crosses:

- Sample sizes: 200 and 800
- Panel sizes: 4, 8, and 16 raters
- Prevalence: .10, .50, and .90
- Threshold standard deviation: 12 values from .05 to .90
- Replications per condition: 400

The analyses use R with a fixed initial random seed of 12345.

## Reproducibility

The R implementation documents the simulation, model-fitting procedures, validation partitions, SDT calibration, bootstrap comparisons, and numerical verification.

The accompanying output files provide the numerical results used to construct the manuscript's main and supplementary analyses.

## Scope

The results are conditional on the specified Gaussian SDT model, the exchangeability assumptions, and the availability of reference classifications. The predictive comparisons concern prospective criterion displacement and do not establish general superiority of one inter-rater reliability statistic across all measurement targets.
