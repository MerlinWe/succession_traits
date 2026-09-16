# Submission-hardening workflow

Run every command from the project root. Analytical code belongs in `scripts/`;
versioned results belong under `tables/`; figures belong under `figures/`.

## Provenance and preflight

```sh
Rscript scripts/00_submission_preflight.R
```

This writes only to `tables/diagnostics/submission_freeze/`.

## Figure 2

`scripts/04_shap_CV_uncertainty.R` is the single maintained SHAP workflow.
`04_shap.R` and `04_shap_with_CV_uncertainty.R` are compatibility entry points.

```sh
Rscript scripts/04_shap_CV_uncertainty.R --smoke-test
Rscript scripts/04_shap_CV_uncertainty.R
Rscript scripts/plot_fig2_shap_cv_uncertainty.R
```

The completed production uncertainty run contains 90 jobs and 900 fold-level
estimates. Its intervals describe stability across grouped CV partitions, not
population confidence intervals. The ratio compares the combined attribution
of five environmental predictors with the attribution of stand age.

## Figure 3: grouped PDP bootstrap

```sh
Rscript scripts/05_pdp.R --smoke-test
Rscript scripts/05_pdp.R --parallel-smoke
PDP_N_CORES=32 Rscript scripts/05_pdp.R
```

Normal outputs are written under
`tables/pdp_grouped_pid/normal_pid_b100_q25_seed42/`. Each permanent plot is
sampled as a cluster, each job is checkpointed, and rerunning skips valid jobs.
The comparison uses an observed common age grid from 10 to 100 years.

## Figure 4: grouped VEcv

```sh
Rscript scripts/06_vecv.R --smoke-test
Rscript scripts/06_vecv.R --parallel-smoke
VECV_N_CORES=32 Rscript scripts/06_vecv.R
```

Normal outputs are written under
`tables/vecv_grouped_pid/normal_pid_f10_r30_seed42/`. All inventories from a
permanent plot remain in one CV fold. Raw results are saved before filtering or
summarisation, and VEcv uses the exact out-of-fold `1 - SSE/SST` definition.

## Review gate

Do not launch the two full grouped runs until the hardened code and smoke-test
outputs have been reviewed. Do not copy grouped outputs over legacy root-level
RDS files. Figure and manuscript integration happens only after the full runs
pass completeness and result-stability audits.
