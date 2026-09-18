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

After pulling the completed production run locally, validate and plot it with:

```sh
Rscript scripts/audit_grouped_pdp.R
Rscript scripts/plot_fig3_grouped_pdp.R
```

The audit reconstructs the raw table and summaries from all 1,800 independent
checkpoints and writes manuscript-facing diagnostics to
`tables/diagnostics/submission_audit_grouped_pdp/`. The audited review figure is
saved as `figures/main/fig3_pdp_grouped_review.{png,pdf}`. Its panels distinguish
changes in environmental separation from the magnitude of slope modification;
a signed slope contrast is not interpreted as an absolute rate contrast.

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

After pulling the completed production run locally, validate and plot it with:

```sh
Rscript scripts/audit_grouped_vecv.R
Rscript scripts/plot_fig4_grouped_vecv.R
```

The audit reconstructs the combined raw table from all 540 independent
checkpoints and writes manuscript-facing summaries to
`tables/diagnostics/submission_audit_grouped_vecv/`. The review figure is saved
as `figures/main/fig4_vecv_grouped_review.{png,pdf}` and does not overwrite the
legacy Figure 4.

## Audited manuscript copy

The submission-facing revision is
`Main_Manuscript_Trait_Succession_audited_revision.docx`. It preserves the
original working manuscript and integrates the audited Figure 2-4 results. The
reproducible replacement logic is in
`scripts/revise_manuscript_submission.py`; its input must be a copy of the
working manuscript in which tracked changes have first been accepted. Figures
are read from `figures/main/`, and the requested output path must remain in the
project root.

## Review gate

The full grouped runs have passed completeness and reconstruction audits. Do not
copy grouped outputs over legacy root-level RDS files. Before submission, rerun
the two audit scripts, regenerate Figures 3 and 4, and then run
`scripts/00_submission_preflight.R`. Any analytical audit failure must be
resolved before manuscript or figure claims are changed.
