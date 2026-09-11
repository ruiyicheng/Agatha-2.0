# Iterative RV diagnosis of all PFS target folders

This extension analyzes the 735 target folders in the audited PFS list, using
**all accepted RV instruments**, without the former >100-measurement cut.
It performs a blind, additive search; catalogue periods and catalogue planet
counts are not used. It does not run MCMC or read photometry.

From the project workspace (the parent of this repository):

```bash
python Agatha-2.0/scripts/pfs_search.py . analysis/pfs_735_vetting \
  --rscript .envs/agatha/bin/Rscript --workers 20
```

The root must contain `combined/`, `rv_statistics/pfs_targets_all_combined.csv`
and `rv_statistics/rv_file_audit.csv`. Python needs numpy, pandas, scipy,
matplotlib and pypdf. R uses the existing Agatha environment. The directory
is resumable: existing PDFs and summaries are skipped. `--retry-failed`
retries failed analyses. `--codes P0001 P0002` selects coded targets.

## Search and stopping rule

1. Select WN+jitter, AR1, MA1 or MA2 by BIC for each instrument. Red-noise
   models require >10 observations. Instruments with <=3 observations retain
   their supplied uncertainties, with an offset and zero fixed jitter.
2. Scan the RV residuals, conditional on selected noise and previously fitted
   components. The grid spans 1 day to max(2*baseline, 10 days), uniformly in
   frequency with at least four points per baseline resolution element. Refine
   the strongest grid peak locally.
3. If its conditional ln BF >=5, add a Keplerian, jointly refit all components,
   select instrument noise again on signal-subtracted RVs, and refit.
4. Repeat until the strongest remaining peak has ln BF <5. There is no arbitrary
   component-count cap. A minimum of ten residual degrees of freedom is
   required after accounting for five parameters per added component and the
   instrument offsets/shared trend. Optimizer and degrees-of-freedom stops
   are explicitly marked as incomplete searches, not threshold non-detections.

The scan statistic is a conditional BIC approximation with two sinusoid
coefficients: gain in log likelihood minus ln N. It is not a marginalized
Bayes factor or a calibrated global false-alarm probability. The fitted
Keplerian drop-component diagnostic separately penalizes five parameters.
The additive search can represent aliases, stellar activity, or inadequately
modeled binary motion with multiple components; these are candidate periodic
RV signals, not declared planets. Noise absorption and local minima remain
limitations. Periods below one day are not searched.

## Reports and audits

- `pdf/Pxxxx.pdf`: one coded PDF per target, including insufficient-data and
  failed-analysis reports so that no target silently disappears.
- `private_identity_key.csv`: code-to-target mapping, kept outside the PDFs.
- `run_summary.csv`: status, retained measurements, number of components,
  stopping evidence, residual quality, measured scan cost and elapsed time.
- `results/Pxxxx/`: complete numerical results, sequential periodograms,
  noise comparisons at each stage, fit components, activity and sensitivity
  diagnostics. `progress.json` records current work; errors use private logs.
- `preprocessing_audit.PRIVATE.csv`: source paths, activity schemas, epoch
  conversions, invalid uncertainties and instrument-family deduplication.
- `batch_status.json`: number of reported targets and status counts.

Spectroscopic diagnostics include all secondary local activity peaks within
5% of fitted periods, seven half-baseline moving windows, instrument-removal
checks, and exploratory +/-10-day RV/activity lag correlations. The lag
screen uses nearest matches within two days; it does not interpolate through
long gaps. Its selected correlations and standardized unit-error proxy
periodograms are exploratory, with no multiple-testing calibration.

All input folders remain distinct, including target aliases. The reports omit
explicit target names, catalogue labels and absolute times, but familiar RV
patterns can still reveal identities. The original five Harbor tasks and
their scores are unchanged.

## Verification

The existing seven non-MCMC Agatha checks cover the R likelihoods and their
Python implementation. `tests/test_pfs_search.py` additionally checks blind
two-signal recovery, noise-only stopping, explicit degrees-of-freedom limits,
noise reselection after accepted signals, and commented activity headers.
It also checks that resuming an accepted component can recover another signal.

```bash
python Agatha-2.0/tests/test_pfs_search.py
```

To rebuild PDFs from cached numerical results (without repeating fits), run
`python Agatha-2.0/scripts/pfs_refresh_reports.py analysis/pfs_735_vetting`.
Add `--watch` to render completed results as an active batch progresses.

Numerical convergence limits can be retried without discarding accepted signals:

```bash
python Agatha-2.0/scripts/pfs_retry_searches.py analysis/pfs_735_vetting \
  --rscript .envs/agatha/bin/Rscript --workers 2 --watch
```

This uses more optimizer evaluations and warm starts in a separate
`numerical_retry/` tree. When both the batch and watchers have finished, run
`python Agatha-2.0/scripts/pfs_finalize.py analysis/pfs_735_vetting` to promote
completed or more advanced retries, archive first-pass reports, and create
`report_index.csv` for all targets. Then refresh PDFs once more. The index
includes total fitting/reporting time across promoted retries. Numerical or
degrees-of-freedom limits remain explicit if the threshold was not reached.

Run `python Agatha-2.0/scripts/pfs_validate.py analysis/pfs_735_vetting` after
finalization to check all PDFs for explicit target names, accepted-peak
thresholds, component counts, scan/history alignment and completion flags.
This checks explicit-name removal; it cannot prevent identifying a familiar
system from its measured periods.

## Analytic derivative option

`fit_planets(..., analytic_jac=True)` computes the exact derivative of the
least-squares projection over the Keplerian amplitudes, including both
projection-derivative terms. It retains the same likelihood, noise transform,
rank cutoff, parameter bounds and convergence tolerances. The default remains
finite differences. Regression checks compare the residuals and derivatives
with the existing implementation and verify blind recovery and resumption.

For expensive targets, an independent accelerated attempt can start from a
saved accepted model, or from zero components if no model checkpoint exists:

```bash
python Agatha-2.0/scripts/pfs_analytic_retry.py analysis/pfs_735_vetting \
  --codes P0002 P0509 --rscript .envs/agatha/bin/Rscript --workers 2
```

The separate `analytic_retry/` tree contains the results and immutable input
snapshots. Periodogram scans retain their existing implementation. Analytic
derivatives can follow a different optimization trajectory from finite
differences; compare diagnostics and retain the run provenance when selecting
final reports.
