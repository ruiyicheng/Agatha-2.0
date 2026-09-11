# Headless RV evaluation

This extension runs instrument noise selection, a fixed-count multi-Keplerian
fit, RV/activity/window/moving diagnostics, all non-MCMC repository tests, and
a coded PDF. It uses the supplied validation/test manifests and excludes
photometry. It does not run chains or produce posterior uncertainties.

## Run

```bash
conda env create -f environment-vetting.yml
conda activate agatha-vetting
Rscript -e 'install.packages(c("JPEN", "ramify"), repos="https://cloud.r-project.org")'
python scripts/headless_vetting.py /path/to/rv_vetting_9plus9 /path/to/output \
    --rscript "$(command -v Rscript)" --workers 6
```

The input directory must contain `host_manifest.csv`, `val.csv`, `test.csv`,
and the RV files referenced by the host manifest. The output directory should
be **outside this repository**. Raw data and identity keys are not published.
The script exits unsuccessfully if software checks fail, while retaining their
logs and a PDF that records those failures. An incomplete target prevents a
complete PDF. Individual stages can also be run:

```bash
python scripts/prepare_vetting.py DATASET OUTPUT
python scripts/run_vetting.py OUTPUT --rscript /path/to/Rscript --estimate-only
python scripts/run_vetting.py OUTPUT --rscript /path/to/Rscript --workers 6
# With test_status.json present:
python scripts/report_vetting.py OUTPUT
```

The estimate stage fits initial instrument noise models if needed. `--codes
T01 V01` restricts the science stage for a targeted rerun. Prepared inputs are
recreated from source files and the source file mapping remains in the audit.

## Statistical scope

* Instrument noise grid: WN plus jitter, AR1, MA1, MA2, selected by minimum BIC
  using `sopt` and three starts. Instruments with <=10 rows use WN only.
* Fit a preliminary model, repeat noise selection after removing its planet
  signal, and fit the final model conditional on the selected noise parameters.
  This is a two-pass approximation, not global joint noise/planet inference.
* `vetting_core.py` reproduces `CircularSig`: MA filters both data and design;
  AR subtracts lagged observed RV. The latter can absorb real signals and make
  amplitudes inconsistent across instruments. The software equivalence tests
  establish implementation agreement, not scientific validity of that choice.
* Every full model contains exactly the number of **selected** confirmed plus
  retracted entries for that system. Catalogue periods, but never truth labels,
  initialize the fits. Known periods vary within +/-10%; e <=0.85. A missing
  period is initialized by a residual scan. This is catalogue-informed vetting,
  not a blind discovery or held-out classification benchmark.
* Simultaneous independent Keplerians share instrument offsets and a linear
  trend. Three starts and a retry after an evaluation-limit exit are used.
  Global optimality is not guaranteed. Strongly interacting systems may require
  dynamical models even if all numerical fits converge.
* Frequency spacing <=1/(4 baseline); range 1 day to max(2 baseline, 1.25 longest
  catalogue period). Scans stream blocks of frequencies to avoid Agatha's dense
  observation-by-frequency result arrays. Conditional sinusoid evidence is
  delta ln L - ln N, with nuisance coefficients refit at each frequency. It is
  not a calibrated FAP or exact marginalized evidence.
* The sequential RV scans subtract the fitted components in descending K.
  Drop-component BIC reoptimizes remaining planets with fixed noise. Moving
  windows span half the baseline, with seven centers and 129 frequencies near
  each fitted period. Instrument-removal tests keep other components fixed.
* S-index and H-alpha are scanned separately per instrument with >10 valid,
  nonconstant values, standardized with unit errors. Missing activity is not
  imputed. Spectral windows, residual/activity Spearman correlations, and
  approximate ordered-residual Ljung–Box tests are included. Exploratory
  p-values are unadjusted; irregular cadence limits their interpretation.
* GP is covered by the software tests but not the science noise-selection grid.

## Outputs and privacy

`agatha_rv_report_anonymized.pdf` contains coded hosts/candidates, relative times,
and instrument names. The generated text is checked against known target and
planet names and source paths. Individual catalogue labels do not appear.
Periods and RV patterns may nevertheless identify familiar systems: this is
pseudonymization, not a guarantee against re-identification.

`private_identity_key.csv`, `private_aliases.csv`, `private_candidates.csv`,
`runtime_by_target_PRIVATE.csv`, and `preprocessing_audit.csv` contain identities
or source paths. Keep them separate from a blinded reviewer. Other top-level
summary CSVs and per-target result files use coded identities. `prepared/`
retains absolute times and is not part of the anonymized report.

Each target has noise comparisons, fitted RVs/components/residuals, candidate
parameters, drop-model evidence, activity and instrument spectra, sampling
window, moving spectra and instrument-removal sensitivity tables. The summary
records convergence and actual timing. Runtime estimates are recorded before
science runs; actual science timing excludes software installation, tests,
preprocessing, queuing and PDF assembly.

Input conversion preserves first rows of headerless files, logs proxy mappings,
uses the dataset's 2,400,000 reduced-JD convention, and removes duplicate full
series and repeated instrument-family epochs. No nightly binning is performed.
Instrument families prefer pre/post-upgrade files over overlapping unsplit
records. RV/error units are assumed to follow the input m/s convention.

## Tests and upstream fixes

The five original deterministic test scripts are run unchanged. The two
MCMC-specific scripts are skipped. Added regression tests cover:

* R/Python likelihood and periodogram agreement for all four noise models.
* Recovery of two injected eccentric planets in a simultaneous model.
* SHO covariance continuity at Q=0.5 and finite long-lag overdamped values.
* Preservation of missing/constant proxies and simultaneous distinct RV rows.

Upstream reader fixes remove synthetic random activity replacement and
restrict its automatic deduplication to identical complete rows. The SHO
critical-damping formula and overdamped numerical stability are corrected in
both independent kernel implementations. No MCMC algorithm is changed.
