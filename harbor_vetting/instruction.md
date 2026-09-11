# Vet the RV fitting report

Read `/workspace/report.pdf` for the anonymized system **__TARGET__**. Assess
all __COUNT__ candidate components: __CANDIDATES__.

Write `/workspace/output/verdict.json`, following
`/workspace/submission.schema.json`. A blank template is in
`/workspace/answer_template.json`. Use only this PDF. You are reviewing an
existing fit; fitting RV data, generating a new report, running MCMC, and using
photometry are outside this task. Do not look up or infer target identities.

This system passed a label-independent selection: every catalogue period has a
nearby qualifying peak in the cached RV scans. A peak may still be activity,
an alias or a nonplanetary companion. Selection does not establish its label.

There are two distinct objectives:

1. Predict each candidate's catalogue label (`confirmed` or `retracted`) and
   supply `probability_confirmed`. This is a prediction against a frozen
   catalogue snapshot, **not** something the PDF necessarily proves. A real RV
   signal can belong to a retracted planet entry. `confirmed` must be the
   prediction at probability >= 0.5, otherwise use `retracted`.
2. Assess the conditional RV evidence, extract the requested report values,
   and cite the PDF. Follow the reproducible review protocol below. This
   protocol is not a physical definition of a planet.

## Review protocol

Use the report's printed precision. Set `poor_model_fit` when reported
chi-squared / N > 2, and set system `fit_adequacy` to `poor` in that case,
otherwise `adequate` (this term refers only to this diagnostic).

Use these candidate flags when the report indicates them:

- `weak_evidence`: drop-component delta BIC <= 10.
- `parameter_boundary`: a fitted parameter is at its optimizer bound.
- `short_baseline`: fewer than two observed orbital cycles.
- `activity_overlap`: an activity peak is reported near the candidate period.
- `time_instability`: fewer than half of tested windows have ln BF > 5.
- `instrument_sensitive`: minimum instrument-removal ln BF <= 5.
- `unknown_period`: catalogue period was missing and initialized from a scan.
- `nonconvergence`: full or reduced fit is reported not to have converged.
- `poor_model_fit`: the system has chi-squared / N > 2.

Set `rv_assessment` to `unsupported` for weak evidence; otherwise
`inconclusive` if any of the other flags apply; otherwise
`conditionally_supported`. This assessment and catalogue prediction may differ.
When reading threshold values near a rounded boundary, follow explicit report
flags. Do not treat missing diagnostics as negative evidence.

For each candidate, include at least two verbatim quotations from **different
physical PDF pages**, each 20–400 characters, and a short evidence-based
rationale explaining the assessment and uncertainty in the catalogue
prediction. Physical pages 1–6 are numbered in the footer: methods, data/noise,
RV scans, fitted components, activity/window/instrument diagnostics, stability.
Quotes must come from pages 2–6; at least one must contain the candidate ID.
Use the candidate's fit and diagnostic pages.
Do not invent numerical precision, posterior uncertainties or false-alarm rates.

## Tools and scoring

`python /workspace/pdf_tool.py text` extracts page-labelled text.
`python /workspace/pdf_tool.py render --page 4` renders a page to PNG; inspect
plots when useful. `pdftotext` and `pdftoppm` are also installed.

The scalar reward equally weights catalogue-label accuracy and a report-review
score. Review scoring checks system facts (15%), instrument noise choices (15%),
candidate numerical facts (25%), flag F1 (20%), protocol assessments (15%), and
valid page quotations (10%). Numeric comparisons allow the PDF's rounding.
Invalid JSON/schema or missing/duplicate candidates receives zero. Confidence
Brier score is reported separately. Free-text reasoning is retained for human
review; this deterministic verifier does not claim to grade its scientific
quality. Catalogue truth and grading references are unavailable to the agent.
