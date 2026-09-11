"""Promote improved numerical retries and rebuild the complete coded report index."""
import argparse
import json
import shutil
from pathlib import Path

import pandas as pd

RETRY_TREES = ['numerical_retry', 'extended_retry', 'checkpoint_retry', 'analytic_retry', 'analytic_extended']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    out = args.output.resolve()
    batch = json.loads((out / 'batch_status.json').read_text())
    if batch['reported'] != batch['requested']:
        raise SystemExit('Wait for all first-pass reports before finalizing.')
    for path in [out / 'report_refresh_status.json', out / 'numerical_retry/retry_status.json']:
        if path.exists() and json.loads(path.read_text()).get('pending', 0):
            raise SystemExit(f'Wait for pending work: {path}')
    promotions = []
    choices = {}
    for tree in RETRY_TREES:
        for source in sorted((out / tree / 'results').glob('P*')):
            path = source / 'summary.json'
            if not path.exists():
                raise SystemExit(f'Retry still running or missing summary: {source}')
            result = json.loads(path.read_text())
            score = (result['status'] == 'complete', result['status'] != 'failed', result['n_signals'])
            if source.name in choices and choices[source.name][0][0]:
                continue  # Keep the first completed search, regardless of component count.
            if source.name not in choices or score > choices[source.name][0]:
                choices[source.name] = (score, source)
    for _, source in choices.values():
        summary_path = source / 'summary.json'
        retry_pdf = source.parent.parent / 'pdf' / (source.name + '.pdf')
        if not summary_path.exists() or not retry_pdf.exists():
            continue
        retry = json.loads(summary_path.read_text())
        dest = out / 'results' / source.name
        original = json.loads((dest / 'summary.json').read_text())
        if original.get('retry_promoted') or original['status'] == 'complete':
            continue
        improved = retry['status'] == 'complete' or (
            retry['status'] == 'search_limited' and retry['n_signals'] > original['n_signals'])
        if not improved:
            continue
        archive = out / 'first_pass_archive'
        (archive / 'results').mkdir(parents=True, exist_ok=True)
        (archive / 'pdf').mkdir(exist_ok=True)
        if (archive / 'results' / source.name).exists():
            raise SystemExit(f'Archive already exists: {source.name}; inspect before proceeding.')
        shutil.move(str(dest), archive / 'results' / source.name)
        shutil.move(str(out / 'pdf' / retry_pdf.name), archive / 'pdf' / retry_pdf.name)
        shutil.copytree(source, dest)
        shutil.copy2(retry_pdf, out / 'pdf' / retry_pdf.name)
        retry.update(retry_promoted=True, pdf=str(out / 'pdf' / retry_pdf.name))
        (dest / 'summary.json').write_text(json.dumps(retry, indent=2) + '\n')
        promotions.append(dict(code=source.name, before=original['status'], after=retry['status'],
                               before_signals=original['n_signals'], after_signals=retry['n_signals']))
    tasks = pd.read_csv(out / 'tasks.csv')
    rows = []
    for code in tasks.code:
        result = json.loads((out / 'results' / code / 'summary.json').read_text())
        if not (out / 'pdf' / (code + '.pdf')).exists():
            raise SystemExit(f'Missing PDF: {code}')
        base = out / 'first_pass_archive/results' / code / 'summary.json'
        seconds = json.loads(base.read_text())['elapsed_seconds'] if base.exists() else result['elapsed_seconds']
        for tree in RETRY_TREES:
            attempt = out / tree / 'results' / code / 'summary.json'
            if attempt.exists():
                seconds += json.loads(attempt.read_text())['elapsed_seconds']
        result['total_elapsed_seconds'] = seconds
        (out / 'results' / code / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
        rows.append(result)
    frame = pd.DataFrame(rows).sort_values('code')
    frame.to_csv(out / 'run_summary.csv', index=False)
    cols = ['code', 'n_raw', 'n', 'pfs_n', 'n_instruments', 'n_signals', 'status',
            'stop_reason', 'last_peak_lnBF', 'total_elapsed_seconds']
    index = frame.reindex(columns=cols).copy()
    index['pdf'] = index.code.map(lambda code: f'pdf/{code}.pdf')
    index.to_csv(out / 'report_index.csv', index=False)
    index[index.status != 'complete'].to_csv(out / 'search_exceptions.csv', index=False)
    batch.update(reported=len(frame), status_counts=frame.status.value_counts().to_dict(),
                 finalized=True, accepted_signals=int(frame.n_signals.sum()))
    (out / 'batch_status.json').write_text(json.dumps(batch, indent=2) + '\n')
    if promotions:
        (out / 'retry_promotions.json').write_text(json.dumps(promotions, indent=2) + '\n')
    print(json.dumps(batch, indent=2))


if __name__ == '__main__':
    main()
