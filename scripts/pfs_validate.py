"""Audit batch coverage, stopping rules and explicit-name removal from PDFs."""
import argparse
import concurrent.futures
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from pypdf import PdfReader


def validate(item):
    out, code, target = item
    dest = out / 'results' / code
    summary = json.loads((dest / 'summary.json').read_text())
    errors = []
    reader = PdfReader(out / 'pdf' / (code + '.pdf'))
    text = '\n'.join(page.extract_text() or '' for page in reader.pages)
    if code not in text:
        errors.append('coded identifier missing')
    normalized = re.sub(r'[^A-Z0-9]', '', text.upper())
    name = re.sub(r'[^A-Z0-9]', '', target.upper())
    if name and name in normalized:
        errors.append('explicit target name found in PDF')
    if target.lower() in str(reader.metadata).lower():
        errors.append('explicit target name found in metadata')
    if summary['status'] in {'complete', 'search_limited'}:
        history = pd.read_csv(dest / 'search_history.csv')
        accepted = history[history.action == 'accepted']
        scans = np.load(dest / 'search_periodograms.npz')
        if len(accepted) != summary['n_signals']:
            errors.append('accepted component count mismatch')
        if (accepted.peak_lnBF < 5).any():
            errors.append('accepted search peak below threshold')
        if scans['power'].shape != (len(history), len(scans['frequency'])):
            errors.append('scan/history shape mismatch')
        if not np.isclose(history.iloc[-1].peak_lnBF, summary['last_peak_lnBF']):
            errors.append('final peak mismatch')
        complete = summary['status'] == 'complete'
        if complete != summary['search_complete']:
            errors.append('completion flag mismatch')
        if complete and not (summary['last_peak_lnBF'] < 5 and
                             summary['stop_reason'] == 'lnBF_below_5'):
            errors.append('invalid threshold completion')
        fit = pd.read_csv(dest / 'fit.csv')
        if len(fit) != summary['n'] or not np.isfinite(fit.residual).all():
            errors.append('invalid fitted data')
    return dict(code=code, pages=len(reader.pages), status='PASS' if not errors else 'FAIL',
                errors=errors)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    out = args.output.resolve()
    key = pd.read_csv(out / 'private_identity_key.csv')
    items = [(out, row.code, row.target) for row in key.itertuples()]
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(validate, items))
    result = dict(targets=len(results), pages=sum(r['pages'] for r in results),
                  passed=sum(r['status'] == 'PASS' for r in results),
                  failed=sum(r['status'] == 'FAIL' for r in results), results=results)
    (out / 'validation.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'results'}))
    if result['failed']:
        raise SystemExit('Report audit failures; inspect validation.json')


if __name__ == '__main__':
    main()
