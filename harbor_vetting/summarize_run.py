"""Aggregate completed Harbor trials by candidate count, not host-system count."""
import argparse,json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('job',type=Path);p.add_argument('--output',type=Path);a=p.parse_args()
paths=sorted(a.job.glob('*/verifier/reward.json'));records=[json.loads(p.read_text()) for p in paths]
if not records:raise SystemExit('No completed verifier rewards found')
metrics=['reward','catalogue_accuracy','evidence_score','assessment_accuracy','brier_score']
n=sum(r['candidate_count'] for r in records)
result={'completed_trials':len(records),'candidate_evaluations':n,'aggregation':'candidate-weighted; repeated attempts count separately',
        **{k:sum(r[k]*r['candidate_count'] for r in records)/n for k in metrics}}
meta=a.job/'result.json'
if meta.exists():
 job=json.loads(meta.read_text());result['scheduled_trials']=job['n_total_trials'];result['complete']=len(records)==job['n_total_trials']
 result['warning']=None if result['complete'] else 'Incomplete job: missing/errored trials are not scored; do not compare this aggregate as a complete benchmark.'
text=json.dumps(result,indent=2);print(text)
if a.output:a.output.write_text(text+'\n')
