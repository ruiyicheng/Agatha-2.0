"""Build one isolated Harbor PDF-review task per fitted host; preserve splits."""
import argparse,csv,hashlib,json,re,shutil,sys
from pathlib import Path
import pymupdf
import numpy as np
from recovery import match_recovery

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'scripts'))
from export_target_report import export_target
from grader import protocol

BASE='rv-pdf-review-base:1.0.0'
FLAG_MAP={'weak conditional evidence':'weak_evidence','parameter boundary':'parameter_boundary',
          'less than 2 cycles':'short_baseline','activity peak nearby':'activity_overlap',
          'limited time stability':'time_instability','instrument sensitivity':'instrument_sensitive',
          'catalogue period missing':'unknown_period','reduced fit not converged':'nonconvergence',
          'full fit not converged':'nonconvergence','poor full-model residuals':'poor_model_fit'}

def rows(path):
    with open(path,newline='') as f:return list(csv.DictReader(f))
def dump(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def rounded(v,digits):return float(format(float(v),f'.{digits}g'))
def citations(texts,cid):
    result=[]
    for page in (4,6):
        text=' '.join(texts[page-1].split());start=text.find(cid)
        if start<0:raise ValueError(f'{cid} missing from page {page}')
        result.append({'page':page,'quote':text[start:start+min(220,len(text)-start)]})
    return result

def build_one(run,destination,row,private):
    code=row['code'];target=destination/row['split']/('rv-pdf-'+code.lower())
    if target.exists():raise FileExistsError(f'{target} already exists; choose a new output directory')
    env=target/'environment';tests=target/'tests';solution=target/'solution'
    for p in (env,tests,solution):p.mkdir(parents=True)
    texts=export_target(run/'agatha_rv_report_anonymized.pdf',code,env/'report.pdf')
    cands=rows(run/'results'/code/'diagnostic_flags.csv')
    summary=json.loads((run/'results'/code/'summary.json').read_text())
    chosen=[r for r in rows(run/'results'/code/'noise_refined.csv') if r['selected'].lower()=='true']
    # Audit the actual extracted target PDF against every private catalogue alias.
    norm=lambda s:re.sub('[^a-z0-9]','',s.lower())
    actual=norm(' '.join(texts));names=[]
    for x in private:names.append(x['planet_name'])
    for x in rows(run/'private_identity_key.csv'):names.extend([x['target'],x['host_id']])
    if (run/'private_aliases.csv').exists():
        for x in rows(run/'private_aliases.csv'):names.extend(x['aliases'].split(';'))
    leaks=[s for s in names if len(norm(s))>=5 and norm(s) in actual]
    if leaks:raise ValueError('PDF identity leak: '+str(leaks))
    for other in rows(run/'tasks.csv'):
        if other['code']!=code and re.search(r'\b'+re.escape(other['code'])+r'\b',' '.join(texts)):
            raise ValueError('PDF includes another target')
    system=dict(n_rv=summary['n'],n_components=summary['n_signals'],n_instruments=summary['n_instruments'],
                chi2_per_measurement=round(summary['conditional_chi2']/summary['n'],2),
                fit_adequacy='poor' if summary['conditional_chi2']/summary['n']>2 else 'adequate')
    reference=dict(target_id=code,system=system,noise_models={x['set_id']:x['model'] for x in chosen},
                   candidates=[],page_text={str(i+1):t for i,t in enumerate(texts)})
    oracle=dict(schema_version='1.0',target_id=code,system={**system,'noise_models':[dict(instrument=x['set_id'],model=x['model']) for x in chosen]},candidates=[])
    labels={x['candidate_id']:x['catalogue_label'].lower() for x in private if x['code']==code}
    for c in cands:
        cid=c['candidate_id'];flags=sorted({FLAG_MAP[x.strip()] for x in c['review_flags'].split(';') if x.strip() in FLAG_MAP})
        facts=dict(period_days=rounded(c['fitted_period_days'],5),K_m_s=rounded(c['K_m_s'],3),eccentricity=round(float(c['eccentricity']),3),delta_BIC=round(float(c['delta_BIC_drop']),1),cycles_observed=round(float(c['cycles_observed']),1))
        tolerance=dict(period_days=max(1e-6,abs(facts['period_days'])*.00006),K_m_s=max(.001,abs(facts['K_m_s'])*.006),eccentricity=.0006,delta_BIC=.051,cycles_observed=.051)
        reference['candidates'].append(dict(candidate_id=cid,catalogue_label=labels[cid],facts=facts,tolerance=tolerance,flags=flags,rv_assessment=protocol(flags)))
        oracle['candidates'].append(dict(candidate_id=cid,catalogue_prediction=labels[cid],probability_confirmed=float(labels[cid]=='confirmed'),rv_assessment=protocol(flags),**facts,flags=flags,citations=citations(texts,cid),rationale='Reference infrastructure answer: apply the published conditional-RV review protocol. Catalogue truth is supplied only to the oracle; it is not proven by RV evidence.'))
    dump(tests/'reference.json',reference);dump(solution/'answer.json',oracle)
    # Public template contains identifiers only, never source-derived answers.
    template=dict(schema_version='1.0',target_id=code,system={k:None for k in system},candidates=[])
    template['system']['noise_models']=[dict(instrument=x['set_id'],model=None) for x in chosen]
    for c in cands:
        template['candidates'].append({k:(c['candidate_id'] if k=='candidate_id' else [] if k in ['flags','citations'] else None) for k in oracle['candidates'][0]})
    dump(env/'answer_template.json',template)
    for p in (env,tests):shutil.copy2(HERE/'submission.schema.json',p/'submission.schema.json')
    shutil.copy2(HERE/'pdf_tool.py',env/'pdf_tool.py');shutil.copy2(HERE/'grader.py',tests/'grader.py')
    instruction=(HERE/'instruction.md').read_text().replace('__TARGET__',code).replace('__COUNT__',str(len(cands))).replace('__CANDIDATES__',', '.join(c['candidate_id'] for c in cands))
    (target/'instruction.md').write_text(instruction)
    (env/'Dockerfile').write_text(f'''FROM {BASE}
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
RUN useradd --create-home --uid 1000 agent && mkdir -p /workspace/output && chown -R agent:agent /workspace
WORKDIR /workspace
COPY report.pdf submission.schema.json answer_template.json pdf_tool.py /workspace/
RUN chmod 444 /workspace/report.pdf /workspace/submission.schema.json /workspace/answer_template.json /workspace/pdf_tool.py
CMD ["sleep", "infinity"]
''')
    (env/'.dockerignore').write_text('*\n!Dockerfile\n!report.pdf\n!submission.schema.json\n!answer_template.json\n!pdf_tool.py\n')
    (tests/'Dockerfile').write_text(f'''FROM {BASE}
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
COPY test.sh grader.py reference.json submission.schema.json /tests/
RUN chmod 555 /tests/test.sh && mkdir -p /logs/verifier /workspace/output
CMD ["sleep", "infinity"]
''')
    (tests/'test.sh').write_text('#!/bin/bash\nset -euo pipefail\nexec python -I /tests/grader.py\n')
    (solution/'solve.sh').write_text('#!/bin/bash\nset -euo pipefail\nmkdir -p /workspace/output\ncp /solution/answer.json /workspace/output/verdict.json\n')
    (target/'task.toml').write_text(f'''schema_version = "1.4"
artifacts = ["/workspace/output/"]

[task]
name = "rv-vetting/rv-pdf-{code.lower()}"
version = "1.0.0"
description = "Predict candidate catalogue labels and assess the evidence in an anonymized RV PDF."
keywords = ["astronomy", "radial-velocity", "pdf", "scientific-reasoning"]

[metadata]
category = "scientific-reasoning"
split = "{row['split']}"
n_candidates = {len(cands)}

[agent]
timeout_sec = 900.0
user = "agent"

[environment]
network_mode = "no-network"
build_timeout_sec = 600.0
cpus = 2
memory_mb = 2048
storage_mb = 4096

[verifier]
timeout_sec = 60.0
environment_mode = "separate"

[verifier.environment]
network_mode = "no-network"
build_timeout_sec = 600.0
cpus = 1
memory_mb = 512
storage_mb = 2048
''')
    return dict(target_id=code,split=row['split'],task_path=str(target.relative_to(destination)),n_candidates=len(cands),pdf_sha256=hashlib.sha256((env/'report.pdf').read_bytes()).hexdigest(),pages=len(texts))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('run',type=Path);p.add_argument('output',type=Path);p.add_argument('--codes',nargs='*');p.add_argument('--period-tolerance',type=float,default=.05);p.add_argument('--min-logbf',default='5',help='Minimum peak conditional ln BF, or none');p.add_argument('--top-peaks',type=int,default=10);a=p.parse_args()
    run=a.run.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
    tasks=rows(run/'tasks.csv');private=rows(run/'private_candidates.csv')
    if a.codes:tasks=[t for t in tasks if t['code'] in a.codes]
    if not tasks:raise ValueError('No matching fitted targets')
    threshold=None if a.min_logbf.lower()=='none' else float(a.min_logbf)
    eligible=[];audit=[];selection=[]
    for t in tasks:
        code=t['code'];cand=[c for c in private if c['code']==code]
        summary=json.loads((run/'results'/code/'summary.json').read_text())
        with np.load(run/'results'/code/'periodograms.npz') as z:
            ok,details=match_recovery(z['frequency'],z['sequence'],z['removal_order'],cand,summary['baseline_days'],a.period_tolerance,threshold,a.top_peaks)
        audit.extend(dict(target_id=code,split=t['split'],**r) for r in details)
        selection.append(dict(target_id=code,split=t['split'],eligible=ok,n_candidates=len(cand),matched_candidates=sum(r['matched'] for r in details)))
        if ok:eligible.append(t)
    dump(out/'recovery_selection.json',dict(period_tolerance=a.period_tolerance,min_logbf=threshold,top_peaks=a.top_peaks,targets=selection))
    with (out/'recovery_audit.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(audit[0]));w.writeheader();w.writerows(audit)
    if not eligible:raise ValueError('No whole systems satisfy all-candidate recovery; see recovery_audit.csv')
    manifest=[build_one(run,out,t,private) for t in eligible]
    dump(out/'manifest.json',manifest)
    dump(out/'provenance.PRIVATE.json',dict(source_run=str(run),source_pdf_sha256=hashlib.sha256((run/'agatha_rv_report_anonymized.pdf').read_bytes()).hexdigest(),source_label_sha256=hashlib.sha256((run/'private_candidates.csv').read_bytes()).hexdigest(),harbor_version='0.22.0',schema_version='1.4'))
    print(json.dumps(dict(tasks=len(manifest),candidates=sum(t['n_candidates'] for t in manifest),output=str(out)),indent=2))
if __name__=='__main__':main()
