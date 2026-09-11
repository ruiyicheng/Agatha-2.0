"""One-command RV batch, non-MCMC tests, runtime estimates and coded PDF."""
import argparse,os,subprocess,sys,time,json,concurrent.futures
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('dataset',type=Path);p.add_argument('output',type=Path);p.add_argument('--rscript',type=Path,required=True);p.add_argument('--workers',type=int,default=6);a=p.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=True);r=a.rscript.resolve();env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
    def command(args):subprocess.run(args,cwd=ROOT,env=env,check=True)
    command([sys.executable,'scripts/prepare_vetting.py',str(a.dataset.resolve()),str(out)])
    logdir=out/'test_logs';logdir.mkdir(exist_ok=True)
    tests=['test_fourier_kepler','test_gp_periodogram','test_moving_periodogram','test_multiset_periodogram','test_plots','test_vetting_regressions','test_vetting_python']
    def test(name):
        cmd=[sys.executable,'tests/test_vetting_python.py',str(r)] if name.endswith('python') else [str(r),'tests/'+name+'.R']
        start=time.perf_counter()
        with (logdir/(name+'.log')).open('w') as f:proc=subprocess.run(cmd,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
        row=dict(test=name,status='PASS' if proc.returncode==0 else 'FAIL',seconds=time.perf_counter()-start,exit_code=proc.returncode)
        (logdir/(name+'.json')).write_text(json.dumps(row));return row
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
        checks=[pool.submit(test,name) for name in tests]
        command([sys.executable,'scripts/run_vetting.py',str(out),'--rscript',str(r),'--estimate-only'])
        command([sys.executable,'scripts/run_vetting.py',str(out),'--rscript',str(r),'--workers',str(a.workers)])
        status=[f.result() for f in checks]
    status += [dict(test=name,status='SKIP (MCMC)',seconds=0) for name in ['test_ptmcmc','test_multiset_mcmc']]
    (out/'test_status.json').write_text(json.dumps(status,indent=2))
    command([sys.executable,'scripts/report_vetting.py',str(out)])
    if any(r['status']=='FAIL' for r in status):raise SystemExit('The report was generated, but one or more software checks failed; inspect test_status.json.')
    print(out/'agatha_rv_report_anonymized.pdf')
if __name__=='__main__':main()
