"""Resume numerically limited searches in a separate output tree; keep the first pass."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse,concurrent.futures,json,time
from pathlib import Path
import pandas as pd
from pfs_search import run_target,dump

RETRY_REASONS={'joint_fit_nonconvergence','noise_refit_nonconvergence','no_joint_fit_improvement'}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path);p.add_argument('--rscript',type=Path,required=True);p.add_argument('--workers',type=int,default=2);p.add_argument('--fit-budget',type=int,default=1000);p.add_argument('--watch',action='store_true');a=p.parse_args()
    out=a.output.resolve();retry=out/'numerical_retry';retry.mkdir(exist_ok=True);tasks=pd.read_csv(out/'tasks.csv').set_index('code');attempted=set();results=[]
    with concurrent.futures.ProcessPoolExecutor(max_workers=a.workers) as pool:
        futures={}
        while True:
            for dest in (out/'results').glob('P*'):
                if dest.name in attempted:continue
                progress=dest/'progress.json';summary=dest/'summary.json'
                if not progress.exists() or not summary.exists() or not (dest/'model_parameters.json').exists():continue
                if json.loads(progress.read_text()).get('stage')!='search_limited':continue
                previous=json.loads(summary.read_text())
                if previous.get('stop_reason') not in RETRY_REASONS:continue
                attempted.add(dest.name);done=retry/'results'/dest.name/'summary.json'
                if done.exists() and (retry/'pdf'/(dest.name+'.pdf')).exists():results.append(json.loads(done.read_text()));continue
                row=tasks.loc[dest.name].to_dict();row['code']=dest.name
                futures[pool.submit(run_target,row,retry,a.rscript.resolve(),dest,a.fit_budget,True)]=dest.name
            for future in list(futures):
                if not future.done():continue
                result=future.result();results.append(result);del futures[future];print(json.dumps({k:result.get(k) for k in ['code','status','n_signals','stop_reason','last_peak_lnBF','elapsed_seconds']}),flush=True)
            dump(retry/'retry_status.json',dict(attempted=len(attempted),finished=len(results),pending=len(futures),complete=sum(r['status']=='complete' for r in results),fit_budget=a.fit_budget,warm_starts=True))
            batch=json.loads((out/'batch_status.json').read_text())
            if not futures and (not a.watch or batch['reported']==batch['requested']):break
            time.sleep(5)
    pd.DataFrame(results).to_csv(retry/'summary.csv',index=False)

if __name__=='__main__':main()
