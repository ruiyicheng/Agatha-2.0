"""Use a validated analytic Jacobian for expensive RV searches and diagnostics."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse
import concurrent.futures
import functools
import json
import shutil
import time
from pathlib import Path

import pandas as pd
import pfs_search
from vetting_core import fit_planets


def run(row,out,rscript,source,budget):
    pfs_search.fit_planets=functools.partial(fit_planets,analytic_jac=True)
    row['optimizer_jacobian']='analytic variable projection'
    return pfs_search.run_target(row,out,rscript,source,budget,True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('--codes',nargs='+',required=True)
    parser.add_argument('--rscript',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--fit-budget',type=int,default=1000)
    args=parser.parse_args();out=args.output.resolve();dest=out/'analytic_retry';dest.mkdir(exist_ok=True)
    tasks=pd.read_csv(out/'tasks.csv').set_index('code');work=[]
    for code in args.codes:
        row=tasks.loc[code].to_dict();row['code']=code;choices=[]
        for tree in ['results','numerical_retry/results','extended_retry/results']:
            source=out/tree/code
            if (source/'model_parameters.json').exists() and (source/'search_history.csv').exists():
                n=len(json.loads((source/'model_parameters.json').read_text())['parameters'])
                choices.append((n,source))
        source=None
        if choices:
            _,original=max(choices,key=lambda item:item[0])
            source=dest/'snapshots'/code
            if source.exists():raise SystemExit(f'Snapshot exists: {source}; use a fresh retry directory.')
            shutil.copytree(original,source)
            if not (source/'summary.json').exists():
                progress=json.loads((original/'progress.json').read_text())
                elapsed=progress['elapsed_seconds']+time.time()-(original/'progress.json').stat().st_mtime
                snapshot=dict(row,n_signals=len(json.loads((source/'model_parameters.json').read_text())['parameters']),elapsed_seconds=elapsed,status='accepted_model_checkpoint')
                pfs_search.dump(source/'summary.json',snapshot)
            pfs_search.dump(source/'checkpoint_provenance.json',dict(source=str(original),copied_at_unix=time.time()))
        work.append((row,source))
    results=[]
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(run,row,dest,args.rscript.resolve(),source,args.fit_budget) for row,source in work]
        for future in concurrent.futures.as_completed(futures):
            result=future.result();results.append(result)
            pfs_search.dump(dest/'status.json',dict(requested=len(work),finished=len(results),complete=sum(r['status']=='complete' for r in results)))
            pd.DataFrame(results).to_csv(dest/'summary.csv',index=False)
            print(json.dumps({k:result.get(k) for k in ['code','status','n_signals','last_peak_lnBF','elapsed_seconds']}),flush=True)


if __name__=='__main__':main()
