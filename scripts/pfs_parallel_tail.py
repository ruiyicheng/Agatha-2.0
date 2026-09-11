"""Finish a saved complete search with parallel component-removal diagnostics."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse
import concurrent.futures
import functools
import json
import math
import shutil
import time
from multiprocessing.pool import ThreadPool
from pathlib import Path

import numpy as np
import pandas as pd
import pfs_search
from vetting_core import ConditionalRV,fit_planets


def initialize(model,parameters):
    global _model,_parameters
    _model=model;_parameters=parameters


def reduced_component(j):
    remaining=[i for i in range(len(_parameters)) if i!=j]
    pars=_parameters[remaining]
    result=fit_planets(_model,pars[:,0],starts=2,max_nfev=250,
                       initial_parameters=pars,polish_nfev=600,analytic_jac=True)
    return dict(component=j,periods=pars[:,0].tolist(),initial_parameters=pars.tolist(),chi2=float(result['chi2']),
                converged=result['converged'],nfev=result['nfev'],parameters=result['parameters'].tolist())


def reduced_fits(model,parameters,workers):
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers,initializer=initialize,
                                                initargs=(model,parameters)) as pool:
        return list(pool.map(reduced_component,range(len(parameters))))


def install_parallel_scans(workers):
    original=ConditionalRV.scan
    def scan(self,frequency,subtract=None,block=128):
        if workers<=1 or len(frequency)<4*block:return original(self,frequency,subtract,block)
        chunk=math.ceil(len(frequency)/(workers*block))*block
        pieces=[frequency[i:i+chunk] for i in range(0,len(frequency),chunk)]
        with ThreadPool(workers) as pool:
            values=pool.map(lambda f:original(self,f,subtract,block),pieces)
        return np.concatenate(values)
    ConditionalRV.scan=scan


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('--code',required=True)
    parser.add_argument('--source-tree',default='analytic_retry/results')
    parser.add_argument('--rscript',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=8)
    parser.add_argument('--scan-workers',type=int,default=4)
    parser.add_argument('--wait',action='store_true')
    args=parser.parse_args();out=args.output.resolve();original=out/args.source_tree/args.code
    while True:
        ready=(original/'model_parameters.json').exists() and (original/'search_history.csv').exists()
        if ready and pd.read_csv(original/'search_history.csv').iloc[-1].action=='stop_below_threshold':break
        if not args.wait or (original/'summary.json').exists():raise SystemExit('Source search has not completed below threshold')
        time.sleep(5)
    dest=out/'parallel_retry';snapshot=dest/'snapshots'/args.code
    if snapshot.exists():raise SystemExit('Snapshot already exists; inspect the previous parallel attempt')
    shutil.copytree(original,snapshot)
    row=pd.read_csv(out/'tasks.csv').set_index('code').loc[args.code].to_dict();row.update(code=args.code,optimizer_jacobian='analytic variable projection',diagnostic_workers=args.workers,scan_workers=args.scan_workers)
    if not (snapshot/'summary.json').exists():
        progress=json.loads((original/'progress.json').read_text());elapsed=progress['elapsed_seconds']+time.time()-(original/'progress.json').stat().st_mtime
        pfs_search.dump(snapshot/'summary.json',dict(row,n_signals=len(json.loads((snapshot/'model_parameters.json').read_text())['parameters']),elapsed_seconds=elapsed,status='accepted_model_checkpoint'))
    pfs_search.dump(snapshot/'checkpoint_provenance.json',dict(source=str(original),copied_at_unix=time.time()))
    pfs_search.fit_planets=functools.partial(fit_planets,analytic_jac=True)
    original_diagnostics=pfs_search.diagnostics
    def diagnostics(d,ins,noise,model,fit,frequency,path,code):
        rows=reduced_fits(model,fit['parameters'],args.workers)
        pfs_search.dump(path/'parallel_reduced_fits.json',rows)
        cache={tuple(np.asarray(r['initial_parameters']).ravel()):r for r in rows};saved=pfs_search.fit_planets
        def cached(model,periods,**kwargs):return cache[tuple(np.asarray(kwargs['initial_parameters']).ravel())]
        pfs_search.fit_planets=cached
        try:return original_diagnostics(d,ins,noise,model,fit,frequency,path,code)
        finally:pfs_search.fit_planets=saved
    pfs_search.diagnostics=diagnostics;install_parallel_scans(args.scan_workers)
    result=pfs_search.run_target(row,dest,args.rscript.resolve(),snapshot,1000,True)
    print(json.dumps({k:result.get(k) for k in ['code','status','n_signals','last_peak_lnBF','elapsed_seconds']}),flush=True)


if __name__=='__main__':main()
