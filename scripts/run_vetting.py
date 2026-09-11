"""Run headless RV vetting. See HEADLESS_VETTING.md for the statistical scope."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse,json,time,subprocess,concurrent.futures,traceback
from pathlib import Path
import numpy as np
import pandas as pd
from vetting_core import ConditionalRV,fit_planets,kepler_basis

REPO=Path(__file__).resolve().parents[1]
def load(task):
    ins=pd.read_csv(Path(task.input_dir)/'instruments.csv'); frames=[]
    for row in ins.itertuples():
        d=pd.read_csv(row.file);d['rv_center']=d.rv.median();d['y']=d.rv-d.rv.median();d['set_id']=row.set_id;frames.append(d)
    d=pd.concat(frames,ignore_index=True);d=d.sort_values(['set_id','time']).reset_index(drop=True)
    d['t']=d.time-d.time.min();return d,ins

def chosen(file):
    noise=pd.read_csv(file);return noise[noise.selected].set_index('set_id')
def model_for(d,noise):return ConditionalRV(d.t.values,d.y.values,d.error.values,d.set_id.values,noise)
def wn(noise):
    w=noise.copy();w['Nma']=0;w['Nar']=0;return w

def noise_run(task,out,R,signal=None):
    cmd=[str(R),'scripts/select_vetting_noise.R',str(task.input_dir),str(out)]
    if signal is not None:cmd.append(str(signal))
    with open(out.with_suffix('.log'),'w') as f:
        p=subprocess.run(cmd,cwd=REPO,stdout=f,stderr=subprocess.STDOUT)
    if p.returncode:raise RuntimeError('Noise fit failed; see '+str(out.with_suffix('.log')))

def grid(d,task):
    span=np.ptp(d.t);pmax=max(2*span,1.25*task.max_catalogue_period,10)
    # Uniform frequency spacing <= 1/(4*baseline); local nonlinear fits refine peaks.
    return np.linspace(1/pmax,1,max(2000,int(np.ceil((1-1/pmax)*span*4))+1))

def estimate(task,R,out):
    d,ins=load(task);p=Path(task.input_dir)
    noise=chosen(p/'noise_initial.csv');model=model_for(d,noise);freq=grid(d,task)
    begin=time.perf_counter();model.scan(freq[np.linspace(0,len(freq)-1,128,dtype=int)])
    scan=(time.perf_counter()-begin)/128*len(freq)
    noise_s=pd.read_csv(p/'noise_initial.csv').seconds.sum()
    # A planning range, including two model fits, reduced fits and diagnostics.
    base=2*noise_s+scan*(task.n_signals+3)+.00035*len(d)*(task.n_signals**2)*12+25
    return dict(code=task.code,measurements=len(d),components=task.n_signals,frequencies=len(freq),noise_initial_seconds=noise_s,full_scan_seconds=scan,estimated_minutes=base/60,low_minutes=base/120,high_minutes=base/30)

def run(task,R,out,estimates):
    start=time.perf_counter();dest=out/'results'/task.code;dest.mkdir(parents=True,exist_ok=True)
    previous=json.loads((dest/'summary.json').read_text()) if (dest/'summary.json').exists() else {}
    previous_seconds=previous.get('elapsed_seconds_including_retries',previous.get('elapsed_seconds',0))
    d,ins=load(task);candidates=pd.read_csv(out/'private_candidates.csv');cand=candidates[candidates.code==task.code].copy().reset_index(drop=True)
    noise_initial=chosen(Path(task.input_dir)/'noise_initial.csv');frequency=grid(d,task)
    prelim=model_for(d,noise_initial);periods=cand.period_days.to_numpy(float)
    if np.any(~np.isfinite(periods)):
        valid=np.isfinite(periods);knownfit=fit_planets(prelim,periods[valid],starts=2,max_nfev=180)
        power=prelim.scan(frequency,knownfit['signal']);periods[~valid]=1/frequency[np.argmax(power)]
    initialfit=fit_planets(prelim,periods,starts=3,max_nfev=250)
    # Second selection on RV minus the preliminary N-component signal, preventing
    # the noise-only jitter from permanently absorbing large known RV signals.
    pd.DataFrame(dict(set_id=d.set_id,signal=initialfit['signal'])).to_csv(dest/'preliminary_signal.csv',index=False)
    noise_run(task,dest/'noise_refined.csv',R,dest/'preliminary_signal.csv')
    noise=chosen(dest/'noise_refined.csv');model=model_for(d,noise);white=model_for(d,wn(noise))
    fit=fit_planets(model,periods,starts=3,max_nfev=300)
    n=len(periods);components=fit['components'];amp=np.linalg.norm(fit['amplitude'].reshape(-1,2),axis=1)
    order=np.argsort(-amp);sequence=[];sub=np.zeros(len(d));stage_peaks=[]
    for stage in range(n+1):
        power=model.scan(frequency,sub);sequence.append(power)
        stage_peaks.append(dict(stage=stage,peak_period_days=1/frequency[np.argmax(power)],peak_logBF_BIC=float(max(power))))
        if stage<n:sub+=components[:,order[stage]]
    whitepower=white.scan(frequency)
    np.savez_compressed(dest/'periodograms.npz',frequency=frequency,sequence=np.array(sequence),white=whitepower,removal_order=order)
    pd.DataFrame(stage_peaks).to_csv(dest/'stage_peaks.csv',index=False)
    # Spectral window is normalized |sum exp(i omega t)|^2/N^2, not a RV fit.
    window=np.empty(len(frequency))
    for a in range(0,len(frequency),128):
        angle=2*np.pi*d.t.to_numpy()[:,None]*frequency[a:a+128]
        window[a:a+128]=(np.mean(np.cos(angle),axis=0)**2+np.mean(np.sin(angle),axis=0)**2)
    np.savez_compressed(dest/'window.npz',frequency=frequency,power=window)
    activity=[];ac_specs={};inst_specs={}
    for row in ins.itertuples():
        part=d[d.set_id==row.set_id].copy();
        if len(part)>10:
            mf=model_for(part,noise);inst_specs[row.set_id]=mf.scan(frequency)
        for proxy in ['sindex','halpha']:
            q=part[np.isfinite(part[proxy])].copy();avail=len(q)>0 and q[proxy].std()>0
            rec=dict(set_id=row.set_id,proxy=proxy,n_valid=len(q),available=bool(avail),scan_available=False,peak_period_days=np.nan,peak_logBF_BIC=np.nan)
            if avail and len(q)>10:
                q['y']=(q[proxy]-q[proxy].mean())/q[proxy].std();q['error']=1.
                no=wn(noise);no['jitter']=0.;mf=model_for(q,no);power=mf.scan(frequency)
                rec.update(scan_available=True,peak_period_days=1/frequency[np.argmax(power)],peak_logBF_BIC=max(power))
                ac_specs[row.set_id+'__'+proxy]=power
            activity.append(rec)
    pd.DataFrame(activity).to_csv(dest/'activity.csv',index=False)
    np.savez_compressed(dest/'activity_periodograms.npz',frequency=frequency,**ac_specs)
    np.savez_compressed(dest/'instrument_periodograms.npz',frequency=frequency,**inst_specs)
    moving=[];candidate_rows=[];cross=[]
    for j,(P,e,phase) in enumerate(fit['parameters']):
        without=[k for k in range(n) if k!=j]
        reduced=fit_planets(model,fit['parameters'][without,0],starts=2,max_nfev=200)
        delta_bic=float(reduced['chi2']-fit['chi2']-5*np.log(len(d)))
        bounds_hit=bool(P<=periods[j]*.9001 or P>=periods[j]*1.0999 or e>=.8499)
        other=fit['signal']-components[:,j]
        at_period=float(model.scan(np.array([1/P]),other)[0])
        width=np.ptp(d.t)*.5;centers=np.linspace(d.t.min()+width/2,d.t.max()-width/2,7)
        moving_grid=np.linspace(1/(P*1.2),1/(P*.8),129)
        for k,center in enumerate(centers):
            keep=(abs(d.t-center)<=width/2).values;q=d[keep].copy()
            if len(q)<max(20,2*q.set_id.nunique()+5):continue
            mf=model_for(q,noise);power=mf.scan(moving_grid,other[keep]);near=float(mf.scan(np.array([1/P]),other[keep])[0])
            for f,pow in zip(moving_grid,power):moving.append(dict(candidate_id=cand.candidate_id[j],window=k,center_days=center,n=len(q),period_days=1/f,logBF_BIC=pow,at_fitted_period=near))
        for sid in ins.set_id:
            keep=(d.set_id!=sid).values;q=d[keep]
            if len(q)<20:continue
            mf=model_for(q,noise)
            cross.append(dict(candidate_id=cand.candidate_id[j],excluded_set=sid,n=len(q),logBF_BIC_at_period=float(mf.scan(np.array([1/P]),other[keep])[0])))
        candidate_rows.append(dict(candidate_id=cand.candidate_id[j],catalogue_period_days=cand.period_days[j],fitted_period_days=P,K_m_s=amp[j],eccentricity=e,phase_radians=phase,delta_BIC_drop=delta_bic,conditional_logBF_at_period=at_period,optimizer_boundary=bounds_hit,cycles_observed=np.ptp(d.t)/P,unseeded=bool(pd.isna(cand.period_days[j])),reduced_fit_converged=reduced['converged']))
    pd.DataFrame(moving).to_csv(dest/'moving_periodograms.csv',index=False)
    pd.DataFrame(candidate_rows).to_csv(dest/'candidate_results.csv',index=False)
    pd.DataFrame(cross).to_csv(dest/'leave_one_instrument_out.csv',index=False)
    # Only coded IDs and relative times are exported in fit tables.
    data=pd.DataFrame(dict(t_days=d.t,set_id=d.set_id,rv_centered=d.y,error=d.error,effective_error=model.error,model=fit['deterministic'],signal=fit['signal'],residual=fit['residual']))
    for j in range(n):data['component_'+str(j+1)]=components[:,j]
    data.to_csv(dest/'fit.csv',index=False)
    summary=dict(code=task.code,split=task.split,n_raw=int(task.raw_rv_points),n=int(len(d)),n_signals=n,n_instruments=len(ins),baseline_days=float(np.ptp(d.t)),frequencies=len(frequency),period_min_days=1.,period_max_days=float(1/frequency[0]),conditional_chi2=float(fit['chi2']),residual_rms=float(np.std(fit['residual'])),optimizer_converged=fit['converged'],optimizer_evaluations=fit['nfev'],elapsed_seconds=time.perf_counter()-start,estimated_minutes=float(estimates.set_index('code').loc[task.code].estimated_minutes),status='completed')
    summary['previous_attempt_seconds']=previous_seconds
    summary['elapsed_seconds_including_retries']=summary['elapsed_seconds']+previous_seconds
    (dest/'summary.json').write_text(json.dumps(summary,indent=2));return summary

def main():
    a=argparse.ArgumentParser();a.add_argument('output',type=Path);a.add_argument('--rscript',type=Path,required=True);a.add_argument('--workers',type=int,default=6);a.add_argument('--estimate-only',action='store_true');a.add_argument('--codes',nargs='*');args=a.parse_args()
    out=args.output.resolve();tasks=pd.read_csv(out/'tasks.csv');R=args.rscript.resolve()
    if args.codes:tasks=tasks[tasks.code.isin(args.codes)]
    if args.estimate_only:
        for t in tasks.itertuples():
            p=Path(t.input_dir)/'noise_initial.csv'
            if not p.exists():noise_run(t,p,R)
        estimates=[estimate(t,R,out) for t in tasks.itertuples()];pd.DataFrame(estimates).to_csv(out/'runtime_estimates.csv',index=False);print(pd.DataFrame(estimates).to_string(index=False));return
    estimates=pd.read_csv(out/'runtime_estimates.csv')
    errors=[]
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as ex:
        # The top-level run function is picklable; errors are recorded by the parent.
        futures={}
        # pandas namedtuples are not picklable; use SimpleNamespace tasks.
        from types import SimpleNamespace
        for row in tasks.to_dict('records'):
            t=SimpleNamespace(**row);futures[ex.submit(run,t,R,out,estimates)]=t.code
        for f in concurrent.futures.as_completed(futures):
            code=futures[f]
            try:print(json.dumps(f.result()),flush=True)
            except Exception:
                errors.append(code)
                dest=out/'results'/code;dest.mkdir(parents=True,exist_ok=True);(dest/'failure.txt').write_text(traceback.format_exc());print(code,traceback.format_exc(),flush=True)
    if errors:raise SystemExit('Incomplete targets: '+', '.join(errors))
if __name__=='__main__':main()
