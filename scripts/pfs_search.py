"""Blind iterative RV search for every audited PFS target; no catalogue priors."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse, concurrent.futures, csv, hashlib, json, math, re, time, traceback
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks
from scipy.stats import spearmanr
from vetting_core import ConditionalRV, fit_planets
from run_vetting import load, chosen, model_for, wn, noise_run

def dump(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(value,indent=2,default=str)+'\n');tmp.replace(path)

def read_table(path):
    with open(path,newline='') as stream:return list(csv.DictReader(stream))

def number(value):
    try:return float(value.replace('D','E'))
    except (ValueError,AttributeError):return math.nan

def parse_source(path):
    """Keep real activity columns and reject invalid RV/error rows."""
    header=[];rows=[]
    for line in path.read_text(errors='replace').splitlines():
        tokens=line.strip().lstrip('#').split()
        if len(tokens)<3:continue
        if not line.lstrip().startswith('#') and all(math.isfinite(number(x)) for x in tokens[:3]):rows.append(tokens)
        elif not rows and any(re.sub('[^a-z0-9]','',t.lower()) in {'rv','vrad','vel','rvel','rvs'} for t in tokens):header=tokens
    if not rows:return None,{'reason':'no numeric RV rows'}
    width=max(map(len,rows));d=np.array([[number(v) for v in row]+[math.nan]*(width-len(row)) for row in rows])
    raw=len(d);d=d[np.all(np.isfinite(d[:,:3]),axis=1)&(d[:,2]>0)]
    if not len(d):return None,{'reason':'no positive finite RV uncertainties','raw_rows':raw}
    ins=path.stem.rsplit('_',1)[-1].upper();names=[re.sub('[^a-z0-9]','',v.lower()) for v in header]
    offset=0.
    if np.median(d[:,0])<2400000:offset=2400000.5 if names and names[0]=='mjd' else 2400000.
    d[:,0]+=offset;sidx=hidx=None;schema='unidentified activity columns omitted'
    if header:
        sidx=next((i for i,s in enumerate(names) if s in {'sindex','smw','shk','s'}),None)
        hidx=next((i for i,s in enumerate(names) if s in {'halpha','ha','iha'}),None);schema='named activity columns'
    elif ins=='PFS' and width==8:
        if np.allclose(d[:,3],1):sidx=4;schema='PFS constant flag: S column 5, H-alpha unidentified'
        else:sidx,hidx=3,4;schema='PFS activity format: S column 4, H-alpha column 5'
    elif ins=='KECK' and width==7:sidx,hidx=3,4;schema='KECK S column 4, H-alpha column 5'
    values={}
    for key,index in [('sindex',sidx),('halpha',hidx)]:
        v=np.full(len(d),np.nan) if index is None or index>=width else d[:,index].copy()
        v[(v<=0)|(v>=99)|~np.isfinite(v)]=np.nan;values[key]=v
    frame=pd.DataFrame(dict(time=d[:,0],rv=d[:,1],error=d[:,2],**values)).sort_values('time').reset_index(drop=True)
    return frame,dict(raw_rows=raw,valid_rows=len(d),jd_offset=offset,activity_schema=schema)

def prepare(root,target_csv,audit_csv,out):
    if (out/'tasks.csv').exists():return pd.read_csv(out/'tasks.csv')
    targets=read_table(target_csv);audit=read_table(audit_csv);by_target=defaultdict(list)
    for row in audit:
        if not row['excluded_reason'] and int(row['rv_points'])>0:by_target[row['target']].append(row)
    identities=[];tasks=[];records=[]
    ordered=sorted(targets,key=lambda r:hashlib.sha256(('pfs-search-v1:'+r['target']).encode()).hexdigest())
    for i,row in enumerate(ordered,1):
        code=f'P{i:04d}';name=row['target'];identities.append(dict(code=code,target=name))
        files=sorted(by_target[name],key=lambda r:('PRE' not in r['instrument'] and 'POST' not in r['instrument'],r['file']))
        seen=defaultdict(set);signatures=set();instruments=[];raw_count=0
        for entry in files:
            src=root/entry['file'];ins=entry['instrument'];frame,record=parse_source(src)
            record.update(code=code,source=str(src),instrument=ins,used_rows=0)
            raw_count+=record.get('raw_rows',0)
            if frame is None:records.append(record);continue
            signature=hashlib.sha256(np.round(frame[['time','rv','error']].to_numpy(),8).tobytes()).hexdigest()
            if signature in signatures:record['reason']='duplicate complete RV series';records.append(record);continue
            signatures.add(signature)
            family={'HIRES':'KECK','KECKPRE':'KECK','KECKPOST':'KECK','HARPSPRE':'HARPS','HARPSPOST':'HARPS','HARPSN':'HARPN','HARPS-N':'HARPN'}.get(ins,ins)
            keep=[]
            for j,t in enumerate(frame.time):
                epoch=round(t,7)
                if epoch not in seen[family]:keep.append(j);seen[family].add(epoch)
            frame=frame.iloc[keep].reset_index(drop=True);record['used_rows']=len(frame)
            record['reason']=f"{record['valid_rows']-len(frame)} duplicate instrument epochs removed"
            records.append(record)
            if not len(frame):continue
            sid=f'I{len(instruments)+1:02d}_{ins}';dest=out/'prepared'/code/(sid+'.csv');dest.parent.mkdir(parents=True,exist_ok=True);frame.to_csv(dest,index=False)
            instruments.append(dict(set_id=sid,instrument=ins,file=str(dest),n=len(frame)))
        directory=out/'prepared'/code;directory.mkdir(parents=True,exist_ok=True)
        pd.DataFrame(instruments,columns=['set_id','instrument','file','n']).to_csv(directory/'instruments.csv',index=False)
        tasks.append(dict(code=code,n_raw=raw_count,n=sum(v['n'] for v in instruments),n_instruments=len(instruments),pfs_n=sum(v['n'] for v in instruments if v['instrument']=='PFS'),input_dir=str(directory)))
    pd.DataFrame(identities).to_csv(out/'private_identity_key.csv',index=False)
    pd.DataFrame(records).to_csv(out/'preprocessing_audit.PRIVATE.csv',index=False)
    pd.DataFrame(tasks).to_csv(out/'tasks.csv',index=False)
    dump(out/'preparation_summary.json',dict(targets=len(tasks),raw_rv_rows=sum(t['n_raw'] for t in tasks),used_rv_rows=sum(t['n'] for t in tasks),instruments=sum(t['n_instruments'] for t in tasks),pfs_rv_rows=sum(t['pfs_n'] for t in tasks),aliases_merged=False))
    return pd.DataFrame(tasks)

def frequency_grid(d):
    span=float(np.ptp(d.t));pmax=max(2*span,10.)
    return np.linspace(1/pmax,1.,max(2000,int(np.ceil((1-1/pmax)*span*4))+1))

def strongest(model,frequency,signal):
    power=model.scan(frequency,signal);finite=np.isfinite(power)
    if not finite.any():raise RuntimeError('No finite periodogram values')
    index=int(np.nanargmax(power));f=float(frequency[index]);value=float(power[index])
    if 0<index<len(frequency)-1:
        opt=minimize_scalar(lambda x:-float(model.scan(np.array([x]),signal)[0]),bounds=(frequency[index-1],frequency[index+1]),method='bounded',options={'xatol':min(1e-10,(frequency[1]-frequency[0])*1e-4)})
        if opt.success and -opt.fun>value:f=float(opt.x);value=float(-opt.fun)
    return power,1/f,value

def search(model,frequency,update_noise=None,threshold=5.,checkpoint=None,initial_fit=None,fit_budget=300,warm_starts=False):
    """Add the strongest residual signal until the periodogram threshold fails.

    Each scan conditions on that stage's noise and jointly fitted prior signals.
    Non-threshold stops are explicit; they are never reported as non-detections.
    """
    fit=fit_planets(model,[]) if initial_fit is None else initial_fit;stages=[];powers=[];total_evals=0
    while True:
        power,period,value=strongest(model,frequency,fit['signal']);powers.append(power)
        k=len(fit['parameters']);row=dict(stage=k,n_existing_signals=k,peak_period_days=period,peak_lnBF=value,action='pending')
        stages.append(row)
        if checkpoint:checkpoint(k,row)
        if value<threshold:row['action']='stop_below_threshold';reason='lnBF_below_5';break
        if len(model.y)-(len(model.levels)+1)-5*(k+1)<10:
            row['action']='stop_insufficient_degrees_of_freedom';reason='insufficient_degrees_of_freedom';break
        previous=fit['parameters'];periods=np.r_[previous[:,0],period]
        initial=np.vstack([previous,[period,.01,0.]])
        proposed=fit_planets(model,periods,starts=3,max_nfev=fit_budget,initial_parameters=initial,polish_nfev=3*fit_budget,warm_starts=warm_starts)
        total_evals+=proposed['nfev']
        if not proposed['converged'] or not np.isfinite(proposed['chi2']):
            row['action']='stop_joint_fit_nonconvergence';reason='joint_fit_nonconvergence';break
        if proposed['chi2']>=fit['chi2']-1e-6:
            row['action']='stop_no_joint_fit_improvement';reason='no_joint_fit_improvement';break
        old_model=model
        if update_noise is not None:
            model=update_noise(proposed['signal'],k+1)
            proposed=fit_planets(model,proposed['parameters'][:,0],starts=2,max_nfev=fit_budget,initial_parameters=proposed['parameters'],polish_nfev=3*fit_budget,warm_starts=warm_starts)
            total_evals+=proposed['nfev']
            if not proposed['converged']:
                model=old_model;row['action']='stop_noise_refit_nonconvergence';reason='noise_refit_nonconvergence';break
        row['action']='accepted';row['accepted_fitted_period_days']=float(proposed['parameters'][-1,0]);fit=proposed
    return model,fit,stages,np.asarray(powers),reason,total_evals

def diagnostics(d,ins,noise,model,fit,frequency,dest,code):
    activity=[];spectra={};instrument_spectra={};matches=[];moving=[];cross=[];candidates=[];lagged=[]
    span=float(np.ptp(d.t));n=len(fit['parameters'])
    for instrument in ins.itertuples():
        part=d[d.set_id==instrument.set_id]
        if len(part)>10:instrument_spectra[instrument.set_id]=model_for(part,noise).scan(frequency)
        for proxy in ['sindex','halpha']:
            q=part[np.isfinite(part[proxy])].copy();available=len(q)>10 and q[proxy].std()>0
            rec=dict(set_id=instrument.set_id,proxy=proxy,n_valid=len(q),scan_available=bool(available),global_peak_period_days=None,global_peak_lnBF=None)
            if available:
                q['y']=(q[proxy]-q[proxy].mean())/q[proxy].std();q['error']=1.;no=wn(noise);no['jitter']=0
                power=model_for(q,no).scan(frequency);spectra[instrument.set_id+'__'+proxy]=power;imax=int(np.argmax(power));rec.update(global_peak_period_days=float(1/frequency[imax]),global_peak_lnBF=float(power[imax]))
                peaks=find_peaks(power)[0]
                for j,par in enumerate(fit['parameters']):
                    near=peaks[abs(1/frequency[peaks]-par[0])/par[0]<=.05]
                    if len(near):
                        idx=near[np.argmax(power[near])];matches.append(dict(candidate_id=f'{code}-S{j+1:02d}',set_id=instrument.set_id,proxy=proxy,peak_period_days=float(1/frequency[idx]),proxy_lnBF=float(power[idx]),relative_period_error=float(abs(1/frequency[idx]-par[0])/par[0])))
                    # Explore +/-10-day lags without interpolating across long gaps.
                    series=(fit['residual']+fit['components'][:,j])[q.index];tt=q.t.to_numpy();proxy_values=q[proxy].to_numpy();lag_options=[]
                    for lag in range(-10,11):
                        query=tt-lag;right=np.clip(np.searchsorted(tt,query),0,len(tt)-1);left=np.maximum(0,right-1)
                        nearest=np.where(abs(tt[left]-query)<abs(tt[right]-query),left,right);valid=abs(tt[nearest]-query)<=2
                        if valid.sum()>10 and np.ptp(series[valid])>0 and np.ptp(proxy_values[nearest[valid]])>0:
                            corr=float(spearmanr(series[valid],proxy_values[nearest[valid]]).statistic)
                            if np.isfinite(corr):lag_options.append((abs(corr),lag,corr,int(valid.sum())))
                    if lag_options:
                        _,lag,corr,count=max(lag_options);lagged.append(dict(candidate_id=f'{code}-S{j+1:02d}',set_id=instrument.set_id,proxy=proxy,best_lag_days=lag,spearman_r=corr,n_pairs=count))
            activity.append(rec)
    for j,(period,e,phase) in enumerate(fit['parameters']):
        cid=f'{code}-S{j+1:02d}';other=fit['signal']-fit['components'][:,j];remaining=[k for k in range(n) if k!=j]
        reduced=fit_planets(model,fit['parameters'][remaining,0],starts=2,max_nfev=250,initial_parameters=fit['parameters'][remaining],polish_nfev=600)
        delta=float(reduced['chi2']-fit['chi2']-5*np.log(len(d)));at=float(model.scan(np.array([1/period]),other)[0])
        width=span*.5;centers=np.linspace(width/2,span-width/2,7);freq=np.linspace(1/(period*1.2),1/(period*.8),129);support=[]
        for window,center in enumerate(centers):
            keep=(abs(d.t-center)<=width/2).to_numpy();part=d[keep]
            if len(part)<max(20,2*part.set_id.nunique()+5):continue
            mf=model_for(part,noise);power=mf.scan(freq,other[keep]);value=float(mf.scan(np.array([1/period]),other[keep])[0]);support.append(value>5)
            for f,y in zip(freq,power):moving.append(dict(candidate_id=cid,window=window,center_days=center,period_days=1/f,lnBF=float(y),at_fitted_period=value))
        drops=[]
        for sid in ins.set_id:
            keep=(d.set_id!=sid).to_numpy();part=d[keep]
            if len(part)<20:continue
            value=float(model_for(part,noise).scan(np.array([1/period]),other[keep])[0]);drops.append(value);cross.append(dict(candidate_id=cid,excluded_set=sid,lnBF=value))
        flags=[]
        if delta<=10:flags.append('weak conditional inclusion evidence')
        if e>=.8499:flags.append('eccentricity boundary')
        if 'bounds' in fit and (period<=fit['bounds'][0][j]*1.0001 or period>=fit['bounds'][1][j]*.9999):flags.append('period boundary')
        elif 'bounds' not in fit:
            if 'period_boundary_flags' in fit:
                if fit['period_boundary_flags'][j]:flags.append('period boundary')
            else:flags.append('period-bound check unavailable in checkpoint')
        if span/period<2:flags.append('fewer than 2 cycles')
        if any(r['candidate_id']==cid and r['proxy_lnBF']>=5 for r in matches):flags.append('nearby secondary activity peak')
        if support and np.mean(support)<.5:flags.append('time instability')
        if drops and min(drops)<=5:flags.append('instrument sensitivity')
        if not reduced['converged']:flags.append('reduced fit nonconvergence')
        if fit['chi2']/len(d)>2:flags.append('poor residual fit')
        amp=float(np.linalg.norm(fit['amplitude'][2*j:2*j+2]))
        assessment='unsupported' if delta<=10 else 'inconclusive' if flags else 'conditionally_supported'
        candidates.append(dict(candidate_id=cid,period_days=period,K_m_s=amp,eccentricity=e,delta_BIC_drop=delta,conditional_lnBF_at_period=at,cycles=span/period,window_support=float(np.mean(support)) if support else None,min_instrument_removal_lnBF=min(drops) if drops else None,flags='; '.join(flags) or 'no threshold flags',rv_assessment=assessment))
    for name,rows in [('activity',activity),('activity_matches',matches),('lagged_activity',lagged),('candidates',candidates),('moving',moving),('instrument_removal',cross)]:
        dump(dest/(name+'.json'),rows)
        if rows:pd.DataFrame(rows).to_csv(dest/(name+'.csv'),index=False)
    np.savez_compressed(dest/'activity_periodograms.npz',frequency=frequency,**spectra)
    np.savez_compressed(dest/'instrument_periodograms.npz',frequency=frequency,**instrument_spectra)
    return candidates

def run_target(row,out,rscript,resume_from=None,fit_budget=300,warm_starts=False):
    from pfs_report import report
    start=time.perf_counter();task=SimpleNamespace(**row);dest=out/'results'/task.code;dest.mkdir(parents=True,exist_ok=True)
    summary=dict(**row,status='running',n_signals=0,search_complete=False)
    try:
        if task.n<8 or task.n<=task.n_instruments+3:summary.update(status='insufficient_data',stop_reason='insufficient_measurements');return finish(summary,dest,out,start,report)
        d,ins=load(task);span=float(np.ptp(d.t));summary['baseline_days']=span
        if span<=0:summary.update(status='insufficient_data',stop_reason='zero_baseline');return finish(summary,dest,out,start,report)
        frequency=frequency_grid(d);summary['frequencies']=len(frequency);summary['period_range_days']=[1.,float(1/frequency[0])]
        previous_fit=None;previous_stages=[];previous_powers=None
        initial=dest/'noise_stage_00.csv'
        if resume_from is not None:
            import shutil
            from vetting_core import kepler_basis
            resume_from=Path(resume_from);prior=json.loads((resume_from/'summary.json').read_text());accepted=int(prior['n_signals'])
            for path in resume_from.glob('noise_stage_*.csv'):shutil.copy2(path,dest/path.name)
            initial=dest/f'noise_stage_{accepted:02d}.csv'
            summary.update(resumed_from=str(resume_from),prior_elapsed_seconds=prior['elapsed_seconds'],retry_fit_budget=fit_budget,warm_starts=warm_starts)
        elif not initial.exists():noise_run(task,initial,rscript)
        noise=chosen(initial);model=model_for(d,noise)
        if resume_from is not None:
            checkpoint_parameters=json.loads((resume_from/'model_parameters.json').read_text())
            pars=np.asarray(checkpoint_parameters['parameters'],float).reshape(-1,3)
            basis=np.column_stack([kepler_basis(model.t,*par) for par in pars]) if len(pars) else np.empty((len(d),0))
            previous_fit=model.fit_basis(basis);previous_fit.update(parameters=pars,components=np.column_stack([basis[:,2*j:2*j+2]@previous_fit['amplitude'][2*j:2*j+2] for j in range(len(pars))]) if len(pars) else np.empty((len(d),0)),converged=True,nfev=0)
            if checkpoint_parameters.get('bounds') is not None:previous_fit['bounds']=tuple(np.asarray(b,float) for b in checkpoint_parameters['bounds'])
            elif checkpoint_parameters.get('period_boundary_flags') is not None:previous_fit['period_boundary_flags']=checkpoint_parameters['period_boundary_flags']
            elif (resume_from/'candidates.json').exists():
                prior_candidates=json.loads((resume_from/'candidates.json').read_text())
                if len(prior_candidates)==len(pars) and not any('period-bound check unavailable' in r['flags'] for r in prior_candidates):previous_fit['period_boundary_flags']=['period boundary' in r['flags'] for r in prior_candidates]
            previous_stages=pd.read_csv(resume_from/'search_history.csv').iloc[:-1].to_dict('records')
            previous_powers=np.load(resume_from/'search_periodograms.npz')['power'][:-1]
        begin=time.perf_counter();model.scan(frequency[np.linspace(0,len(frequency)-1,min(128,len(frequency)),dtype=int)])
        scan_seconds=(time.perf_counter()-begin)*len(frequency)/min(128,len(frequency))
        summary['estimated_scan_seconds']=scan_seconds
        dump(dest/'progress.json',dict(code=task.code,stage='noise_selected',estimated_scan_seconds=scan_seconds,elapsed_seconds=time.perf_counter()-start))
        def update(signal,k):
            nonlocal noise
            path=dest/f'noise_stage_{k:02d}.csv';signal_path=dest/'noise_subtraction.csv';pd.DataFrame(dict(set_id=d.set_id,signal=signal)).to_csv(signal_path,index=False)
            noise_run(task,path,rscript,signal_path);noise=chosen(path);return model_for(d,noise)
        def checkpoint(k,stage):
            dump(dest/'progress.json',dict(code=task.code,stage='residual_search',n_signals=k,peak_period_days=stage['peak_period_days'],peak_lnBF=stage['peak_lnBF'],elapsed_seconds=time.perf_counter()-start))
        model,fit,stages,powers,reason,evals=search(model,frequency,update,checkpoint=checkpoint,initial_fit=previous_fit,fit_budget=fit_budget,warm_starts=warm_starts)
        if previous_stages:
            stages=previous_stages+stages;powers=np.concatenate([previous_powers,powers])
        # If a new noise refit failed, restore the noise of the last accepted model.
        accepted=len(fit['parameters']);noise=chosen(dest/f'noise_stage_{accepted:02d}.csv')
        noise.to_csv(dest/'selected_noise.csv');pd.DataFrame(stages).to_csv(dest/'search_history.csv',index=False)
        np.savez_compressed(dest/'search_periodograms.npz',frequency=frequency,power=powers)
        data=pd.DataFrame(dict(t_days=d.t,set_id=d.set_id,rv=d.y,error=d.error,effective_error=model.error,model=fit['deterministic'],residual=fit['residual']))
        for j in range(accepted):data[f'component_{j+1}']=fit['components'][:,j]
        data.to_csv(dest/'fit.csv',index=False)
        dump(dest/'model_parameters.json',dict(parameters=fit['parameters'].tolist(),amplitudes=fit['amplitude'].tolist(),bounds=[b.tolist() for b in fit['bounds']] if 'bounds' in fit else None,period_boundary_flags=fit.get('period_boundary_flags')))
        summary.update(n_signals=accepted,stop_reason=reason,search_complete=reason=='lnBF_below_5',last_peak_lnBF=stages[-1]['peak_lnBF'],last_peak_period_days=stages[-1]['peak_period_days'],conditional_chi2=float(fit['chi2']),chi2_per_measurement=float(fit['chi2']/len(d)),optimizer_evaluations=evals)
        dump(dest/'progress.json',dict(code=task.code,stage='diagnostics',n_signals=accepted,elapsed_seconds=time.perf_counter()-start))
        diagnostics(d,ins,noise,model,fit,frequency,dest,task.code)
        summary['status']='complete' if summary['search_complete'] else 'search_limited'
        return finish(summary,dest,out,start,report)
    except Exception as exc:
        (dest/'failure.PRIVATE.txt').write_text(traceback.format_exc())
        summary.update(status='failed',error_type=type(exc).__name__,stop_reason='analysis_error')
        return finish(summary,dest,out,start,report)

def finish(summary,dest,out,start,report):
    summary['elapsed_seconds']=time.perf_counter()-start;dump(dest/'summary.json',summary)
    pdf=out/'pdf'/f"{summary['code']}.pdf";report(dest,pdf)
    summary['pdf']=str(pdf);summary['elapsed_seconds']=time.perf_counter()-start;dump(dest/'summary.json',summary)
    dump(dest/'progress.json',dict(code=summary['code'],stage=summary['status'],n_signals=summary['n_signals'],elapsed_seconds=summary['elapsed_seconds']))
    return summary

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('output',type=Path);p.add_argument('--rscript',type=Path,required=True);p.add_argument('--workers',type=int,default=24);p.add_argument('--codes',nargs='*');p.add_argument('--prepare-only',action='store_true');p.add_argument('--retry-failed',action='store_true');a=p.parse_args()
    root=a.root.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
    tasks=prepare(root,root/'rv_statistics/pfs_targets_all_combined.csv',root/'rv_statistics/rv_file_audit.csv',out)
    dump(out/'run_config.json',dict(target_count=len(tasks),threshold_lnBF=5,period_min_days=1,period_max='max(2*baseline,10 days)',frequency_oversampling=4,noise_models=['WN+jitter','AR1','MA1','MA2'],noise_refit='after every accepted signal',catalogue_priors=False,photometry=False,MCMC=False,workers=a.workers,activity_period_tolerance=.05,component_limit='available degrees of freedom; no arbitrary component cap'))
    if a.prepare_only:print((out/'preparation_summary.json').read_text());return
    if a.codes:tasks=tasks[tasks.code.isin(a.codes)]
    pending=[];completed=[]
    for row in tasks.to_dict('records'):
        path=out/'results'/row['code']/'summary.json';pdf=out/'pdf'/(row['code']+'.pdf')
        if path.exists() and pdf.exists():
            previous=json.loads(path.read_text())
            if not(a.retry_failed and previous['status']=='failed'):completed.append(previous);continue
        pending.append(row)
    def record():
        pd.DataFrame(completed).to_csv(out/'run_summary.csv',index=False)
        dump(out/'batch_status.json',dict(requested=len(tasks),reported=len(completed),status_counts=dict(pd.Series([r['status'] for r in completed],dtype=str).value_counts().items()),pdf_directory=str(out/'pdf')))
    # Start expensive targets early so that a large system does not become the last straggler.
    pending.sort(key=lambda r:r['n'],reverse=True);record()
    with concurrent.futures.ProcessPoolExecutor(max_workers=a.workers) as executor:
        futures={executor.submit(run_target,row,out,a.rscript.resolve()):row['code'] for row in pending}
        for future in concurrent.futures.as_completed(futures):
            code=futures[future]
            try:result=future.result()
            except Exception:
                dest=out/'results'/code;dest.mkdir(parents=True,exist_ok=True);(dest/'worker_failure.PRIVATE.txt').write_text(traceback.format_exc());print(code,'worker failure',flush=True);continue
            completed.append(result);record();print(json.dumps({k:result.get(k) for k in ['code','status','n','n_signals','last_peak_lnBF','elapsed_seconds']}),flush=True)
    record()
    if len(completed)!=len(tasks):raise SystemExit('Missing target reports; inspect worker failures')

if __name__=='__main__':main()
