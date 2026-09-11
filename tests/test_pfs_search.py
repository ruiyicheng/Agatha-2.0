"""Scientific regression checks for blind discovery and explicit stopping states."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from vetting_core import ConditionalRV,kepler_basis,fit_planets,kepler_objective
from pfs_search import search,parse_source

def model(t,y,error=1.):
    noise=pd.DataFrame([dict(set_id='A',jitter=0,Nma=0,Nar=0,tau=1,tauAR=1)]).set_index('set_id')
    return ConditionalRV(t,y,np.full(len(t),error),np.repeat('A',len(t)),noise)

def test_blind_two_signals_and_threshold_stop():
    rng=np.random.default_rng(746);t=np.sort(rng.uniform(0,500,240))
    y=kepler_basis(t,17,.10,.3)@np.array([7.,2.])+kepler_basis(t,43,.08,1.2)@np.array([4.,1.])+rng.normal(0,.4,len(t))
    _,fit,history,_,reason,_=search(model(t,y,.4),np.linspace(1/1000,1,2001))
    assert reason=='lnBF_below_5'
    assert len(fit['parameters'])==2
    assert np.max(np.abs(np.sort(fit['parameters'][:,0])-[17,43]))<.1
    assert history[-1]['peak_lnBF']<5
    assert all(row['peak_lnBF']>=5 for row in history[:-1])

def test_zero_signal_is_a_valid_result():
    rng=np.random.default_rng(8);t=np.sort(rng.uniform(0,400,120));y=rng.normal(0,1,len(t))
    _,fit,history,_,reason,_=search(model(t,y),np.linspace(1/800,1,1601))
    assert reason=='lnBF_below_5' and len(fit['parameters'])==0

def test_insufficient_dof_not_misreported_as_nondetection():
    t=np.linspace(0,100,12);y=20*np.cos(2*np.pi*t/17)
    _,fit,history,_,reason,_=search(model(t,y,.1),np.linspace(1/200,1,2000))
    assert reason=='insufficient_degrees_of_freedom'
    assert history[-1]['peak_lnBF']>=5 and len(fit['parameters'])==0

def test_noise_is_updated_after_each_accepted_signal():
    rng=np.random.default_rng(43);t=np.sort(rng.uniform(0,300,140));y=5*np.cos(2*np.pi*t/23)+rng.normal(0,.3,len(t));initial=model(t,y,.3);calls=[]
    def update(signal,k):calls.append((k,signal.copy()));return model(t,y,.3)
    _,fit,_,_,reason,_=search(initial,np.linspace(1/600,1,2000),update)
    assert reason=='lnBF_below_5' and len(fit['parameters'])==1
    assert len(calls)==1 and np.std(calls[0][1])>1

def test_commented_activity_header_and_invalid_errors(tmp_path):
    path=tmp_path/'target_PFS.vels';path.write_text('# BJD RV error Sindex Halpha\n2450001 3 1 .2 .3\n2450002 4 -1 .4 .5\n')
    frame,meta=parse_source(path)
    assert len(frame)==1 and frame.sindex.iloc[0]==.2 and frame.halpha.iloc[0]==.3
    assert meta['raw_rows']==2

def test_resume_search_preserves_prior_component():
    rng=np.random.default_rng(52);t=np.sort(rng.uniform(0,400,180))
    y=8*np.cos(2*np.pi*t/17)+3*np.cos(2*np.pi*t/43)+rng.normal(0,.3,len(t));m=model(t,y,.3)
    prior=fit_planets(m,[17],starts=2)
    _,fit,history,_,reason,_=search(m,np.linspace(1/800,1,2000),initial_fit=prior,warm_starts=True)
    assert history[0]['n_existing_signals']==1
    assert reason=='lnBF_below_5' and len(fit['parameters'])==2
    assert np.max(abs(np.sort(fit['parameters'][:,0])-[17,43]))<.1

def test_analytic_jacobian_with_arma_noise():
    rng=np.random.default_rng(71);t=np.sort(rng.uniform(0,300,180));ids=np.where(np.arange(len(t))%2,'A','B')
    noise=pd.DataFrame([dict(set_id=s,jitter=.2,Nma=2,Nar=1,tau=4,tauAR=3,m1=.3,m2=-.1,l1=.2) for s in ['A','B']]).set_index('set_id')
    pars=np.array([[17.,.2,.4],[43.,.55,1.1]])
    B=np.column_stack([kepler_basis(t,*p) for p in pars]);m=ConditionalRV(t,B@np.array([5.,1.,3.,2.])+rng.normal(0,.3,len(t)),np.full(len(t),.3),ids,noise)
    v=pars.copy();v[:,0]=np.log(v[:,0]);v=v.ravel();fun,jac=kepler_objective(m);h=1e-6
    numeric=np.column_stack([(fun(v+np.eye(len(v))[k]*h)-fun(v-np.eye(len(v))[k]*h))/(2*h) for k in range(len(v))])
    assert np.linalg.norm(jac(v)-numeric)/np.linalg.norm(numeric)<1e-6
    assert np.max(abs(fun(v)-m.fit_basis(B)['residual']/m.error))<1e-10

def test_analytic_blind_recovery_and_resume():
    import functools,pfs_search
    original=pfs_search.fit_planets
    try:
        pfs_search.fit_planets=functools.partial(fit_planets,analytic_jac=True)
        test_blind_two_signals_and_threshold_stop()
        test_resume_search_preserves_prior_component()
        test_zero_signal_is_a_valid_result()
        test_noise_is_updated_after_each_accepted_signal()
        test_insufficient_dof_not_misreported_as_nondetection()
    finally:pfs_search.fit_planets=original

def test_missing_checkpoint_bounds_are_not_reported_as_clear():
    import json,tempfile
    from pfs_finalize import checkpoint_bound_metadata
    with tempfile.TemporaryDirectory() as directory:
        out=Path(directory);dest=out/'results/P0001';source=out/'snapshot';dest.mkdir(parents=True);source.mkdir()
        record={'parameters':[[17,.2,.4]],'amplitudes':[3,1]}
        for path in [dest,source]:(path/'model_parameters.json').write_text(json.dumps(record))
        (dest/'candidates.json').write_text(json.dumps([{'candidate_id':'P0001-S01','flags':'no threshold flags','rv_assessment':'conditionally_supported'}]))
        checkpoint_bound_metadata(out,dict(code='P0001',n_signals=1,resumed_from=str(source)))
        candidate=json.loads((dest/'candidates.json').read_text())[0]
        assert 'period-bound check unavailable' in candidate['flags']
        assert candidate['rv_assessment']=='inconclusive'

if __name__=='__main__':
    import tempfile
    for name,function in list(globals().items()):
        if not name.startswith('test_'):continue
        with tempfile.TemporaryDirectory() as directory:
            function(Path(directory)) if name=='test_commented_activity_header_and_invalid_errors' else function()
        print(name,'PASS',flush=True)
