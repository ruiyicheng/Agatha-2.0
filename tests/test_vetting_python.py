import sys,tempfile,subprocess,os
from pathlib import Path
import numpy as np,pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from vetting_core import ConditionalRV,kepler_basis,fit_planets
root=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory() as td:
 subprocess.run([sys.argv[1],'tests/vetting_likelihood_fixture.R',td],cwd=root,check=True)
 data=pd.read_csv(Path(td)/'data.csv');noise=pd.read_csv(Path(td)/'noise.csv')
 for row in noise.itertuples():
  model=ConditionalRV(data.t,data.y,data.error,np.repeat('A',len(data)),noise[noise.model==row.model].set_index('set_id'))
  fit=model.fit_basis(np.column_stack((np.cos(2*np.pi*data.t/17),np.sin(2*np.pi*data.t/17))))
  ref=pd.read_csv(Path(td)/(row.model+'.csv')).residual
  assert np.max(np.abs(fit['residual']-ref))<2e-6,(row.model,np.max(np.abs(fit['residual']-ref)))
  ll=-.5*fit['chi2']-np.sum(np.log(np.sqrt(2*np.pi)*model.error))
  assert abs(ll-row.logL)<1e-5,(row.model,ll,row.logL)
  scan=model.scan(np.array([1/17]))[0]
  assert abs(scan-((model.chi0-fit['chi2'])/2-np.log(len(data))))<1e-6
 print('Python likelihood and scan match Agatha R for WN, AR1, MA1 and MA2')
# Non-circular two-planet recovery checks the joint Keplerian model, not just sinusoids.
rng=np.random.default_rng(812);t=np.sort(rng.uniform(0,500,180));true=[(17,.25,.7),(43,.4,1.6)]
y=kepler_basis(t,*true[0])@np.array([5.,2.])+kepler_basis(t,*true[1])@np.array([3.,-1.])+rng.normal(0,.15,len(t))
no=pd.DataFrame([dict(set_id='A',jitter=0,Nma=0,Nar=0,tau=1,tauAR=1)]).set_index('set_id')
m=ConditionalRV(t,y,np.full(len(t),.15),np.repeat('A',len(t)),no);f=fit_planets(m,[17,43],starts=3)
assert np.max(abs(f['parameters'][:,0]-[17,43]))<.05
assert np.max(abs(f['parameters'][:,1]-[.25,.4]))<.05
assert f['components'].shape==(len(t),2)
print('Joint eccentric two-planet injection recovered')
