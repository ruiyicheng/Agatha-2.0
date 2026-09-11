"""Conditional Agatha likelihood, streaming scans and variable-projection Keplerians.

MA transforms both the data and design; AR subtracts lagged observed RV, exactly
as CircularSig in periodoframe.R. Hyperparameters are fixed within each scan.
No posterior samples, dynamical integration, photometry or label-based choices.
"""
import numpy as np
from scipy.linalg import qr
from scipy.optimize import least_squares, minimize_scalar

class ConditionalRV:
    def __init__(self,t,y,error,ids,noise):
        self.t=np.asarray(t,float);self.y=np.asarray(y,float); self.ids=np.asarray(ids)
        self.ma=[];self.ar=[];self.error=np.array(error,float).copy()
        self.levels=np.unique(self.ids)
        for sid in self.levels:
            r=np.where(self.ids==sid)[0];n=noise.loc[sid]
            self.error[r]=np.hypot(self.error[r],n.jitter)
            for prefix,count,tau,out in [('m',int(n.Nma),n.tau,self.ma),('l',int(n.Nar),n.tauAR,self.ar)]:
                for lag in range(1,count+1):
                    dest,src=r[lag:],r[:-lag]
                    out.append((dest,src,float(n[prefix+str(lag)])*np.exp(-np.abs(self.t[dest]-self.t[src])/tau)))
        self.X=np.column_stack([self.ids==i for i in self.levels]+[self.t/max(np.ptp(self.t),1)])*1.
        self.A=self.transform(self.X)/self.error[:,None]
        self.Q,self.R=qr(self.A,mode='economic')
        self.Y=self.response(self.y)/self.error
        self.z=self.Y-self.Q@(self.Q.T@self.Y)
        self.chi0=self.z@self.z
    def transform(self,x):
        out=np.array(x,float).copy()
        for dest,src,c in self.ma:
            out[dest]-=c[:,None]*x[src] if out.ndim==2 else c*x[src]
        return out
    def response(self,y):
        out=self.transform(y)
        # Agatha's AR uses the observed RV, not MA innovations.
        for dest,src,c in self.ar:out[dest]-=c*self.y[src]
        return out
    def scan(self,freq,subtract=None,block=128):
        z=self.z if subtract is None else self.Y-self.transform(subtract)/self.error
        if subtract is not None:z=z-self.Q@(self.Q.T@z)
        power=np.empty(len(freq)); chi=z@z
        for a in range(0,len(freq),block):
            f=freq[a:a+block];ang=2*np.pi*self.t[:,None]*f
            C=self.transform(np.cos(ang))/self.error[:,None];S=self.transform(np.sin(ang))/self.error[:,None]
            C-=self.Q@(self.Q.T@C);S-=self.Q@(self.Q.T@S)
            cc=np.sum(C*C,axis=0);ss=np.sum(S*S,axis=0);cs=np.sum(C*S,axis=0)
            yc=z@C;ys=z@S;det=cc*ss-cs*cs
            gain=np.divide(ss*yc*yc+cc*ys*ys-2*cs*yc*ys,det,out=np.zeros_like(det),where=det>1e-20)
            power[a:a+len(f)]=np.maximum(0,gain)/2-np.log(len(self.y))
        return power
    def fit_basis(self,B):
        if B.shape[1]:
            C=self.transform(B)/self.error[:,None];Cr=C-self.Q@(self.Q.T@C)
            amp=np.linalg.lstsq(Cr,self.z,rcond=1e-10)[0]
            gamma=np.linalg.lstsq(self.A,self.Y-C@amp,rcond=1e-10)[0]
            signal=B@amp
        else:
            amp=np.array([]);gamma=np.linalg.lstsq(self.A,self.Y,rcond=1e-10)[0];signal=np.zeros(len(self.y))
        deterministic=signal+self.X@gamma
        resid=(self.response(self.y)-self.transform(deterministic))
        return dict(amplitude=amp,gamma=gamma,signal=signal,deterministic=deterministic,residual=resid,chi2=np.sum((resid/self.error)**2))

def kepler_basis(t,period,e,phase):
    M=(2*np.pi*t/period+phase)%(2*np.pi); E=M.copy()
    for _ in range(30):
        step=(E-e*np.sin(E)-M)/(1-e*np.cos(E)); E-=step
        if np.max(abs(step))<1e-11:break
    nu=2*np.arctan2(np.sqrt(1+e)*np.sin(E/2),np.sqrt(1-e)*np.cos(E/2))
    return np.column_stack((np.cos(nu)+e,np.sin(nu)))

def kepler_objective(model):
    """Residual and exact variable-projection Jacobian in log(P), e, phase.

    Both terms in the derivative of the least-squares projection are retained.
    Instrument offsets/trend and fixed MA noise are projected identically to
    fit_basis. The SVD uses the same relative rank cutoff as fit_basis.
    """
    from scipy.linalg import svd
    cache={}
    def calculate(v):
        if 'v' in cache and np.array_equal(v,cache['v']):return cache
        pars=v.reshape(-1,3);blocks=[];derivatives=[]
        for logp,e,phase in pars:
            period=np.exp(logp);B=kepler_basis(model.t,period,e,phase)
            c=B[:,0]-e;s=B[:,1];dnu_dM=(1+e*c)**2/(1-e*e)**1.5
            dnu_de=s*(2+e*c)/(1-e*e)
            angular=np.column_stack((-s,c))
            derivatives.extend([angular*(-2*np.pi*model.t/period*dnu_dM)[:,None],
                                angular*dnu_de[:,None]+np.array([1.,0.]),
                                angular*dnu_dM[:,None]])
            blocks.append(B)
        B=np.column_stack(blocks);C=model.transform(B)/model.error[:,None]
        C-=model.Q@(model.Q.T@C)
        U,sigma,Vt=svd(C,full_matrices=False,check_finite=False)
        keep=sigma>sigma[0]*1e-10;U=U[:,keep];sigma=sigma[keep];Vt=Vt[keep]
        amp=Vt.T@((U.T@model.z)/sigma);r=model.z-C@amp
        J=np.empty((len(model.t),len(v)))
        for k,D in enumerate(derivatives):
            D=model.transform(D)/model.error[:,None];D-=model.Q@(model.Q.T@D)
            j=k//3;Da=D@amp[2*j:2*j+2]
            J[:,k]=-Da+U@(U.T@Da)-(U/sigma)@(Vt[:,2*j:2*j+2]@(D.T@r))
        cache.update(v=v.copy(),residual=r,jacobian=J)
        return cache
    return lambda v:calculate(v)['residual'],lambda v:calculate(v)['jacobian']


def fit_planets(model,periods,starts=3,max_nfev=250,bounds=None,initial_parameters=None,polish_nfev=1500,warm_starts=False,analytic_jac=False):
    periods=np.asarray(periods);n=len(periods)
    if n==0:return {**model.fit_basis(np.empty((len(model.t),0))), 'parameters':np.empty((0,3)),'components':np.empty((len(model.t),0)),'converged':True,'nfev':0}
    lo=periods*.9 if bounds is None else bounds[0];hi=periods*1.1 if bounds is None else bounds[1]
    lower=np.column_stack((np.log(lo),np.zeros(n),np.full(n,-4*np.pi))).ravel()
    upper=np.column_stack((np.log(hi),np.full(n,.85),np.full(n,4*np.pi))).ravel()
    def basis(v):return np.column_stack([kepler_basis(model.t,np.exp(p),e,phase) for p,e,phase in v.reshape(-1,3)])
    def residual(v):
        fit=model.fit_basis(basis(v));return fit['residual']/model.error
    jac='2-point'
    if analytic_jac:residual,jac=kepler_objective(model)
    best=None
    for j in range(starts):
        v=np.column_stack((np.log(periods),np.full(n,[.01,.25,.55][j%3]),np.full(n,j*np.pi/2))).ravel()
        if initial_parameters is not None and (j==0 or warm_starts):
            initial=np.asarray(initial_parameters,float).copy()
            if initial.shape!=(n,3):raise ValueError('Initial parameters must be (n_signals, 3)')
            if j>0:
                initial[-1,1]=[.01,.25,.55][j%3];initial[-1,2]+=j*np.pi/2
            initial[:,0]=np.log(initial[:,0]);v=np.clip(initial.ravel(),lower+1e-9,upper-1e-9)
        opt=least_squares(residual,v,jac=jac,bounds=(lower,upper),max_nfev=max_nfev,ftol=1e-7,xtol=1e-7,gtol=1e-6,x_scale='jac')
        if best is None or np.sum(opt.fun**2)<np.sum(best.fun**2):best=opt
    # Retry the best solution before accepting an evaluation-limit exit.
    if not best.success:
        retry=least_squares(residual,best.x,jac=jac,bounds=(lower,upper),max_nfev=polish_nfev,ftol=1e-7,xtol=1e-7,gtol=1e-6,x_scale='jac')
        if np.sum(retry.fun**2)<=np.sum(best.fun**2):
            retry.nfev+=best.nfev;best=retry
    B=basis(best.x);fit=model.fit_basis(B);par=best.x.reshape(-1,3);par[:,0]=np.exp(par[:,0])
    components=np.column_stack([B[:,2*i:2*i+2]@fit['amplitude'][2*i:2*i+2] for i in range(n)])
    return {**fit,'parameters':par,'components':components,'converged':bool(best.success),'nfev':int(best.nfev),'bounds':(lo,hi)}
