"""Label-independent eligibility: every catalogue period needs a distinct RV peak.

Uses cached, catalogue-informed sequential RV scans before the candidate was
subtracted. Never treats a period-constrained orbital fit as signal recovery.
"""
import math
import numpy as np
from scipy.signal import find_peaks

def match_recovery(frequency,sequence,removal_order,candidates,baseline_days,
                   period_tolerance=.05,min_logbf=5.,top_peaks=10):
    if not 0<period_tolerance<1 or top_peaks<1 or baseline_days<=0:
        raise ValueError('Invalid recovery settings')
    f=np.asarray(frequency);seq=np.asarray(sequence);order=np.asarray(removal_order,dtype=int)
    n=len(candidates)
    if f.ndim!=1 or seq.shape!=(n+1,len(f)) or sorted(order.tolist())!=list(range(n)):
        raise ValueError('Scan shape/removal order does not match catalogue component count')
    if not np.all(np.isfinite(f)) or np.any(f<=0):raise ValueError('Invalid frequency grid')
    peaks=[]
    for stage,y in enumerate(seq[:-1]):
        local=find_peaks(y)[0];local=local[np.argsort(-y[local],kind='stable')][:top_peaks]
        for rank,index in enumerate(local,1):
            if np.isfinite(y[index]) and (min_logbf is None or y[index]>=min_logbf):
                peaks.append(dict(stage=stage,rank=rank,index=int(index),frequency=float(f[index]),period=float(1/f[index]),logbf=float(y[index])))
    # Peaks closer than the Rayleigh resolution are one unresolved signal,
    # even if they reappear in several subtraction stages.
    peaks.sort(key=lambda p:p['frequency']);cluster=-1;previous=None
    for peak in peaks:
        if previous is None or peak['frequency']-previous>1/baseline_days:cluster+=1
        peak['cluster']=cluster;previous=peak['frequency']
    options=[];audit=[]
    for j,c in enumerate(candidates):
        try:period=float(c['period_days'])
        except (ValueError,TypeError):period=math.nan
        before=int(np.where(order==j)[0][0]);eligible={};nearest=None
        valid=math.isfinite(period) and period>0
        if valid:
            for peak in peaks:
                if peak['stage']>before:continue
                error=abs(peak['period']-period)/period;p={**peak,'relative_error':error}
                if nearest is None or error<nearest['relative_error']:nearest=p
                if error<=period_tolerance:
                    old=eligible.get(peak['cluster'])
                    if old is None or (error,-peak['logbf'])<(old['relative_error'],-old['logbf']):eligible[peak['cluster']]=p
        options.append(dict(sorted(eligible.items(),key=lambda v:(v[1]['relative_error'],-v[1]['logbf']))))
        audit.append(dict(candidate_id=c['candidate_id'],catalogue_period_days=period if valid else None,
                          matched=False,reason='missing_catalogue_period' if not valid else 'no_matching_top_peak',nearest=nearest))
    # Maximum bipartite matching avoids both reusing one peak for two planets
    # and a greedy match needlessly blocking another candidate.
    owner={}
    def assign(j,visited):
        for group in options[j]:
            if group in visited:continue
            visited.add(group)
            if group not in owner or assign(owner[group],visited):owner[group]=j;return True
        return False
    for j in sorted(range(n),key=lambda j:len(options[j])):assign(j,set())
    assigned={j:group for group,j in owner.items()}
    for j,row in enumerate(audit):
        peak=options[j][assigned[j]] if j in assigned else row['nearest']
        row.pop('nearest')
        if j in assigned:row.update(matched=True,reason='matched')
        elif options[j]:row['reason']='distinct_peak_conflict'
        row.update(peak_period_days=None if peak is None else peak['period'],relative_period_error=None if peak is None else peak['relative_error'],
                   peak_logbf=None if peak is None else peak['logbf'],stage=None if peak is None else peak['stage'],rank=None if peak is None else peak['rank'])
    return all(x['matched'] for x in audit),audit
