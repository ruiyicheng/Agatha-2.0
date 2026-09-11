"""Prepare auditable RV-only inputs and a private anonymization key."""
import argparse,csv,json,hashlib,math,re
from pathlib import Path
from collections import defaultdict
import numpy as np

def read(p):
    with p.open(newline='') as f:return list(csv.DictReader(f))
def write(p,rows):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def num(v):
    try:return float(v)
    except:return math.nan
def cn(v):return re.sub('[^a-z0-9]','',v.lower())

def main():
    a=argparse.ArgumentParser();a.add_argument('dataset',type=Path);a.add_argument('output',type=Path);args=a.parse_args()
    ds=args.dataset.resolve();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    hosts=read(ds/'host_manifest.csv');planets=read(ds/'val.csv')+read(ds/'test.csv')
    hosts=sorted(hosts,key=lambda h:(h['split'],hashlib.sha256(('agatha-blind-v1:'+h['host_id']).encode()).hexdigest()))
    key=[];aliases=[];tasks=[];cand=[];audit=[];split_counts=defaultdict(int)
    for h in hosts:
        split=h['split'];split_counts[split]+=1;code=('V' if split=='val' else 'T')+f'{split_counts[split]:02d}'
        rows=[p for p in planets if p['target']==h['target']]
        rows.sort(key=lambda p:hashlib.sha256(p['planet_name'].encode()).hexdigest())
        key.append(dict(code=code,split=split,target=h['target'],host_id=h['host_id']))
        aliases.append(dict(code=code,aliases=h.get('reserved_host_aliases','')))
        for i,p in enumerate(rows,1):
            cand.append(dict(code=code,candidate_id=f'{code}-C{i:02}',period_days=p['orbital_period_days'],
                catalogue_label=p['planet_status'],planet_name=p['planet_name']))
        signatures=set();seen=defaultdict(set);instrows=[]
        files=sorted(h['rv_files'].split(';'),key=lambda x:('pre' not in x.lower() and 'post' not in x.lower(),x))
        for path in files:
            src=ds/path;ls=[s.split() for s in src.read_text().splitlines() if s.strip() and not s.lstrip().startswith('#')]
            header=not math.isfinite(num(ls[0][0]));names=ls[0] if header else []
            rs=ls[int(header):];width=max(map(len,rs));d=np.array([[num(v) for v in r]+[math.nan]*(width-len(r)) for r in rs])
            d=d[np.all(np.isfinite(d[:,:3]),axis=1)&(d[:,2]>0)]
            ins=src.stem.split('_')[-1].upper();offset=2400000 if np.nanmedian(d[:,0])<2400000 else 0
            d[:,0]+=offset
            sidx=hidx=None;schema='named header'
            if header:
                normalized=[cn(s) for s in names]
                sidx=next((i for i,s in enumerate(normalized) if s in {'sindex','smw','shk','s'}),None)
                hidx=next((i for i,s in enumerate(normalized) if s in {'halpha','ha','iha'}),None)
            elif ins=='PFS' and width==8:
                if np.allclose(d[:,3],1):sidx=4;schema='PFS constant-column format: S=5; no identified H-alpha'
                else:sidx,hidx=3,4;schema='PFS activity format: S=4, H-alpha=5'
            elif ins=='KECK' and width==7:sidx,hidx=3,4;schema='KECK seven-column format: S=4, H-alpha=5'
            else:schema='time/RV/error only; unidentified extra columns unused'
            signature=hashlib.sha256(np.round(d[:,:3],8).tobytes()).hexdigest()
            reason='';raw=len(d)
            if signature in signatures:reason='duplicate full RV series'
            signatures.add(signature)
            if reason:
                audit.append(dict(code=code,source=str(src),instrument=ins,raw_rows=raw,used_rows=0,
                    jd_offset=offset,activity_schema=schema,notes=reason));continue
            family={'HIRES':'KECK','HARPSPRE':'HARPS','HARPSPOST':'HARPS','HARPSN':'HARPN','HARPS-N':'HARPN'}.get(ins,ins)
            keep=[]
            for j,t in enumerate(d[:,0]):
                ident=round(t,7)
                if ident not in seen[family]:keep.append(j);seen[family].add(ident)
            d=d[keep];d=d[np.argsort(d[:,0])]
            if not len(d):continue
            si=np.full(len(d),np.nan) if sidx is None else d[:,sidx].copy()
            ha=np.full(len(d),np.nan) if hidx is None else d[:,hidx].copy()
            for v in [si,ha]:v[(v<=0)|(v>=99)|~np.isfinite(v)]=np.nan
            # Instrument-labelled sets are kept distinct; no nightly binning.
            setid=f'I{len(instrows)+1:02d}_{ins}'
            dest=out/'prepared'/code/f'{setid}.csv'
            write(dest,[dict(time=t,rv=y,error=e,sindex=s,halpha=a) for (t,y,e),s,a in zip(d[:,:3],si,ha)])
            instrows.append(dict(set_id=setid,instrument=ins,file=str(dest),n=len(d)))
            audit.append(dict(code=code,source=str(src),instrument=ins,raw_rows=raw,used_rows=len(d),
                jd_offset=offset,activity_schema=schema,notes=f'{raw-len(d)} duplicate instrument epochs removed'))
        write(out/'prepared'/code/'instruments.csv',instrows)
        known=[float(p['orbital_period_days']) for p in rows if p['orbital_period_days']]
        tasks.append(dict(code=code,split=split,n_signals=len(rows),raw_rv_points=int(h['total_rv_points']),
            used_rv_points=sum(i['n'] for i in instrows),n_instruments=len(instrows),
            max_catalogue_period=max(known,default=0),input_dir=str(out/'prepared'/code)))
    write(out/'private_aliases.csv',aliases);write(out/'private_identity_key.csv',key);write(out/'private_candidates.csv',cand)
    write(out/'preprocessing_audit.csv',audit);write(out/'tasks.csv',tasks)
    print(json.dumps(dict(targets=len(tasks),signals=sum(t['n_signals'] for t in tasks),
        raw_rows=sum(t['raw_rv_points'] for t in tasks),used_rows=sum(t['used_rv_points'] for t in tasks)),indent=2))

if __name__=='__main__':main()
