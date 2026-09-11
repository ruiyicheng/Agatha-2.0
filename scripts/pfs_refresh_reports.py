"""Refresh completed PDFs from cached results, optionally watching an active batch."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse,concurrent.futures,hashlib,json,time,traceback
from pathlib import Path
from pfs_report import report

def render(args):
    out,code,version=args;dest=out/'results'/code;start=time.perf_counter()
    try:
        report(dest,out/'pdf'/(code+'.pdf'))
        record=dict(code=code,status='complete',renderer_sha256=version,pdf_seconds=time.perf_counter()-start)
    except Exception:
        (dest/'report_failure.PRIVATE.txt').write_text(traceback.format_exc())
        record=dict(code=code,status='failed',renderer_sha256=version,pdf_seconds=time.perf_counter()-start)
    (dest/'report_status.json').write_text(json.dumps(record,indent=2)+'\n');return record

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path);p.add_argument('--workers',type=int,default=4);p.add_argument('--watch',action='store_true');a=p.parse_args();out=a.output.resolve()
    version=hashlib.sha256(Path(__file__).with_name('pfs_report.py').read_bytes()).hexdigest();attempted=set();records=[]
    with concurrent.futures.ProcessPoolExecutor(max_workers=a.workers) as pool:
        futures={}
        while True:
            for dest in (out/'results').glob('P*'):
                if dest.name in attempted:continue
                progress=dest/'progress.json';pdf=out/'pdf'/(dest.name+'.pdf')
                if not progress.exists() or not pdf.exists():continue
                if json.loads(progress.read_text()).get('stage') not in {'complete','failed','search_limited','insufficient_data'}:continue
                status=dest/'report_status.json'
                if status.exists():
                    previous=json.loads(status.read_text())
                    if previous['status']=='complete' and previous['renderer_sha256']==version:
                        attempted.add(dest.name);records.append(previous);continue
                attempted.add(dest.name);futures[pool.submit(render,(out,dest.name,version))]=dest.name
            for future in list(futures):
                if not future.done():continue
                result=future.result();records.append(result);del futures[future];print(json.dumps(result),flush=True)
            summary=dict(refreshed=len(records),complete=sum(r['status']=='complete' for r in records),failed=sum(r['status']=='failed' for r in records),pending=len(futures),renderer_sha256=version)
            (out/'report_refresh_status.json').write_text(json.dumps(summary,indent=2)+'\n')
            batch=json.loads((out/'batch_status.json').read_text())
            if not futures and (not a.watch or (batch['reported']==batch['requested'] and len(attempted)>=batch['requested'])):break
            time.sleep(3)
    if any(r['status']!='complete' for r in records):raise SystemExit('Some PDF refreshes failed')

if __name__=='__main__':main()
