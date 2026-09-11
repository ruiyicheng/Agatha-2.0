"""Read or render the task's PDF using 1-based physical page numbers."""
import argparse
from pathlib import Path
import pymupdf
p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['text','render']);p.add_argument('--page',type=int);p.add_argument('--output',default='/workspace/pages');a=p.parse_args()
with pymupdf.open('/workspace/report.pdf') as doc:
    indices=[a.page-1] if a.page else range(len(doc))
    for i in indices:
        if not 0<=i<len(doc):raise SystemExit('Invalid page number')
        if a.mode=='text':print(f'--- PAGE {i+1} ---\n{doc[i].get_text()}')
        else:
            folder=Path(a.output);folder.mkdir(parents=True,exist_ok=True)
            file=folder/f'page-{i+1}.png';doc[i].get_pixmap(matrix=pymupdf.Matrix(1.6,1.6)).save(file);print(file)
