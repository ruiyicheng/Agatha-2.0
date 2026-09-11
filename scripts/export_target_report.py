"""Export a methods page and one target's five pages as a self-contained PDF.

Page resources are cleaned so unrelated targets' shared image objects are not
copied into the task input. Requires PyMuPDF; no RV fitting is repeated.
"""
from pathlib import Path
import argparse
import pymupdf as fitz

SECTIONS = ('Data and noise', 'RV periodograms', 'Full planet model',
            'Activity, sampling and instruments', 'Stability and review flags')

def export_target(source, target, destination):
    source, destination = Path(source), Path(destination)
    with fitz.open(source) as original:
        texts = [p.get_text() for p in original]
        method = [i for i,t in enumerate(texts) if 'Methods and interpretation' in t]
        indices = []
        for section in SECTIONS:
            matches = [i for i,t in enumerate(texts) if f'{target}  |  {section}' in t or f'{target} | {section}' in t]
            if len(matches) != 1:
                # Font text extraction may expand/collapse runs of spaces.
                matches = [i for i,t in enumerate(texts) if ' '.join(f'{target} | {section}'.split()) in ' '.join(t.split())]
            if len(matches) != 1: raise ValueError(f'Expected one {target}/{section} page, got {matches}')
            indices.append(matches[0])
        if len(method) != 1: raise ValueError('Expected one methods page')
        with fitz.open() as report:
            for n,i in enumerate(method+indices,1):
                report.insert_pdf(original, from_page=i, to_page=i, links=False, annots=False)
                p=report[-1]; w,h=p.rect.width,p.rect.height
                # Replace aggregate-page footers with physical, local page numbers.
                p.add_redact_annot(fitz.Rect(.88*w,.96*h,w,h),fill=(1,1,1))
                p.apply_redactions(images=0,graphics=0)
                p.insert_text((.92*w,.975*h),str(n),fontsize=8,color=(.4,.4,.4))
                p.clean_contents(sanitize=True)
            report.set_metadata({'title':f'RV fitting report {target}', 'author':'', 'subject':'Anonymized RV diagnostics'})
            destination.parent.mkdir(parents=True,exist_ok=True)
            report.save(destination,garbage=4,deflate=True,clean=True)
    with fitz.open(destination) as report:
        return [p.get_text() for p in report]

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('source',type=Path);p.add_argument('target');p.add_argument('output',type=Path);a=p.parse_args()
    pages=export_target(a.source,a.target,a.output);print(f'{a.output}: {len(pages)} pages')
if __name__=='__main__':main()
