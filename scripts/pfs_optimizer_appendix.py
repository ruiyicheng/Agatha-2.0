"""Append a coded comparison of complete alternative optimization runs."""
import argparse
import json
from pathlib import Path

import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from pypdf import PdfReader, PdfWriter
from pfs_report import page, paragraph


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('--codes',nargs='+',required=True)
    args=parser.parse_args();out=args.output.resolve()
    for code in args.codes:
        primary=out/'results'/code;alternate=out/'analytic_retry/results'/code
        marker=primary/'optimizer_comparison.json'
        target=out/'pdf'/(code+'.pdf')
        if marker.exists() and 'Optimization sensitivity' in (PdfReader(target).pages[-1].extract_text() or ''):
            print(code,'comparison already present');continue
        a=json.loads((primary/'summary.json').read_text());b=json.loads((alternate/'summary.json').read_text())
        if a['status']!='complete' or b['status']!='complete':continue
        if a.get('optimizer_jacobian') and not a.get('resumed_from'):continue
        ap=np.array(json.loads((primary/'model_parameters.json').read_text())['parameters'])
        bp=np.array(json.loads((alternate/'model_parameters.json').read_text())['parameters'])
        if ap.shape==bp.shape and np.allclose(ap,bp,rtol=1e-9,atol=1e-9):continue
        record=dict(code=code,primary_components=len(ap),alternative_components=len(bp),
                    primary_last_lnBF=a['last_peak_lnBF'],alternative_last_lnBF=b['last_peak_lnBF'],
                    primary_periods_days=ap[:,0].tolist(),alternative_periods_days=bp[:,0].tolist(),
                    interpretation='Both searches meet the conditional threshold, but the fitted decomposition depends on optimization and noise reselection.')
        appendix=primary/'optimizer_appendix.pdf'
        with PdfPages(appendix,metadata={'Title':code+' | Optimization sensitivity','Author':''}) as pdf:
            fig=page(code+' | Optimization sensitivity');y=.90
            y=paragraph(fig,y,'An independent optimization run on the same RV measurements produced a different complete model. It used analytic projection derivatives and warm starts; the primary search used finite differences and its original starting strategy. Both use the same noise-model candidates, period grid and conditional ln BF threshold. These are alternative fits to the same data, not independent observations.',size=10,width=90)
            y=paragraph(fig,y,f"Primary model: {len(ap)} components; strongest remaining peak ln BF = {a['last_peak_lnBF']:.3f}. Alternative model: {len(bp)} components; strongest remaining peak ln BF = {b['last_peak_lnBF']:.3f}.")
            rows=[]
            for i in range(max(len(ap),len(bp))):
                rows.append([i+1,f'{ap[i,0]:.7g}' if i<len(ap) else '',f'{bp[i,0]:.7g}' if i<len(bp) else ''])
            ax=fig.add_axes([.10,.25,.80,max(.20,y-.28)]);ax.axis('off')
            table=ax.table(cellText=rows,colLabels=['Addition order','Primary period (d)','Alternative period (d)'],cellLoc='center',loc='upper center');table.auto_set_font_size(False);table.set_fontsize(9);table.scale(1,1.6)
            paragraph(fig,.22,'Rows follow addition order and do not imply a physical one-to-one match. Changes can arise from local optima, evolving nuisance/noise estimates, or alternative descriptions of activity, aliases and long-term motion. Reaching the residual threshold does not establish a unique orbital solution. Treat component-level interpretations affected by these differences as uncertain.',size=10,width=90)
            fig.text(.93,.025,str(len(PdfReader(target).pages)+1),ha='right',fontsize=8)
            pdf.savefig(fig);plt.close(fig)
        writer=PdfWriter();reader=PdfReader(target)
        for p in reader.pages:writer.add_page(p)
        for p in PdfReader(appendix).pages:writer.add_page(p)
        writer.add_metadata({'/Title':code+': RV signal vetting','/Author':'','/Subject':'RV evidence and optimization sensitivity'})
        temp=target.with_suffix('.tmp.pdf')
        with temp.open('wb') as stream:writer.write(stream)
        temp.replace(target);marker.write_text(json.dumps(record,indent=2)+'\n');print(code,'comparison appended')


if __name__=='__main__':main()
