"""Coded target reports for blind additive RV searches, including non-detections."""
import json, textwrap
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pypdf import PdfReader
from vetting_core import kepler_basis

plt.rcParams.update({'font.size':9,'axes.titlesize':10,'axes.labelsize':9,'pdf.fonttype':42,'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False})

def page(title):
    fig=plt.figure(figsize=(8.27,11.69));fig.suptitle(title,x=.08,y=.96,ha='left',fontsize=16,weight='bold');return fig

def paragraph(fig,y,text,size=9,width=100):
    lines=textwrap.wrap(text,width=width);fig.text(.08,y,'\n'.join(lines),fontsize=size,va='top');return y-len(lines)*size/840-.025

def read(path,default=None):return json.loads(path.read_text()) if path.exists() else default

def draw_scan(ax,frequency,power,label=None):
    order=np.argsort(1/frequency);ax.plot(1/frequency[order],power[order],lw=.65,label=label,rasterized=True);ax.set_xscale('log');ax.axhline(5,color='gray',ls='--',lw=.7);ax.set_ylabel('Conditional ln BF')

def report(dest,path):
    dest=Path(dest);path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);summary=read(dest/'summary.json');code=summary['code'];count=0
    tmp=path.with_suffix('.tmp.pdf')
    with PdfPages(tmp,metadata={'Title':f'{code}: RV signal vetting','Author':'','Subject':'RV and spectroscopic activity only; no catalogue priors or MCMC'}) as pdf:
        def save(fig):
            nonlocal count
            count+=1;fig.text(.93,.025,str(count),ha='right',fontsize=8);pdf.savefig(fig);plt.close(fig)
        fig=page(f'{code} | RV signal search');y=.90
        y=paragraph(fig,y,f"Status: {summary['status']}. Measurements: {summary['n']} used / {summary['n_raw']} numeric source rows; {summary['n_instruments']} instrument sets; {summary['pfs_n']} retained PFS RVs. All available audited RV instruments are used. No photometry is used.",size=11,width=87)
        y=paragraph(fig,y,f"Accepted periodic components: {summary['n_signals']}. Stopping reason: {summary.get('stop_reason','pending')}. Search complete to the requested threshold: {summary.get('search_complete',False)}.")
        if 'last_peak_lnBF' in summary:y=paragraph(fig,y,f"Strongest final residual peak: P = {summary['last_peak_period_days']:.6g} days, conditional ln BF = {summary['last_peak_lnBF']:.3f}. Threshold stop requires ln BF < 5. A degrees-of-freedom or optimizer limit is not a non-detection.")
        if 'chi2_per_measurement' in summary:y=paragraph(fig,y,f"Full-model chi-squared/N = {summary['chi2_per_measurement']:.3f}. Baseline = {summary['baseline_days']:.1f} days. Search range = 1 to {summary['period_range_days'][1]:.1f} days. Processing time before PDF generation: {summary['elapsed_seconds']/60:.2f} minutes.")
        y=paragraph(fig,y,'Method: choose WN+jitter, AR1, MA1 or MA2 per instrument by BIC. Scan residual RVs on a uniform frequency grid with four samples per baseline resolution element; locally refine the strongest peak. When its conditional ln BF is at least 5, add one Keplerian and jointly refit all components, then reselect noise on signal-subtracted RVs and refit. Repeat until the strongest residual peak is below 5. No catalogue periods or planet counts initialize the search.')
        y=paragraph(fig,y,'Evidence is a conditional BIC approximation: periodogram ln BF = likelihood gain minus ln N for two sinusoid coefficients, with noise fixed within each scan. It is neither a marginalized Bayes factor nor a calibrated false-alarm probability. Five orbital parameters per component are penalized separately in drop-component delta BIC. Local optimizers, aliases and noise absorption can affect results. Eccentricity is bounded at 0.85; period fits start within +/-10% of the detected peak. Signals may be stellar activity or companions of nonplanetary mass.')
        y=paragraph(fig,y,'Review: activity diagnostics inspect secondary local peaks within 5% of fitted periods. Standardized unit-error activity scans and exploratory lag correlations have different calibration from RV evidence. Activity coincidences warrant review, not automatic rejection. No posterior errors, MCMC, photometry, astrometric mass classification or catalogue labels are supplied.')
        y=paragraph(fig,y,'Data handling: invalid RV uncertainties are excluded, matching instrument-family epochs are deduplicated, and pre/post subsets take precedence over combined copies. All retained sets contribute RVs; red-noise selection and proxy scans require more than 10 points. Instruments with at most 3 RVs use supplied errors with zero fixed jitter. Relative times and coded IDs hide explicit target identities; periods can remain recognizable.')
        if summary['status'] in {'failed','insufficient_data'}:
            paragraph(fig,y,'No reliable complete search result is available. This PDF records coverage and the failure/insufficient-data status. Detailed source and failure records are stored separately from this coded report.');save(fig)
        else:
            save(fig)
            noise=pd.read_csv(dest/'selected_noise.csv');data=pd.read_csv(dest/'fit.csv');candidates=read(dest/'candidates.json',[])
            for start in range(0,len(noise),24):
                fig=page(f'{code} | Data and selected noise');chunk=noise.iloc[start:start+24]
                ax=fig.add_axes([.08,.43,.85,.45]);ax.axis('off');rows=[[r.set_id,int(r.n),r.model,f'{r.jitter:.3g}',f'{r.BIC:.2f}'] for r in chunk.itertuples()]
                table=ax.table(cellText=rows,colLabels=['Instrument set','N','Noise','Jitter m/s','BIC'],cellLoc='center',loc='upper center');table.auto_set_font_size(False);table.set_fontsize(8);table.scale(1,1.45)
                ax=fig.add_axes([.11,.12,.81,.24])
                for sid,part in data.groupby('set_id'):ax.scatter(part.t_days/365.25,part.rv,s=5,alpha=.6,label=sid,rasterized=True)
                ax.set(xlabel='Relative time (years)',ylabel='Centered RV (m/s)');save(fig)
            scan=np.load(dest/'search_periodograms.npz');history=pd.read_csv(dest/'search_history.csv')
            for start in range(0,len(history),4):
                fig=page(f'{code} | Sequential residual periodograms');axes=fig.subplots(min(4,len(history)-start),1,squeeze=False).ravel();fig.subplots_adjust(left=.12,right=.94,top=.90,bottom=.08,hspace=.55)
                for ax,(idx,row) in zip(axes,history.iloc[start:start+4].iterrows()):
                    draw_scan(ax,scan['frequency'],scan['power'][idx]);ax.axvline(row.peak_period_days,color='#b91c1c',ls=':',lw=.8);ax.set_title(f"{int(row.n_existing_signals)} components removed: peak {row.peak_period_days:.5g} d, ln BF {row.peak_lnBF:.2f}; {row.action.replace('_',' ')}",fontsize=8);ax.set_xlabel('Period (days)')
                save(fig)
            for start in range(0,len(candidates),4):
                fig=page(f'{code} | Fitted signals and review');chunk=candidates[start:start+4];y=.91
                for c in chunk:y=paragraph(fig,y,f"{c['candidate_id']}: P={c['period_days']:.7g} d; K={c['K_m_s']:.4g} m/s; e={c['eccentricity']:.4f}; drop delta BIC={c['delta_BIC_drop']:.2f}; {c['cycles']:.1f} cycles. Assessment: {c['rv_assessment']}. Flags: {c['flags']}.",size=8,width=112)
                cols=1 if len(chunk)==1 else 2;rows=int(np.ceil(len(chunk)/cols))
                axes=fig.subplots(rows,cols,squeeze=False).ravel();fig.subplots_adjust(left=.11,right=.94,top=min(.78,y-.025),bottom=.09,hspace=.43,wspace=.30)
                parameters=read(dest/'model_parameters.json')
                for ax,c in zip(axes,chunk):
                    j=int(c['candidate_id'].rsplit('S',1)[1]);phase=(data.t_days/c['period_days'])%1;component=data['component_'+str(j)];order=np.argsort(phase)
                    ax.scatter(phase,component+data.residual,s=4,alpha=.5,rasterized=True)
                    if parameters:
                        grid=np.linspace(0,1,500);par=parameters['parameters'][j-1];amplitude=np.array(parameters['amplitudes'][2*(j-1):2*j]);curve=kepler_basis(grid*par[0],*par)@amplitude;ax.plot(grid,curve,color='black',lw=.8)
                    else:ax.plot(phase.iloc[order],component.iloc[order],color='black',lw=.8)
                    ax.set(title=c['candidate_id'],xlabel='Phase',ylabel='Conditional RV m/s')
                for ax in axes[len(chunk):]:ax.set_axis_off()
                save(fig)
            activity=np.load(dest/'activity_periodograms.npz');fig=page(f'{code} | Spectroscopic activity');axes=fig.subplots(2,1);fig.subplots_adjust(left=.12,right=.94,top=.90,bottom=.22,hspace=.35)
            for ax,proxy,title in zip(axes,['sindex','halpha'],['S-index','H-alpha']):
                keys=[k for k in activity.files if k.endswith('__'+proxy)]
                for key in keys:draw_scan(ax,activity['frequency'],activity[key],key.split('__')[0])
                for c in candidates:ax.axvline(c['period_days'],color='black',lw=.5,ls=':')
                ax.set_title(title,loc='left');ax.set_xlabel('Period (days)')
                if keys:ax.legend(fontsize=5,ncol=3)
                else:ax.text(.5,.5,'No usable series with >10 values',ha='center',transform=ax.transAxes)
            matches=read(dest/'activity_matches.json',[]);strong=[r for r in matches if r['proxy_lnBF']>=5]
            text='Nearby local activity peaks (5% period tolerance): '+('; '.join(f"{r['candidate_id']} {r['set_id']} {r['proxy']}: {r['peak_period_days']:.4g} d, proxy ln BF {r['proxy_lnBF']:.1f}" for r in strong[:8]) or 'none meeting the exploratory threshold')
            paragraph(fig,.16,text,size=8,width=115);save(fig)
            if len(strong)>8:
                for start in range(8,len(strong),20):
                    fig=page(f'{code} | Additional activity matches');y=.9
                    for r in strong[start:start+20]:y=paragraph(fig,y,f"{r['candidate_id']} {r['set_id']} {r['proxy']}: {r['peak_period_days']:.6g} d; proxy ln BF {r['proxy_lnBF']:.2f}; period difference {r['relative_period_error']:.2%}",size=8,width=115)
                    save(fig)
            instruments=np.load(dest/'instrument_periodograms.npz');frequency=instruments['frequency']
            window_path=dest/'sampling_window.npz'
            if window_path.exists():window=np.load(window_path)['power']
            else:
                window=np.empty(len(frequency));times=data.t_days.to_numpy()
                for offset in range(0,len(frequency),128):
                    angle=2*np.pi*times[:,None]*frequency[offset:offset+128]
                    window[offset:offset+128]=np.mean(np.cos(angle),axis=0)**2+np.mean(np.sin(angle),axis=0)**2
                np.savez_compressed(window_path,frequency=frequency,power=window)
            fig=page(f'{code} | Instruments, sampling and residuals');axes=fig.subplots(3,1);fig.subplots_adjust(left=.12,right=.94,top=.90,bottom=.09,hspace=.45)
            for key in instruments.files:
                if key!='frequency':draw_scan(axes[0],frequency,instruments[key],key)
            axes[0].set_title('Individual instruments with >10 RV measurements',loc='left');axes[0].set_xlabel('Period (days)')
            if len(instruments.files)>1:axes[0].legend(fontsize=5,ncol=3)
            order=np.argsort(1/frequency);axes[1].plot(1/frequency[order],window[order],lw=.65,rasterized=True);axes[1].set_xscale('log');axes[1].set(xlabel='Period (days)',ylabel='Sampling-window power',title='Normalized sampling spectrum; aliases require review')
            for sid,part in data.groupby('set_id'):axes[2].scatter(part.t_days/365.25,part.residual/part.effective_error,s=4,label=sid,alpha=.6,rasterized=True)
            axes[2].axhline(0,color='gray',lw=.7);axes[2].set(xlabel='Relative time (years)',ylabel='Residual / effective error',title='Residual time series after conditional noise and signal model');save(fig)
            moving=read(dest/'moving.json',[])
            if moving:
                moving=pd.DataFrame(moving)
                for start in range(0,len(candidates),4):
                    fig=page(f'{code} | Moving periodograms');chunk=candidates[start:start+4];cols=1 if len(chunk)==1 else 2;rows=int(np.ceil(len(chunk)/cols));axes=fig.subplots(rows,cols,squeeze=False).ravel();fig.subplots_adjust(left=.12,right=.88,top=.90,bottom=.12,hspace=.40,wspace=.45)
                    for ax,candidate in zip(axes,chunk):
                        subset=moving[moving.candidate_id==candidate['candidate_id']]
                        if subset.empty:ax.text(.5,.5,'Insufficient window coverage',ha='center',transform=ax.transAxes);continue
                        values=subset.pivot(index='center_days',columns='period_days',values='lnBF').sort_index(axis=1)
                        im=ax.pcolormesh(values.columns,values.index/365.25,values.values,shading='nearest',cmap='viridis',rasterized=True);fig.colorbar(im,ax=ax,label='Conditional ln BF',fraction=.05,pad=.03);ax.axvline(candidate['period_days'],color='white',ls='--',lw=.7);ax.set(title=candidate['candidate_id'],xlabel='Period (days)',ylabel='Window center (relative years)')
                    for ax in axes[len(chunk):]:ax.set_axis_off()
                    fig.text(.08,.055,'Seven half-baseline windows; other fitted components subtracted. Noise is fixed at its final full-data estimate.',fontsize=8);save(fig)
            for start in range(0,max(1,len(candidates)),8):
                fig=page(f'{code} | Stability and evidence limits');y=.90
                for c in candidates[start:start+8]:
                    support='not available' if c['window_support'] is None else f"{c['window_support']:.0%}";drop='not available' if c['min_instrument_removal_lnBF'] is None else f"{c['min_instrument_removal_lnBF']:.3g}"
                    y=paragraph(fig,y,f"{c['candidate_id']}: support across seven half-baseline windows = {support}; minimum leave-one-instrument-out ln BF = {drop}. {c['flags']}.")
                lags=read(dest/'lagged_activity.json',[]);selected=[r for r in lags if r['candidate_id'] in {c['candidate_id'] for c in candidates[start:start+8]}];selected.sort(key=lambda r:abs(r['spearman_r']),reverse=True)
                y=paragraph(fig,y,'Exploratory RV/activity lag screen: candidate-plus-residual RV against each available proxy, scanning -10 to +10 days with nearest matches within 2 days. The strongest associations below are selected over many tests and are not calibrated significances.')
                for r in selected[:8]:y=paragraph(fig,y,f"{r['candidate_id']} {r['set_id']} {r['proxy']}: r={r['spearman_r']:.3f}, lag={r['best_lag_days']} d, N={r['n_pairs']}.",size=8)
                paragraph(fig,max(.10,y),'These fits establish conditional periodic RV evidence. They do not establish planetary mass, independent discovery, or confirmed/retracted catalogue status. Missing activity indicators are not evidence against stellar activity. A search stopped by an optimizer or degrees-of-freedom limit remains incomplete.');save(fig)
    # Structural check: readable pages and only this target's coded ID.
    doc=PdfReader(tmp);text='\n'.join(p.extract_text() or '' for p in doc.pages)
    if code not in text:raise ValueError('PDF audit failed: missing code')
    tmp.replace(path)
