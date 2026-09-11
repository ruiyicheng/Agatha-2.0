"""Write the coded, RV-only PDF and audit its extracted text for identity leaks."""
import argparse,json,textwrap,re
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import spearmanr,chi2
from pypdf import PdfReader

plt.rcParams.update({'font.size':9,'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42,'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False})
COLORS=plt.cm.tab20.colors

def page(title,subtitle=''):
    fig=plt.figure(figsize=(8.27,11.69));fig.suptitle(title,x=.07,y=.97,ha='left',fontsize=17,fontweight='bold')
    if subtitle:fig.text(.07,.938,subtitle,fontsize=9,color='#444444')
    return fig

def paragraph(fig,y,text,size=10,width=105):
    lines=[]
    for p in text.split('\n'):
        lines+=textwrap.wrap(p,width=width) if p else ['']
    fig.text(.07,y,'\n'.join(lines),va='top',fontsize=size,linespacing=1.45)
    return y-len(lines)*size/72/11.69*1.5-.014

def table(ax,rows,columns,font=8,widths=None):
    ax.axis('off');t=ax.table(cellText=rows,colLabels=columns,cellLoc='left',colLoc='left',loc='upper center',colWidths=widths)
    t.auto_set_font_size(False);t.set_fontsize(font);t.scale(1,1.45)
    for (r,c),cell in t.get_celld().items():
        cell.set_linewidth(.3);cell.set_edgecolor('#cccccc')
        if r==0:cell.set_facecolor('#e9eef4');cell.set_text_props(weight='bold')
    return t

def save(pdf,fig,page_no):
    fig.text(.07,.025,'AGATHA  |  RV evaluation  |  Coded identities',fontsize=8,color='#666666')
    fig.text(.93,.025,str(page_no),ha='right',fontsize=8);pdf.savefig(fig);plt.close(fig)

def plotper(ax,f,p,label=None,color=None):
    order=np.argsort(1/f);ax.plot(1/f[order],p[order],lw=.65,label=label,color=color,rasterized=True);ax.set_xscale('log');ax.set_xlim(min(1/f),max(1/f));ax.set_ylabel('conditional ln BF');ax.axhline(5,ls='--',lw=.7,color='#999999')

def main():
    a=argparse.ArgumentParser();a.add_argument('output',type=Path);args=a.parse_args();out=args.output.resolve()
    tasks=pd.read_csv(out/'tasks.csv');keys=pd.read_csv(out/'private_identity_key.csv');private=pd.read_csv(out/'private_candidates.csv');est=pd.read_csv(out/'runtime_estimates.csv')
    summaries=[];allc=[];allnoise=[];correlations=[];residualtests=[]
    for task in tasks.itertuples():
        dest=out/'results'/task.code;summary=json.loads((dest/'summary.json').read_text());summary['total_attempt_seconds']=summary.get('elapsed_seconds_including_retries',summary['elapsed_seconds']);summaries.append(summary)
        c=pd.read_csv(dest/'candidate_results.csv');c.insert(0,'code',task.code)
        activity=pd.read_csv(dest/'activity.csv');moving=pd.read_csv(dest/'moving_periodograms.csv');cross=pd.read_csv(dest/'leave_one_instrument_out.csv') if (dest/'leave_one_instrument_out.csv').stat().st_size>2 else pd.DataFrame()
        for j,row in c.iterrows():
            resolution=1/summary['baseline_days'];matches=activity[(activity.peak_logBF_BIC>5)&(abs(1/activity.peak_period_days-1/row.fitted_period_days)<resolution)]
            windows=moving[moving.candidate_id==row.candidate_id].groupby('window').at_fitted_period.first()
            c.loc[j,'activity_peak_match']=len(matches)>0;c.loc[j,'moving_windows']=len(windows);c.loc[j,'moving_fraction_lnBF_gt5']=np.mean(windows>5) if len(windows) else np.nan
            drops=cross[cross.candidate_id==row.candidate_id].logBF_BIC_at_period if len(cross) else []
            c.loc[j,'min_leave_one_out_lnBF']=min(drops) if len(drops) else np.nan
            flags=[]
            if not summary['optimizer_converged']:flags.append('full fit not converged')
            if summary['conditional_chi2']/summary['n']>2:flags.append('poor full-model residuals')
            if row.delta_BIC_drop<=10:flags.append('weak conditional evidence')
            if row.optimizer_boundary:flags.append('parameter boundary')
            if row.cycles_observed<2:flags.append('less than 2 cycles')
            if len(matches):flags.append('activity peak nearby')
            if len(windows) and np.mean(windows>5)<.5:flags.append('limited time stability')
            if len(drops) and min(drops)<=5:flags.append('instrument sensitivity')
            if row.unseeded:flags.append('catalogue period missing')
            if not row.reduced_fit_converged:flags.append('reduced fit not converged')
            c.loc[j,'review_flags']='; '.join(flags) or 'no threshold flags'
        c.to_csv(dest/'diagnostic_flags.csv',index=False);allc.append(c)
        noise=pd.read_csv(dest/'noise_refined.csv');noise.insert(0,'code',task.code);allnoise.append(noise)
        fit=pd.read_csv(dest/'fit.csv');ins=pd.read_csv(Path(task.input_dir)/'instruments.csv')
        for row in ins.itertuples():
            raw=pd.read_csv(row.file);part=fit[fit.set_id==row.set_id];z=part.residual.to_numpy()/part.effective_error.to_numpy();lag=min(10,max(1,len(z)//5))
            if len(z)>10:
                ac=np.array([np.corrcoef(z[k:],z[:-k])[0,1] for k in range(1,lag+1)]);Q=len(z)*(len(z)+2)*np.sum(ac**2/(len(z)-np.arange(1,lag+1)))
                residualtests.append(dict(code=task.code,set_id=row.set_id,n=len(z),lag1=ac[0],ljung_box_Q=Q,lags=lag,p_approx=chi2.sf(Q,lag)))
            for proxy in ['sindex','halpha']:
                valid=np.isfinite(raw[proxy]);nn=sum(valid)
                if nn>10 and raw.loc[valid,proxy].std()>0:
                    r,p=spearmanr(part.residual.to_numpy()[valid],raw.loc[valid,proxy])
                    correlations.append(dict(code=task.code,set_id=row.set_id,proxy=proxy,n=nn,spearman_r=r,p_unadjusted=p))
    sums=pd.DataFrame(summaries);cands=pd.concat(allc,ignore_index=True);noises=pd.concat(allnoise,ignore_index=True);corr=pd.DataFrame(correlations);residtests=pd.DataFrame(residualtests)
    noises['delta_BIC_from_best']=noises.BIC-noises.groupby(['code','set_id']).BIC.transform('min')
    gaps=[]
    for _,group in noises.groupby(['code','set_id']):
        scores=group.BIC.dropna().sort_values()
        if len(scores)>1:gaps.append(scores.iloc[1]-scores.iloc[0])
    sums.to_csv(out/'run_summary_anonymized.csv',index=False);cands.to_csv(out/'candidate_diagnostics_anonymized.csv',index=False);noises.to_csv(out/'noise_comparison_anonymized.csv',index=False);corr.to_csv(out/'activity_correlations_anonymized.csv',index=False);residtests.to_csv(out/'residual_tests_anonymized.csv',index=False)
    runtime=keys.merge(est,on='code').merge(sums[['code','elapsed_seconds','total_attempt_seconds']],on='code');runtime['actual_science_minutes']=(runtime.total_attempt_seconds+runtime.noise_initial_seconds)/60;runtime.to_csv(out/'runtime_by_target_PRIVATE.csv',index=False)
    rt=runtime.drop(columns=['target','host_id']);rt.to_csv(out/'runtime_anonymized.csv',index=False)
    status=json.loads((out/'test_status.json').read_text())
    pdfpath=out/'agatha_rv_report_anonymized.pdf';page_no=0
    with PdfPages(pdfpath,metadata={'Title':'Coded radial-velocity evaluation report','Author':'','Subject':'RV-only diagnostics; no MCMC','Keywords':'RV, periodogram, noise models'}) as pdf:
        fig=page('RV evaluation report',f'{len(sums)} coded systems  |  {len(cands)} fitted planet components  |  Validation + test sets')
        y=paragraph(fig,.89,f'All {len(sums)} systems completed the headless analysis. The validation and test sets each contain 9 confirmed and 9 retracted catalogue entries, with disjoint host systems. Only RV measurements and spectroscopic activity indicators were used. No photometry or MCMC was used.',size=11,width=91)
        y=paragraph(fig,y,f'The source manifests contain {int(sums.n_raw.sum()):,} RV rows; {int(sums.n.sum()):,} remain after duplicate removal. Each system retains more than 100 measurements. Every full model contains exactly the number of selected catalogue candidates for that system, including weak or unsupported components. A fitted component is not a planet confirmation. {", ".join(sums.loc[sums.conditional_chi2/sums.n>2,"code"])} have poor full-model residuals despite numerical convergence; their orbital elements and evidence scores require particular caution.')
        y=paragraph(fig,y,'Identities are coded throughout this PDF. Target names, sky coordinates, absolute observing epochs, source paths and individual catalogue truth labels are omitted. The identity key and labelled candidate table are separate private files. Numerical periods and RV patterns can still be recognizable to a specialist; this is identity masking, not a guarantee against re-identification.')
        selected=noises[noises.selected];counts=selected.model.value_counts().to_dict()
        y=paragraph(fig,y,'Final instrument noise choices: '+', '.join(f'{k}: {v}' for k,v in counts.items())+f'. {sum(np.array(gaps)<2)}/{len(gaps)} comparisons have a runner-up within 2 BIC units, so these choices are weakly distinguished. {int((cands.delta_BIC_drop>10).sum())}/{len(cands)} components have conditional drop-component delta BIC > 10. This threshold alone does not establish planetary origin. A retracted catalogue planet can retain a real RV signal, for example if its companion is stellar or substellar; these diagnostics do not determine companion mass class.')
        y=paragraph(fig,y,'The following pages provide per-target timing, the test results, assumptions, and five diagnostic pages for each coded system: data/noise, sequential RV scans, planet fits, activity/window/instrument scans, and time stability plus sensitivity checks.')
        paragraph(fig,y,'Source: github.com/phillippro/Agatha-2.0, upstream commit fef88699411e3afbe2ab3ba711f1fec8a173ee85. The batch extension is in github.com/ruiyicheng/Agatha-2.0 on branch rv-vetting-headless-report. Results reflect the supplied RV dataset and selected catalogue snapshot.',size=9)
        page_no+=1;save(pdf,fig,page_no)
        fig=page('Runtime for every target','CPU execution; six concurrent target workers; one numerical thread per worker')
        rows=[]
        for r in runtime.itertuples():rows.append([r.code,r.split,str(r.measurements),str(r.components),f'{r.estimated_minutes:.1f}',f'{r.low_minutes:.1f}–{r.high_minutes:.1f}',f'{r.actual_science_minutes:.2f}'])
        table(fig.add_axes([.07,.38,.86,.52]),rows,['Code','Split','RV N','Planets','Est. min','Range min','Actual min'],font=9)
        paragraph(fig,.34,'Estimates were recorded before the batch using a 128-frequency timing sample for each system, the measured first noise-selection cost, and allowances for fitting and diagnostics. The range is half to twice the central estimate; it is a planning allowance, not a statistical interval.')
        paragraph(fig,.23,'Actual time includes both noise selections and all per-target science calculations, including recorded retries. It excludes dependency installation, repository tests, preprocessing, report assembly and time waiting for a worker. Timings are specific to this machine and workload; sums are CPU-job elapsed times, not batch wall time.')
        page_no+=1;save(pdf,fig,page_no)
        fig=page('Methods and interpretation','Read before using the conditional evidence or fitted orbital elements')
        paragraphs=[
          'Noise selection. For each instrument, compare white noise plus fitted jitter (WN), AR(1), MA(1), and MA(2) with exponential lag decay, an offset and a trend. Select minimum BIC using Agatha’s likelihood and three optimizer restarts. Sets with at most 10 measurements receive WN only. Repeat selection on RV minus a preliminary full planet model, then hold hyperparameters fixed for the final calculations. This two-pass local procedure is not a joint global noise-and-planet search. Agatha’s AR term uses lagged observed RV, whereas its MA term uses lagged deterministic-model residuals; AR can absorb true signals and compromise shared amplitudes. GP noise is tested in software but is not in the science model grid.',
          'Signal model. Sum independent Keplerians with one offset per instrument and a shared linear trend. Known catalogue periods initialize fits within +/-10%; eccentricity is bounded at 0.85. Three starts are tried. The one candidate without a catalogue period is initialized from the RV residual scan. The full component count is fixed by the selected entries; no truth label enters model selection or fitting. Interacting systems are not integrated dynamically. Point estimates and local convergence are reported; no posterior uncertainties are claimed.',
          'RV scans. Uniform frequency spacing is at most 1/(4 x baseline), from 1 day to max(2 x baseline, 1.25 x longest catalogue period). Selected jitter and lag coefficients/timescales are fixed; harmonic amplitudes, instrument offsets and trend are refitted at every frequency. The plotted conditional ln BF is delta log-likelihood minus ln N (two sinusoid coefficients). It is a BIC approximation conditional on estimated noise, not a marginalized Bayes factor or a calibrated false-alarm probability. Selected components are removed in descending fitted amplitude; this is a catalogue-informed residual diagnosis, not a blind discovery search.',
          'Candidate evidence. Drop each planet and reoptimize the remaining orbital models; delta BIC = chi-squared(reduced) - chi-squared(full) - 5 ln N. Noise is fixed in both fits. Larger positive values favor inclusion. Fixed noise, search bounds and local minima can make this optimistic. Boundary fits, fewer than two observed cycles, nearby activity peaks, weak time stability and sensitivity to instrument removal are flagged.',
          'Other diagnostics. S-index and H-alpha are analyzed separately per instrument when more than 10 valid, nonconstant values exist. They are standardized and use unit errors because uniform proxy uncertainties are unavailable. No missing values are invented. The spectral window is normalized squared complex sampling amplitude. Seven windows span half the time baseline; period neighborhoods are 0.8–1.2 of each fitted period. Instrument-removal scans keep other fitted components fixed. Residual correlations and Ljung–Box p-values are exploratory and unadjusted for multiple testing.'
        ]
        y=.895
        for txt in paragraphs:y=paragraph(fig,y,txt,size=9,width=116)
        page_no+=1;save(pdf,fig,page_no)
        fig=page('Verification and data handling','Repository checks plus independent likelihood and input regressions')
        table(fig.add_axes([.07,.62,.86,.27]),[[r['test'],r['status'],('—' if r.get('seconds') is None else f"{r.get('seconds',0):.1f}")] for r in status],['Check','Status','Seconds'],font=9,widths=[.62,.22,.16])
        y=paragraph(fig,.55,'The two MCMC test scripts are intentionally skipped. The GP tests exercise covariance and deterministic fitting without running chains. The Python conditional likelihood is checked against R for WN, AR1, MA1 and MA2; an eccentric two-planet injection verifies the joint signal fit. Added R regressions check the critical-damping SHO limit and preservation of real activity values and simultaneous distinct measurements.')
        y=paragraph(fig,y,'Input handling: first three numeric fields are time, RV and RV error. Headerless files retain their first row. Activity columns use explicit per-format mappings recorded in the preprocessing audit. Reduced Julian dates receive the dataset’s documented 2,400,000-day offset. Full duplicate RV series and repeated epochs within an instrument family are removed and counted; instruments remain distinct and observations are not nightly binned. Instrument families with pre/post-upgrade files prefer those over overlapping unsplit records.')
        paragraph(fig,y,'Limits: RV values and errors use the source files’ m/s convention; no external recalibration is applied. The tests cannot establish that every source file has correct provenance or units. The largest system has dense sampling, so temporal correlation and alias structure deserve particular care. The PDF does not use withheld per-candidate catalogue labels to score model performance.')
        page_no+=1;save(pdf,fig,page_no)
        for summary in summaries:
            code=summary['code'];dest=out/'results'/code;fit=pd.read_csv(dest/'fit.csv');c=cands[cands.code==code].reset_index(drop=True);noise=noises[noises.code==code];sel=noise[noise.selected];n=len(c);ids=fit.set_id.unique();color={sid:COLORS[i%20] for i,sid in enumerate(ids)}
            fig=page(f'{code}  |  Data and noise',f"{summary['split']}  ·  {summary['n']:,} RV measurements  ·  {summary['n_signals']} Keplerians  ·  {summary['baseline_days']/365.25:.1f} years")
            rows=[]
            for r in sel.itertuples():
                alternatives=noise[(noise.set_id==r.set_id)&np.isfinite(noise.BIC)].sort_values('BIC');gap=alternatives.BIC.iloc[1]-r.BIC if len(alternatives)>1 else np.nan
                tau=r.tau if r.Nma else r.tauAR if r.Nar else np.nan
                rows.append([r.set_id,str(r.n),r.model,f'{r.jitter:.3g}',f'{tau:.3g}' if np.isfinite(tau) else '—',f'{gap:.2f}' if np.isfinite(gap) else 'WN only'])
            table(fig.add_axes([.07,.54,.86,.36]),rows,['Instrument','N','Noise','Jitter m/s','Decay d','BIC gap'],font=8)
            ax=fig.add_axes([.1,.315,.83,.20]);ar=fig.add_axes([.1,.13,.83,.14],sharex=ax)
            for sid in ids:
                q=fit[fit.set_id==sid];ax.scatter(q.t_days,q.rv_centered,s=4,color=color[sid],rasterized=True);ax.plot(q.t_days,q.model,lw=.6,color=color[sid]);ar.scatter(q.t_days,q.residual,s=4,color=color[sid],rasterized=True)
            ax.set_ylabel('Centered RV (m/s)');ax.set_title('RV and deterministic full model; instrument median removed');ar.set_ylabel('Innovation (m/s)');ar.set_xlabel('Days from first observation');ar.axhline(0,lw=.7,color='gray')
            fig.text(.1,.535,'Noise-table BIC gap: runner-up minus winner; a gap below 2 weakly distinguishes the models.',fontsize=7)
            fig.text(.1,.065,f"RMS {summary['residual_rms']:.3g} m/s | chi-squared / N = {summary['conditional_chi2']/summary['n']:.2f} | converged: {summary['optimizer_converged']}",fontsize=8)
            page_no+=1;save(pdf,fig,page_no)
            fig=page(f'{code}  |  RV periodograms','Conditional BIC evidence; grey dashed line is ln BF = 5, not a calibrated false-alarm level')
            spec=np.load(dest/'periodograms.npz');f=spec['frequency'];order=spec['removal_order'];axes=fig.subplots(n+1,1,squeeze=False).ravel();fig.subplots_adjust(left=.1,right=.95,top=.90,bottom=.07,hspace=.42)
            for k,ax in enumerate(axes):
                plotper(ax,f,spec['sequence'][k],color='#245a81')
                if k==0:plotper(ax,f,spec['white'],label='WN, same jitter',color='#c77d37');ax.legend(fontsize=7,loc='upper right')
                for row in c.itertuples():ax.axvline(row.fitted_period_days,color='#b6b6b6',lw=.5)
                title='Full RV scan' if k==0 else f"After removing {', '.join(c.candidate_id.iloc[order[:k]])}"
                ax.set_title(title,loc='left',fontsize=8)
            axes[-1].set_xlabel('Period (days)');page_no+=1;save(pdf,fig,page_no)
            fig=page(f'{code}  |  Full planet model','All requested components are retained, including weak or boundary solutions')
            rows=[[r.candidate_id,f'{r.fitted_period_days:.5g}',f'{r.K_m_s:.3g}',f'{r.eccentricity:.3f}',f'{r.delta_BIC_drop:.1f}',f'{r.cycles_observed:.1f}'] for r in c.itertuples()]
            table(fig.add_axes([.07,.72,.86,.18]),rows,['Candidate','Period d','K m/s','e','Drop ΔBIC','Cycles'],font=8)
            nrows=int(np.ceil(n/2));gs=fig.add_gridspec(nrows,2,left=.1,right=.95,top=.68,bottom=.12,hspace=.4,wspace=.3)
            for j,r in enumerate(c.itertuples()):
                ax=fig.add_subplot(gs[j//2,j%2]);phase=(fit.t_days/r.fitted_period_days)%1;comp=fit['component_'+str(j+1)];o=np.argsort(phase)
                for sid in ids:
                    keep=fit.set_id==sid;ax.scatter(phase[keep],(comp+fit.residual)[keep],s=5,color=color[sid],alpha=.6,rasterized=True)
                ax.plot(phase.iloc[o],comp.iloc[o],color='black',lw=1);ax.set_title(r.candidate_id+(' · boundary' if r.optimizer_boundary else ''));ax.set_xlabel('Orbital phase');ax.set_ylabel('Conditional RV (m/s)')
            fig.text(.08,.065,'Phase panels subtract other planets and the fitted trend, offsets and conditional red-noise prediction.\nPeriod bounds: ±10% around initialization; eccentricity ≤ 0.85. No posterior error bars.',fontsize=8)
            page_no+=1;save(pdf,fig,page_no)
            fig=page(f'{code}  |  Activity, sampling and instruments','Proxy periodograms are standardized per instrument; missing or constant series are omitted')
            acts=np.load(dest/'activity_periodograms.npz');inst=np.load(dest/'instrument_periodograms.npz');window=np.load(dest/'window.npz');gs=fig.add_gridspec(4,1,left=.1,right=.94,top=.89,bottom=.10,hspace=.47)
            for k,proxy in enumerate(['sindex','halpha']):
                ax=fig.add_subplot(gs[k]);available=[name for name in acts.files if name.endswith('__'+proxy)]
                for name in available:plotper(ax,acts['frequency'],acts[name],label=name.split('__')[0],color=color[name.split('__')[0]])
                ax.set_title('S-index' if proxy=='sindex' else 'H-alpha',loc='left')
                if available:ax.legend(fontsize=6,ncol=4,loc='upper right')
                else:ax.text(.5,.5,'No usable series with >10 values',ha='center',transform=ax.transAxes);ax.set_axis_off()
            ax=fig.add_subplot(gs[2]);o=np.argsort(1/window['frequency']);ax.plot(1/window['frequency'][o],window['power'][o],lw=.6,rasterized=True);ax.set_xscale('log');ax.set_ylabel('Window power');ax.set_title('Sampling window',loc='left')
            ax=fig.add_subplot(gs[3]);
            for name in inst.files:
                if name!='frequency':plotper(ax,inst['frequency'],inst[name],label=name,color=color[name])
            ax.set_title('RV scans by instrument (>10 measurements)',loc='left');ax.legend(fontsize=6,ncol=4,loc='upper right');ax.set_xlabel('Period (days)')
            fig.text(.08,.055,'Activity coincidences are review flags. Unit-error proxy scans and RV conditional evidence have different scales.',fontsize=8)
            page_no+=1;save(pdf,fig,page_no)
            fig=page(f'{code}  |  Stability and review flags','Seven half-baseline windows; other planet components subtracted before each scan')
            moving=pd.read_csv(dest/'moving_periodograms.csv');nr=int(np.ceil(n/2));gs=fig.add_gridspec(nr,2,left=.10,right=.94,top=.89,bottom=.48,hspace=.55,wspace=.4)
            for j,r in enumerate(c.itertuples()):
                ax=fig.add_subplot(gs[j//2,j%2]);q=moving[moving.candidate_id==r.candidate_id]
                if len(q):
                    pivot=q.pivot(index='center_days',columns='period_days',values='logBF_BIC');im=ax.pcolormesh(pivot.columns,pivot.index/365.25,pivot.values,shading='nearest',cmap='viridis',rasterized=True);fig.colorbar(im,ax=ax,pad=.02,fraction=.04);ax.axvline(r.fitted_period_days,color='white',lw=.8,ls='--')
                else:ax.text(.5,.5,'Insufficient window coverage',ha='center',transform=ax.transAxes)
                ax.set_title(r.candidate_id,fontsize=9);ax.set_xlabel('Period d' if j//2==nr-1 else '');ax.set_ylabel('Window center (yr)')
            y=.43
            for r in c.itertuples():y=paragraph(fig,y,f'{r.candidate_id}: {r.review_flags}. Window support: {r.moving_fraction_lnBF_gt5:.0%}; minimum instrument-removal ln BF: {r.min_leave_one_out_lnBF:.2g}.',size=8,width=123)
            rr=residtests[residtests.code==code];cc=corr[corr.code==code] if len(corr) else corr
            txt=f'Residual serial-correlation screen: {sum(rr.p_approx<.01)}/{len(rr)} tested instruments have approximate Ljung–Box p < 0.01.'
            if len(cc):txt+=f' Residual/activity correlation: {sum(cc.p_unadjusted<.01)}/{len(cc)} tests have unadjusted p < 0.01.'
            paragraph(fig,min(y,.13),txt+' These exploratory tests ignore irregular cadence and multiple testing; inspect the data before interpreting significance.',size=8,width=123)
            page_no+=1;save(pdf,fig,page_no)
    extracted='\n'.join(p.extract_text() or '' for p in PdfReader(pdfpath).pages)
    names=set(keys.target)|set(keys.host_id)|set(private.planet_name)
    if (out/'private_aliases.csv').exists():
        for aliases in pd.read_csv(out/'private_aliases.csv').aliases.fillna(''):
            names.update(aliases.split(';'))
    # Normalize punctuation/spacing to catch typical catalogue name variants.
    normalize=lambda s:re.sub('[^a-z0-9]','',s.lower())
    textnorm=normalize(extracted);leaks=[name for name in names if len(normalize(name))>=5 and normalize(name) in textnorm]
    if leaks:raise RuntimeError('Identity audit failed: '+str(leaks))
    assert '/home/' not in extracted  # absolute paths are forbidden
    audit=dict(pages=page_no,targets=len(sums),components=len(cands),name_leaks=leaks,absolute_paths_present='/home/' in extracted,pdf=str(pdfpath),test_checks=status)
    (out/'report_audit.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit,indent=2))
if __name__=='__main__':main()
