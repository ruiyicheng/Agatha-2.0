import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from recovery import match_recovery

F=np.linspace(.01,.2,1001)
def peak(period,height=20):return height*np.exp(-((F-1/period)/.0002)**2)
def cat(period,name='A'):return {'candidate_id':name,'period_days':period}

def test_period_constrained_fit_is_not_a_recovery_signal():
    ok,a=match_recovery(F,np.zeros((2,len(F))),[0],[cat(20)],1000)
    assert not ok and not a[0]['matched']

def test_all_candidates_required_and_missing_period_excluded():
    seq=np.array([peak(20),np.zeros(len(F)),np.zeros(len(F))])
    ok,a=match_recovery(F,seq,[0,1],[cat(20),cat('', 'B')],1000)
    assert not ok and a[1]['reason']=='missing_catalogue_period'

def test_one_peak_cannot_count_as_two_signals():
    y=peak(20);seq=np.array([y,y,np.zeros(len(F))])
    ok,a=match_recovery(F,seq,[0,1],[cat(20),cat(20.01,'B')],1000)
    assert not ok and sum(x['matched'] for x in a)==1

def test_signal_after_subtraction_does_not_count_for_removed_candidate():
    seq=np.array([np.zeros(len(F)),peak(20),peak(30)])
    ok,a=match_recovery(F,seq,[0,1],[cat(20),cat(30,'B')],1000)
    assert not ok and not a[0]['matched']

def test_distinct_recovered_signals_and_label_independence():
    y=peak(20)+peak(35);seq=np.array([y,peak(35),np.zeros(len(F))]);c=[cat(20),cat(35,'B')]
    ok,a=match_recovery(F,seq,[0,1],c,1000)
    assert ok
    c[0]['catalogue_label']='retracted';c[1]['catalogue_label']='confirmed'
    assert match_recovery(F,seq,[0,1],c,1000)==(ok,a)

def test_significance_tolerance_and_rank_are_enforced():
    seq=np.array([peak(20,4),np.zeros(len(F))])
    assert not match_recovery(F,seq,[0],[cat(20)],1000)[0]
    assert match_recovery(F,seq,[0],[cat(20)],1000,min_logbf=None)[0]
    seq=np.array([peak(20),np.zeros(len(F))])
    assert not match_recovery(F,seq,[0],[cat(22)],1000,period_tolerance=.01)[0]
    seq=np.array([peak(20)+peak(35,30),np.zeros(len(F))])
    assert not match_recovery(F,seq,[0],[cat(20)],1000,top_peaks=1)[0]
