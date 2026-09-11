"""Deterministic Harbor verifier. Run only after the agent, preferably isolated."""
import argparse,json,math,re,unicodedata
from pathlib import Path
from jsonschema import Draft202012Validator

def normalize(text):return ' '.join(unicodedata.normalize('NFKC',text).split())
def strict_json(path):
    def pairs(items):
        out={}
        for k,v in items:
            if k in out:raise ValueError('duplicate JSON key: '+k)
            out[k]=v
        return out
    def bad_constant(s):raise ValueError('Non-finite JSON number: '+s)
    if path.stat().st_size>1_000_000:raise ValueError('Submission exceeds 1 MB')
    result=json.loads(path.read_text(),object_pairs_hook=pairs,parse_constant=bad_constant)
    def finite(v):
        if isinstance(v,float) and not math.isfinite(v):raise ValueError('Non-finite number')
        if isinstance(v,dict):
            for x in v.values():finite(x)
        if isinstance(v,list):
            for x in v:finite(x)
    finite(result);return result

def flag_f1(actual,expected):
    actual,expected=set(actual),set(expected)
    return 1. if not actual and not expected else 2*len(actual&expected)/(len(actual)+len(expected))

def protocol(flags):
    return 'unsupported' if 'weak_evidence' in flags else 'inconclusive' if flags else 'conditionally_supported'

def score(submission,reference,schema):
    errors=sorted(Draft202012Validator(schema).iter_errors(submission),key=lambda e:str(e.path))
    if errors:raise ValueError('Schema: '+errors[0].message)
    if submission['target_id']!=reference['target_id']:raise ValueError('Wrong target ID')
    candidates=submission['candidates'];expected={x['candidate_id']:x for x in reference['candidates']}
    ids=[x['candidate_id'] for x in candidates]
    if len(ids)!=len(set(ids)) or set(ids)!=set(expected):raise ValueError('Candidate set must match exactly')
    noise=submission['system']['noise_models'];noise_ids=[x['instrument'] for x in noise]
    if len(noise_ids)!=len(set(noise_ids)) or set(noise_ids)!=set(reference['noise_models']):raise ValueError('Instrument set must match exactly')
    labels=[];facts=[];flags=[];assessments=[];citations=[];briers=[];detail=[]
    for actual in candidates:
        expected_row=expected[actual['candidate_id']]
        if actual['catalogue_prediction']!=('confirmed' if actual['probability_confirmed']>=.5 else 'retracted'):
            raise ValueError('Probability and catalogue prediction disagree')
        labels.append(float(actual['catalogue_prediction']==expected_row['catalogue_label']))
        briers.append((actual['probability_confirmed']-(expected_row['catalogue_label']=='confirmed'))**2)
        matches={}
        for field in ('period_days','K_m_s','eccentricity','delta_BIC','cycles_observed'):
            matches[field]=float(abs(actual[field]-expected_row['facts'][field])<=expected_row['tolerance'][field])
        facts.append(sum(matches.values())/len(matches))
        flags.append(flag_f1(actual['flags'],expected_row['flags']))
        assessments.append(float(actual['rv_assessment']==expected_row['rv_assessment']))
        valid=[];valid_quotes=[]
        for cite in actual['citations']:
            page=cite['page'];quote=normalize(cite['quote']);text=normalize(reference['page_text'][str(page)])
            if len(quote)>=20 and quote in text:valid.append(page);valid_quotes.append(quote)
        has_candidate_quote=any(actual['candidate_id'] in quote for quote in valid_quotes)
        citations.append(float(len(set(valid))>=2 and has_candidate_quote))
        detail.append(dict(candidate_id=actual['candidate_id'],catalogue_correct=labels[-1],numeric_facts=matches,
                           flag_f1=flags[-1],assessment_correct=assessments[-1],citations_valid=citations[-1],
                           rationale=actual['rationale']))
    sys=submission['system'];truth=reference['system']
    sys_checks=[sys['n_rv']==truth['n_rv'],sys['n_components']==truth['n_components'],sys['n_instruments']==truth['n_instruments'],
                abs(sys['chi2_per_measurement']-truth['chi2_per_measurement'])<=.011,sys['fit_adequacy']==truth['fit_adequacy']]
    system_score=sum(sys_checks)/len(sys_checks)
    noise_score=sum(x['model']==reference['noise_models'][x['instrument']] for x in noise)/len(noise)
    mean=lambda x:sum(x)/len(x)
    evidence_score=.15*system_score+.15*noise_score+.25*mean(facts)+.20*mean(flags)+.15*mean(assessments)+.10*mean(citations)
    rewards=dict(reward=.5*mean(labels)+.5*evidence_score,catalogue_accuracy=mean(labels),evidence_score=evidence_score,
                 assessment_accuracy=mean(assessments),brier_score=mean(briers),valid_submission=1.,candidate_count=len(candidates))
    return rewards,dict(system_score=system_score,noise_score=noise_score,n_candidates=len(candidates),
                        candidates=detail,note='Rationale is retained, not semantically graded. Protocol is not planetary ground truth.')

def main():
    p=argparse.ArgumentParser();p.add_argument('--submission',type=Path,default=Path('/workspace/output/verdict.json'));p.add_argument('--reference',type=Path,default=Path('/tests/reference.json'));p.add_argument('--schema',type=Path,default=Path('/tests/submission.schema.json'));p.add_argument('--output',type=Path,default=Path('/logs/verifier'));a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True)
    reference=strict_json(a.reference)
    try:
        rewards,details=score(strict_json(a.submission),reference,strict_json(a.schema))
    except Exception as exc:
        rewards=dict(reward=0.,catalogue_accuracy=0.,evidence_score=0.,assessment_accuracy=0.,brier_score=1.,valid_submission=0.,candidate_count=len(reference['candidates']))
        details={'error':type(exc).__name__+': '+str(exc)}
    (a.output/'reward.json').write_text(json.dumps(rewards,allow_nan=False)+'\n')
    (a.output/'details.json').write_text(json.dumps(details,indent=2,allow_nan=False)+'\n')
    print(json.dumps(rewards))
if __name__=='__main__':main()
