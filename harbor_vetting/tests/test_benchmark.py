import copy,json,os,sys
from pathlib import Path
import pytest
import pymupdf
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from grader import score,strict_json,protocol
TASKS=Path(os.environ.get('RV_HARBOR_TASKS',Path(__file__).resolve().parents[3]/'datasets/harbor_rv_recovered'))
CASES=sorted(TASKS.glob('*/*/task.toml'))
assert CASES,'Build the task dataset first, or set RV_HARBOR_TASKS'

def case(task):
 p=task.parent
 return json.loads((p/'solution/answer.json').read_text()),json.loads((p/'tests/reference.json').read_text()),json.loads((p/'tests/submission.schema.json').read_text())

@pytest.mark.parametrize('task',CASES,ids=lambda p:p.parent.name)
def test_oracle_and_pdf_isolation(task):
 a,r,s=case(task);reward,_=score(a,r,s)
 assert reward['reward']==pytest.approx(1)
 assert reward['brier_score']==0
 from harbor.models.task.task import Task
 assert Task(task.parent).config.verifier.environment_mode.value=='separate'
 with pymupdf.open(task.parent/'environment/report.pdf') as d:
  assert len(d)==6
  text=' '.join(p.get_text() for p in d)
  for other in CASES:
   if other!=task:
    other_id=other.parent.name.split('-')[-1].upper()
    assert other_id not in text
  for page in d:
   assert all(page.get_image_rects(x[0]) for x in page.get_images()),'Unrelated image resource retained'
 assert not any((task.parent/'environment').glob('*reference*'))
 assert not any((task.parent/'environment').glob('*answer.json'))

@pytest.mark.parametrize('task',CASES,ids=lambda p:p.parent.name)
def test_flipping_labels_does_not_change_evidence(task):
 a,r,s=case(task)
 for c in a['candidates']:
  c['probability_confirmed']=1-c['probability_confirmed'];c['catalogue_prediction']='confirmed' if c['probability_confirmed'] else 'retracted'
 reward,_=score(a,r,s)
 assert reward['catalogue_accuracy']==0
 assert reward['evidence_score']==pytest.approx(1)
 assert reward['reward']==pytest.approx(.5)
 assert reward['brier_score']==1

@pytest.mark.parametrize('mutation', ['missing','duplicate','wrong_target','bad_probability','bool_number','unknown_flag'])
def test_invalid_output_is_rejected(mutation):
 a,r,s=case(CASES[0])
 if mutation=='missing':a['candidates']=[]
 if mutation=='duplicate':a['candidates'].append(copy.deepcopy(a['candidates'][0]))
 if mutation=='wrong_target':a['target_id']='fake'
 if mutation=='bad_probability':a['candidates'][0]['probability_confirmed']=1-a['candidates'][0]['probability_confirmed']
 if mutation=='bool_number':a['system']['n_rv']=True
 if mutation=='unknown_flag':a['candidates'][0]['flags'].append('planet_definitely_real')
 with pytest.raises(ValueError):score(a,r,s)

def test_fake_quotes_and_numerical_hallucinations_lose_points():
 a,r,s=case(CASES[0])
 for c in a['candidates']:
  c['period_days']=987654321
  for q in c['citations']:q['quote']='This is a completely invented report quotation.'
 reward,d=score(a,r,s)
 assert reward['evidence_score']<1
 assert all(c['citations_valid']==0 and c['numeric_facts']['period_days']==0 for c in d['candidates'])

@pytest.mark.parametrize('text',['{"x":NaN}','{"x":Infinity}','{"x":1e999}','{"x":1,"x":2}'])
def test_nonfinite_and_duplicate_json_rejected(tmp_path,text):
 p=tmp_path/'bad.json';p.write_text(text)
 with pytest.raises(ValueError):strict_json(p)

def test_protocol_does_not_equate_signal_and_catalogue_status():
 assert protocol([])=='conditionally_supported'
 assert protocol(['weak_evidence'])=='unsupported'
 assert protocol(['poor_model_fit'])=='inconclusive'
 assert protocol(['instrument_sensitive'])=='inconclusive'


def test_filtered_membership_and_host_separation():
 m=json.loads((TASKS/'manifest.json').read_text())
 assert {x['target_id'] for x in m if x['split']=='val'}.isdisjoint(x['target_id'] for x in m if x['split']=='test')
 selection=json.loads((TASKS/'recovery_selection.json').read_text())
 expected={x['target_id'] for x in selection['targets'] if x['eligible']}
 assert {x['target_id'] for x in m}==expected
 assert all(x['matched_candidates']==x['n_candidates'] for x in selection['targets'] if x['eligible'])
 assert sum(x['n_candidates'] for x in m)==sum(x['n_candidates'] for x in selection['targets'] if x['eligible'])

def test_fake_candidate_quote_cannot_piggyback_on_valid_same_page():
 a,r,s=case(CASES[0])
 for c in a['candidates']:
  quotes=[]
  for page in (2,4):
   text=' '.join(r['page_text'][str(page)].split())
   # The end-of-page methods/footer is real but not candidate-specific.
   quote=text[-160:]
   assert c['candidate_id'] not in quote
   quotes.append(dict(page=page,quote=quote))
  quotes.append(dict(page=4,quote=c['candidate_id']+' completely invented candidate quotation with false evidence'))
  c['citations']=quotes
 _,details=score(a,r,s)
 assert all(c['citations_valid']==0 for c in details['candidates'])
