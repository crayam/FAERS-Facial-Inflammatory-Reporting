"""Read-only source verification and independent TTO-only recalculation."""
from pathlib import Path
from datetime import date
from collections import Counter
import csv, gzip, hashlib, json, math, re, sqlite3, sys
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
OLD=ROOT/'outputs/primary11_tto_v1_0'
OUT=ROOT/'outputs/tto_correction_v1_1'
WORK=Path(__file__).parent
DRUGS=['dupilumab','tralokinumab','upadacitinib']
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(p,rows):
    with p.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def equal(a,b):
    if a in [None,''] or b in [None,'']:assert a in [None,''] and b in [None,''],(a,b)
    elif isinstance(a,(int,float,np.number)):assert math.isclose(float(a),float(b),rel_tol=1e-13,abs_tol=1e-13),(a,b)
    else:assert str(a)==str(b),(a,b)
oldhash={line.split('  ',1)[1]:line.split('  ',1)[0] for line in (OLD/'checksums_v1.1.sha256').read_text().splitlines()}
for name in ['case_TTO_audit.csv','retained_date_evidence.json','analysis_results.json','TTO_summary.csv','TTO_distribution.csv','exclusion_flow.csv','extraction_manifest_v1.1.json']:
    assert sha(OLD/name)==oldhash[name],name
evidence=json.loads((OLD/'retained_date_evidence.json').read_text())
historical=read(OLD/'case_TTO_audit.csv');h={r['caseid']:r for r in historical}
meta=json.loads((OLD/'extraction_manifest_v1.1.json').read_text())
oldresults=json.loads((OLD/'analysis_results.json').read_text())
protected=[OLD/n for n in oldhash if (OLD/n).exists()]+[ROOT/'work/tto_primary11/analyse.py',ROOT/'outputs/final_closeout_v1_1/Results_Evidence_Pack.md',ROOT/'outputs/final_closeout_v1_1/Methods_Fact_Sheet_v1.1.md']
protected+=list((ROOT/'outputs/supplementary_results_v1_0').glob('*'))
for folder in ['firth_primary11_main_v1_0','firth_primary_narrow_broad_v1_0','final_statistical_results_v1_0/models','restriction_sensitivity_v1_0/models']:
    protected+=list((ROOT/'outputs'/folder).rglob('fit_results.json'))
protected += [ROOT/'outputs/final_statistical_results_v1_0/crude_comparisons.csv',ROOT/'outputs/prr_pt_composition_v1_0/prr_results.csv']
snapshot={str(p.relative_to(ROOT)):sha(p) for p in protected if p.is_file()}
(WORK/'protected_hashes.json').write_text(json.dumps(snapshot,indent=2))
cohort_path=ROOT/'outputs/all_outcomes_v1_0/case_outcome_flags.csv.gz'
links_path=ROOT/'outputs/primary11_v1_0/eligible_drug_indication_links.csv.gz'
for p in [cohort_path,links_path]:assert sha(p)==next(x['sha256'] for x in meta['sources'] if x['path']==str(p))
with gzip.open(cohort_path,'rt',encoding='utf-8-sig') as f:frozen={r['caseid']:r for r in csv.DictReader(f) if r['Primary11']=='1'}
with gzip.open(links_path,'rt',encoding='utf-8-sig') as f:
    frozen_links=[r for r in csv.DictReader(f) if r['caseid'] in frozen]
assert len(historical)==len(h)==len(evidence)==len(frozen)==1108
assert set(h)==set(frozen)=={a['cohort']['caseid'] for a in evidence}
db=sqlite3.connect(f'file:{ROOT}/work/primary11/latest.sqlite?mode=ro',uri=True)
def parse(raw):
    if not raw:return 'missing',None
    if re.fullmatch(r'[0-9]{4}|[0-9]{6}',raw):return 'incomplete',None
    if not re.fullmatch(r'[0-9]{8}',raw):return 'invalid',None
    y,m,d=int(raw[:4]),int(raw[4:6]),int(raw[6:8])
    if m==0 or d==0:return 'incomplete',None
    try:return 'complete',date(y,m,d)
    except ValueError:return 'invalid',None
verified=[];valid=[];statuses=Counter()
for a in evidence:
    c=a['cohort'];cid=c['caseid'];old=h[cid];demo=a['demo'];drug=c['drug'];q=c['source_quarter']
    for k,v in c.items():assert frozen[cid][k]==v==old[k],(cid,k)
    assert db.execute('SELECT pid,ver,fda,q FROM latest WHERE cid=?',(int(cid),)).fetchone()==(int(c['primaryid']),int(c['caseversion']),int(c['fda_dt']),q)
    for k in ['caseid','primaryid','caseversion','fda_dt','init_fda_dt']:assert demo[k]==c[k]
    fl=[x for x in frozen_links if x['caseid']==cid]
    assert sorted(json.dumps(x,sort_keys=True) for x in fl)==sorted(json.dumps(x,sort_keys=True) for x in a['links'])
    for l in a['links']:
        assert l['primaryid']==c['primaryid'] and l['caseid']==cid and l['drug']==drug and l['source_quarter']==q
        assert l['role_cod']=='PS' and l['indi_pt']=='Dermatitis atopic'
    seq={x['drug_seq'] for x in a['links']}
    for t in a['ther']:
        assert t['primaryid']==c['primaryid'] and t['caseid']==cid and t.get('dsg_drug_seq',t.get('drug_seq')) in seq
    starts=[t['start_dt'] for t in a['ther']];ev=demo['event_dt']
    assert all(x==x.strip() for x in starts+[ev])
    assert starts==json.loads(old['start_dt_raw_values']) and ev==old['event_dt_raw']
    distinct=sorted(set(x for x in starts if x));es,ed=parse(ev)
    ss,sd=parse(distinct[0]) if len(distinct)==1 else ('missing' if not distinct else 'ambiguous',None)
    if es=='missing':status='01 Missing EVENT_DT'
    elif not distinct:status='02 Missing START_DT'
    elif es=='incomplete':status='03 Incomplete EVENT_DT'
    elif es=='invalid':status='04 Invalid EVENT_DT'
    elif len(distinct)>1:status='07 Ambiguous distinct nonblank START_DT'
    elif ss=='incomplete':status='05 Unique START_DT incomplete'
    elif ss=='invalid':status='06 Unique START_DT invalid'
    elif (ed-sd).days<0:status='08 Negative TTO'
    else:status='09 Valid TTO'
    days=(ed-sd).days if ed and sd else None
    assert status==old['exclusion_status'] and int(status=='09 Valid TTO')==int(old['valid'])
    equal(days,old['tto_days']);assert (sd.isoformat() if sd else '')==old['selected_start_date']
    assert (ed.isoformat() if ed else '')==old['event_date']
    statuses[(drug,status)]+=1
    verified.append(dict(drug=drug,caseid=cid,status=status,tto_days=days))
    if status=='09 Valid TTO':
        assert days>=0 and len(distinct)==1
        archive=next(x for x in meta['source_archives'] if x['quarter']==q)
        dm=next(x for x in meta['source_members'] if x['quarter']==q and x['table']=='DEMO')
        tm=next(x for x in meta['source_members'] if x['quarter']==q and x['table']=='THER')
        valid.append(dict(drug=drug,CASEID=cid,PRIMARYID=c['primaryid'],CASEVERSION=c['caseversion'],source_quarter=q,DRUG_SEQ=';'.join(sorted(seq)),START_DT=distinct[0],EVENT_DT=ev,TTO_days=days,source_archive=archive['filename'],archive_sha256=archive['sha256'],DEMO_member=dm['member'],DEMO_line=demo['source_line'],THER_member=tm['member'],THER_lines=';'.join(str(t['source_line']) for t in a['ther']),evidence_file='outputs/primary11_tto_v1_0/retained_date_evidence.json'))
db.close()
assert len(valid)==len({x['CASEID'] for x in valid})==119
summary=[];distribution=[];flow=[]
steps=['01 Missing EVENT_DT','02 Missing START_DT','03 Incomplete EVENT_DT','04 Invalid EVENT_DT','05 Unique START_DT incomplete','06 Unique START_DT invalid','07 Ambiguous distinct nonblank START_DT','08 Negative TTO']
for drug in DRUGS:
    vs=[r['TTO_days'] for r in valid if r['drug']==drug];n=len(vs);N=sum(x['drug']==drug for x in verified)
    qs=np.quantile(vs,[.25,.5,.75],method='linear').tolist() if n else [None]*3
    sr=dict(drug=drug,Primary11_total_n=N,valid_TTO_n=n,completeness=n/N,median_days=qs[1],Q1_days=qs[0],Q3_days=qs[2],IQR_width_days=qs[2]-qs[0] if n else None,min_days=min(vs) if n else None,max_days=max(vs) if n else None)
    if n:assert sr['min_days']<=qs[0]<=qs[1]<=qs[2]<=sr['max_days']
    summary.append(sr)
    masks=[lambda v:0<=v<=30,lambda v:31<=v<=90,lambda v:91<=v<=180,lambda v:v>180]
    for label,mask in zip(['0–30 d','31–90 d','91–180 d','＞180 d'],masks):
        k=sum(map(mask,vs));distribution.append(dict(drug=drug,TTO_bin=label,n=k,denominator=n,fraction=k/n if n else None))
    assert sum(r['n'] for r in distribution if r['drug']==drug)==n
    remain=N;flow.append(dict(drug=drug,step='00 Primary11 positive',excluded_n=0,remaining_n=remain))
    for s in steps:
        k=statuses[(drug,s)];remain-=k;flow.append(dict(drug=drug,step=s,excluded_n=k,remaining_n=remain))
    assert remain==n
    flow.append(dict(drug=drug,step='09 Valid TTO',excluded_n=0,remaining_n=n))
for rows,name,key in [(summary,'TTO_summary.csv','summary'),(distribution,'TTO_distribution.csv','distribution'),(flow,'exclusion_flow.csv','flow')]:
    prior=read(OLD/name)
    assert len(rows)==len(prior)
    for a,b,j in zip(rows,prior,oldresults[key]):
        for k,v in a.items():equal(v,b[k]);equal(v,j[k])
write(OUT/'TTO_summary_v1.1.csv',summary)
write(OUT/'valid_TTO_verified_v1.1.csv',valid)
write(WORK/'totals.csv',[dict(drug=r['drug'],Primary11_total_n=r['Primary11_total_n']) for r in summary])
result=dict(summary=summary,distribution=distribution,flow=flow,valid=valid,checks=dict(positive_cases=1108,valid_cases=119,unique_valid_CASEID=119,case_classification_differences=0,date_link_differences=0,summary_differences=0,distribution_differences=0,exclusion_flow_differences=0),versions=dict(Python=sys.version.split()[0],NumPy=np.__version__),protected_files=snapshot)
(WORK/'verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
print(json.dumps(dict(summary=summary,checks=result['checks']),ensure_ascii=False,indent=2))
