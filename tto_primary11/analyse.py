from pathlib import Path
from collections import Counter
from datetime import date
import argparse,csv,json,re,hashlib
import numpy as np

R=Path(__file__).resolve().parents[2];O=R/'outputs/primary11_tto_v1_0'
p=argparse.ArgumentParser();p.add_argument('--multiple-start-rule',required=True,choices=['strict_distinct_nonblank']);args=p.parse_args()
raw=json.loads((O/'retained_date_evidence.json').read_text());meta=json.loads((O/'extraction_manifest.json').read_text())
for x in meta['sources']:assert hashlib.sha256(Path(x['path']).read_bytes()).hexdigest()==x['sha256']
auditdir=R/'outputs/primary11_tto_date_audit_v1_0'
checks={x.split('  ',1)[1]:x.split('  ',1)[0] for x in (auditdir/'checksums.sha256').read_text().splitlines()}
for name in ['audit_results.json','delivery_validation.json','raw_source_records.csv']:
 assert hashlib.sha256((auditdir/name).read_bytes()).hexdigest()==checks[name]
prior=json.loads((auditdir/'audit_results.json').read_text())
for x in prior['multi_cases']:
 a=next(a for a in raw if a['cohort']['caseid']==x['CASEID'])
 assert json.loads(x['all_raw_START_DT'])==[t['start_dt'] for t in a['ther']]
for x in prior['tralokinumab']:
 a=next(a for a in raw if a['cohort']['caseid']==x['CASEID'])
 assert x['raw_DEMO_EVENT_DT']==a['demo']['event_dt']==''
drugs=['dupilumab','tralokinumab','upadacitinib']
def parse(s):
 s=s.strip()
 if not s:return 'missing',None
 if re.fullmatch(r'\d{4}|\d{6}',s):return 'incomplete',None
 if re.fullmatch(r'\d{8}',s):
  if s[4:6]=='00' or s[6:8]=='00':return 'incomplete',None
  try:return 'complete',date(int(s[:4]),int(s[4:6]),int(s[6:8]))
  except ValueError:return 'invalid',None
 return 'invalid',None
cases=[]
for a in raw:
 c=a['cohort'];event=a['demo']['event_dt'];es,ed=parse(event)
 vals=[t['start_dt'] for t in a['ther']];assert all(x==x.strip() for x in vals+[event])
 nonblank=sorted(set(x for x in vals if x));parts=[parse(x) for x in nonblank]
 complete=sorted({v for s,v in parts if s=='complete'})
 missing_start=not nonblank
 incomplete_start=any(s=='incomplete' for s,v in parts)
 invalid_start=any(s=='invalid' for s,v in parts)
 multiple=len(complete)>1;ambiguous=len(nonblank)>=2
 sd=complete[0] if len(nonblank)==1 and len(complete)==1 else None
 tto=(ed-sd).days if ed and sd else None
 if es=='missing':status='01 Missing EVENT_DT'
 elif missing_start:status='02 Missing START_DT'
 elif es=='incomplete':status='03 Incomplete EVENT_DT'
 elif es=='invalid':status='04 Invalid EVENT_DT'
 elif ambiguous:status='07 Ambiguous distinct nonblank START_DT'
 elif not complete:
  status='05 Unique START_DT incomplete' if incomplete_start else '06 Unique START_DT invalid'
 elif tto<0:status='08 Negative TTO'
 else:status='09 Valid TTO'
 cases.append(dict(**c,event_dt_raw=event,event_status=es,qualifying_drug_seq=';'.join(sorted({x['drug_seq'] for x in a['links']})),ther_rows=len(vals),start_dt_raw_values=json.dumps(vals),distinct_nonblank_START_DT=json.dumps(nonblank),distinct_nonblank_n=len(nonblank),exposure_start_ambiguous=int(ambiguous),complete_start_dates=';'.join(x.isoformat() for x in complete),selected_start_date=sd.isoformat() if sd else '',event_date=ed.isoformat() if ed else '',missing_start=int(missing_start),missing_event=int(es=='missing'),any_incomplete_start=int(incomplete_start),any_incomplete_event=int(es=='incomplete'),any_invalid_start=int(invalid_start),any_invalid_event=int(es=='invalid'),multiple_complete_starts=int(multiple),tto_days=tto,valid=int(status=='09 Valid TTO'),exclusion_status=status,demo_source_line=a['demo']['source_line'],ther_source_lines=';'.join(str(x['source_line']) for x in a['ther'])))
summary=[];distribution=[];flow=[];qc=[]
statuses=['01 Missing EVENT_DT','02 Missing START_DT','03 Incomplete EVENT_DT','04 Invalid EVENT_DT','05 Unique START_DT incomplete','06 Unique START_DT invalid','07 Ambiguous distinct nonblank START_DT','08 Negative TTO']
for drug in drugs:
 rows=[a for a in cases if a['drug']==drug];valid=[a['tto_days'] for a in rows if a['valid']];N=len(rows);n=len(valid)
 assert len({a['caseid'] for a in rows})==N
 q=np.quantile(valid,[.25,.5,.75],method='linear') if n else [None,None,None]
 summary.append(dict(drug=drug,Primary11_total_n=N,valid_TTO_n=n,completeness=n/N,median_days=q[1],Q1_days=q[0],Q3_days=q[2],IQR_width_days=q[2]-q[0] if n else None,min_days=min(valid) if n else None,max_days=max(valid) if n else None))
 for label,fn in [('0–30 d',lambda v:0<=v<=30),('31–90 d',lambda v:31<=v<=90),('91–180 d',lambda v:91<=v<=180),('＞180 d',lambda v:v>180)]:
  k=sum(fn(v) for v in valid);distribution.append(dict(drug=drug,TTO_bin=label,n=k,denominator=n,fraction=k/n if n else None))
 remain=N;flow.append(dict(drug=drug,step='00 Primary11 positive',excluded_n=0,remaining_n=remain))
 for status in statuses:
  k=sum(x['exclusion_status']==status for x in rows);remain-=k
  flow.append(dict(drug=drug,step=status,excluded_n=k,remaining_n=remain))
 assert remain==n;flow.append(dict(drug=drug,step='09 Valid TTO',excluded_n=0,remaining_n=n))
 qc.append(dict(drug=drug,Primary11_total=N,missing_START_DT=sum(x['missing_start'] for x in rows),no_matched_THER=sum(x['ther_rows']==0 for x in rows),matched_THER_all_START_blank=sum(x['ther_rows']>0 and x['missing_start'] for x in rows),missing_EVENT_DT=sum(x['missing_event'] for x in rows),both_dates_missing=sum(x['missing_event'] and x['missing_start'] for x in rows),any_incomplete_START_DT=sum(x['any_incomplete_start'] for x in rows),any_incomplete_EVENT_DT=sum(x['any_incomplete_event'] for x in rows),any_incomplete_date=sum(x['any_incomplete_start'] or x['any_incomplete_event'] for x in rows),any_invalid_date=sum(x['any_invalid_start'] or x['any_invalid_event'] for x in rows),multiple_complete_START_DT=sum(x['multiple_complete_starts'] for x in rows),multiple_THER_rows=sum(x['ther_rows']>1 for x in rows),complete_date_pair=sum(x['tto_days'] is not None for x in rows),negative_TTO_excluded=sum(x['exclusion_status']=='08 Negative TTO' for x in rows),same_day_TTO=sum(x['valid'] and x['tto_days']==0 for x in rows),valid_TTO=n))
 qc[-1].update(ambiguous_distinct_nonblank_START_DT=sum(x['exposure_start_ambiguous'] for x in rows),exactly_180_days=sum(v==180 for v in valid),unique_incomplete_START_DT=sum(x['distinct_nonblank_n']==1 and x['any_incomplete_start'] for x in rows))
def write(name,rows):
 with (O/name).open('w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
for name,rows in [('case_TTO_audit.csv',cases),('TTO_summary.csv',summary),('TTO_distribution.csv',distribution),('exclusion_flow.csv',flow),('date_QC.csv',qc)]:write(name,rows)
assert Counter(x['drug'] for x in cases)=={'dupilumab':1068,'tralokinumab':19,'upadacitinib':21}
assert len(cases)==len({x['caseid'] for x in cases})==1108
assert sum(x['exposure_start_ambiguous'] for x in cases)==24
data=dict(summary=summary,distribution=distribution,flow=flow,qc=qc,multiple_start_rule=args.multiple_start_rule,meta=meta,statistics='Exploratory descriptive only; no between-drug comparison, significance test, Weibull or Kaplan–Meier analysis',quantile='Linear interpolation, type 7 / NumPy linear',case_count=len(cases),prior_date_audit_sha256=checks['audit_results.json'],boundary_note='Mutually exclusive bins: 0–30 d, 31–90 d, 91–180 d, >180 d')
(O/'analysis_results.json').write_text(json.dumps(data,ensure_ascii=False,indent=2))
print(json.dumps(dict(summary=summary,qc=qc),ensure_ascii=False,indent=2))
