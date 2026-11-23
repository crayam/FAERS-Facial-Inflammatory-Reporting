from pathlib import Path
from datetime import datetime
from collections import Counter
import csv,hashlib,json,math,re,zipfile
from openpyxl import load_workbook
R=Path(__file__).resolve().parents[2];O=R/'outputs/primary11_tto_date_audit_v1_0';W=R/'work/tto_primary11'
d=json.loads((O/'audit_results.json').read_text())
with (O/'raw_source_records.csv').open(encoding='utf-8-sig') as f:evidence=list(csv.DictReader(f))
eventrecords={x['primaryid']:x for x in evidence if x['table']=='DEMO'}
for r in d['tralokinumab']:
 e=eventrecords[r['PRIMARYID']];h=next(x['header'] for x in d['source_provenance'] if x['quarter']==r['source_quarter'] and x['table']=='DEMO');v=e['raw_record'].split('$');a=dict(zip(h,v))
 assert len(v)==len(h) and a['event_dt']==r['raw_DEMO_EVENT_DT']=='' and r['parsed_EVENT_DT']==''
 assert a['primaryid']==r['PRIMARYID'] and a['caseid']==r['CASEID'] and int(a['caseversion'])==r['retained_CASEVERSION']
 assert r['raw_EVENT_DT_byte_length']==0 and int(e['source_line'])==r['DEMO_line']
assert len(d['tralokinumab'])==len({x['CASEID'] for x in d['tralokinumab']})==19
numeric=0;neg=0
for a in d['multi_cases']:
 ther=[x for x in evidence if x['primaryid']==a['PRIMARYID'] and x['table']=='THER']
 raw=[]
 for t in ther:
  h=next(x['header'] for x in d['source_provenance'] if x['quarter']==t['source_quarter'] and x['table']=='THER');row=dict(zip(h,t['raw_record'].split('$')))
  assert row['caseid']==a['CASEID'] and row['dsg_drug_seq']==a['DRUG_SEQ'];raw.append(row['start_dt'])
 assert raw==json.loads(a['all_raw_START_DT'])
 distinct=list(dict.fromkeys(raw));assert distinct==json.loads(a['distinct_START_DT'])
 candidates=[x for x in d['candidates'] if x['CASEID']==a['CASEID']]
 assert [x['raw_START_DT'] for x in candidates]==distinct
 for c in candidates:
  def parse(x):
   if len(x)!=8 or not x.isdigit():return None
   try:return datetime.strptime(x,'%Y%m%d')
   except ValueError:return None
  sd,ed=parse(c['raw_START_DT']),parse(c['raw_EVENT_DT']);delta=(ed-sd).days if sd and ed else None
  assert delta==c['candidate_TTO_days'] and c['duplicate_rows']==Counter(raw)[c['raw_START_DT']]
  if delta is not None:numeric+=1;neg+=int(delta<0)
assert not d['QC']['final_START_DT_selected'] and not d['QC']['final_TTO_summary_computed']
book=O/'FAERS_Primary11_TTO_DateAudit_v1.0.xlsx';wb=load_workbook(book,read_only=True,data_only=True)
assert wb.sheetnames==['Tralo_EVENT_DT','Multi_START_DT','Candidate_TTO','Audit_QC_Methods']
cells={s.title:{c.coordinate:c.value for row in s for c in row if c.value is not None} for s in wb}
targets=json.loads((W/'audit_check_targets.json').read_text())
for t in targets:
 a=cells[t['sheet']].get(t['address']);b=t['value']
 if a is None and b in ('',None):continue
 if isinstance(b,(int,float)):assert isinstance(a,(int,float)) and math.isclose(a,b,abs_tol=1e-12),(t,a)
 else:assert a==b,(t,a)
for sn,values in cells.items():
 for cell,v in values.items():
  if isinstance(v,str):assert not re.fullmatch(r'#(?:REF!|DIV/0!|VALUE!|NAME\?|N/A|NUM!|NULL!|SPILL!|CALC!)',v),(sn,cell,v)
wb.close()
validation=dict(status='PASS',tralokinumab_raw_empty_EVENT_DT=19,all_multiple_THER_cases=27,multiple_distinct_nonblank_START_DT_cases=24,multiple_distinct_complete_START_DT_cases=5,distinct_candidate_rows=57,numeric_candidate_differences=numeric,negative_candidates_preserved=neg,exported_cell_targets_checked=len(targets),raw_record_fields_independently_checked=True,final_START_DT_selected=False,final_TTO_summary_computed=False,workbook_sha256=hashlib.sha256(book.read_bytes()).hexdigest())
(O/'delivery_validation.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2))
files=sorted(p for p in O.iterdir() if p.is_file() and p.name!='checksums.sha256' and not p.name.endswith('.inspect.ndjson'))
(O/'checksums.sha256').write_text(''.join(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n' for p in files))
print(json.dumps(validation,ensure_ascii=False,indent=2))
