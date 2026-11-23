from pathlib import Path
from collections import Counter,defaultdict
from datetime import date,datetime
import csv,gzip,hashlib,json,re,sqlite3,zipfile
R=Path(__file__).resolve().parents[2];BASE=R/'outputs/primary11_tto_v1_0';O=R/'outputs/primary11_tto_date_audit_v1_0';O.mkdir(parents=True,exist_ok=True)
SRC=Path('inputs/faers_ascii')
raw=json.loads((BASE/'retained_date_evidence.json').read_text());oldmeta=json.loads((BASE/'extraction_manifest.json').read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def read(p):
 with (gzip.open(p,'rt',encoding='utf-8-sig') if p.suffix=='.gz' else p.open(encoding='utf-8-sig')) as f:return list(csv.DictReader(f))
for s in oldmeta['sources']:assert sha(Path(s['path']))==s['sha256']
frozen=read(R/'outputs/all_outcomes_v1_0/case_outcome_flags.csv.gz');cohort={a['caseid']:a for a in frozen if a['Primary11']=='1'}
assert len(cohort)==len(raw)==1108
assert Counter(a['drug'] for a in cohort.values())=={'dupilumab':1068,'tralokinumab':19,'upadacitinib':21}
for a in raw:
 c=a['cohort'];assert all(cohort[c['caseid']][k]==v for k,v in c.items())
multi=[a for a in raw if len(a['ther'])>1]
tralo=[a for a in raw if a['cohort']['drug']=='tralokinumab'];assert len(tralo)==19
targets={a['cohort']['primaryid']:a for a in tralo+multi}
db=sqlite3.connect(f'file:{R}/work/primary11/latest.sqlite?mode=ro',uri=True)
for a in targets.values():
 c=a['cohort'];assert db.execute('SELECT pid,ver,fda,q FROM latest WHERE cid=?',(int(c['caseid']),)).fetchone()==(int(c['primaryid']),int(c['caseversion']),int(c['fda_dt']),c['source_quarter'])
db.close()
fresh={};provenance=[];archivechecks=[];evidence=[]
for src in oldmeta['source_archives']:
 q=src['quarter'];sel={p:a for p,a in targets.items() if a['cohort']['source_quarter']==q}
 if not sel:continue
 path=SRC/src['filename'];assert sha(path)==src['sha256'];archivechecks.append(src)
 pids={x.encode() for x in sel};demo={};ther=defaultdict(list)
 with zipfile.ZipFile(path) as z:
  for table in ['DEMO','THER']:
   names=[n for n in z.namelist() if re.search(r'/'+table+r'\d{2}Q[1-4]\.txt$',n,re.I)];assert len(names)==1
   name=names[0];info=z.getinfo(name)
   with z.open(name) as f:
    header=f.readline().decode('utf-8-sig').rstrip('\r\n').lower().split('$')
    ix=header.index('event_dt' if table=='DEMO' else 'start_dt')
    provenance.append(dict(quarter=q,table=table,zip_filename=src['filename'],member=name,crc32=f'{info.CRC:08x}',header=header,date_column_1based=ix+1))
    for line,b in enumerate(f,2):
     if b.split(b'$',1)[0] not in pids:continue
     byte_fields=b.rstrip(b'\r\n').split(b'$');assert len(byte_fields)==len(header)
     vals=[x.decode('utf-8') for x in byte_fields];a=dict(zip(header,vals));pid=a['primaryid'];c=sel[pid]['cohort'];assert a['caseid']==c['caseid']
     if table=='THER' and a['dsg_drug_seq'] not in {x['drug_seq'] for x in sel[pid]['links']}:continue
     a['source_line']=line
     if table=='DEMO':
      assert pid not in demo;assert all(a[k]==c[k] for k in ['caseversion','fda_dt','init_fda_dt'])
      assert a==sel[pid]['demo'];demo[pid]=a
     else:ther[pid].append(a)
     evidence.append(dict(drug=c['drug'],caseid=c['caseid'],primaryid=pid,retained_caseversion=c['caseversion'],source_quarter=q,table=table,zip_filename=src['filename'],member=name,source_line=line,raw_date_field=vals[ix],raw_date_json=json.dumps(vals[ix]),raw_date_byte_length=len(byte_fields[ix]),raw_record_sha256=hashlib.sha256(b).hexdigest(),raw_record=b.decode('utf-8').rstrip('\r\n')))
  assert set(demo)==set(sel)
  for pid,a in sel.items():
   assert ther[pid]==a['ther'],('THER differs',pid)
   fresh[pid]=dict(cohort=a['cohort'],links=a['links'],demo=demo[pid],ther=ther[pid])
 print(q,'raw DEMO/THER verified',len(sel),flush=True)
def parse(s):
 if not s.strip():return None,'Missing'
 if not re.fullmatch('[0-9]{8}',s):return None,'Incomplete' if re.fullmatch('[0-9]{4}|[0-9]{6}',s) else 'Invalid format'
 try:
  x=date(int(s[:4]),int(s[4:6]),int(s[6:8]));assert x==datetime.strptime(s,'%Y%m%d').date();return x,'Complete'
 except ValueError:return None,'Invalid calendar date'
def member(q,t):return next(x['member'] for x in provenance if x['quarter']==q and x['table']==t)
eventrows=[]
for old in tralo:
 a=fresh[old['cohort']['primaryid']];c=a['cohort'];e=a['demo']['event_dt'];parsed,status=parse(e)
 ev=next(x for x in evidence if x['primaryid']==c['primaryid'] and x['table']=='DEMO')
 eventrows.append(dict(CASEID=c['caseid'],PRIMARYID=c['primaryid'],retained_CASEVERSION=int(c['caseversion']),source_quarter=c['source_quarter'],raw_DEMO_EVENT_DT=e,parsed_EVENT_DT=parsed.isoformat() if parsed else '',raw_EVENT_DT_json=json.dumps(e),raw_EVENT_DT_byte_length=ev['raw_date_byte_length'],parse_status=status,retained_version_QC='PASS',quarter_key_QC='PASS',raw_reread_QC='PASS',DEMO_member=ev['member'],DEMO_line=ev['source_line']))
multirows=[];candidates=[]
for old in multi:
 a=fresh[old['cohort']['primaryid']];c=a['cohort'];seqs=sorted({x['drug_seq'] for x in a['links']});assert len(seqs)==1
 values=[x['start_dt'] for x in a['ther']];distinct=list(dict.fromkeys(values));nonblank={x for x in distinct if x.strip()};dates=[parse(x)[0] for x in distinct];full={x for x in dates if x is not None}
 e=a['demo']['event_dt'];ed,es=parse(e)
 kind='Multiple distinct complete dates' if len(full)>1 else 'Multiple distinct raw dates (includes incomplete)' if len(nonblank)>1 else 'Repeated same date / blank records'
 multirows.append(dict(drug=c['drug'],CASEID=c['caseid'],PRIMARYID=c['primaryid'],DRUG_SEQ=seqs[0],retained_CASEVERSION=int(c['caseversion']),source_quarter=c['source_quarter'],all_raw_START_DT=json.dumps(values),distinct_START_DT=json.dumps(distinct),raw_EVENT_DT=e,raw_EVENT_DT_json=json.dumps(e),THER_rows=len(values),distinct_nonblank_dates=len(nonblank),distinct_complete_dates=len(full),audit_class=kind,THER_member=member(c['source_quarter'],'THER'),THER_lines=';'.join(str(t['source_line']) for t in a['ther'])))
 for start in distinct:
  sd,ss=parse(start);matching=[t for t in a['ther'] if t['start_dt']==start];delta=(ed-sd).days if ed and sd else None
  reason='Candidate only; negative retained' if delta is not None and delta<0 else 'Candidate only; no selection' if delta is not None else '; '.join(([f'EVENT_DT {es}'] if ed is None else [])+([f'START_DT {ss}'] if sd is None else []))
  candidates.append(dict(drug=c['drug'],CASEID=c['caseid'],PRIMARYID=c['primaryid'],DRUG_SEQ=seqs[0],raw_START_DT=start,raw_START_DT_json=json.dumps(start),START_parse_status=ss,raw_EVENT_DT=e,EVENT_parse_status=es,candidate_TTO_days=delta,status=reason,duplicate_rows=len(matching),source_quarter=c['source_quarter'],THER_lines=';'.join(str(x['source_line']) for x in matching)))
assert len(eventrows)==19 and all(x['raw_EVENT_DT_byte_length']==0 and x['raw_DEMO_EVENT_DT']=='' and x['parse_status']=='Missing' for x in eventrows)
assert len(multirows)==27 and sum(x['distinct_nonblank_dates']>1 for x in multirows)==24 and sum(x['distinct_complete_dates']>1 for x in multirows)==5
def write(name,rows):
 with (O/name).open('w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
for name,rows in [('tralokinumab_EVENT_DT_audit.csv',eventrows),('multiple_START_DT_case_audit.csv',multirows),('candidate_TTO_audit.csv',candidates),('raw_source_records.csv',evidence)]:write(name,rows)
result=dict(tralokinumab=eventrows,multi_cases=multirows,candidates=candidates,QC=dict(tralokinumab_cases=19,raw_EVENT_DT_empty=19,version_mismatch=0,quarter_or_key_mismatch=0,parse_failure_nonempty_EVENT_DT=0,multiple_THER_cases=27,multiple_distinct_nonblank_START_cases=24,multiple_distinct_complete_START_cases=5,candidate_rows=len(candidates),final_START_DT_selected=False,final_TTO_summary_computed=False),source_provenance=provenance,archives_reverified=archivechecks,frozen_sources=oldmeta['sources'])
(O/'audit_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2));print(json.dumps(result['QC'],ensure_ascii=False,indent=2))
