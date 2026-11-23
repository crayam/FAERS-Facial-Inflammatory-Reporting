from pathlib import Path
import csv,gzip,hashlib,json,re,sqlite3,zipfile
from collections import Counter,defaultdict
R=Path(__file__).resolve().parents[2];O=R/'outputs/primary11_tto_v1_0';O.mkdir(parents=True,exist_ok=True)
SRC=Path('inputs/faers_ascii')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def read(p):
 with (gzip.open(p,'rt',encoding='utf-8-sig') if p.suffix=='.gz' else p.open(encoding='utf-8-sig')) as f:return list(csv.DictReader(f))
sources=[]
def check(folder,name):
 p=R/'outputs'/folder/name
 hashes={x.split('  ',1)[1]:x.split('  ',1)[0] for x in (p.parent/'checksums.sha256').read_text().splitlines()}
 assert sha(p)==hashes[name],name
 sources.append(dict(path=str(p),sha256=hashes[name]));return p
cohort=read(check('all_outcomes_v1_0','case_outcome_flags.csv.gz'))
assert Counter(a['drug'] for a in cohort)=={'dupilumab':164515,'tralokinumab':3637,'upadacitinib':6254}
assert len(cohort)==len({a['caseid'] for a in cohort})==174406
pos=[a for a in cohort if a['Primary11']=='1'];assert Counter(a['drug'] for a in pos)=={'dupilumab':1068,'tralokinumab':19,'upadacitinib':21}
byid={a['primaryid']:a for a in pos};assert len(byid)==1108
links=read(check('primary11_v1_0','eligible_drug_indication_links.csv.gz'))
selected=defaultdict(list)
for a in links:
 if a['primaryid'] not in byid:continue
 b=byid[a['primaryid']]
 assert all(a[k]==b[k] for k in ['drug','caseid','source_quarter']) and a['role_cod']=='PS' and a['indi_pt']=='Dermatitis atopic'
 selected[a['primaryid']].append(a)
assert set(selected)==set(byid)
db=sqlite3.connect(f'file:{R}/work/primary11/latest.sqlite?mode=ro',uri=True)
for a in pos:assert db.execute('SELECT pid,ver,fda,q FROM latest WHERE cid=?',(int(a['caseid']),)).fetchone()==(int(a['primaryid']),int(a['caseversion']),int(a['fda_dt']),a['source_quarter'])
db.close()
manifest=read(check('all_outcomes_v1_0','input_manifest.csv'));raw=[];members=[];quarter_qc=[]
for src in manifest:
 q=src['quarter'];path=SRC/src['filename'];assert sha(path)==src['sha256'],q
 eligible={p:a for p,a in byid.items() if a['source_quarter']==q};pids={x.encode() for x in eligible};demos={};ther=defaultdict(list);qc=Counter()
 with zipfile.ZipFile(path) as z:
  for prefix in ['DEMO','THER']:
   names=[n for n in z.namelist() if re.search(r'/'+prefix+r'\d{2}Q[1-4]\.txt$',n,re.I)];assert len(names)==1,(q,prefix,names)
   name=names[0];info=z.getinfo(name)
   with z.open(name) as f:
    h=f.readline().decode('utf-8-sig').rstrip('\r\n').lower().split('$')
    members.append(dict(quarter=q,table=prefix,member=name,crc32=f'{info.CRC:08x}',header=h))
    for ln,b in enumerate(f,2):
     qc[prefix+'_scanned_rows']+=1
     if b.split(b'$',1)[0] not in pids:continue
     vals=b.decode('utf-8').rstrip('\r\n').split('$');assert len(vals)==len(h),(q,prefix,ln)
     a=dict(zip(h,vals));pid=a['primaryid'];c=eligible[pid];assert a['caseid']==c['caseid']
     a['source_line']=ln
     if prefix=='DEMO':
      assert pid not in demos
      assert all(a[k]==c[k] for k in ['caseversion','fda_dt','init_fda_dt'])
      demos[pid]=a
     else:
      seqkey='dsg_drug_seq' if 'dsg_drug_seq' in a else 'drug_seq'
      if a[seqkey] in {x['drug_seq'] for x in selected[pid]}:ther[pid].append(a)
   print(q,prefix,'selected',len(demos) if prefix=='DEMO' else sum(map(len,ther.values())),flush=True)
  assert set(demos)==set(eligible)
  for pid,c in eligible.items():
   raw.append(dict(cohort={k:c[k] for k in ['drug','caseid','primaryid','caseversion','source_quarter','fda_dt','init_fda_dt']},links=selected[pid],demo=demos[pid],ther=ther.get(pid,[])))
 quarter_qc.append(dict(quarter=q,positive_cases=len(eligible),matched_DEMO=len(demos),matched_THER_cases=sum(bool(records) for records in ther.values()),matched_THER_rows=sum(map(len,ther.values())),**qc))
raw.sort(key=lambda a:(['dupilumab','tralokinumab','upadacitinib'].index(a['cohort']['drug']),int(a['cohort']['caseid'])))
(O/'retained_date_evidence.json').write_text(json.dumps(raw,ensure_ascii=False,indent=2))
audit=[]
for a in raw:
 starts=sorted({x['start_dt'].strip() for x in a['ther'] if x['start_dt'].strip()})
 if len(starts)>1:audit.append(dict(**a['cohort'],event_dt=a['demo']['event_dt'],start_dates=starts,ther_rows=len(a['ther'])))
meta=dict(sources=sources,source_archives=manifest,source_members=members,quarter_qc=quarter_qc,formal_N=dict(Counter(a['drug'] for a in cohort)),Primary11_N=dict(Counter(a['drug'] for a in pos)),retained_versions_verified=True,multiple_nonblank_start_dates=audit,THER_header=next(x['header'] for x in members if x['table']=='THER'))
(O/'extraction_manifest.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2))
print(json.dumps(dict(complete=True,positive_cases=len(raw),multiple_start_dates=audit,ther_count_distribution=dict(Counter(len(a['ther']) for a in raw)),unique_seq_distribution=dict(Counter(len({x['drug_seq'] for x in a['links']}) for a in raw))),ensure_ascii=False),flush=True)
