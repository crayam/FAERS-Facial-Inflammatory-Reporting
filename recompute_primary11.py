"""Recompute the frozen Primary11 union from the 17 user-supplied FAERS ZIPs.

Python standard library except python-docx, used to extract the frozen dictionary.
No original input is modified. Run with --source-dir and --dictionary if relocated.
"""
from pathlib import Path
import argparse, collections, csv, gzip, hashlib, json, re, sqlite3, time, zipfile
from docx import Document

P = argparse.ArgumentParser()
P.add_argument('--source-dir', default='inputs/faers_ascii')
P.add_argument('--dictionary', default='inputs/FAERS_Phenotype_Dictionary_v1.0.docx')
P.add_argument('--output-dir', default='outputs/primary11_v1_0')
P.add_argument('--work-dir', default='work/primary11')
A = P.parse_args()
SRC, OUT, WORK = Path(A.source_dir), Path(A.output_dir), Path(A.work_dir)
OUT.mkdir(parents=True, exist_ok=True); WORK.mkdir(parents=True, exist_ok=True)
DRUGS = ['dupilumab', 'tralokinumab', 'upadacitinib']
ROLES = ['index drug', 'primary active comparator', 'secondary active comparator']
ALIASES = {'dupilumab':['dupilumab','dupixent'], 'tralokinumab':['tralokinumab','adbry','adtralza'], 'upadacitinib':['upadacitinib','rinvoq']}
PAT = {d: re.compile(r'\b(?:'+'|'.join(a)+r')\b', re.I) for d,a in ALIASES.items()}
PERIODS = [f'{y}Q{q}' for y in range(2022,2027) for q in range(1,5) if (y,q)>=(2022,2) and (y,q)<=(2026,2)]
FILES = {}
for f in SRC.glob('*.zip'):
    m = re.fullmatch(r'faers_ascii_(20\d{2})q([1-4])\.zip',f.name,re.I)
    if m and (q:=f'{m[1]}Q{m[2]}') in PERIODS:
        assert q not in FILES, ('duplicate quarter',q)
        FILES[q]=f
assert sorted(FILES)==PERIODS, ('missing quarters',set(PERIODS)-set(FILES))

def log(*args): print(time.strftime('%H:%M:%S'), *args, flush=True)
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for block in iter(lambda:f.read(2**20),b''): h.update(block)
    return h.hexdigest()
def dump(p,obj): p.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')
def write_csv(name,rows,fields=None):
    rows=list(rows)
    fields=fields or list(rows[0])
    p=OUT/name
    op=gzip.open if name.endswith('.gz') else open
    with op(p,'wt',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fields);w.writeheader();w.writerows(rows)
    return p
def member(z,prefix):
    a=[n for n in z.namelist() if re.search(r'/(?:'+prefix+r')\d{2}Q[1-4]\.txt$',n,re.I)]
    assert len(a)==1,(prefix,a)
    return a[0]
STATS={}
def read_table(z,prefix,q):
    name=member(z,prefix); num=0
    with z.open(name) as f:
        hdr=f.readline().decode('utf-8-sig').rstrip('\r\n').lower().split('$')
        for num,b in enumerate(f,1):
            s=b.decode('utf-8').rstrip('\r\n').split('$')
            if len(s)!=len(hdr): raise ValueError((q,name,num+1,'field count',len(s),len(hdr)))
            yield num+1,s
    STATS[f'{q}/{prefix}']={'rows':num,'header':hdr,'member':name,'crc32':f'{z.getinfo(name).CRC:08x}','uncompressed_bytes':z.getinfo(name).file_size}
def read_deletes(z,q):
    name=member(z,'DELETE'); vals=set(); blanks=0
    with z.open(name) as f:
        for line in f:
            s=line.decode('utf-8-sig').strip()
            if not s: blanks+=1;continue
            if not s.isdigit(): raise ValueError((q,'DELETE nonnumeric',s))
            vals.add(int(s))
    STATS[f'{q}/DELETE']={'unique_caseids':len(vals),'blank_lines':blanks,'member':name,'crc32':f'{z.getinfo(name).CRC:08x}'}
    return vals

doc=Document(A.dictionary)
primary_table=next(t for t in doc.tables if [c.text.strip() for c in t.rows[0].cells]==['分区','原始 REAC.PT','纳入理由'])
primary=[{'region':r.cells[0].text.strip(),'pt':r.cells[1].text.strip(),'rationale':r.cells[2].text.strip()} for r in primary_table.rows[1:]]
PTS=[x['pt'] for x in primary]
assert len(PTS)==len(set(PTS))==11
assert set(PTS)=={'Erythema of eyelid','Eczema eyelids','Eyelid rash','Periorbital dermatitis','Periorbital eczema','Periorbital inflammation','Perioral dermatitis','Eczematous cheilitis','Cheilitis','Angular cheilitis','Lip erythema'}
PTSET=set(PTS)
dump(OUT/'primary11_definition.json',{'version':'1.0','frozen_date':'2026-09-20','dictionary_sha256':sha(A.dictionary),'primary11':primary,'drug_aliases':ALIASES,'strict_ad_indi_pt':'Dermatitis atopic','periods':PERIODS})

db=sqlite3.connect(WORK/'latest.sqlite')
db.execute('PRAGMA journal_mode=OFF');db.execute('PRAGMA synchronous=OFF');db.execute('PRAGMA temp_store=MEMORY')
db.execute('DROP TABLE IF EXISTS latest')
db.execute('CREATE TABLE latest (cid INTEGER PRIMARY KEY, pid INTEGER, ver INTEGER, fda INTEGER, q TEXT, init TEXT)')
sql='''INSERT INTO latest VALUES(?,?,?,?,?,?) ON CONFLICT(cid) DO UPDATE SET pid=excluded.pid,ver=excluded.ver,fda=excluded.fda,q=excluded.q,init=excluded.init WHERE (excluded.ver,excluded.fda,excluded.pid,excluded.q)>(latest.ver,latest.fda,latest.pid,latest.q)'''
deleted=set();manifest=[];quarter_rows={}
for q in PERIODS:
    f=FILES[q];manifest.append({'quarter':q,'filename':f.name,'bytes':f.stat().st_size,'sha256':sha(f)})
    with zipfile.ZipFile(f) as z:
        batch=[];seen_pid=set()
        for line,r in read_table(z,'DEMO',q):
            pid,cid,ver=int(r[0]),int(r[1]),int(r[2]); fda=int(r[7])
            assert pid not in seen_pid, ('duplicate DEMO PRIMARYID',q,pid)
            seen_pid.add(pid)
            batch.append((cid,pid,ver,fda,q,r[6]))
            if len(batch)==20000: db.executemany(sql,batch);batch=[]
        if batch: db.executemany(sql,batch)
        quarter_rows[q]=len(seen_pid)
        deleted.update(read_deletes(z,q))
    db.commit();log('DEMO',q,quarter_rows[q])
db.execute('CREATE INDEX latest_q ON latest(q)')
raw_unique=db.execute('SELECT COUNT(*) FROM latest').fetchone()[0]
db.execute('CREATE TEMP TABLE deleted (cid INTEGER PRIMARY KEY)')
db.executemany('INSERT INTO deleted VALUES(?)',((x,) for x in deleted))
deleted_intersection=db.execute('SELECT COUNT(*) FROM latest JOIN deleted USING(cid)').fetchone()[0]
db.execute('DELETE FROM latest WHERE cid IN (SELECT cid FROM deleted)');db.commit()
retained_count=db.execute('SELECT COUNT(*) FROM latest').fetchone()[0]
outside_dates=db.execute('SELECT COUNT(*) FROM latest WHERE fda<20220401 OR fda>20260630').fetchone()[0]
dump(WORK/'stage1.json',{'stats':STATS,'manifest':manifest,'raw_unique':raw_unique,'deleted_union':len(deleted),'deleted_intersection':deleted_intersection,'retained_count':retained_count,'outside_dates':outside_dates})
log('Global dedup complete',retained_count,'cases; removed',deleted_intersection)

cohort={};casepts={};case_meta={};eligible_links=[]
drug_name_audit=collections.Counter();indi_audit=collections.Counter();flow={d:collections.Counter() for d in DRUGS}
ambiguous=[];linked_mismatch=0;empty_reac=[];seq_collisions=[]
for q in PERIODS:
    retained={str(r[1]):r for r in db.execute('SELECT cid,pid,ver,fda,q,init FROM latest WHERE q=? AND fda BETWEEN 20220401 AND 20260630',(q,))}
    with zipfile.ZipFile(FILES[q]) as z:
        demomap={r[0]:r[1] for line,r in read_table(z,'DEMO',q)}
        selected={};all_ps={d:set() for d in DRUGS};any_role={d:set() for d in DRUGS}
        for line,r in read_table(z,'DRUG',q):
            pid,cid,seq,role,name,ai=r[:6]
            assert demomap.get(pid)==cid, ('DRUG key mismatch',q,line,pid,cid)
            if pid not in retained:continue
            ds=[d for d in DRUGS if PAT[d].search(name) or PAT[d].search(ai)]
            if len(ds)>1:
                ambiguous.append({'quarter':q,'primaryid':pid,'caseid':cid,'role':role,'drugname':name,'prod_ai':ai})
                assert role!='PS',('ambiguous PS drug',q,line)
                continue
            if not ds:continue
            d=ds[0];any_role[d].add(cid)
            drug_name_audit[(d,name,ai,role)]+=1
            if role!='PS':continue
            all_ps[d].add(cid)
            selected.setdefault((pid,seq),[]).append({'drug':d,'caseid':cid,'primaryid':pid,'drug_seq':seq,'drugname':name,'prod_ai':ai,'role_cod':role,'source_quarter':q,'drug_source_line':line})
        good={};any_ad_reports=set()
        for line,r in read_table(z,'INDI',q):
            pid,cid,seq,pt=r
            assert demomap.get(pid)==cid, ('INDI key mismatch',q,line,pid,cid)
            if pid in retained and pt=='Dermatitis atopic': any_ad_reports.add(pid)
            if (pid,seq) not in selected:continue
            for drugrow in selected[(pid,seq)]:
                indi_audit[(drugrow['drug'],pt)]+=1
            if pt!='Dermatitis atopic':continue
            for drugrow in selected[(pid,seq)]:
                d=drugrow['drug'];key=(d,cid)
                good[key]=pid
                eligible_links.append({**drugrow,'indi_pt':pt,'indi_source_line':line})
        for d in DRUGS:
            flow[d]['any_role_cases']+=len(any_role[d]);flow[d]['ps_cases']+=len(all_ps[d])
            flow[d]['ps_case_level_ad_any_drug']+=len({r['caseid'] for records in selected.values() for r in records if r['drug']==d and r['primaryid'] in any_ad_reports})
        eligible_pids=set(good.values());pts_by_pid={p:set() for p in eligible_pids};matched_rows=[]
        for line,r in read_table(z,'REAC',q):
            pid,cid,pt=r[:3]
            assert demomap.get(pid)==cid, ('REAC key mismatch',q,line,pid,cid)
            if pid not in eligible_pids: continue
            pts_by_pid[pid].add(pt)
            if pt in PTSET: matched_rows.append({'source_quarter':q,'reac_source_line':line,'primaryid':pid,'caseid':cid,'pt':pt})
        for key,pid in good.items():
            assert key not in cohort, ('duplicate drug CASEID',key)
            cid,_,ver,fda,_,init=retained[pid]
            pts=pts_by_pid[pid]
            if not pts:empty_reac.append(key)
            hits=[pt for pt in PTS if pt in pts]
            cohort[key]={'drug':key[0],'comparison_role':ROLES[DRUGS.index(key[0])],'caseid':str(cid),'primaryid':pid,'caseversion':ver,'fda_dt':str(fda),'init_fda_dt':init,'source_quarter':q,'Primary11':int(bool(hits)),'matched_pt_count':len(hits),'matched_pts':'; '.join(hits),'all_reac_unique_pt_count':len(pts)}
            casepts[key]=pts
        with gzip.open(WORK/f'{q}_matched_reac_rows.json.gz','wt',encoding='utf-8') as f:json.dump(matched_rows,f)
        log('DRUG INDI REAC',q,'eligible',len(good),'Primary11',sum(bool(pts_by_pid[p]&PTSET) for p in good.values()))

counts=collections.Counter(c for d,c in cohort)
for (d,c),row in cohort.items(): row['multi_drug_cohort']=int(counts[c]>1)
ordered=sorted(cohort.values(),key=lambda r:(DRUGS.index(r['drug']),int(r['caseid'])))
summary=[];pt_counts=[];patterns=collections.Counter();quarters=[]
for d,role in zip(DRUGS,ROLES):
    rows=[r for r in ordered if r['drug']==d];positive=[r for r in rows if r['Primary11']]
    n=len(rows);a=len(positive);s=sum(r['matched_pt_count'] for r in rows)
    summary.append({'drug':d,'comparison_role':role,**dict(flow[d]),'ps_strict_ad_cases':n,'Primary11_cases':a,'Primary11_negative_cases':n-a,'Primary11_proportion':a/n,'sum_individual_pt_cases':s,'duplicate_count_excess':s-a,'cases_with_multiple_Primary11_pts':sum(r['matched_pt_count']>1 for r in rows),'unique_reac_pts':len(set().union(*(casepts[(d,r['caseid'])] for r in rows))),'initial_report_before_window':sum(bool(r['init_fda_dt']) and int(r['init_fda_dt'])<20220401 for r in rows)})
    for r in positive:patterns[(d,r['matched_pts'],r['matched_pt_count'])]+=1
    for q in PERIODS:
        qr=[r for r in rows if r['source_quarter']==q]
        quarters.append({'drug':d,'source_quarter':q,'ps_strict_ad_cases':len(qr),'Primary11_cases':sum(r['Primary11'] for r in qr)})
for item in primary:
    pt_counts.append({**item,**{d:sum(item['pt'] in pts for (drug,cid),pts in casepts.items() if drug==d) for d in DRUGS}})
write_csv('union_summary.csv',summary)
write_csv('case_Primary11_flags.csv.gz',ordered)
write_csv('positive_case_evidence.csv',[r for r in ordered if r['Primary11']])
write_csv('primary11_pt_counts.csv',pt_counts)
write_csv('quarterly_counts.csv',quarters)
write_csv('matched_pt_patterns.csv',[{'drug':d,'matched_pts':pts,'matched_pt_count':k,'cases':n} for (d,pts,k),n in sorted(patterns.items())])
write_csv('eligible_drug_indication_links.csv.gz',eligible_links)
write_csv('cohort_case_PT.csv.gz',({'drug':d,'caseid':c,'primaryid':cohort[(d,c)]['primaryid'],'source_quarter':cohort[(d,c)]['source_quarter'],'pt':pt} for (d,c),pts in sorted(casepts.items()) for pt in sorted(pts)))
write_csv('drug_name_audit.csv',[{'drug':d,'drugname':name,'prod_ai':ai,'role_cod':role,'retained_drug_rows':n} for (d,name,ai,role),n in sorted(drug_name_audit.items())])
write_csv('index_drug_indication_audit.csv',[{'drug':d,'indi_pt':pt,'linked_ps_drug_indi_rows':n} for (d,pt),n in sorted(indi_audit.items())])
write_csv('input_manifest.csv',manifest)
matched_evidence=[]
for q in PERIODS:
    with gzip.open(WORK/f'{q}_matched_reac_rows.json.gz','rt',encoding='utf-8') as f:matched_evidence.extend(json.load(f))
write_csv('primary11_raw_reac_evidence.csv',matched_evidence)
qc={'source_periods':PERIODS,'source_stats':STATS,'DEMO_rows':sum(quarter_rows.values()),'unique_CASEID_before_delete':raw_unique,'DELETE_CASEID_union':len(deleted),'deleted_CASEID_in_DEMO':deleted_intersection,'retained_CASEID_after_delete':retained_count,'retained_cases_outside_date_window':outside_dates,'cohort_drug_case_rows':len(cohort),'cohort_unique_CASEID':len(counts),'cross_drug_overlap_cases':sum(v>1 for v in counts.values()),'empty_REAC_cohort_cases':empty_reac,'ambiguous_non_PS_drug_rows':ambiguous,'decode_errors':0,'field_count_errors':0,'key_mismatches':0,'duplicate_drug_case_rows':len(ordered)-len(cohort),'cohort_case_PT_pairs':sum(len(v) for v in casepts.values()),'cohort_unique_PT':len(set().union(*casepts.values())),'PT_match_mode':'exact source text; no case folding, synonyms or MedDRA remapping','dedup_order':['CASEVERSION numeric descending','FDA_DT numeric descending','PRIMARYID numeric descending','source_quarter descending'],'input_sha256':{m['filename']:m['sha256'] for m in manifest},'dictionary_sha256':sha(A.dictionary),'drug_aliases':ALIASES,'strict_AD':'Drug-specific exact Dermatitis atopic linked via PRIMARYID + CASEID + DRUG_SEQ','date_filter':'Retained FDA_DT 20220401 through 20260630 inclusive','summary':summary}
assert not empty_reac
assert qc['duplicate_drug_case_rows']==0
dump(OUT/'analysis_manifest.json',qc)
dump(WORK/'workbook_data.json',{'summary':summary,'pt_counts':pt_counts,'quarters':quarters,'qc':qc,'patterns':[{'drug':d,'matched_pts':pts,'matched_pt_count':k,'cases':n} for (d,pts,k),n in sorted(patterns.items())]})
db.close()
log('COMPLETE',json.dumps(summary,ensure_ascii=False))
