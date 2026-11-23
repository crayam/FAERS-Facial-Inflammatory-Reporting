"""Extract only the frozen retained report versions; no cohort selection changes."""
from pathlib import Path
from collections import Counter, defaultdict
import csv, gzip, hashlib, json, re, sqlite3, time, zipfile

SRC=Path('inputs/faers_ascii')
FROZEN=Path('outputs/all_outcomes_v1_0')
WORK=Path('work/table1')
EXPECTED={'dupilumab':164515,'tralokinumab':3637,'upadacitinib':6254}
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()
def dump(name,obj):
    (WORK/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2))
def log(*args):print(time.strftime('%H:%M:%S'),*args,flush=True)

frozen_checks={line.split('  ',1)[1]:line.split('  ',1)[0] for line in (FROZEN/'checksums.sha256').read_text().splitlines()}
for name in ['case_outcome_flags.csv.gz','input_manifest.csv','analysis_manifest.json']:
    assert sha(FROZEN/name)==frozen_checks[name],('source checksum changed',name)
with gzip.open(FROZEN/'case_outcome_flags.csv.gz','rt',encoding='utf-8-sig') as f:
    cohort=list(csv.DictReader(f))
counts=dict(Counter(r['drug'] for r in cohort))
preflight={'expected_denominators':EXPECTED,'observed_denominators':counts,'rows':len(cohort),'unique_caseids':len({r['caseid'] for r in cohort}),'unique_primaryids':len({r['primaryid'] for r in cohort}),'frozen_source_sha256':sha(FROZEN/'case_outcome_flags.csv.gz')}
dump('preflight.json',preflight)
assert counts==EXPECTED,preflight
assert len(cohort)==len({r['caseid'] for r in cohort})==len({r['primaryid'] for r in cohort})==174406
log('PREFLIGHT PASSED',counts)
byq=defaultdict(dict)
for r in cohort:byq[r['source_quarter']][r['primaryid']]=r
db=sqlite3.connect('file:work/primary11/latest.sqlite?mode=ro',uri=True)
for r in cohort:
    winning=db.execute('SELECT pid,ver,fda,q FROM latest WHERE cid=?',(int(r['caseid']),)).fetchone()
    assert winning==(int(r['primaryid']),int(r['caseversion']),int(r['fda_dt']),r['source_quarter']),('retained version mismatch',r['caseid'])
db.close()
with (FROZEN/'input_manifest.csv').open(encoding='utf-8-sig') as f:manifest=list(csv.DictReader(f))
rows=[];quarter_qc=[];members=[]
for source in manifest:
    q=source['quarter'];path=SRC/source['filename'];eligible=byq[q]
    assert sha(path)==source['sha256'],('ZIP checksum changed',q)
    qc=Counter();selected={};outccodes=defaultdict(set);outcrows=Counter();outclines=defaultdict(list)
    with zipfile.ZipFile(path) as z:
        def reader(prefix):
            names=[n for n in z.namelist() if re.search(r'/(?:'+prefix+r')\d{2}Q[1-4]\.txt$',n,re.I)]
            assert len(names)==1,(q,prefix,names)
            name=names[0];info=z.getinfo(name)
            with z.open(name) as f:
                header=f.readline().decode('utf-8-sig').rstrip('\r\n').lower().split('$')
                members.append({'quarter':q,'table':prefix,'member':name,'crc32':f'{info.CRC:08x}','uncompressed_bytes':info.file_size,'header':header})
                for ln,b in enumerate(f,2):
                    vals=b.decode('utf-8').rstrip('\r\n').split('$')
                    assert len(vals)==len(header),(q,prefix,ln,'field count')
                    qc[prefix+'_rows']+=1
                    yield ln,header,vals
        demomap={}
        for ln,h,v in reader('DEMO'):
            pid,cid=v[:2]
            assert pid not in demomap,(q,'duplicate DEMO primaryid',pid)
            demomap[pid]=cid
            if pid not in eligible:continue
            d=dict(zip(h,v));r=eligible[pid]
            for key in ['primaryid','caseid','caseversion','fda_dt','init_fda_dt']:
                assert d[key]==r[key],(q,pid,key,d[key],r[key])
            selected[pid]={k:r[k] for k in ['drug','caseid','primaryid','caseversion','source_quarter','fda_dt','init_fda_dt']}
            selected[pid].update({k+'_raw':d[k] for k in ['age','age_cod','age_grp','sex','occp_cod','reporter_country','occr_country']})
            selected[pid]['demo_source_line']=ln
        assert set(selected)==set(eligible),(q,'unmatched DEMO',set(eligible)-set(selected))
        for ln,h,v in reader('OUTC'):
            assert h==['primaryid','caseid','outc_cod'],(q,h)
            pid,cid,code=v
            if demomap.get(pid)!=cid:
                qc['OUTC_key_mismatches']+=1
                if pid in eligible or cid in {r['caseid'] for r in eligible.values()}:
                    raise ValueError(('cohort OUTC linkage mismatch',q,pid,cid))
            if pid not in eligible:continue
            assert eligible[pid]['caseid']==cid
            qc['OUTC_matched_rows']+=1
            outcrows[pid]+=1;outclines[pid].append(ln)
            if code in outccodes[pid]:qc['OUTC_duplicate_same_code_rows']+=1
            outccodes[pid].add(code)
        for pid,r in selected.items():
            codes=sorted(outccodes[pid])
            r.update(outc_raw_codes=';'.join(codes),outc_row_count=outcrows[pid],outc_source_lines=';'.join(map(str,outclines[pid])))
            rows.append(r)
        qc['cohort_DEMO_matched_cases']=len(selected)
        qc['cohort_OUTC_matched_cases']=len(outcrows)
        qc['cohort_no_OUTC_rows']=len(selected)-len(outcrows)
        assert qc['OUTC_key_mismatches']==0,(q,'OUTC key mismatches',qc['OUTC_key_mismatches'])
    quarter_qc.append({'quarter':q,**dict(qc)})
    log(q,dict(qc))
assert len(rows)==174406
rows.sort(key=lambda r:(list(EXPECTED).index(r['drug']),int(r['caseid'])))
with gzip.open(WORK/'retained_characteristics_raw.csv.gz','wt',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
value_audit={field:{drug:dict(Counter(r[field] for r in rows if r['drug']==drug)) for drug in EXPECTED} for field in ['age_cod_raw','sex_raw','occp_cod_raw','reporter_country_raw','occr_country_raw','outc_raw_codes']}
dump('raw_value_audit.json',value_audit)
dump('extraction_qc.json',{'preflight':preflight,'final_row_count':len(rows),'retained_versions_equal_global_dedup_database':True,'quarter_qc':quarter_qc,'source_members':members,'source_manifest':manifest,'raw_extract_sha256':sha(WORK/'retained_characteristics_raw.csv.gz')})
log('COMPLETE')
