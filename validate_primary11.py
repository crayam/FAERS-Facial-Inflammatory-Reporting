"""Independent CSV-based reconstruction from ZIP sources. Does not import production code."""
from pathlib import Path
from collections import defaultdict, Counter
import csv,gzip,io,json,re,zipfile,hashlib,time
import openpyxl

OUT=Path('outputs/primary11_v1_0'); BASE=Path('inputs/faers_ascii')
definition=json.loads((OUT/'primary11_definition.json').read_text())
quarters=definition['periods']; pts=set(x['pt'] for x in definition['primary11'])
names={'dupilumab':r'\b(dupilumab|dupixent)\b','tralokinumab':r'\b(tralokinumab|adbry|adtralza)\b','upadacitinib':r'\b(upadacitinib|rinvoq)\b'}
regex={d:re.compile(p,re.I) for d,p in names.items()}
files={re.search(r'(20\d{2}q[1-4])',p.name,re.I)[0].upper():p for p in BASE.glob('*.zip')}
def rows(z,prefix):
    n=next(n for n in z.namelist() if Path(n).name.upper().startswith(prefix) and n.lower().endswith('.txt'))
    with z.open(n) as f:
        yield from csv.DictReader(io.TextIOWrapper(f,encoding='utf-8-sig',newline=''),delimiter='$',quoting=csv.QUOTE_NONE)
def log(*s):print(time.strftime('%H:%M:%S'),*s,flush=True)

# Build eligibility for every source report first, then remove superseded versions globally.
candidate={}; synonym=Counter(); deleted=set()
for q in quarters:
    with zipfile.ZipFile(files[q]) as z:
        drug_map=defaultdict(set)
        for r in rows(z,'DRUG'):
            if r['role_cod']!='PS':continue
            ds={d for d,pat in regex.items() if pat.search(r['drugname']) or pat.search(r['prod_ai'])}
            assert len(ds)<=1
            if ds:drug_map[(r['primaryid'],r['caseid'],r['drug_seq'])].update(ds)
        for r in rows(z,'INDI'):
            if r['indi_pt'] in {'Atopic dermatitis','Atopic eczema'}:synonym[r['indi_pt']]+=1
            if r['indi_pt']!='Dermatitis atopic':continue
            key=(r['primaryid'],r['caseid'],r['indi_drug_seq'])
            for d in drug_map.get(key,set()):candidate[(q,r['primaryid'],r['caseid'],d)]=True
        n=next(n for n in z.namelist() if Path(n).name.upper().startswith('DELETE') and n.lower().endswith('.txt'))
        with z.open(n) as f:
            for line in f:
                s=line.strip()
                if s:deleted.add(s.decode())
    log('Validation eligibility',q)
cids={k[2] for k in candidate};latest={}
for q in quarters:
    with zipfile.ZipFile(files[q]) as z:
        for r in rows(z,'DEMO'):
            c=r['caseid']
            if c not in cids or c in deleted:continue
            rank=(int(r['caseversion']),int(r['fda_dt']),int(r['primaryid']),q)
            if c not in latest or rank>latest[c]:latest[c]=rank
    log('Validation latest',q)
expected={}
for q,p,c,d in candidate:
    if c not in latest:continue
    v,date,pid,source=latest[c]
    if int(p)==pid and q==source and 20220401<=date<=20260630:expected[(d,c)]=(q,p)
with gzip.open(OUT/'case_Primary11_flags.csv.gz','rt',encoding='utf-8-sig',newline='') as f:actual={(r['drug'],r['caseid']):r for r in csv.DictReader(f)}
assert set(actual)==set(expected),('cohort mismatch',len(set(actual)-set(expected)),len(set(expected)-set(actual)))
allpts=defaultdict(set); hits=defaultdict(set)
for q in quarters:
    pm=defaultdict(list)
    for key,(qr,pid) in expected.items():
        if qr==q:pm[pid].append(key)
    with zipfile.ZipFile(files[q]) as z:
        for r in rows(z,'REAC'):
            for key in pm.get(r['primaryid'],[]):
                assert key[1]==r['caseid']
                allpts[key].add(r['pt'])
                if r['pt'] in pts:hits[key].add(r['pt'])
    log('Validation REAC',q)
for k,r in actual.items():
    assert int(r['Primary11'])==int(bool(hits[k]))
    assert int(r['matched_pt_count'])==len(hits[k])
    assert (set(r['matched_pts'].split('; ')) if r['matched_pts'] else set())==hits[k]
    assert int(r['all_reac_unique_pt_count'])==len(allpts[k])
    assert (r['source_quarter'],r['primaryid'])==expected[k]
independent=Counter(d for (d,c) in expected if hits[(d,c)])
pt_counts={d:Counter() for d in names}
for (d,c),pp in allpts.items():pt_counts[d].update(pp)
prior='inputs/historical_full_PT_matrix.xlsx'
wb=openpyxl.load_workbook(prior,read_only=True,data_only=True)
s=wb['全量PT索引']; mismatches=[]; cells=0;prior_pts=set()
for row in s.iter_rows(min_row=6,values_only=True):
    if not row[1]:continue
    pt=row[1];prior_pts.add(pt)
    for i,d in enumerate(names,2):
        cells+=1
        if row[i]!=pt_counts[d][pt]:mismatches.append({'pt':pt,'drug':d,'prior':row[i],'new':pt_counts[d][pt]})
observed_pts=set().union(*(set(c) for c in pt_counts.values()))
assert not mismatches,('prior PT matrix mismatch',mismatches[:10])
assert prior_pts==observed_pts
assert {d:sum(k[0]==d for k in expected) for d in names}=={'dupilumab':164515,'tralokinumab':3637,'upadacitinib':6254}
report={'method':'Independent csv.DictReader reconstruction from all 17 raw ZIPs: qualify all PS drug-specific exact AD reports first, then select each candidate CASEID latest record across all DEMO and apply DELETE, then rebuild all PT sets from raw REAC. Production code and its intermediate database were not imported.','cohort_rows_validated':len(actual),'Primary11_union':dict(independent),'every_case_flag_and_match_list_equal':True,'source_lineage_equal':True,'prior_full_PT_matrix_cells_compared':cells,'prior_full_PT_matrix_mismatches':mismatches,'prior_unique_PT_set_equal':True,'unique_PT_count':len(observed_pts),'synonym_audit_all_raw_INDI':{'Atopic dermatitis':synonym['Atopic dermatitis'],'Atopic eczema':synonym['Atopic eczema']},'prior_workbook_sha256':hashlib.sha256(Path(prior).read_bytes()).hexdigest()}
(OUT/'independent_validation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
log('VALIDATION PASSED',report)
