"""Frozen Dictionary v1.0 Narrow6 / BroadSX23 union on the verified three-drug cohort."""
from pathlib import Path
from collections import defaultdict,Counter
import csv,gzip,hashlib,json
from docx import Document

SRC=Path('outputs/primary11_v1_0');OUT=Path('outputs/narrow_broad_v1_0');OUT.mkdir(parents=True,exist_ok=True)
DICT=Path('inputs/FAERS_Phenotype_Dictionary_v1.0.docx')
DRUGS=['dupilumab','tralokinumab','upadacitinib'];OUTCOMES=['Narrow6','Primary11','BroadSX23']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(name):
    with gzip.open(SRC/name,'rt',encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(name,rows):
    rows=list(rows);op=gzip.open if name.endswith('.gz') else open
    with op(OUT/name,'wt',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def dump(name,x): (OUT/name).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
checks={s.split('  ',1)[1]:s.split('  ',1)[0] for s in (SRC/'checksums.sha256').read_text().splitlines()}
for name in ['case_Primary11_flags.csv.gz','cohort_case_PT.csv.gz','analysis_manifest.json','independent_validation.json']:
    assert sha(SRC/name)==checks[name],('prior source changed',name)
prior=json.loads((SRC/'analysis_manifest.json').read_text())
assert sha(DICT)==prior['dictionary_sha256'], 'Dictionary changed since verified Primary11 run'
doc=Document(DICT)
def extract(hdr):
    tables=[t for t in doc.tables if [c.text.strip() for c in t.rows[0].cells]==hdr];assert len(tables)==1
    return [{'region':r.cells[0].text.strip(),'pt':r.cells[1].text.strip(),'rationale':r.cells[2].text.strip()} for r in tables[0].rows[1:]]
primary=extract(['分区','原始 REAC.PT','纳入理由']);narrow=extract(['分区','原始 REAC.PT','理由']);added=extract(['分区','新增 REAC.PT','用途/边界'])
definitions={'Narrow6':[r['pt'] for r in narrow],'Primary11':[r['pt'] for r in primary],'BroadSX23':[r['pt'] for r in primary+added]}
sets={k:set(v) for k,v in definitions.items()}
assert [len(sets[o]) for o in OUTCOMES]==[6,11,23]
assert sets['Narrow6']<sets['Primary11']<sets['BroadSX23']
assert len(added)==12 and not(sets['Primary11']&{r['pt'] for r in added})
old=read('case_Primary11_flags.csv.gz');pairs=read('cohort_case_PT.csv.gz')
bykey={(r['drug'],r['caseid']):r for r in old};assert len(bykey)==len(old)==174406
pts=defaultdict(set);seen=set()
for r in pairs:
    key=(r['drug'],r['caseid']);assert key in bykey
    assert r['primaryid']==bykey[key]['primaryid'] and r['source_quarter']==bykey[key]['source_quarter']
    pk=key+(r['pt'],);assert pk not in seen;seen.add(pk);pts[key].add(r['pt'])
assert len(pairs)==455922 and set(pts)==set(bykey)
flags=[];patterns=Counter();pt_counter={d:Counter() for d in DRUGS};outcome_sets={d:{o:set() for o in OUTCOMES} for d in DRUGS}
for oldr in old:
    key=(oldr['drug'],oldr['caseid']);sourcepts=pts[key]
    r={k:v for k,v in oldr.items() if k not in ['Primary11','matched_pt_count','matched_pts']}
    for o in OUTCOMES:
        hits=[pt for pt in definitions[o] if pt in sourcepts]
        r[o]=int(bool(hits));r[o+'_matched_pt_count']=len(hits);r[o+'_matched_pts']='; '.join(hits)
        if hits:
            outcome_sets[r['drug']][o].add(r['caseid'])
            patterns[(r['drug'],o,r[o+'_matched_pts'],len(hits))]+=1
    assert r['Primary11']==int(oldr['Primary11'])
    assert r['Primary11_matched_pt_count']==int(oldr['matched_pt_count']) and r['Primary11_matched_pts']==oldr['matched_pts']
    assert r['Narrow6']<=r['Primary11']<=r['BroadSX23']
    r['Primary11_not_Narrow6']=int(r['Primary11']==1 and r['Narrow6']==0)
    r['BroadSX23_not_Primary11']=int(r['BroadSX23']==1 and r['Primary11']==0)
    flags.append(r);pt_counter[r['drug']].update(sourcepts)
summary=[];overlap=[];quarterly=[]
for d in DRUGS:
    rows=[r for r in flags if r['drug']==d];n=len(rows)
    row={'drug':d,'comparison_role':rows[0]['comparison_role'],'ps_strict_ad_cases':n}
    for o in OUTCOMES:
        positives=len(outcome_sets[d][o]);ptsum=sum(pt_counter[d][pt] for pt in sets[o])
        row.update({o+'_cases':positives,o+'_negative_cases':n-positives,o+'_proportion':positives/n})
        overlap.append({'drug':d,'outcome':o,'dictionary_pt_count':len(sets[o]),'sum_individual_pt_cases':ptsum,'case_union':positives,'duplicate_count_excess':ptsum-positives,'cases_with_multiple_pts':sum(r[o+'_matched_pt_count']>1 for r in rows),'max_pts_per_case':max(r[o+'_matched_pt_count'] for r in rows)})
    row['Primary11_not_Narrow6_cases']=len(outcome_sets[d]['Primary11']-outcome_sets[d]['Narrow6'])
    row['BroadSX23_not_Primary11_cases']=len(outcome_sets[d]['BroadSX23']-outcome_sets[d]['Primary11'])
    summary.append(row)
    for q in prior['source_periods']:
        qr=[r for r in rows if r['source_quarter']==q]
        quarterly.append({'drug':d,'source_quarter':q,'ps_strict_ad_cases':len(qr),**{o+'_cases':sum(r[o] for r in qr) for o in OUTCOMES}})
dictionary=[]
for r in primary+added:
    dictionary.append({**r,**{o:int(r['pt'] in sets[o]) for o in OUTCOMES},**{d:pt_counter[d][r['pt']] for d in DRUGS}})
write('union_summary.csv',summary);write('case_outcome_flags.csv.gz',flags)
write('positive_case_evidence.csv',[r for r in flags if r['BroadSX23']])
write('primary11_additional_to_narrow6.csv',[r for r in flags if r['Primary11_not_Narrow6']])
write('broadsx23_additional_to_primary11.csv',[r for r in flags if r['BroadSX23_not_Primary11']])
write('dictionary_pt_counts.csv',dictionary);write('overlap_audit.csv',overlap);write('quarterly_counts.csv',quarterly)
write('matched_pt_patterns.csv',[{'drug':d,'outcome':o,'matched_pts':pp,'matched_pt_count':k,'case_count':v} for (d,o,pp,k),v in sorted(patterns.items())])
dump('phenotype_definition.json',{'version':'1.0','dictionary_sha256':sha(DICT),'definitions':definitions,'BroadSX23_rule':'Primary11 union 12 additional symptom PTs; 23 total unique PTs'})
manifest={'source':'Verified 2022Q2–2026Q2 PS + drug-specific exact Dermatitis atopic cohort from previous Primary11 run; no change to enrollment, deduplication, dates, DELETE, drug names or role rules.','source_sha256':{k:sha(SRC/k) for k in ['case_Primary11_flags.csv.gz','cohort_case_PT.csv.gz','analysis_manifest.json','independent_validation.json']},'dictionary_sha256':sha(DICT),'raw_zip_sha256':prior['input_sha256'],'source_periods':prior['source_periods'],'cohort_rows':len(flags),'unique_caseid':len({r['caseid'] for r in flags}),'cohort_case_pt_pairs':len(pairs),'dictionary_counts':{o:len(sets[o]) for o in OUTCOMES},'all_primary11_flags_counts_and_match_lists_unchanged':True,'Narrow6_subset_Primary11_subset_BroadSX23_every_case':True,'nesting_violations':0,'duplicate_case_pt_pairs':0,'summary':summary,'overlap':overlap}
dump('analysis_manifest.json',manifest)
Path('work/narrow_broad_workbook_data.json').write_text(json.dumps({'summary':summary,'overlap':overlap,'dictionary':dictionary,'manifest':manifest},ensure_ascii=False,indent=2))
print(json.dumps({'summary':summary,'overlap':overlap},ensure_ascii=False,indent=2))
