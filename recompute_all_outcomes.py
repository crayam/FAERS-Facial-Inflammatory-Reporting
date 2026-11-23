"""Generate all seven prespecified outcomes on exactly the frozen formal cohort."""
from pathlib import Path
from collections import defaultdict,Counter
import csv,gzip,hashlib,json
from docx import Document
BASE=Path('outputs/primary11_v1_0');PREV=Path('outputs/narrow_broad_v1_0');OUT=Path('outputs/all_outcomes_v1_0');OUT.mkdir(parents=True,exist_ok=True)
DICT=Path('inputs/FAERS_Phenotype_Dictionary_v1.0.docx')
DRUGS=['dupilumab','tralokinumab','upadacitinib']
EXISTING=['Primary11','Narrow6','BroadSX23']
NEW=['Rosacea_family','Demodex_family','Blepharitis_family','Seborrhoeic_dermatitis']
OUTCOMES=EXISTING+NEW
TIERS={'Primary11':'Primary','Narrow6':'Sensitivity: high specificity','BroadSX23':'Sensitivity: symptom enriched',**{o:'Secondary' for o in NEW}}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(folder,name):
    op=gzip.open if name.endswith('.gz') else open
    with op(folder/name,'rt',encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(name,rows):
    rows=list(rows);op=gzip.open if name.endswith('.gz') else open
    with op(OUT/name,'wt',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def dump(name,obj):(OUT/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2))
source_sha={}
for folder,names in [(BASE,['cohort_case_PT.csv.gz','analysis_manifest.json']),(PREV,['case_outcome_flags.csv.gz','phenotype_definition.json','independent_validation.json'])]:
    checks={line.split('  ',1)[1]:line.split('  ',1)[0] for line in (folder/'checksums.sha256').read_text().splitlines()}
    for name in names:
        source_sha[f'{folder.name}/{name}']=sha(folder/name)
        assert source_sha[f'{folder.name}/{name}']==checks[name],('Changed input',name)
prior=json.loads((BASE/'analysis_manifest.json').read_text());assert sha(DICT)==prior['dictionary_sha256']
doc=Document(DICT)
def table(hdr):
    t=[t for t in doc.tables if [c.text.strip() for c in t.rows[0].cells]==hdr];assert len(t)==1
    return [r.cells[1].text.strip() for r in t[0].rows[1:]]
primary=table(['分区','原始 REAC.PT','纳入理由']);narrow=table(['分区','原始 REAC.PT','理由']);extra=table(['分区','新增 REAC.PT','用途/边界'])
def family(first):
    ts=[t for t in doc.tables if len(t.rows)>1 and t.rows[0].cells[0].text=='PT' and t.rows[1].cells[0].text==first];assert len(ts)==1
    return [r.cells[0].text.strip() for r in ts[0].rows[1:]]
seb=next(p.text for p in doc.paragraphs if p.text.startswith('正式 secondary outcome 仅使用精确 PT：')).split('：',1)[1].rstrip('。').strip()
definitions={'Primary11':primary,'Narrow6':narrow,'BroadSX23':primary+extra,'Rosacea_family':family('Rosacea'),'Demodex_family':family('Demodicidosis'),'Blepharitis_family':family('Blepharitis'),'Seborrhoeic_dermatitis':[seb]}
sets={o:set(v) for o,v in definitions.items()}
assert [len(sets[o]) for o in OUTCOMES]==[11,6,23,3,3,6,1]
assert sets['Narrow6']<sets['Primary11']<sets['BroadSX23']
assert sets['Demodex_family']&sets['Blepharitis_family']=={'Demodex blepharitis'}
frozen=json.loads((PREV/'phenotype_definition.json').read_text())['definitions'];assert all(definitions[o]==frozen[o] for o in EXISTING)
old=read(PREV,'case_outcome_flags.csv.gz');pairs=read(BASE,'cohort_case_PT.csv.gz')
bykey={(r['drug'],r['caseid']):r for r in old};assert len(bykey)==len(old)==174406
denominators=Counter(r['drug'] for r in old);assert denominators==dict(zip(DRUGS,[164515,3637,6254]))
pts=defaultdict(set);seen=set()
for r in pairs:
    key=(r['drug'],r['caseid']);assert key in bykey
    assert r['primaryid']==bykey[key]['primaryid'] and r['source_quarter']==bykey[key]['source_quarter']
    assert (*key,r['pt']) not in seen;seen.add((*key,r['pt']));pts[key].add(r['pt'])
assert set(pts)==set(bykey) and len(pairs)==455922
flags=[];patterns=Counter();ptcounts={d:Counter() for d in DRUGS}
for oldr in old:
    r=dict(oldr);key=(r['drug'],r['caseid']);ptcounts[r['drug']].update(pts[key])
    for o in OUTCOMES:
        hits=[pt for pt in definitions[o] if pt in pts[key]]
        flag=int(bool(hits));hitcount=len(hits);hitlist='; '.join(hits)
        if o in EXISTING:
            assert flag==int(oldr[o]) and hitcount==int(oldr[o+'_matched_pt_count']) and hitlist==oldr[o+'_matched_pts']
        else:r.update({o:flag,o+'_matched_pt_count':hitcount,o+'_matched_pts':hitlist})
        if flag:patterns[(r['drug'],o,hitlist,hitcount)]+=1
    assert all(r[k]==v for k,v in oldr.items())
    assert int(r['Narrow6'])<=int(r['Primary11'])<=int(r['BroadSX23'])
    flags.append(r)
summary=[];wide=[];overlap=[];dictionary=[];quarters=[]
for o in OUTCOMES:
    w={'outcome':o,'analysis_tier':TIERS[o],'dictionary_pt_count':len(sets[o])}
    for d in DRUGS:
        rows=[r for r in flags if r['drug']==d];n=denominators[d];a=sum(int(r[o]) for r in rows)
        summary.append({'outcome':o,'analysis_tier':TIERS[o],'dictionary_pt_count':len(sets[o]),'drug':d,'comparison_role':rows[0]['comparison_role'],'denominator_N':n,'positive_n':a,'negative_n':n-a,'proportion':a/n})
        w.update({d+'_N':n,d+'_n':a,d+'_proportion':a/n})
        total=sum(ptcounts[d][pt] for pt in sets[o]);overlap.append({'drug':d,'outcome':o,'sum_individual_pt_cases':total,'case_union':a,'duplicate_count_excess':total-a,'cases_with_multiple_pts':sum(int(r[o+'_matched_pt_count'])>1 for r in rows),'max_pts_per_case':max(int(r[o+'_matched_pt_count']) for r in rows)})
    wide.append(w)
    for pt in definitions[o]:dictionary.append({'outcome':o,'pt':pt,**{d:ptcounts[d][pt] for d in DRUGS}})
for d in DRUGS:
    for q in prior['source_periods']:
        rows=[r for r in flags if r['drug']==d and r['source_quarter']==q]
        quarters.append({'drug':d,'source_quarter':q,'denominator_N':len(rows),**{o:sum(int(r[o]) for r in rows) for o in OUTCOMES}})
write('case_outcome_flags.csv.gz',flags)
write('outcome_summary.csv',summary);write('outcome_summary_wide.csv',wide)
write('positive_case_evidence.csv',[r for r in flags if any(int(r[o]) for o in OUTCOMES)])
write('secondary_positive_case_evidence.csv',[r for r in flags if any(int(r[o]) for o in NEW)])
write('dictionary_pt_counts.csv',dictionary);write('within_outcome_overlap_audit.csv',overlap);write('quarterly_counts.csv',quarters)
write('matched_pt_patterns.csv',[{'drug':d,'outcome':o,'matched_pts':pp,'matched_pt_count':k,'case_count':v} for (d,o,pp,k),v in sorted(patterns.items())])
dump('phenotype_definition.json',{'version':'1.0','frozen_date':'2026-09-20','dictionary_sha256':sha(DICT),'definitions':definitions,'count_unit':'drug + CASEID; one retained PRIMARYID and source quarter per CASEID','match_mode':'Exact REAC.PT text. No synonym expansion or MedDRA remapping.'})
manifest={'cohort_source':'Same verified 2022Q2–2026Q2 formal PS + drug-specific exact Dermatitis atopic cohort. Enrollment not rerun or modified.','source_sha256':source_sha,'dictionary_sha256':sha(DICT),'raw_zip_sha256':prior['input_sha256'],'source_periods':prior['source_periods'],'cohort_rows':len(flags),'cohort_unique_caseid':len({r['caseid'] for r in flags}),'denominators':dict(denominators),'all_existing_columns_unchanged':True,'all_existing_three_flags_and_hit_lists_regenerated_equal':True,'new_outcomes_only':NEW,'all_seven_outcomes':OUTCOMES,'nesting_violations':0,'source_case_PT_pairs':len(pairs),'cohort_window':'20220401 <= retained FDA_DT <= 20260630','strict_AD':'ROLE_COD=PS and INDI_PT exactly Dermatitis atopic linked on PRIMARYID + CASEID + DRUG_SEQ/INDI_DRUG_SEQ','dedup_order':prior['dedup_order'],'delete_rule':'Same previously applied 17-quarter DELETE union','summary':summary}
dump('analysis_manifest.json',manifest)
Path('work/all_outcomes_workbook_data.json').write_text(json.dumps({'summary':summary,'wide':wide,'dictionary':dictionary,'overlap':overlap,'definitions':definitions,'manifest':manifest},ensure_ascii=False,indent=2))
print(json.dumps(wide,ensure_ascii=False,indent=2))
