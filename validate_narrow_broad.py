"""Independent raw-REAC verification on the previously validated, unchanged cohort."""
from pathlib import Path
from collections import defaultdict,Counter
import csv,gzip,hashlib,json,zipfile,time
OUT=Path('outputs/narrow_broad_v1_0');OLD=Path('outputs/primary11_v1_0');RAW=Path('inputs/faers_ascii')
N={'Erythema of eyelid','Eczema eyelids','Periorbital dermatitis','Periorbital eczema','Perioral dermatitis','Eczematous cheilitis'}
P=N|{'Eyelid rash','Periorbital inflammation','Cheilitis','Angular cheilitis','Lip erythema'}
B=P|{'Eyelid erosion','Eyelid exfoliation','Eyelid irritation','Eyelid margin crusting','Eyelid skin dryness','Eyelid thickening','Eyelids pruritus','Periorbital irritation','Chapped lips','Lip dry','Lip exfoliation','Lip pruritus'}
D={'Narrow6':N,'Primary11':P,'BroadSX23':B}
definition=json.loads((OUT/'phenotype_definition.json').read_text())
assert all(set(definition['definitions'][o])==s for o,s in D.items())
def read(folder,name):
    op=gzip.open if name.endswith('.gz') else open
    with op(folder/name,'rt',encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
old=read(OLD,'case_Primary11_flags.csv.gz');new=read(OUT,'case_outcome_flags.csv.gz')
assert [(r['drug'],r['caseid']) for r in new]==[(r['drug'],r['caseid']) for r in old]
bykey={(r['drug'],r['caseid']):r for r in new};rawpts=defaultdict(set);evidence=[];raw_rows=0
manifest=read(OLD,'input_manifest.csv')
for m in manifest:
    q=m['quarter'];file=RAW/m['filename']
    h=hashlib.sha256()
    with file.open('rb') as f:
        for chunk in iter(lambda:f.read(2**20),b''):h.update(chunk)
    assert h.hexdigest()==m['sha256']
    selected=defaultdict(list)
    for r in old:
        if r['source_quarter']==q:selected[r['primaryid']].append((r['drug'],r['caseid']))
    with zipfile.ZipFile(file) as z:
        member=next(n for n in z.namelist() if Path(n).name.upper().startswith('REAC') and n.lower().endswith('.txt'))
        with z.open(member) as f:
            header=f.readline().decode('utf-8-sig').rstrip('\r\n').split('$')
            assert header[:3]==['primaryid','caseid','pt']
            for line,chunk in enumerate(f,2):
                parts=chunk.decode('utf-8').rstrip('\r\n').split('$');assert len(parts)==len(header)
                raw_rows+=1;pid,cid,pt=parts[:3]
                for k in selected.get(pid,[]):
                    assert cid==k[1];rawpts[k].add(pt)
                    if pt in B:evidence.append({'drug':k[0],'caseid':cid,'primaryid':pid,'source_quarter':q,'reac_source_line':line,'pt':pt,'Narrow6_pt':int(pt in N),'Primary11_pt':int(pt in P),'BroadSX23_pt':1})
    print(time.strftime('%H:%M:%S'),'Verified raw REAC',q,flush=True)
assert len(rawpts)==len(old)
for r,oldr in zip(new,old):
    key=(r['drug'],r['caseid'])
    for c,v in oldr.items():
        if c not in {'matched_pt_count','matched_pts'}:assert r[c]==v
    for o,s in D.items():
        matched=rawpts[key]&s
        assert int(r[o])==bool(matched)
        assert int(r[o+'_matched_pt_count'])==len(matched)
        assert (set(r[o+'_matched_pts'].split('; ')) if r[o+'_matched_pts'] else set())==matched
    assert int(r['Narrow6'])<=int(r['Primary11'])<=int(r['BroadSX23'])
    assert int(r['Primary11_not_Narrow6'])==int(r['Primary11'])-int(r['Narrow6'])
    assert int(r['BroadSX23_not_Primary11'])==int(r['BroadSX23'])-int(r['Primary11'])
    assert len(rawpts[key])==int(r['all_reac_unique_pt_count'])
positive=read(OUT,'positive_case_evidence.csv');assert positive==[r for r in new if r['BroadSX23']=='1']
for name,col in [('primary11_additional_to_narrow6.csv','Primary11_not_Narrow6'),('broadsx23_additional_to_primary11.csv','BroadSX23_not_Primary11')]:
    assert read(OUT,name)==[r for r in new if r[col]=='1']
summary=read(OUT,'union_summary.csv')
for r in summary:
    drug_rows=[x for x in new if x['drug']==r['drug']];assert len(drug_rows)==int(r['ps_strict_ad_cases'])
    for o in D:assert sum(int(x[o]) for x in drug_rows)==int(r[o+'_cases'])
counter={d:Counter() for d in ['dupilumab','tralokinumab','upadacitinib']}
for (d,c),pts in rawpts.items():counter[d].update(pts)
for r in read(OUT,'dictionary_pt_counts.csv'):
    for d in counter:assert int(r[d])==counter[d][r['pt']]
for r in read(OUT,'overlap_audit.csv'):
    d,o=r['drug'],r['outcome'];values=[len(pp&D[o]) for (drug,c),pp in rawpts.items() if drug==d]
    assert sum(values)==int(r['sum_individual_pt_cases'])
    assert sum(v>0 for v in values)==int(r['case_union'])
    assert sum(values)-sum(v>0 for v in values)==int(r['duplicate_count_excess'])
    assert sum(v>1 for v in values)==int(r['cases_with_multiple_pts'])
    assert max(values)==int(r['max_pts_per_case'])
with (OUT/'raw_REAC_evidence.csv').open('w',encoding='utf-8-sig',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(evidence[0]));writer.writeheader();writer.writerows(evidence)
validation={'method':'Independent fixed PT definitions and re-reading all 17 original raw REAC files on the unchanged, previously independently validated PS + drug-specific strict AD cohort. No rebuild of enrollment was needed.','raw_zip_sha256_all_match_previous_run':True,'raw_REAC_rows_read':raw_rows,'case_rows_verified':len(new),'all_case_metadata_unchanged':True,'all_three_flags_and_matched_PT_sets_equal_raw_REAC':True,'Primary11_unchanged':True,'nesting_violations':0,'summary_dictionary_and_overlap_equal_raw_REAC':True,'positive_and_incremental_case_files_equal_flags':True,'positive_broad_cases':len(positive),'raw_broad_evidence_rows':len(evidence),'summary':{r['drug']:{o:int(r[o+'_cases']) for o in D} for r in summary}}
(OUT/'independent_validation.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2))
print(json.dumps(validation,ensure_ascii=False,indent=2))
