"""Independent raw REAC validation; frozen enrollment is inherited and hash-checked."""
from pathlib import Path
from collections import defaultdict,Counter
import csv,gzip,hashlib,json,zipfile,time
OUT=Path('outputs/all_outcomes_v1_0');PREV=Path('outputs/narrow_broad_v1_0');BASE=Path('outputs/primary11_v1_0');RAW=Path('inputs/faers_ascii')
N={'Erythema of eyelid','Eczema eyelids','Periorbital dermatitis','Periorbital eczema','Perioral dermatitis','Eczematous cheilitis'}
P=N|{'Eyelid rash','Periorbital inflammation','Cheilitis','Angular cheilitis','Lip erythema'}
B=P|{'Eyelid erosion','Eyelid exfoliation','Eyelid irritation','Eyelid margin crusting','Eyelid skin dryness','Eyelid thickening','Eyelids pruritus','Periorbital irritation','Chapped lips','Lip dry','Lip exfoliation','Lip pruritus'}
D={'Primary11':P,'Narrow6':N,'BroadSX23':B,'Rosacea_family':{'Rosacea','Ocular rosacea','Papulopustular rosacea'},'Demodex_family':{'Demodicidosis','Ocular demodicosis','Demodex blepharitis'},'Blepharitis_family':{'Blepharitis','Blepharitis allergic','Demodex blepharitis','Meibomianitis','Blepharokeratoconjunctivitis','Staphylococcal blepharitis'},'Seborrhoeic_dermatitis':{'Seborrhoeic dermatitis'}}
def read(folder,name):
    op=gzip.open if name.endswith('.gz') else open
    with op(folder/name,'rt',encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
defs=json.loads((OUT/'phenotype_definition.json').read_text())['definitions'];assert all(set(defs[o])==v for o,v in D.items())
old=read(PREV,'case_outcome_flags.csv.gz');new=read(OUT,'case_outcome_flags.csv.gz')
assert len(old)==len(new)==174406
assert [(r['drug'],r['caseid']) for r in new]==[(r['drug'],r['caseid']) for r in old]
for r,prev in zip(new,old):assert all(r[k]==v for k,v in prev.items())
rawpts=defaultdict(set);evidence=[];rawrows=0;allpts=set().union(*D.values())
for m in read(BASE,'input_manifest.csv'):
    file=RAW/m['filename'];q=m['quarter'];h=hashlib.sha256()
    with file.open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    assert h.hexdigest()==m['sha256']
    selected=defaultdict(list)
    for r in old:
        if r['source_quarter']==q:selected[r['primaryid']].append((r['drug'],r['caseid']))
    with zipfile.ZipFile(file) as z:
        member=next(n for n in z.namelist() if Path(n).name.upper().startswith('REAC') and n.lower().endswith('.txt'))
        with z.open(member) as f:
            header=f.readline().decode('utf-8-sig').rstrip('\r\n').split('$');assert header[:3]==['primaryid','caseid','pt']
            for line,b in enumerate(f,2):
                r=b.decode('utf-8').rstrip('\r\n').split('$');assert len(r)==len(header)
                rawrows+=1;pid,cid,pt=r[:3]
                for k in selected.get(pid,[]):
                    assert cid==k[1];rawpts[k].add(pt)
                    if pt in allpts:evidence.append({'drug':k[0],'caseid':cid,'primaryid':pid,'source_quarter':q,'reac_source_line':line,'pt':pt,'matching_outcomes':'; '.join(o for o in D if pt in D[o])})
    print(time.strftime('%H:%M:%S'),'Verified',q,flush=True)
assert len(rawpts)==len(old)
shared=[]
for r in new:
    key=(r['drug'],r['caseid']);pp=rawpts[key]
    for o,s in D.items():
        hits=pp&s
        assert int(r[o])==bool(hits) and int(r[o+'_matched_pt_count'])==len(hits)
        assert (set(r[o+'_matched_pts'].split('; ')) if r[o+'_matched_pts'] else set())==hits
    assert int(r['Narrow6'])<=int(r['Primary11'])<=int(r['BroadSX23'])
    if 'Demodex blepharitis' in pp:
        assert r['Demodex_family']==r['Blepharitis_family']=='1'
        shared.append({'drug':r['drug'],'caseid':r['caseid'],'primaryid':r['primaryid']})
ptcounts={d:Counter() for d in ['dupilumab','tralokinumab','upadacitinib']}
for (d,c),pp in rawpts.items():ptcounts[d].update(pp)
for r in read(OUT,'outcome_summary.csv'):
    d,o=r['drug'],r['outcome'];n=sum(x['drug']==d for x in new);a=sum(bool(pp&D[o]) for (drug,c),pp in rawpts.items() if drug==d)
    assert int(r['denominator_N'])==n and int(r['positive_n'])==a and int(r['negative_n'])==n-a
    assert abs(float(r['proportion'])-a/n)<1e-12
for r in read(OUT,'dictionary_pt_counts.csv'):
    for d in ptcounts:assert int(r[d])==ptcounts[d][r['pt']]
for r in read(OUT,'within_outcome_overlap_audit.csv'):
    d,o=r['drug'],r['outcome'];v=[len(pp&D[o]) for (drug,c),pp in rawpts.items() if drug==d]
    assert [int(r[k]) for k in ['sum_individual_pt_cases','case_union','duplicate_count_excess','cases_with_multiple_pts','max_pts_per_case']]==[sum(v),sum(x>0 for x in v),sum(v)-sum(x>0 for x in v),sum(x>1 for x in v),max(v)]
assert read(OUT,'positive_case_evidence.csv')==[r for r in new if any(r[o]=='1' for o in D)]
assert read(OUT,'secondary_positive_case_evidence.csv')==[r for r in new if any(r[o]=='1' for o in list(D)[3:])]
with (OUT/'raw_REAC_evidence.csv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(evidence[0]));w.writeheader();w.writerows(evidence)
report={'method':'Independent fixed PT sets and rereading all 17 raw REAC files on the same frozen formal cohort. Enrollment inherited from previously independently validated cohort.','raw_zip_sha256_all_unchanged':True,'raw_REAC_rows_read':rawrows,'case_rows_verified':len(new),'all_existing_columns_unchanged':True,'all_seven_outcome_flags_and_hit_lists_equal_raw_REAC':True,'all_summary_and_dictionary_and_within_outcome_overlap_counts_equal_raw_REAC':True,'positive_evidence_files_equal_case_flags':True,'nesting_violations':0,'Demodex_blepharitis_in_both_families_verified_cases':shared,'all_outcome_positive_cases':len(read(OUT,'positive_case_evidence.csv')),'any_secondary_positive_cases':len(read(OUT,'secondary_positive_case_evidence.csv')),'raw_evidence_rows':len(evidence)}
(OUT/'independent_validation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
print(json.dumps(report,ensure_ascii=False,indent=2))
