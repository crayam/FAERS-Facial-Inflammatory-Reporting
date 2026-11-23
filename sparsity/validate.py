from pathlib import Path
import csv,gzip,json,hashlib
from collections import Counter
R=Path(__file__).resolve().parents[2];O=R/'outputs/outcome_covariate_sparsity_v1_0';W=Path(__file__).parent
data=json.loads((W/'audit_data.json').read_text())
def read(p):
    op=gzip.open if str(p).endswith('.gz') else open
    with op(p,'rt',encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
cs=read(Path(data['sources'][0]['path']));os=read(Path(data['sources'][1]['path']))
cc={r['caseid']:r for r in cs};oc={r['caseid']:r for r in os}
assert len(cc)==len(cs)==len(oc)==len(os)==174406
assert cc.keys()==oc.keys()
counts=Counter();joint=Counter();base=Counter()
keys=['drug','primaryid','caseversion','source_quarter','fda_dt','init_fda_dt']
for cid,c in cc.items():
    o=oc[cid]
    assert all(c[k]==o[k] for k in keys)
    assert '20220401'<=c['fda_dt']<='20260630'
    assert c['calendar_year']==c['fda_dt'][:4]
    c['age_status']='Missing' if not c['age_years'] else 'Observed'
    base[c['drug']]+=1
    for y in data['outcomes']:
        for v in data['levels']:counts[y,v,c['drug'],c[v],int(o[y])]+=1
    for model in ['Main','S1_geography','S2_year_factor','S3_age_observed_base','S3_age_observed_adjusted']:
        if model.startswith('S3') and c['age_status']=='Missing':continue
        vs=['sex','hcp']+(['us_reporter'] if model=='S1_geography' else [])+['calendar_year']
        pat='; '.join(f'{v}={c[v]}' for v in vs)
        for y in data['outcomes']:joint[model,y,pat,c['drug'],int(o[y])]+=1
assert base==data['N']
for r in read(O/'stratum_counts.csv'):
    k=(r['outcome'],r['variable'],r['drug'],r['level'])
    assert int(r['positive'])==counts[k+(1,)] and int(r['negative'])==counts[k+(0,)]
    assert int(r['total'])==int(r['positive'])+int(r['negative'])
    for v in ['positive','negative']:
        n=int(r[v]);assert r[v+'_flag']==('ZERO' if n==0 else 'TINY 1–4' if n<5 else '')
jrows=read(O/'candidate_joint_counts.csv')
for r in jrows:
    k=(r['model'],r['outcome'],r['pattern'],r['drug'])
    assert int(r['positive'])==joint[k+(1,)] and int(r['negative'])==joint[k+(0,)]
v={'status':'PASS','source_key_pairs_verified':len(cc),'strata_verified':144,'positive_negative_cells_verified':288,'candidate_joint_rows_verified':len(jrows),'denominators':dict(base),'regressions_fitted':0,'independent_method':'Python csv, CASEID dictionaries and Counter; independent of pandas aggregation.'}
(O/'independent_validation.json').write_text(json.dumps(v,indent=2))
print(json.dumps(v,indent=2))
