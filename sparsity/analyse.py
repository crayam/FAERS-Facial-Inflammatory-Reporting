from pathlib import Path
import pandas as pd
import numpy as np
import hashlib, json, csv

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/outcome_covariate_sparsity_v1_0'
OUT.mkdir(parents=True,exist_ok=True)
WORK=Path(__file__).parent
paths=[ROOT/'outputs/table1_missingness_v1_0/case_characteristics.csv.gz',ROOT/'outputs/all_outcomes_v1_0/case_outcome_flags.csv.gz']
sources=[]
for p in paths:
    h=hashlib.sha256(p.read_bytes()).hexdigest()
    checks=(p.parent/'checksums.sha256').read_text()
    assert any(line.split()[0]==h and line.rstrip().endswith(p.name) for line in checks.splitlines()),p
    sources.append({'path':str(p),'sha256':h})
c,o=[pd.read_csv(p,dtype=str,keep_default_na=False,encoding='utf-8-sig') for p in paths]
KEY=['drug','caseid','primaryid','caseversion','source_quarter','fda_dt','init_fda_dt']
outcomes=['Primary11','Narrow6','BroadSX23']
drugs=['dupilumab','tralokinumab','upadacitinib']
Ns=dict(zip(drugs,[164515,3637,6254]))
expected={'Primary11':[1068,19,21],'Narrow6':[754,13,14],'BroadSX23':[2347,51,30]}
assert len(c)==len(o)==174406
assert c.caseid.is_unique and o.caseid.is_unique and c.primaryid.is_unique and o.primaryid.is_unique
assert set(map(tuple,c[KEY].values))==set(map(tuple,o[KEY].values))
d=c.merge(o[KEY+outcomes],on=KEY,validate='one_to_one')
assert d.groupby('drug').size().to_dict()==Ns
for y in outcomes:
    assert set(d[y])=={'0','1'}
    d[y]=d[y].astype(int)
    assert d.groupby('drug')[y].sum().reindex(drugs).tolist()==expected[y]
assert ((d.Narrow6<=d.Primary11)&(d.Primary11<=d.BroadSX23)).all()
d['age_status']=np.where(d.age_missing=='1','Missing','Observed')
assert (d.age_years.eq('')==d.age_missing.eq('1')).all()
M='Missing / unknown'
levels={'sex':['F','M',M],'hcp':['HCP','Non-HCP',M],'us_reporter':['US','Non-US',M],'calendar_year':[str(i) for i in range(2022,2027)],'age_status':['Observed','Missing']}
def flag(n):return 'ZERO' if n==0 else 'TINY 1–4' if n<5 else ''
rows=[]
for y in outcomes:
    for v,ls in levels.items():
        assert set(d[v])<=set(ls)
        for drug in drugs:
            x=d[d.drug==drug]
            for lev in ls:
                z=x[x[v]==lev]; p=int(z[y].sum());n=len(z)-p
                rows.append(dict(outcome=y,variable=v,drug=drug,level=lev,positive=p,negative=n,total=len(z),positive_flag=flag(p),negative_flag=flag(n),empty_stratum=int(len(z)==0)))
a=pd.DataFrame(rows)
a.to_csv(OUT/'stratum_counts.csv',index=False,encoding='utf-8-sig')
assert len(a)==144
for (y,v,drug),z in a.groupby(['outcome','variable','drug']):
    assert z.total.sum()==Ns[drug]
    assert z.positive.sum()==expected[y][drugs.index(drug)]
summary=[]
for (y,v,drug),z in a.groupby(['outcome','variable','drug'],sort=False):
    summary.append(dict(outcome=y,variable=v,drug=drug,positive_total=int(z.positive.sum()),zero_positive_strata=int(z.positive.eq(0).sum()),tiny_positive_strata=int(z.positive.between(1,4).sum()),zero_negative_strata=int(z.negative.eq(0).sum()),tiny_negative_strata=int(z.negative.between(1,4).sum()),flagged_levels='; '.join(f'{r.level}: {r.positive}/{r.negative}' for r in z.itertuples() if r.positive<5 or r.negative<5)))
pd.DataFrame(summary).to_csv(OUT/'sparsity_summary.csv',index=False,encoding='utf-8-sig')
data={'drugs':drugs,'N':Ns,'outcomes':outcomes,'levels':levels,'counts':rows,'summary':summary,'sources':sources,'expected':expected}
(WORK/'audit_data.json').write_text(json.dumps(data,ensure_ascii=False,indent=2))
# A linked analytical snapshot preserves exactly one retained record per CASEID.
d[KEY+outcomes+list(levels)+['age_years']].to_csv(OUT/'audit_case_snapshot.csv.gz',index=False,encoding='utf-8-sig',compression='gzip')
print(a.pivot(index=['variable','level'],columns=['outcome','drug'],values='positive').to_string())
print('\nZero / tiny cells:',int(a.positive.eq(0).sum()+a.negative.eq(0).sum()),int(a.positive.between(1,4).sum()+a.negative.between(1,4).sum()))
