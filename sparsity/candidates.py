from analyse import d, OUT, WORK, outcomes, drugs, levels, M, data
import numpy as np
import pandas as pd
import json

# Assess design matrices and observed joint strata only; no outcome model is fitted.
models={
 'Main':dict(cats=['sex','hcp'],year='linear',age=False),
 'S1_geography':dict(cats=['sex','hcp','us_reporter'],year='linear',age=False),
 'S2_year_factor':dict(cats=['sex','hcp'],year='factor',age=False),
 'S3_age_observed_base':dict(cats=['sex','hcp'],year='linear',age=False,subset_age=True),
 'S3_age_observed_adjusted':dict(cats=['sex','hcp'],year='linear',age=True,subset_age=True),
}
design=[]; joint=[]; retention=[]
for comp in drugs[1:]:
    for name,spec in models.items():
        z=d[d.drug.isin(['dupilumab',comp])].copy()
        if spec.get('subset_age'):z=z[z.age_status=='Observed']
        cols=[np.ones(len(z)),z.drug.eq('dupilumab').to_numpy(dtype=float)]
        for v in spec['cats']:
            cols.extend([z[v].eq(lev).to_numpy(dtype=float) for lev in levels[v][1:]])
        if spec['year']=='linear':cols.append(z.calendar_year.astype(float).to_numpy()-2024)
        else:cols.extend([z.calendar_year.eq(yr).to_numpy(dtype=float) for yr in levels['calendar_year'][1:]])
        if spec['age']:cols.append((z.age_years.astype(float).to_numpy()-50)/10)
        X=np.column_stack(cols);rank=int(np.linalg.matrix_rank(X))
        assert rank==X.shape[1]
        design.append(dict(comparator=comp,model=name,N=len(z),slope_df=X.shape[1]-1,total_parameters=X.shape[1],rank=rank,full_rank=True,regression_fitted=False))
        # Group exact categorical + year combinations, not continuous age values.
        keys=spec['cats']+['calendar_year']
        combos=list(z.groupby(keys,sort=True,dropna=False))
        for y in outcomes:
            overlap_cases=overlap_events=0;both_support=0;zero_by_drug=0;any_zero=0;tiny_by_drug=0
            for combo,g in combos:
                support=set(g.drug)=={'dupilumab',comp}
                if support:both_support+=1;overlap_cases+=len(g);overlap_events+=int(g[y].sum())
                zp=False
                for drug in ['dupilumab',comp]:
                    h=g[g.drug==drug];p=int(h[y].sum());n=len(h)-p
                    if len(h):
                        zero_by_drug+=int(p==0);tiny_by_drug+=int(1<=p<=4);zp=zp or p==0
                    joint.append(dict(comparator=comp,model=name,outcome=y,pattern='; '.join(f'{k}={v}' for k,v in zip(keys,combo)),drug=drug,positive=p,negative=n,total=len(h),empty_stratum=int(len(h)==0),both_drugs_present=int(support)))
                any_zero+=int(zp)
            for drug in ['dupilumab',comp]:
                h=z[z.drug==drug]
                retention.append(dict(comparator=comp,model=name,outcome=y,drug=drug,N=len(h),positive=int(h[y].sum()),negative=len(h)-int(h[y].sum())))
            design[-1][f'{y}_joint_patterns']=len(combos)
            design[-1][f'{y}_both_drugs_patterns']=both_support
            design[-1][f'{y}_overlap_cases']=overlap_cases
            design[-1][f'{y}_overlap_events']=overlap_events
            design[-1][f'{y}_zero_positive_nonempty_drug_strata']=zero_by_drug
            design[-1][f'{y}_tiny_positive_nonempty_drug_strata']=tiny_by_drug
pd.DataFrame(design).to_csv(OUT/'candidate_design_audit.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(joint).to_csv(OUT/'candidate_joint_counts.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(retention).to_csv(OUT/'candidate_case_retention.csv',index=False,encoding='utf-8-sig')
data.update(design=design,retention=retention,joint=joint)
(WORK/'audit_data.json').write_text(json.dumps(data,ensure_ascii=False,indent=2))
print(pd.DataFrame(design).drop(columns=[k for k in design[0] if k.startswith(('Narrow6_','BroadSX23_'))]).to_string(index=False))
print(pd.DataFrame(retention).query("model=='S3_age_observed_adjusted'").to_string(index=False))
