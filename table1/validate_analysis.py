"""Independent raw-file linkage, transformations, counts, quantiles and SMD verification."""
import csv,gzip,io,json,math,statistics,zipfile
from pathlib import Path
from collections import Counter,defaultdict
from datetime import datetime
import numpy as np
BASE=Path('outputs/table1_missingness_v1_0');WORK=Path('work/table1')
data=json.loads((WORK/'workbook_data.json').read_text());N=data['N'];missing='Missing / unknown'
features=list(csv.DictReader(gzip.open(BASE/'case_characteristics.csv.gz','rt',encoding='utf-8-sig')))
idx={(r['source_quarter'],r['primaryid'],r['caseid']):r for r in features}
assert len(features)==len(idx)==len({r['caseid'] for r in features})==174406
assert dict(Counter(r['drug'] for r in features))==N
demo_seen=set();outc=defaultdict(list);raw_rows=Counter()
for item in data['extraction_qc']['source_manifest']:
    q=item['quarter']
    with zipfile.ZipFile(Path('inputs/faers_ascii')/item['filename']) as z:
        for prefix in ['DEMO','OUTC']:
            name=next(n for n in z.namelist() if n.upper().endswith(f'{prefix}{q[2:]}.TXT'))
            with z.open(name) as b:
                reader=csv.DictReader(io.TextIOWrapper(b,encoding='utf-8-sig',newline=''),delimiter='$',quoting=csv.QUOTE_NONE)
                for raw in reader:
                    raw_rows[prefix]+=1
                    key=(q,raw['primaryid'],raw['caseid'])
                    if key not in idx:continue
                    r=idx[key]
                    if prefix=='DEMO':
                        assert key not in demo_seen;demo_seen.add(key)
                        for f in ['caseversion','fda_dt','init_fda_dt']:assert r[f]==raw[f]
                        for f in ['age','age_cod','age_grp','sex','occp_cod','reporter_country','occr_country']:assert r[f+'_raw']==raw[f]
                    else:outc[key].append(raw['outc_cod'])
assert demo_seen==set(idx)
for key,r in idx.items():
    assert r['outc_raw_codes']==';'.join(sorted(set(outc[key])))
    assert int(r['outc_row_count'])==len(outc[key])
    assert int(r['serious_recorded_flag'])==bool(set(outc[key])&{'DE','LT','HO','DS','CA','RI','OT'})
    if r['age_raw']:
        value=float(r['age_raw']);unit=r['age_cod_raw']
        if unit=='DEC':expected=value*10
        elif unit=='YR':expected=value
        elif unit=='MON':expected=value/12
        elif unit=='WK':expected=value/365.25*7
        elif unit=='DY':expected=value/365.25
        elif unit=='HR':expected=value/365.25/24
        else:raise ValueError(unit)
        assert math.isclose(float(r['age_years']),expected,abs_tol=1e-12)
    else:assert r['age_years']=='' and r['age_missing']=='1'
    assert r['sex']==(r['sex_raw'] if r['sex_raw'] in ['F','M'] else missing)
    occ=r['occp_cod_raw']
    assert r['reporter_type']==(occ if occ else missing)
    assert r['hcp']==('HCP' if occ in ['MD','PH','HP','OT'] else 'Non-HCP' if occ in ['CN','LW'] else missing)
    assert r['us_reporter']==('US' if r['reporter_country_raw']=='US' else missing if r['reporter_country_raw'] in ['', 'COUNTRY NOT SPECIFIED'] else 'Non-US')
    dt=datetime.strptime(r['fda_dt'],'%Y%m%d')
    assert r['calendar_year']==str(dt.year) and r['calendar_quarter']==f'{dt.year}Q{(dt.month+2)//3}'

ages={d:sorted(float(r['age_years']) for r in features if r['drug']==d and r['age_years']!='') for d in N}
def quantile(a,p):
    loc=(len(a)-1)*p;lo=int(loc);frac=loc-lo
    return a[lo]*(1-frac)+a[min(lo+1,len(a)-1)]*frac
for d,a in ages.items():
    s=data['age_stats'][d]
    for k,e in [('median',statistics.median(a)),('q1',quantile(a,.25)),('q3',quantile(a,.75)),('mean',statistics.mean(a)),('sd',statistics.stdev(a))]:assert math.isclose(s[k],e,abs_tol=1e-10)
    assert s['observed_n']==len(a) and s['missing_n']==N[d]-len(a)
cnt={v:{d:Counter(r[v] for r in features if r['drug']==d) for d in N} for v in data['levels']}
for v,bydrug in cnt.items():
    for d,c in bydrug.items():
        assert sum(c.values())==N[d]
        assert c==Counter(data['counts'][v][d])
for m in data['missingness']:
    v,d=m['variable'],m['drug'];n=N[d]-len(ages[d]) if v=='age_years' else cnt[v][d][missing]
    assert m['missing_n']==n and m['observed_n']+n==N[d]
    assert math.isclose(m['missing_fraction'],n/N[d],abs_tol=1e-15)

for row in data['smd']:
    v,d,basis,metric=row['variable'],row['comparator'],row['basis'],row['metric']
    if row['method']=='continuous_mean_pooled_SD':
        e=(statistics.mean(ages['dupilumab'])-statistics.mean(ages[d]))/math.sqrt((statistics.variance(ages['dupilumab'])+statistics.variance(ages[d]))/2)
    elif row['method']=='multinomial_Mahalanobis':
        levels=[x for x in data['levels'][v] if (metric!='observed_only' or x!=missing)]
        levels=[x for x in levels if cnt[v]['dupilumab'][x]+cnt[v][d][x]>0]
        p=np.array([cnt[v]['dupilumab'][x] for x in levels],dtype=float);p/=p.sum()
        q=np.array([cnt[v][d][x] for x in levels],dtype=float);q/=q.sum()
        S=(np.diag(p)-np.outer(p,p)+np.diag(q)-np.outer(q,q))/2
        delta=p-q
        # Full-k singular covariance with Moore-Penrose inverse, independent of dropped reference.
        e=math.sqrt(max(0,float(delta@np.linalg.pinv(S,hermitian=True)@delta)))
    else:
        if metric=='missingness':
            c1=N['dupilumab']-len(ages['dupilumab']) if v=='age_years' else cnt[v]['dupilumab'][missing]
            c0=N[d]-len(ages[d]) if v=='age_years' else cnt[v][d][missing]
        else:c1=cnt[v]['dupilumab'][row['level']];c0=cnt[v][d][row['level']]
        p,q=c1/N['dupilumab'],c0/N[d];variance=(p-p*p+q-q*q)/2
        e=(p-q)/math.sqrt(variance) if variance>0 else 0.
    assert math.isclose(abs(e),row['absolute_smd'],abs_tol=1e-10,rel_tol=1e-10),(v,d,metric,e,row)
    if row['signed_smd'] is not None:assert math.isclose(e,row['signed_smd'],abs_tol=1e-10)
report={'full_cohort_denominators':N,'rows_verified':len(features),'unique_CASEID':len(idx),'raw_DEMO_rows_rescanned':raw_rows['DEMO'],'raw_OUTC_rows_rescanned':raw_rows['OUTC'],'all_retained_demographics_equal_raw_source':True,'all_OUTC_code_sets_and_multiplicities_equal_raw_source':True,'all_transformed_values_verified':True,'all_categorical_counts_and_missingness_verified':True,'age_quantiles_mean_SD_independently_verified':True,'SMD_rows_independently_verified':len(data['smd']),'multinomial_SMD_matches_full_matrix_pseudoinverse':True,'regressions_run':0,'baseline_P_values_computed':0}
(BASE/'independent_validation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
print(json.dumps(report,ensure_ascii=False,indent=2))
