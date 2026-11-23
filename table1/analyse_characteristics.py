"""Table 1 / missingness / SMD for the frozen full cohort; no regressions or P values."""
from pathlib import Path
from collections import Counter,defaultdict
from datetime import datetime
import csv,gzip,hashlib,json,math,statistics
import numpy as np

WORK=Path('work/table1');OUT=Path('outputs/table1_missingness_v1_0');OUT.mkdir(exist_ok=True)
DRUGS=['dupilumab','tralokinumab','upadacitinib']
N=dict(zip(DRUGS,[164515,3637,6254]));MISSING='Missing / unknown'
SCALE={'YR':1.,'MON':1/12,'WK':7/365.25,'DY':1/365.25,'HR':1/(24*365.25),'DEC':10.}
HCP={'MD','PH','HP','OT'};NONHCP={'CN','LW'};SERIOUS={'DE','LT','HO','DS','CA','RI','OT'}
LABELS={'age_years':'Age (years)','sex':'Sex','reporter_type':'Reporter type','hcp':'HCP vs non-HCP','us_reporter':'US vs non-US (reporter)','calendar_year':'Calendar year','calendar_quarter':'Calendar quarter','serious_recorded':'Recorded serious outcome'}
LEVELS={'sex':['F','M',MISSING],'reporter_type':['MD','PH','HP','OT','CN','LW',MISSING],'hcp':['HCP','Non-HCP',MISSING],'us_reporter':['US','Non-US',MISSING],'calendar_year':[str(y) for y in range(2022,2027)]+[MISSING],'calendar_quarter':[f'{y}Q{q}' for y in range(2022,2027) for q in range(1,5) if (2022,2)<=(y,q)<=(2026,2)]+[MISSING],'serious_recorded':['Yes: recorded','No serious OUTC recorded',MISSING]}
UNSPECIFIED={'','COUNTRY NOT SPECIFIED','UNKNOWN','UNK','NOT SPECIFIED','NULL','UN','ZZ'}
def save_json(name,obj): (OUT/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False))
def writecsv(name,rows):
    rows=list(rows);op=gzip.open if name.endswith('.gz') else open
    with op(OUT/name,'wt',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
def geo(code):
    c=code.strip().upper()
    if c in UNSPECIFIED:return MISSING
    if c=='US':return 'US'
    if len(c)==2 and c.isalpha():return 'Non-US'
    raise ValueError(('Unmapped country',c))

raw=list(csv.DictReader(gzip.open(WORK/'retained_characteristics_raw.csv.gz','rt',encoding='utf-8-sig')))
assert dict(Counter(r['drug'] for r in raw))==N
features=[];quality=Counter();reason=Counter()
for r in raw:
    x=dict(r);d=r['drug'];a=r['age_raw'].strip();unit=r['age_cod_raw'].strip().upper()
    if not a:
        age=None;age_reason='Blank AGE and AGE_COD' if not unit else 'Blank AGE'
    elif unit not in SCALE:age=None;age_reason='Missing or unrecognized AGE_COD'
    else:
        try:age=float(a)*SCALE[unit]
        except ValueError:age=None
        if age is None or not math.isfinite(age) or age<0:age=None;age_reason='Invalid age'
        else:age_reason='Observed'
    if age is not None and age>120:quality['age_over_120']+=1
    if age==0:quality['zero_age_retained']+=1
    # No upper-age exclusion is applied. All observed ages in this extract are 0–104 years.
    s=r['sex_raw'].strip().upper();sex=s if s in {'M','F'} else MISSING
    if s not in {'M','F','','UNK'}:quality['unmapped_sex']+=1
    o=r['occp_cod_raw'].strip().upper();reporter=o if o in HCP|NONHCP else MISSING
    if o and o not in HCP|NONHCP:quality['unmapped_reporter']+=1
    hcp='HCP' if o in HCP else 'Non-HCP' if o in NONHCP else MISSING
    country=geo(r['reporter_country_raw']);occur=geo(r['occr_country_raw'])
    dt=datetime.strptime(r['fda_dt'],'%Y%m%d');year=str(dt.year);q=f'{dt.year}Q{(dt.month-1)//3+1}'
    if q!=r['source_quarter']:quality['fda_quarter_source_disagreement']+=1
    codes=set(r['outc_raw_codes'].split(';')) if r['outc_raw_codes'] else set()
    invalid=codes-SERIOUS
    if invalid:quality['unknown_OUTC_cases']+=1
    if codes&SERIOUS:serious='Yes: recorded'
    elif int(r['outc_row_count'])==0:serious='No serious OUTC recorded'
    else:serious=MISSING
    x.update(age_years=age,age_missing=int(age is None),age_missing_reason=age_reason,sex=sex,reporter_type=reporter,hcp=hcp,us_reporter=country,us_occurrence_audit=occur,calendar_year=year,calendar_quarter=q,serious_recorded=serious,serious_recorded_flag=1 if serious=='Yes: recorded' else 0 if serious=='No serious OUTC recorded' else None,outc_no_record=int(r['outc_row_count'])==0)
    for v in LEVELS: x[v+'_missing']=int(x[v]==MISSING)
    for v in ['age_years',*LEVELS]:
        if v=='age_years':reason[(v,d,age_reason)]+=1
        elif x[v]==MISSING:
            rawval=x[{'sex':'sex_raw','reporter_type':'occp_cod_raw','hcp':'occp_cod_raw','us_reporter':'reporter_country_raw'}.get(v,'outc_raw_codes')]
            reason[(v,d,rawval or '(blank)')]+=1
    features.append(x)
assert len(features)==len({r['caseid'] for r in features})==174406
assert quality['unmapped_reporter']==quality['unmapped_sex']==quality['unknown_OUTC_cases']==0,quality
writecsv('case_characteristics.csv.gz',features)
counts={v:{d:Counter(r[v] for r in features if r['drug']==d) for d in DRUGS} for v in LEVELS}
age_stats={}
for d in DRUGS:
    vals=np.array([r['age_years'] for r in features if r['drug']==d and r['age_years'] is not None])
    q1,med,q3=np.quantile(vals,[.25,.5,.75],method='linear')
    age_stats[d]={'N':N[d],'observed_n':len(vals),'missing_n':N[d]-len(vals),'median':float(med),'q1':float(q1),'q3':float(q3),'mean':float(vals.mean()),'sd':float(vals.std(ddof=1)),'min':float(vals.min()),'max':float(vals.max())}
missing=[];cat_summary=[]
for v in LABELS:
    for d in DRUGS:
        n=age_stats[d]['missing_n'] if v=='age_years' else counts[v][d][MISSING]
        missing.append({'variable':v,'drug':d,'N':N[d],'observed_n':N[d]-n,'missing_n':n,'missing_fraction':n/N[d]})
for v,levels in LEVELS.items():
    for lev in levels:
        for d in DRUGS:
            n=counts[v][d][lev]
            cat_summary.append({'variable':v,'level':lev,'drug':d,'n':n,'N':N[d],'fraction':n/N[d]})
writecsv('missingness.csv',missing);writecsv('categorical_summary.csv',cat_summary)
writecsv('age_summary.csv',[{'drug':d,**age_stats[d]} for d in DRUGS])
writecsv('missingness_reasons.csv',[{'variable':v,'drug':d,'reason':r,'cases':n} for (v,d,r),n in reason.items()])
writecsv('calendar_year.csv',[r for r in cat_summary if r['variable']=='calendar_year'])
writecsv('calendar_quarter.csv',[r for r in cat_summary if r['variable']=='calendar_quarter'])

def binary_smd(p1,p0):
    den=math.sqrt((p1*(1-p1)+p0*(1-p0))/2)
    return (p1-p0)/den if den>0 else 0. if p1==p0 else 'Inf'
def overall_smd(c1,c0,levels):
    n1=sum(c1[x] for x in levels);n0=sum(c0[x] for x in levels)
    lev=[x for x in levels if c1[x]+c0[x]>0]
    if not n1 or not n0:return None
    if len(lev)<2:return 0.
    p=np.array([c1[x]/n1 for x in lev[:-1]]);q=np.array([c0[x]/n0 for x in lev[:-1]])
    cov=(np.diag(p)-np.outer(p,p)+np.diag(q)-np.outer(q,q))/2
    delta=p-q
    # Remove the last nonempty level so this is a nonsingular k-1 covariance in this cohort.
    return float(np.sqrt(delta@np.linalg.solve(cov,delta)))

smd=[]
def add_smd(v,metric,level,basis,comp,n1,n0,stat1,stat0,value,method):
    smd.append({'variable':v,'metric':metric,'level':level,'basis':basis,'index_drug':'dupilumab','comparator':comp,'index_n':n1,'comparator_n':n0,'index_statistic':stat1,'comparator_statistic':stat0,'signed_smd':value if method!='multinomial_Mahalanobis' else None,'absolute_smd':abs(value) if isinstance(value,(float,int)) else value,'method':method})
for comp in DRUGS[1:]:
    a,b=age_stats['dupilumab'],age_stats[comp]
    value=(a['mean']-b['mean'])/math.sqrt((a['sd']**2+b['sd']**2)/2)
    add_smd('age_years','overall','Mean age','Observed age only',comp,a['observed_n'],b['observed_n'],a['mean'],b['mean'],value,'continuous_mean_pooled_SD')
    for v,levels in LEVELS.items():
        c1,c0=counts[v]['dupilumab'],counts[v][comp]
        val=overall_smd(c1,c0,levels)
        add_smd(v,'overall','All categories','Full cohort; missing is a category',comp,N['dupilumab'],N[comp],None,None,val,'multinomial_Mahalanobis')
        # Observed-only results separate recorded distributions from unequal missingness.
        observed=[lev for lev in levels if lev!=MISSING]
        val=overall_smd(c1,c0,observed)
        add_smd(v,'observed_only','All observed categories','Observed values only',comp,N['dupilumab']-c1[MISSING],N[comp]-c0[MISSING],None,None,val,'multinomial_Mahalanobis')
        for lev in levels:
            p1,p0=c1[lev]/N['dupilumab'],c0[lev]/N[comp]
            add_smd(v,'level',lev,'Full cohort',comp,N['dupilumab'],N[comp],p1,p0,binary_smd(p1,p0),'binary_proportion_pooled_SD')
    for v in LABELS:
        m1=next(x for x in missing if x['variable']==v and x['drug']=='dupilumab')['missing_n']
        m0=next(x for x in missing if x['variable']==v and x['drug']==comp)['missing_n']
        p1,p0=m1/N['dupilumab'],m0/N[comp]
        add_smd(v,'missingness','Missing / unknown','Full cohort',comp,N['dupilumab'],N[comp],p1,p0,binary_smd(p1,p0),'binary_proportion_pooled_SD')
writecsv('pairwise_smd.csv',smd)

sets=[('Sex + year',['sex','calendar_year']),('Sex + HCP + country + year',['sex','hcp','us_reporter','calendar_year']),('Age + sex + HCP + country + year',['age_years','sex','hcp','us_reporter','calendar_year'])]
complete=[]
for label,vs in sets:
    for d in DRUGS:
        n=sum(all(r[v] is not None and r[v]!=MISSING for v in vs) for r in features if r['drug']==d)
        complete.append({'variable_set':label,'drug':d,'complete_n':n,'N':N[d],'complete_fraction':n/N[d]})
patterns=Counter((r['drug'],int(r['age_years'] is None),r['sex_missing'],r['hcp_missing'],r['us_reporter_missing']) for r in features)
pattern_rows=[{'drug':d,'age_missing':a,'sex_missing':s,'hcp_missing':h,'country_missing':c,'cases':n,'fraction':n/N[d]} for (d,a,s,h,c),n in sorted(patterns.items())]
writecsv('complete_case_availability.csv',complete);writecsv('missingness_patterns.csv',pattern_rows)

geo_audit=[];outcome_detail=[];code_audit=[]
for d in DRUGS:
    rs=[r for r in features if r['drug']==d]
    geo_audit.append({'drug':d,'reporter_US':sum(r['us_reporter']=='US' for r in rs),'reporter_missing':sum(r['us_reporter']==MISSING for r in rs),'occurrence_US':sum(r['us_occurrence_audit']=='US' for r in rs),'occurrence_missing':sum(r['us_occurrence_audit']==MISSING for r in rs),'both_observed_US_nonUS_disagreement':sum(r['us_reporter']!=r['us_occurrence_audit'] and MISSING not in {r['us_reporter'],r['us_occurrence_audit']} for r in rs),'outc_no_record':sum(r['outc_no_record'] for r in rs),'outc_recorded_serious':sum(r['serious_recorded']=='Yes: recorded' for r in rs)})
    for code in sorted(SERIOUS):
        n=sum(code in r['outc_raw_codes'].split(';') for r in rs)
        outcome_detail.append({'drug':d,'OUTC_code':code,'cases':n,'N':N[d],'fraction':n/N[d]})
    for field in ['age_cod_raw','sex_raw','occp_cod_raw','reporter_country_raw','occr_country_raw']:
        for code,n in sorted(Counter(r[field] for r in rs).items()):code_audit.append({'drug':d,'source_field':field,'raw_code':code,'cases':n})
writecsv('geography_OUTC_audit.csv',geo_audit);writecsv('OUTC_code_counts.csv',outcome_detail);writecsv('raw_code_audit.csv',code_audit)
quarter_exceptions=[{k:r[k] for k in ['drug','caseid','primaryid','caseversion','source_quarter','fda_dt','calendar_quarter']} for r in features if r['source_quarter']!=r['calendar_quarter']]
if quarter_exceptions:writecsv('calendar_quarter_source_exceptions.csv',quarter_exceptions)

extract=json.loads((WORK/'extraction_qc.json').read_text())
manifest={'analysis':'FAERS Table 1 + Missingness v1.0','input_case_flags_sha256':extract['preflight']['frozen_source_sha256'],'input_raw_extract_sha256':extract['raw_extract_sha256'],'denominators':N,'unique_CASEID':174406,'characteristics_from_retained_version_only':True,'age_units_to_years':SCALE,'age_quantile_method':'Hyndman-Fan type 7; linear interpolation / Excel PERCENTILE.INC','age_zero_retained':True,'age_upper_cutoff_applied':False,'reporter_HCP_codes':sorted(HCP),'reporter_nonHCP_codes':sorted(NONHCP),'geography_primary_source':'DEMO.REPORTER_COUNTRY; US exact code, other specified locations non-US; no OCCR fallback','calendar_source':'retained DEMO.FDA_DT','serious_definition':'At least one DE/LT/HO/DS/CA/RI/OT in retained OUTC; absent OUTC = no recorded serious outcome, not proven non-serious','SMD_continuous':'(mean_D-mean_C)/sqrt((sample_variance_D+sample_variance_C)/2)','SMD_binary':'(p_D-p_C)/sqrt((p_D*(1-p_D)+p_C*(1-p_C))/2)','SMD_multilevel':'sqrt(delta.T inv(S) delta), S average multinomial covariance, k-1 nonempty levels','SMD_missing_policy':'Primary categorical SMD includes missing category; observed-only sensitivity and missing-indicator SMD separately supplied. Age uses observed ages. No imputation.','baseline_P_values':False,'regression_run':False,'quality_counts':dict(quality),'source_manifest':extract['source_manifest'],'references':{'FAERS_definitions':'User-supplied 2026Q2 ZIP / ASCII/ASC_NTS.pdf, DEMO and OUTC sections','HP_mapping':'https://formative.jmir.org/2025/1/e71519/','SMD_continuous_binary':'https://pmc.ncbi.nlm.nih.gov/articles/PMC3472075/','SMD_multilevel':'https://support.sas.com/resources/papers/proceedings12/335-2012.pdf'}}
save_json('analysis_manifest.json',manifest)
data={'drugs':DRUGS,'N':N,'labels':LABELS,'levels':LEVELS,'age_stats':age_stats,'counts':counts,'missingness':missing,'smd':smd,'complete_case':complete,'missing_patterns':pattern_rows,'geo_audit':geo_audit,'outcome_details':outcome_detail,'quarter_exceptions':quarter_exceptions,'quality':dict(quality),'extraction_qc':extract,'manifest':manifest}
(WORK/'workbook_data.json').write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False))
print(json.dumps({'age':age_stats,'missing':missing,'overall_SMD':[x for x in smd if x['metric']=='overall'],'completeness':complete,'geography_OUTC':geo_audit,'year':counts['calendar_year']},ensure_ascii=False,indent=2))
