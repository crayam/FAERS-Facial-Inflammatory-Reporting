"""Frozen FAERS case-level cROR: raw sample OR, Wald CI or zero-cell exact CI."""
from pathlib import Path
from collections import Counter
from statistics import NormalDist
import csv,gzip,hashlib,json,math,platform

SRC=Path('outputs/all_outcomes_v1_0');OUT=Path('outputs/crude_comparative_v1_0');OUT.mkdir(parents=True,exist_ok=True)
OUTCOMES=['Primary11','Narrow6','BroadSX23','Rosacea_family','Demodex_family','Blepharitis_family','Seborrhoeic_dermatitis']
COMPS=['tralokinumab','upadacitinib'];Z=NormalDist().inv_cdf(.975)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(name):
    op=gzip.open if name.endswith('.gz') else open
    with op(SRC/name,'rt',encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def dump(name,x):(OUT/name).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False))
def write(name,rows):
    with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def logsumexp(v):
    m=max(v);return m+math.log(sum(math.exp(x-m) for x in v))
def exact_conditional_ci(a,b,c,d,alpha=.05):
    """Central exact Fisher noncentral-hypergeometric confidence limits.

    X = outcome-positive cases in index arm, conditioned on both margins.
    w(x;theta) = choose(a+b,x)*choose(c+d,a+c-x)*theta**x.
    Lower inverts Pr_theta(X >= a)=alpha/2, upper Pr_theta(X <= a)=alpha/2.
    Boundary upper is +Inf if the observation is at the maximum support.
    """
    n1,n0,k=a+b,c+d,a+c
    lo,hi=max(0,k-n0),min(n1,k)
    support=list(range(lo,hi+1))
    def lc(n,x):return math.lgamma(n+1)-math.lgamma(x+1)-math.lgamma(n-x+1)
    base=[lc(n1,x)+lc(n0,k-x) for x in support]
    def logtail(t,upper):
        w=[q+x*t for q,x in zip(base,support)]
        sel=[q for q,x in zip(w,support) if (x>=a if upper else x<=a)]
        return logsumexp(sel)-logsumexp(w)
    target=math.log(alpha/2)
    def root(upper):
        left,right=-64.,64.
        for _ in range(160):
            mid=(left+right)/2;v=logtail(mid,upper)
            if (v<target if upper else v>target):left=mid
            else:right=mid
        root=(left+right)/2
        assert abs(logtail(root,upper)-target)<1e-8
        return math.exp(root)
    return (0. if a==lo else root(True), math.inf if a==hi else root(False))

checks={l.split('  ',1)[1]:l.split('  ',1)[0] for l in (SRC/'checksums.sha256').read_text().splitlines()}
source_names=['case_outcome_flags.csv.gz','outcome_summary.csv','phenotype_definition.json','analysis_manifest.json','independent_validation.json']
source_hash={n:sha(SRC/n) for n in source_names}
assert all(source_hash[n]==checks[n] for n in source_names)
flags=read('case_outcome_flags.csv.gz');assert len(flags)==174406
bydrug={d:[r for r in flags if r['drug']==d] for d in ['dupilumab']+COMPS}
denominators={d:len(rr) for d,rr in bydrug.items()};assert denominators=={'dupilumab':164515,'tralokinumab':3637,'upadacitinib':6254}
assert len({(r['drug'],r['caseid']) for r in flags})==len(flags)
for comp in COMPS:assert not({r['caseid'] for r in bydrug['dupilumab']}&{r['caseid'] for r in bydrug[comp]})
counts={(d,o):sum(int(r[o]) for r in rr) for d,rr in bydrug.items() for o in OUTCOMES}
for r in read('outcome_summary.csv'):
    assert counts[(r['drug'],r['outcome'])]==int(r['positive_n']) and denominators[r['drug']]==int(r['denominator_N'])
results=[];cells=[]
for comp in COMPS:
    for o in OUTCOMES:
        a=counts['dupilumab',o];b=denominators['dupilumab']-a;c=counts[comp,o];d=denominators[comp]-c
        table=[a,b,c,d];zeros=[k for k,v in zip('abcd',table) if v==0]
        assert all(isinstance(x,int) and x>=0 for x in table)
        if zeros:
            assert zeros==['c'] and a>0 and b>0 and d>0
            estimate=math.inf;lower,upper=exact_conditional_ci(a,b,c,d);se=None;method='conditional_exact_95_CI';status='infinite_sample_cROR_due_to_c_zero'
        else:
            estimate=a*d/(b*c);se=math.sqrt(sum(1/x for x in table));lower=math.exp(math.log(estimate)-Z*se);upper=math.exp(math.log(estimate)+Z*se);method='log_Wald_95_CI';status='finite_sample_cROR'
        def serial(x):return 'Inf' if x==math.inf else x
        r={'analysis_version':'Crude comparative analysis v1.0','outcome':o,'index_drug':'dupilumab','comparator':comp,'comparator_role':'primary active comparator' if comp=='tralokinumab' else 'secondary active comparator','a_index_positive':a,'b_index_negative':b,'c_comparator_positive':c,'d_comparator_negative':d,'N_index':a+b,'N_comparator':c+d,'zero_cell':int(bool(zeros)),'zero_cell_names':';'.join(zeros),'zero_cell_count':len(zeros),'continuity_correction':0,'cROR':serial(estimate),'CI95_lower':serial(lower),'CI95_upper':serial(upper),'CI_method':method,'SE_log_cROR':se,'estimate_status':status,'CI_confidence_level':.95}
        results.append(r)
        for drug,pos,neg in [('dupilumab',a,b),(comp,c,d)]:
            for flag,count in [(1,pos),(0,neg)]:cells.append({'outcome':o,'index_drug':'dupilumab','comparator':comp,'row_drug':drug,'outcome_flag':flag,'cell_label':('a' if flag else 'b') if drug=='dupilumab' else ('c' if flag else 'd'),'case_count':count})
assert len(results)==14 and sum(r['zero_cell'] for r in results)==3
write('crude_comparative_results.csv',results);write('case_level_2x2_cells.csv',cells)
write('zero_cell_audit.csv',[r for r in results if r['zero_cell']])
display=[]
for o in OUTCOMES:
    r={'outcome':o}
    for comp in COMPS:
        x=next(t for t in results if t['outcome']==o and t['comparator']==comp)
        def fmt(v):return '∞' if v=='Inf' else f'{v:.3f}'
        r[comp+'_cROR_95CI']=f"{fmt(x['cROR'])} ({fmt(x['CI95_lower'])}–{fmt(x['CI95_upper'])})"
        r[comp+'_zero_cell']=x['zero_cell']
    display.append(r)
write('comparison_summary.csv',display)
dump('methodology.json',{'version':'Crude comparative analysis v1.0','count_unit':'Unique CASEID within drug; fixed formal cohort; same seven prespecified binary outcomes','table_layout':{'rows':['dupilumab','comparator'],'columns':['outcome=1','outcome=0'],'cells':[['a','b'],['c','d']]},'orientation':'dupilumab relative to comparator','point_estimate':'Uncorrected sample cROR = (a/b)/(c/d) = a*d/(b*c)','nonzero_CI':'exp(log(cROR) +/- z_0.975 * sqrt(1/a+1/b+1/c+1/d))','z_0.975':Z,'zero_cell_policy':'Keep original integer counts. For observed c=0 and a,b,d>0: uncorrected sample cROR=Inf. Use two-sided central conditional exact 95% confidence limits; upper bound=Inf. No continuity correction. Log-Wald CI/SE is unavailable for these rows.','exact_algorithm':'Invert Fisher noncentral hypergeometric tail probabilities at 0.025 with fixed row and column margins; log-space bisection.','mixed_CI_methods_explicitly_labeled':True,'confidence_level':.95,'continuity_correction':0,'CI_scope':'Individual nominal 95% intervals; no simultaneous or multiplicity adjustment. Outcomes are correlated.','covariate_adjustment':'None; crude reporting odds ratio, not a causal or incidence risk ratio.','references':[{'title':'SciPy OddsRatioResult.confidence_interval','url':'https://docs.scipy.org/doc/scipy-1.15.0/reference/generated/scipy.stats._result_classes.OddsRatioResult.confidence_interval.html'},{'title':'SciPy contingency.odds_ratio','url':'https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.contingency.odds_ratio.html'}]})
dump('analysis_manifest.json',{'version':'Crude comparative analysis v1.0','source_sha256':source_hash,'source_directory_label':SRC.name,'cohort_rows':len(flags),'denominators':denominators,'outcomes':OUTCOMES,'comparators':COMPS,'comparisons':len(results),'source_window':'2022Q2–2026Q2; retained FDA_DT 20220401–20260630','cohort_unchanged':True,'case_level_2x2_counts_rebuilt_from_flags':True,'counts_equal_frozen_outcome_summary':True,'pairwise_cohort_CASEID_overlap':0,'zero_cell_rows':sum(r['zero_cell'] for r in results),'zero_cell_comparisons':[{'outcome':r['outcome'],'comparator':r['comparator'],'zero_cell_names':r['zero_cell_names']} for r in results if r['zero_cell']],'python_version':platform.python_version(),'dependencies':'Python standard library only; exact interval implemented explicitly in analysis script.'})
Path('work/crude_workbook_data.json').write_text(json.dumps({'results':results,'display':display,'denominators':denominators},ensure_ascii=False,indent=2,allow_nan=False))
print(json.dumps(display,ensure_ascii=False,indent=2));print('All 14 tables and estimates generated; zero-cell rows: 3.')
