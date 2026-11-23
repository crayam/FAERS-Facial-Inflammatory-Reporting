from pathlib import Path
import csv,gzip,json,math,hashlib
from collections import Counter
import numpy as np
R=Path(__file__).resolve().parents[2];O=R/'outputs/firth_primary11_main_v1_0'
fit=json.loads((O/'fit_results.json').read_text());qc=json.loads((O/'cohort_qc.json').read_text())
with open(O/'grouped_design_counts.csv',encoding='utf-8-sig') as f:g=list(csv.DictReader(f))
X=np.array([[float(r[c]) for c in qc['design_columns']] for r in g]);s=np.array([float(r['positive']) for r in g]);n=np.array([float(r['N']) for r in g])
def calc(b,X=X,s=s,n=n):
    eta=X@b;p=1/(1+np.exp(-eta));w=n*p*(1-p);I=X.T@(X*w[:,None]);sg,ld=np.linalg.slogdet(I)
    assert sg==1
    inv=np.linalg.inv(I);h=w*np.einsum('ij,jk,ik->i',X,inv,X)
    score=X.T@(s-n*p+h*(.5-p));ll=float(np.sum(s*eta-n*np.logaddexp(0,eta))+.5*ld)
    return ll,score,I
def solve(start,j=None,v=None):
    b=np.array(start,dtype=float);free=np.arange(9) if j is None else np.delete(np.arange(9),j)
    if j is not None:b[j]=v
    for k in range(500):
        ll,score,I=calc(b)
        if max(abs(score[free]))<1e-7:break
        step=np.linalg.solve(I[np.ix_(free,free)],score[free]);fac=1
        for q in range(50):
            new=b.copy();new[free]+=fac*step;nll,ng,_=calc(new)
            if nll>=ll-1e-10 and (nll>ll or max(abs(ng[free]))<max(abs(score[free]))):break
            fac/=2
        b=new
    assert max(abs(calc(b)[1][free]))<1e-6
    return b,k+1
beta=np.array(fit['full']['beta']);starts=[[-5]+[0]*8,[-3,.3,-.3,.2,-.2,.2,-.2,.1,-.1],[-7,-.5,.5,-.3,.3,-.3,.3,-.2,.1]]
checks=[]
for start in starts:
    b,it=solve(start);diff=float(max(abs(b-beta)));assert diff<1e-7
    checks.append(dict(type='Full model scoring validation',beta_max_difference=diff,iterations=it))
for result in fit['results']:
    j=result['beta_index']-1
    for key in ['null_fit','lower_profile_fit','upper_profile_fit']:
        z=result[key];b,it=solve(starts[0],j,z['beta'][j]);dif=float(max(abs(b-z['beta'])));ld=abs(calc(b)[0]-z['ll'])
        assert dif<1e-6 and ld<1e-8
        checks.append(dict(type=result['comparator']+' '+key,beta_max_difference=dif,ll_difference=ld,iterations=it))
    pv=math.erfc(math.sqrt(result['PLR_statistic']/2));assert abs(pv-result['P_value'])<1e-12
    assert max(abs(result['lower_LR_residual']),abs(result['upper_LR_residual']))<1e-7
# Independently reconstruct individual X and source membership from the two frozen masters.
with gzip.open(qc['source_hashes'][0]['path'],'rt',encoding='utf-8-sig') as f:chars={r['caseid']:r for r in csv.DictReader(f)}
with gzip.open(qc['source_hashes'][1]['path'],'rt',encoding='utf-8-sig') as f:outcomes={r['caseid']:r for r in csv.DictReader(f)}
with gzip.open(O/'analysis_cases.csv.gz','rt',encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
assert len({r['caseid'] for r in rows})==174341
XI=[];ys=[];groupcheck=Counter()
for r in rows:
    c=chars[r['caseid']];y=int(outcomes[r['caseid']]['Primary11'])
    assert c['us_reporter'] in ['US','Non-US']
    q=(int(c['fda_dt'][:4])-2022)*4+(int(c['fda_dt'][4:6])+2)//3-2
    x=[1,int(c['drug']=='tralokinumab'),int(c['drug']=='upadacitinib'),int(c['sex']=='M'),int(c['sex']=='Missing / unknown'),int(c['hcp']=='Non-HCP'),int(c['hcp']=='Missing / unknown'),int(c['us_reporter']=='Non-US'),q]
    assert int(r['quarter_index'])==q and int(r['Primary11'])==y
    XI.append(x);ys.append(y);groupcheck[tuple(x),y]+=1
assert {r['caseid'] for r in rows}=={cid for cid,c in chars.items() if c['us_reporter'] in ['US','Non-US']}
for i,x in enumerate(X):assert groupcheck[tuple(x),1]==s[i] and groupcheck[tuple(x),0]==n[i]-s[i]
ill,iscore,iinfo=calc(beta,np.array(XI),np.array(ys),np.ones(len(rows)))
assert abs(ill-fit['full']['ll'])<1e-8 and max(abs(iscore))<1e-6
# Numerical gradient check off the optimum validates the implemented penalty derivative.
bt=beta+np.array([.01,.02,-.02,.03,-.03,.01,-.01,.02,-.001]);grad=[]
for j in range(9):
    step=np.zeros(9);step[j]=1e-5;grad.append((calc(bt+step)[0]-calc(bt-step)[0])/2e-5)
grad_error=float(max(abs(np.array(grad)-calc(bt)[1])));assert grad_error<1e-4
crude=[]
for scope,den in [('Analysis cohort',qc['analysis_N']),('Formal cohort',qc['formal_N'])]:
    for comp in ['tralokinumab','upadacitinib']:
        a=qc['drug_positive']['dupilumab'];b=den['dupilumab']-a;c=qc['drug_positive'][comp];dd=den[comp]-c
        lor=math.log(a*dd/(b*c));se=math.sqrt(1/a+1/b+1/c+1/dd)
        crude.append(dict(scope=scope,comparator=comp,a=a,b=b,c=c,d=dd,pair_N=a+b+c+dd,crude_ROR=math.exp(lor),CI_lower=math.exp(lor-1.959963984540054*se),CI_upper=math.exp(lor+1.959963984540054*se),P_value=math.erfc(abs(lor/se)/math.sqrt(2)),zero_cell=False))
with open(O/'crude_comparisons.csv','w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,fieldnames=list(crude[0]));w.writeheader();w.writerows(crude)
v=dict(status='PASS',individual_cases_verified=len(rows),group_counts_verified=len(g),full_fit_individual_ll_difference=abs(ill-fit['full']['ll']),full_fit_individual_max_score=float(max(abs(iscore))),off_optimum_finite_difference_gradient_max_error=grad_error,independent_Fisher_scoring_checks=checks,profile_roots_verified=4,contrast_null_fits_verified=2,independent_algorithm='NumPy adjusted Fisher scoring with step halving; R fit uses BFGS plus numerical-Hessian Newton refinement',other_outcomes_fitted=False)
(O/'independent_validation.json').write_text(json.dumps(v,indent=2))
bundle=dict(qc=qc,fit=fit,validation=v,crude=crude)
(Path(__file__).parent/'workbook_data.json').write_text(json.dumps(bundle,ensure_ascii=False,indent=2))
print(json.dumps(v,indent=2));print(json.dumps(crude,indent=2))
