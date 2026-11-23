from pathlib import Path
import csv,gzip,json,math,ast
from collections import Counter
import numpy as np
R=Path(__file__).resolve().parents[2];O=R/'outputs/firth_primary_narrow_broad_v1_0';W=Path(__file__).parent;OLD=R/'outputs/firth_primary11_main_v1_0'
m=json.loads((O/'analysis_manifest.json').read_text());spec=m['model_specification'];cols=json.loads((OLD/'cohort_qc.json').read_text())['design_columns']
def read(p):
 op=gzip.open if str(p).endswith('.gz') else open
 with op(p,'rt',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
cases=read(O/'analysis_cases_three_outcomes.csv.gz');master={a['caseid']:a for a in read(Path(spec['sources'][1]['path']))}
prior=read(OLD/'analysis_cases.csv.gz');assert len(cases)==174341
for a,b in zip(cases,prior):assert all(a[k]==v for k,v in b.items())
XI=np.array([[1,int(a['drug']=='tralokinumab'),int(a['drug']=='upadacitinib'),int(a['sex']=='Male'),int(a['sex']=='Missing'),int(a['HCP']=='Non-HCP'),int(a['HCP']=='Missing'),int(a['geography']=='Non-US'),int(a['quarter_index'])] for a in prior])
assert np.linalg.matrix_rank(XI)==9
assert all(int(a['Narrow6'])<=int(a['Primary11'])<=int(a['BroadSX23']) for a in cases)
# Reuse the independently verified Prompt-1 Fisher-scoring code unchanged.
tree=ast.parse((R/'work/firth_primary11/validate.py').read_text())
functions=ast.Module(body=[x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name in ['calc','solve']],type_ignores=[])
fitall={'Primary11':json.loads((OLD/'fit_results.json').read_text())};validation={};crude=[];results=[]
for y in ['Narrow6','BroadSX23']:
 fit=json.loads((O/y/'fit_results.json').read_text());fitall[y]=fit
 g=read(O/y/'grouped_design_counts.csv');X=np.array([[float(r[c]) for c in cols] for r in g]);s=np.array([float(r['positive']) for r in g]);n=np.array([float(r['N']) for r in g]);ns=dict(np=np,X=X,s=s,n=n)
 exec(compile(functions,'Prompt1_verified_Fisher_scoring','exec'),ns)
 calc,solve=ns['calc'],ns['solve'];beta=np.array(fit['full']['beta']);checks=[]
 for start in [[-5]+[0]*8,[-3,.3,-.3,.2,-.2,.2,-.2,.1,-.1],[-7,-.5,.5,-.3,.3,-.3,.3,-.2,.1]]:
  b,it=solve(start);dif=float(max(abs(b-beta)));assert dif<1e-7;checks.append(dict(check='Full model from independent start',beta_difference=dif,iterations=it))
 for r in fit['results']:
  j=r['beta_index']-1
  for key in ['null_fit','lower_profile_fit','upper_profile_fit']:
   z=r[key];b,it=solve([-5]+[0]*8,j,z['beta'][j]);dif=float(max(abs(b-z['beta'])));ld=abs(calc(b)[0]-z['ll']);assert dif<1e-6 and ld<1e-8
   assert z['converged'] and z['optim_code']==0 and z['score']<1e-6
   checks.append(dict(check=r['comparator']+' '+key,beta_difference=dif,ll_difference=ld,iterations=it))
  assert abs(math.erfc(math.sqrt(r['PLR_statistic']/2))-r['P_value'])<1e-12
  assert max(abs(r['lower_LR_residual']),abs(r['upper_LR_residual']))<1e-7
 ys=np.array([int(master[a['caseid']][y]) for a in cases]);assert ys.tolist()==[int(a[y]) for a in cases]
 counts=Counter((tuple(x),int(yv)) for x,yv in zip(XI,ys))
 for i,x in enumerate(X):assert counts[tuple(x),1]==s[i] and counts[tuple(x),0]==n[i]-s[i]
 ill,sc,_=calc(beta,XI,ys,np.ones(len(ys)));assert abs(ill-fit['full']['ll'])<1e-8 and max(abs(sc))<1e-6
 assert min(fit['information_eigenvalues'])>0 and min(fit['penalized_negative_Hessian_eigenvalues'])>0
 v=dict(status='PASS',outcome=y,N=len(ys),positive=int(ys.sum()),individual_ll_difference=abs(ill-fit['full']['ll']),individual_max_score=float(max(abs(sc))),checks=checks,profile_endpoints_verified=4,null_constraints_verified=2,same_design_as_Primary11=True)
 validation[y]=v;(O/y/'independent_validation.json').write_text(json.dumps(v,indent=2))
for y in ['Primary11','Narrow6','BroadSX23']:
 fit=fitall[y];audit=m['cohort_audits'][y]
 for scope,den in [('Analysis cohort',spec['analysis_N']),('Formal cohort',spec['formal_N'])]:
  for comp in ['tralokinumab','upadacitinib']:
   a=audit['positive']['dupilumab'];b=den['dupilumab']-a;c=audit['positive'][comp];dd=den[comp]-c;lor=math.log(a*dd/b/c);se=math.sqrt(1/a+1/b+1/c+1/dd)
   crude.append(dict(outcome=y,scope=scope,comparator=comp,a=a,b=b,c=c,d=dd,pair_N=a+b+c+dd,crude_ROR=math.exp(lor),CI_lower=math.exp(lor-1.959963984540054*se),CI_upper=math.exp(lor+1.959963984540054*se),zero_cell=False))
 for x in fit['results']:
  c=next(v for v in crude if v['outcome']==y and v['scope']=='Analysis cohort' and v['comparator']==x['comparator']);cf=next(v for v in crude if v['outcome']==y and v['scope']=='Formal cohort' and v['comparator']==x['comparator'])
  results.append(dict(outcome=y,comparison=x['comparison'],comparator=x['comparator'],crude_ROR=c['crude_ROR'],crude_CI_lower=c['CI_lower'],crude_CI_upper=c['CI_upper'],adjusted_OR=x['adjusted_OR'],CI_lower=x['CI_lower'],CI_upper=x['CI_upper'],P_value=x['P_value'],analysis_N=174341,convergence='Converged',separation='Absent',QC='PASS',formal_crude_ROR=cf['crude_ROR'],geography_excluded=65,source='Prompt 1 verified fit reused' if y=='Primary11' else 'New outcome; identical locked model'))
# Reconcile all three formal crude comparisons with the prior seven-outcome table.
old=read(R/'outputs/crude_comparative_v1_0/crude_comparative_results.csv')
for c in [x for x in crude if x['scope']=='Formal cohort']:
 a=next(x for x in old if x['outcome']==c['outcome'] and x['comparator']==c['comparator'])
 assert math.isclose(float(a['cROR']),c['crude_ROR'],rel_tol=1e-12)
def write(p,rs):
 with open(p,'w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rs[0]));w.writeheader();w.writerows(rs)
write(O/'combined_results.csv',results);write(O/'crude_comparisons.csv',crude)
bundle=dict(manifest=m,fits=fitall,validation=validation,results=results,crude=crude)
(W/'workbook_data.json').write_text(json.dumps(bundle,ensure_ascii=False,indent=2));print(json.dumps(results,indent=2))
