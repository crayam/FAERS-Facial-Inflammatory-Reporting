from pathlib import Path
import json,csv,gzip,hashlib,math
from collections import Counter
import numpy as np
from openpyxl import load_workbook
R=Path(__file__).resolve().parents[2];O=R/'outputs/final_statistical_results_v1_0';W=Path(__file__).parent
m=json.loads((O/'analysis_manifest.json').read_text())
def read(p):
 op=gzip.open if str(p).endswith('.gz') else open
 with op(p,'rt',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
chars=read(Path(m['sources'][0]['path']));outcomes={a['caseid']:a for a in read(Path(m['sources'][1]['path']))};M='Missing / unknown'
def source_design(meta):
 sec=meta['analysis_type']=='Secondary_minimal'
 rows=[a for a in chars if sec or a['us_reporter'] in ['US','Non-US']]
 if meta['analysis_type'].startswith('S2'):rows=[a for a in rows if a['age_years']!='']
 assert len({a['caseid'] for a in rows})==len(rows)==meta['N']
 def row(a):
  yr=int(a['fda_dt'][:4]);q=(yr-2022)*4+(int(a['fda_dt'][4:6])+2)//3-2
  v={'Intercept':1,'drug_tralokinumab':int(a['drug']=='tralokinumab'),'drug_upadacitinib':int(a['drug']=='upadacitinib'),'sex_Male':int(a['sex']=='M'),'sex_Missing':int(a['sex']==M),'HCP_NonHCP':int(a['hcp']=='Non-HCP'),'HCP_Missing':int(a['hcp']==M),'geo_NonUS':int(a['us_reporter']=='Non-US'),'quarter_index':q,**{f'year_{y}':int(yr==y) for y in range(2023,2027)}}
  if 'age10' in meta['columns']:v['age10']=(float(a['age_years'])-50)/10
  return [v[k] for k in meta['columns']]
 return np.array([row(a) for a in rows],dtype=float),np.array([int(outcomes[a['caseid']][meta['outcome']]) for a in rows]),rows
def calc(b,X,s,n):
 eta=X@b;p=1/(1+np.exp(-eta));w=n*p*(1-p);I=X.T@(X*w[:,None]);sg,ld=np.linalg.slogdet(I);assert sg==1
 inv=np.linalg.inv(I);h=w*np.einsum('ij,jk,ik->i',X,inv,X)
 return float(np.sum(s*eta-n*np.logaddexp(0,eta))+.5*ld),X.T@(s-n*p+h*(.5-p)),I
def solve(start,X,s,n,j=None,v=None):
 b=np.array(start,dtype=float);free=np.arange(len(b)) if j is None else np.delete(np.arange(len(b)),j)
 if j is not None:b[j]=v
 for it in range(500):
  ll,score,I=calc(b,X,s,n)
  if max(abs(score[free]))<1e-8:break
  step=np.linalg.solve(I[np.ix_(free,free)],score[free]);fac=1
  for k in range(55):
   nb=b.copy();nb[free]+=fac*step;nll,ng,_=calc(nb,X,s,n)
   if nll>=ll-1e-9 and (nll>ll or max(abs(ng[free]))<max(abs(score[free]))):break
   fac/=2
  b=nb
 assert max(abs(calc(b,X,s,n)[1][free]))<1e-6
 return b
fits={};qcs=[];est=[];validations=[]
for meta in m['new_models']:
 mid=meta['model_id'];p=O/'models'/mid;status=json.loads((p/'run_status.json').read_text())
 if status['status']!='PASS':qcs.append({**meta,'convergence':'FAILED','profile_CI':'FAILED','validation':'Not performed','warning':status['reason']});continue
 f=json.loads((p/'fit_results.json').read_text());fits[mid]=f;g=read(p/'grouped_design_counts.csv');X=np.array([[float(a[k]) for k in meta['columns']] for a in g]);s=np.array([int(a['positive']) for a in g]);n=np.array([int(a['N']) for a in g]);P=X.shape[1];b=np.array(f['full']['beta']);checks=[]
 init=np.zeros(P);init[0]=math.log(s.sum()/(n.sum()-s.sum()))
 for start in [init,init+np.linspace(-.15,.15,P)]:
  bb=solve(start,X,s,n);dif=float(max(abs(bb-b)));assert dif<1e-7;checks.append(dif)
 for a in f['results']:
  for key in ['null_fit','lower_profile_fit','upper_profile_fit']:
   z=a[key];j=a['beta_index']-1;bb=solve(init,X,s,n,j,z['beta'][j]);dif=float(max(abs(bb-z['beta'])));assert dif<1e-6 and abs(calc(bb,X,s,n)[0]-z['ll'])<1e-7;checks.append(dif)
   assert z['converged'] and z['score']<1e-6
  assert abs(math.erfc(math.sqrt(a['PLR_statistic']/2))-a['P_value'])<1e-12
  assert max(abs(a['lower_LR_residual']),abs(a['upper_LR_residual']))<1e-7
 XI,yi,rows=source_design(meta);ct=Counter((tuple(x),int(y)) for x,y in zip(XI,yi))
 for i,x in enumerate(X):assert ct[tuple(x),1]==s[i] and ct[tuple(x),0]==n[i]-s[i]
 ll,sc,_=calc(b,XI,yi,np.ones(len(yi)));assert abs(ll-f['full']['ll'])<1e-7 and max(abs(sc))<1e-6
 assert np.linalg.matrix_rank(XI)==P and min(f['information_eigenvalues'])>0 and min(f['penalized_negative_Hessian_eigenvalues'])>0
 warn=''
 if meta['separation']=='Quasi-complete':warn='Zero upadacitinib positives; quasi-complete separation; very wide profile CI; exploratory reporting pattern only.'
 elif meta['analysis_type']=='Secondary_minimal':warn='Sparse comparator positives (3–9); exploratory reporting pattern only.'
 v=dict(model_id=mid,status='PASS',source_cases_verified=len(yi),group_counts_verified=len(g),max_independent_beta_difference=max(checks),individual_ll_difference=abs(ll-f['full']['ll']),individual_max_score=float(max(abs(sc))),profile_endpoints_verified=4)
 (p/'independent_validation.json').write_text(json.dumps(v,indent=2));validations.append(v)
 qcs.append({**meta,'convergence':'Converged','profile_CI':'4/4 success','validation':'PASS','warning':warn,'max_score':f['full']['score']})
 print(mid,'independent validation PASS',flush=True)
for a in m['main_sources']:
 assert hashlib.sha256(Path(a['path']).read_bytes()).hexdigest()==a['sha256']
 y=a['outcome'];f=json.loads(Path(a['path']).read_text());mid=y+'__Main';fits[mid]=f
 qcs.append(dict(model_id=mid,outcome=y,analysis_type='Main',N=174341,positive=f['positive'],drug_N=m['geography_known_N'],drug_positive={v['drug']:v['positive'] for v in m['secondary_counts'] if v['outcome']==y},parameters=9,rank=9,separation='Absent',convergence='Converged',profile_CI='4/4 success',validation='Prior verified fit reused',warning='',max_score=f['full']['score']))
for qc in qcs:
 if qc['model_id'] not in fits:continue
 for a in fits[qc['model_id']]['results']:
  est.append(dict(model_id=qc['model_id'],outcome=qc['outcome'],analysis_type=qc['analysis_type'],comparison=a['comparison'],comparator=a['comparator'],OR=a['adjusted_OR'],lower_CI=a['CI_lower'],upper_CI=a['CI_upper'],P_value=a['P_value'],analysis_N=qc['N'],positive_cases=qc['positive'],convergence=qc['convergence'],separation=qc['separation'],warning=qc['warning']))
crude=[]
for y in ['Primary11','Narrow6','BroadSX23','Rosacea_family','Demodex_family','Blepharitis_family','Seborrhoeic_dermatitis']:
 for scope,den in [('Formal',m['formal_N'])]+([('Main analysis',m['geography_known_N'])] if y in ['Primary11','Narrow6','BroadSX23'] else []):
  counts={a['drug']:a['positive'] for a in m['secondary_counts'] if a['outcome']==y}
  for comp in ['tralokinumab','upadacitinib']:
   a=counts['dupilumab'];b=den['dupilumab']-a;c=counts[comp];dd=den[comp]-c;zero=min(a,b,c,dd)==0
   if zero:rr=lo=hi=None;pv=None;method='Not estimated: zero cell; no correction'
   else:
    lg=math.log(a*dd/(b*c));se=math.sqrt(1/a+1/b+1/c+1/dd);rr=math.exp(lg);lo=math.exp(lg-1.959963984540054*se);hi=math.exp(lg+1.959963984540054*se);pv=math.erfc(abs(lg/se)/math.sqrt(2));method='Uncorrected crude ROR; log-Wald CI'
   crude.append(dict(outcome=y,scope=scope,comparison='dupilumab vs '+comp,comparator=comp,a=a,b=b,c=c,d=dd,OR=rr,lower_CI=lo,upper_CI=hi,P_value=pv,zero_cell=zero,method=method,N=a+b+c+dd))
qcs.append(dict(model_id='Demodex_family__Descriptive',outcome='Demodex_family',analysis_type='Descriptive only',N=174406,positive=6,parameters=None,rank=None,separation='Two zero-positive comparator groups',convergence='Not fitted by design',profile_CI='Not attempted',validation='Counts verified',warning='6/0/0; descriptive reporting pattern only; no drug OR or P value estimated.'))
# Paired age-subset comparison separates sample selection from adding age10.
stability=[]
for y in ['Primary11','Narrow6','BroadSX23']:
 for comp in ['tralokinumab','upadacitinib']:
  z={a['analysis_type']:a for a in est if a['outcome']==y and a['comparator']==comp}
  main=z['Main']
  for typ in ['Main','S1_year','S2_age_observed_base','S2_age_adjusted']:
   a=z.get(typ)
   if not a:continue
   stability.append({**a,'OR_ratio_vs_main':a['OR']/main['OR'],'OR_ratio_vs_S2_base':a['OR']/z['S2_age_observed_base']['OR'] if typ=='S2_age_adjusted' and 'S2_age_observed_base' in z else None,'same_direction_as_main':(a['OR']>1)==(main['OR']>1),'CI_position':'Above 1' if a['lower_CI']>1 else 'Below 1' if a['upper_CI']<1 else 'Includes 1'})
def write(p,rows):
 fields=list(dict.fromkeys(k for a in rows for k in a))
 with open(p,'w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
write(O/'all_model_results.csv',est);write(O/'sensitivity_comparison.csv',stability);write(O/'crude_comparisons.csv',crude);write(O/'model_qc.csv',qcs)
forest=[dict(outcome=a['outcome'],comparison=a['comparison'],OR=a['OR'],lower_CI=a['lower_CI'],upper_CI=a['upper_CI'],analysis_type=a['analysis_type'],analysis_N=a['analysis_N'],P_value=a['P_value']) for a in est]
forest += [dict(outcome=a['outcome'],comparison=a['comparison'],OR=a['OR'],lower_CI=a['lower_CI'],upper_CI=a['upper_CI'],analysis_type='Crude '+a['scope'],analysis_N=a['N'],P_value=a['P_value']) for a in crude if a['OR'] is not None and a['outcome']!='Demodex_family']
write(O/'forest_plot_data.csv',forest)
tablefile=R/'outputs/table1_missingness_v1_0/FAERS_Table1_Missingness_v1.0.xlsx';wb=load_workbook(tablefile,data_only=True,read_only=True)
table1=[[c.value for c in row] for row in wb['Table1'].iter_rows(min_row=1,max_row=60,max_col=5)];wb.close()
priorcrude=read(R/'outputs/crude_comparative_v1_0/crude_comparative_results.csv')
bundle=dict(manifest=m,models=est,stability=stability,crude=crude,qc=qcs,forest=forest,validation=validations,table1=table1,missingness=read(R/'outputs/table1_missingness_v1_0/missingness.csv'),calendar=read(R/'outputs/table1_missingness_v1_0/calendar_quarter.csv'),smd=read(R/'outputs/table1_missingness_v1_0/pairwise_smd.csv'),priorcrude=priorcrude)
(W/'workbook_data.json').write_text(json.dumps(bundle,ensure_ascii=False,indent=2));print('INTEGRATION COMPLETE',len(est),'estimates',len(qcs),'QC rows',flush=True)
