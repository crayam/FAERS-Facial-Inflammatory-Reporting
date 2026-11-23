from pathlib import Path
import json,csv,gzip,hashlib,shutil
from collections import Counter
import numpy as np
R=Path(__file__).resolve().parents[2];W=Path(__file__).parent;O=R/'outputs/restriction_sensitivity_v1_0';M='Missing / unknown'
def read(p):
 with gzip.open(p,'rt',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def write(p,rows):
 with open(p,'w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
prior=json.loads((R/'outputs/final_statistical_results_v1_0/analysis_manifest.json').read_text())
for a in prior['sources']+prior['main_sources']:assert hashlib.sha256(Path(a['path']).read_bytes()).hexdigest()==a['sha256']
chars=read(prior['sources'][0]['path']);flags={a['caseid']:a for a in read(prior['sources'][1]['path'])};assert len(chars)==len(flags)==len({a['caseid'] for a in chars})==174406
ys=['Primary11','Narrow6','BroadSX23'];drugs=['dupilumab','tralokinumab','upadacitinib']
for a in chars:
 b=flags[a['caseid']]
 for k in ['drug','primaryid','caseversion','fda_dt','source_quarter']:assert a[k]==b[k]
 a.update({y:int(b[y]) for y in ys});a['quarter_index']=(int(a['fda_dt'][:4])-2022)*4+(int(a['fda_dt'][4:6])-1)//3-1
 assert 0<=a['quarter_index']<=16 and a['sex'] in ['F','M',M] and a['hcp'] in ['HCP','Non-HCP',M] and a['us_reporter'] in ['US','Non-US',M]
assert dict(Counter(a['drug'] for a in chars))=={'dupilumab':164515,'tralokinumab':3637,'upadacitinib':6254}
meta=[];flow=[];counts=[];members=[]
for typ in ['Main','HCP_only','US_only']:
 selected=[a for a in chars if typ=='Main' or (typ=='HCP_only' and a['hcp']=='HCP') or (typ=='US_only' and a['us_reporter']=='US')]
 rows=[a for a in selected if a['us_reporter'] in ['US','Non-US']]
 if typ=='Main':assert {a['caseid'] for a in rows}=={a['caseid'] for a in read(R/'outputs/firth_primary11_main_v1_0/analysis_cases.csv.gz')}
 for drug in drugs:
  full=sum(a['drug']==drug for a in chars);pre=sum(a['drug']==drug for a in selected);n=sum(a['drug']==drug for a in rows)
  flow.append(dict(analysis_type=typ,drug=drug,formal_N=full,restriction_excluded=full-pre,restriction_selected=pre,geography_missing_excluded=pre-n,analysis_N=n))
 members.extend(dict(analysis_type=typ,drug=a['drug'],caseid=a['caseid'],primaryid=a['primaryid']) for a in rows)
 for y in ys:
  mid=y+'__'+typ;cols=['Intercept','drug_tralokinumab','drug_upadacitinib','sex_Male','sex_Missing']
  if typ!='HCP_only':cols+=['HCP_NonHCP','HCP_Missing']
  if typ!='US_only':cols+=['geo_NonUS']
  cols+=['quarter_index'];g=Counter();Ns=Counter();pos=Counter()
  for a in rows:
   v=dict(Intercept=1,drug_tralokinumab=int(a['drug']=='tralokinumab'),drug_upadacitinib=int(a['drug']=='upadacitinib'),sex_Male=int(a['sex']=='M'),sex_Missing=int(a['sex']==M),HCP_NonHCP=int(a['hcp']=='Non-HCP'),HCP_Missing=int(a['hcp']==M),geo_NonUS=int(a['us_reporter']=='Non-US'),quarter_index=a['quarter_index'])
   g[tuple(v[k] for k in cols),a[y]]+=1;Ns[a['drug']]+=1;pos[a['drug']]+=a[y]
  xs=sorted({x for x,event in g});mixed=[x for x in xs if g[x,0] and g[x,1]];rank=int(np.linalg.matrix_rank(xs));mr=int(np.linalg.matrix_rank(mixed))
  sep='Absent' if mr==len(cols) else 'Unresolved';assert rank==len(cols)
  q=dict(model_id=mid,outcome=y,analysis_type=typ,N=len(rows),positive=sum(pos.values()),drug_N=dict(Ns),drug_positive=dict(pos),columns=cols,parameters=len(cols),rank=rank,mixed_pattern_rank=mr,separation=sep,separation_proof='Full-rank mixed-outcome design patterns exclude a nonzero separating direction.' if sep=='Absent' else 'Requires further diagnostic',geography_excluded=len(selected)-len(rows),patterns=len(xs))
  meta.append(q)
  for drug in drugs:counts.append(dict(analysis_type=typ,outcome=y,drug=drug,N=Ns[drug],positive=pos[drug],negative=Ns[drug]-pos[drug],zero_cell=pos[drug]==0 or Ns[drug]==pos[drug]))
  if typ!='Main':
   p=O/'models'/mid;p.mkdir(parents=True,exist_ok=True);write(p/'grouped_design_counts.csv',[dict(zip(cols,x),positive=g[x,1],negative=g[x,0],N=g[x,1]+g[x,0]) for x in xs]);(p/'design_qc.json').write_text(json.dumps(q,indent=2))
shutil.copyfile(R/'work/final_stats/fit.R',W/'fit.R')
manifest=dict(sources=prior['sources'],main_sources=prior['main_sources'],new_models=[a for a in meta if a['analysis_type']!='Main'],main_models=[a for a in meta if a['analysis_type']=='Main'],flow=flow,counts=counts,fitter_sha256=hashlib.sha256((W/'fit.R').read_bytes()).hexdigest(),locked_design_changes=False)
(O/'analysis_manifest.json').write_text(json.dumps(manifest,indent=2));write(O/'cohort_flow.csv',flow);write(O/'cohort_counts.csv',counts)
with gzip.open(O/'analysis_membership.csv.gz','wt',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(members[0]));w.writeheader();w.writerows(members)
print(json.dumps(dict(flow=flow,models=[{k:a[k] for k in ['model_id','N','drug_positive','rank','mixed_pattern_rank','separation']} for a in meta]),indent=2))
s=(R/'work/final_stats/run_models.py').read_text().replace('outputs/final_statistical_results_v1_0','outputs/restriction_sensitivity_v1_0');(W/'run_models.py').write_text(s)
