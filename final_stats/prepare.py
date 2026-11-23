from pathlib import Path
import json,csv,gzip,hashlib
from collections import Counter
import numpy as np
R=Path(__file__).resolve().parents[2];W=Path(__file__).parent;O=R/'outputs/final_statistical_results_v1_0';O.mkdir(parents=True,exist_ok=True)
def read(p):
 op=gzip.open if str(p).endswith('.gz') else open
 with op(p,'rt',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def verify(p):
 h=hashlib.sha256(p.read_bytes()).hexdigest();assert any(x.split()[0]==h and x.rstrip().endswith(p.name) for x in (p.parent/'checksums.sha256').read_text().splitlines()),p
 return dict(path=str(p),sha256=h)
cp=R/'outputs/table1_missingness_v1_0/case_characteristics.csv.gz';op=R/'outputs/all_outcomes_v1_0/case_outcome_flags.csv.gz';sources=[verify(cp),verify(op)]
c=read(cp);flags={x['caseid']:x for x in read(op)};assert len(c)==len(flags)==174406
for a in c:
 b=flags[a['caseid']]
 for key in ['drug','primaryid','caseversion','fda_dt','source_quarter']:assert a[key]==b[key]
 a.update({y:int(b[y]) for y in ['Primary11','Narrow6','BroadSX23','Rosacea_family','Demodex_family','Blepharitis_family','Seborrhoeic_dermatitis']})
 a['quarter_index']=(int(a['fda_dt'][:4])-2022)*4+(int(a['fda_dt'][4:6])-1)//3-1
drugs=['dupilumab','tralokinumab','upadacitinib'];outcomes=['Primary11','Narrow6','BroadSX23'];secondary=['Rosacea_family','Demodex_family','Blepharitis_family','Seborrhoeic_dermatitis'];M='Missing / unknown'
geo=[a for a in c if a['us_reporter'] in ['US','Non-US']];age=[a for a in geo if a['age_years']!=''];assert len(geo)==174341
old=read(R/'outputs/firth_primary11_main_v1_0/analysis_cases.csv.gz');assert {a['caseid'] for a in geo}=={a['caseid'] for a in old}
def write(p,rows):
 with open(p,'w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
models=[]
for y in outcomes:
 for typ,rows in [('S1_year',geo),('S2_age_observed_base',age),('S2_age_adjusted',age)]:models.append((y,typ,rows))
# All non-Demodex secondary outcomes are exploratory and use the same prechosen minimal design.
for y in secondary:
 if y!='Demodex_family':models.append((y,'Secondary_minimal',c))
meta=[]
for y,typ,rows in models:
 modelid=y+'__'+typ;p=O/'models'/modelid;p.mkdir(parents=True,exist_ok=True)
 cols=['Intercept','drug_tralokinumab','drug_upadacitinib']
 if typ!='Secondary_minimal':cols+=['sex_Male','sex_Missing','HCP_NonHCP','HCP_Missing','geo_NonUS']
 cols+=['year_2023','year_2024','year_2025','year_2026'] if typ=='S1_year' else ['quarter_index']
 if typ=='S2_age_adjusted':cols+=['age10']
 g=Counter();Ns=Counter();pos=Counter()
 for a in rows:
  x=[1,int(a['drug']=='tralokinumab'),int(a['drug']=='upadacitinib')]
  if typ!='Secondary_minimal':x += [int(a['sex']=='M'),int(a['sex']==M),int(a['hcp']=='Non-HCP'),int(a['hcp']==M),int(a['us_reporter']=='Non-US')]
  x+= [int(a['calendar_year']==str(yr)) for yr in range(2023,2027)] if typ=='S1_year' else [a['quarter_index']]
  if typ=='S2_age_adjusted':x += [(float(a['age_years'])-50)/10]
  g[tuple(x),a[y]]+=1;Ns[a['drug']]+=1;pos[a['drug']]+=a[y]
 xs=sorted(set(x for x,v in g));X=np.array(xs);rank=int(np.linalg.matrix_rank(X));mixed=[x for x in xs if g[x,0] and g[x,1]];mr=int(np.linalg.matrix_rank(mixed))
 if mr==len(cols):sep='Absent';proof='Full-rank mixed-outcome patterns exclude any nonzero separating direction.'
 elif y=='Blepharitis_family' and typ=='Secondary_minimal':
  v=np.zeros(len(cols));v[2]=-1
  assert all((np.dot(x,v)>=0 if event else np.dot(x,v)<=0) for (x,event) in g)
  sep='Quasi-complete';proof='All upadacitinib cases negative; direction beta_upadacitinib=-1, others 0 is a nonzero weak separator. Mixed outcomes exclude complete separation.'
 else:sep='Unresolved';proof='Mixed-pattern rank alone insufficient; diagnostic pending.'
 assert rank==len(cols)
 write(p/'grouped_design_counts.csv',[dict(zip(cols,x),positive=g[x,1],negative=g[x,0],N=g[x,1]+g[x,0]) for x in xs])
 a=dict(model_id=modelid,outcome=y,analysis_type=typ,N=len(rows),positive=sum(pos.values()),drug_N=dict(Ns),drug_positive=dict(pos),columns=cols,parameters=len(cols),rank=rank,mixed_pattern_rank=mr,separation=sep,separation_proof=proof,patterns=len(xs),geography_excluded=0 if typ=='Secondary_minimal' else 65,age_excluded=len(geo)-len(age) if typ.startswith('S2') else 0)
 (p/'design_qc.json').write_text(json.dumps(a,indent=2));meta.append(a)
main_sources=[]
for y in outcomes:
 p=R/'outputs/firth_primary11_main_v1_0/fit_results.json' if y=='Primary11' else R/f'outputs/firth_primary_narrow_broad_v1_0/{y}/fit_results.json'
 h=hashlib.sha256(p.read_bytes()).hexdigest();main_sources.append(dict(outcome=y,path=str(p),sha256=h))
 counts=[dict(outcome=y,drug=drug,N=sum(a['drug']==drug for a in c),positive=sum(a[y] for a in c if a['drug']==drug)) for y in outcomes+secondary for drug in drugs]
manifest=dict(sources=sources,main_sources=main_sources,new_models=meta,secondary_counts=counts,secondary_plan='Formal cohort descriptive counts and crude ROR when all cells nonzero. Rosacea, Blepharitis and Seborrhoeic: exploratory three-drug Firth + continuous quarter only, no HCP addition. Demodex 6/0/0 descriptive only, no model or continuity-corrected inferential claim.',age_subset_N=dict(Counter(a['drug'] for a in age)),formal_N=dict(Counter(a['drug'] for a in c)),geography_known_N=dict(Counter(a['drug'] for a in geo)),locked_design_changes=False)
(O/'analysis_manifest.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps({'age_subset':manifest['age_subset_N'],'models':meta},indent=2))
# Generalize only the design-column count in the already verified fitter. Optimization and inference are unchanged.
src=(R/'work/firth_primary11/fit.R').read_text();src=src.replace("X <- as.matrix(a[,1:9]); y <- a$positive; n <- a$N\ncolnames(X) <- c('Intercept','drug_tralokinumab','drug_upadacitinib','sex_Male','sex_Missing','HCP_NonHCP','HCP_Missing','geo_NonUS','quarter_index')", "X <- as.matrix(a[,setdiff(names(a),c('positive','negative','N')),drop=FALSE]); y <- a$positive; n <- a$N")
src=src.replace('9-column','full-design-column');(W/'fit.R').write_text(src)
