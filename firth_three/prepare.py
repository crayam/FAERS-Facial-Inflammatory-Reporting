from pathlib import Path
import csv,gzip,json,hashlib
from collections import Counter
from fractions import Fraction
import numpy as np
R=Path(__file__).resolve().parents[2];O=R/'outputs/firth_primary_narrow_broad_v1_0';W=Path(__file__).parent
O.mkdir(parents=True,exist_ok=True)
OLD=R/'outputs/firth_primary11_main_v1_0'
def verify(p):
 h=hashlib.sha256(p.read_bytes()).hexdigest()
 assert any(x.split()[0]==h and x.rstrip().endswith(p.name) for x in (p.parent/'checksums.sha256').read_text().splitlines()),str(p)
 return {'path':str(p),'sha256':h}
def read(p):
 op=gzip.open if str(p).endswith('.gz') else open
 with op(p,'rt',encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
sources=[verify(OLD/x) for x in ['analysis_cases.csv.gz','Statistical_Model_v1.0.json','fit_results.json','independent_validation.json','cohort_qc.json']]
spec=json.loads((OLD/'Statistical_Model_v1.0.json').read_text());qc=json.loads((OLD/'cohort_qc.json').read_text())
for a in spec['sources']:assert verify(Path(a['path']))['sha256']==a['sha256']
master={a['caseid']:a for a in read(Path(spec['sources'][1]['path']))}
cases=read(OLD/'analysis_cases.csv.gz');assert len(cases)==len({a['caseid'] for a in cases})==174341
cols=qc['design_columns'];groups=Counter();expected={'Primary11':[1068,19,21],'Narrow6':[754,13,14],'BroadSX23':[2347,51,30]};drugs=spec['drug_levels'];Ns=Counter()
for a in cases:
 b=master[a['caseid']]
 for k in ['drug','primaryid','caseversion','source_quarter','fda_dt']:assert a[k]==b[k]
 assert a['Primary11']==b['Primary11'];a.update(Narrow6=b['Narrow6'],BroadSX23=b['BroadSX23']);Ns[a['drug']]+=1
 x=(1,int(a['drug']=='tralokinumab'),int(a['drug']=='upadacitinib'),int(a['sex']=='Male'),int(a['sex']=='Missing'),int(a['HCP']=='Non-HCP'),int(a['HCP']=='Missing'),int(a['geography']=='Non-US'),int(a['quarter_index']))
 assert a['geography'] in ['US','Non-US']
 for y in expected:groups[y,x,int(a[y])]+=1
assert dict(Ns)==spec['analysis_N']
def write(p,rs):
 op=gzip.open if str(p).endswith('.gz') else open
 with op(p,'wt',newline='',encoding='utf-8-sig') as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]));w.writeheader();w.writerows(rs)
write(O/'analysis_cases_three_outcomes.csv.gz',cases)
prior=read(OLD/'grouped_design_counts.csv');xs=[tuple(int(a[c]) for c in cols) for a in prior]
audits={}
for y in expected:
 out=O/y
 if y!='Primary11':out.mkdir(exist_ok=True)
 pos=[sum(int(a[y]) for a in cases if a['drug']==drug) for drug in drugs];assert pos==expected[y]
 # No positive case is excluded for any of the three phenotypes.
 assert pos==[sum(int(a[y]) for a in master.values() if a['drug']==drug) for drug in drugs]
 agg=[dict(zip(cols,x),positive=groups[y,x,1],negative=groups[y,x,0],N=groups[y,x,1]+groups[y,x,0]) for x in xs]
 assert [a['N'] for a in agg]==[int(a['N']) for a in prior]
 mixed=[x for x in xs if groups[y,x,1]>0 and groups[y,x,0]>0];rank=int(np.linalg.matrix_rank(mixed));assert rank==9
 basis=[]
 for x in mixed:
  if np.linalg.matrix_rank(basis+[x])>len(basis):basis.append(x)
  if len(basis)==9:break
 B=[[Fraction(z) for z in x] for x in basis];det=Fraction(1)
 for i in range(9):
  j=next(j for j in range(i,9) if B[j][i])
  if j!=i:B[i],B[j]=B[j],B[i];det=-det
  p=B[i][i];det*=p
  for k in range(i+1,9):
   f=B[k][i]/p
   for j in range(i,9):B[k][j]-=f*B[i][j]
 assert det!=0
 audit=dict(outcome=y,analysis_N=dict(Ns),total_N=174341,positive=dict(zip(drugs,pos)),total_positive=sum(pos),geography_excluded=spec['excluded'],excluded_positive=0,design_rank=9,mixed_pattern_rank=rank,mixed_patterns=len(mixed),mixed_basis_determinant=str(det),separation=False,design_identical_to_Primary11=True)
 audits[y]=audit
 if y!='Primary11':
  write(out/'grouped_design_counts.csv',agg);write(out/'no_separation_certificate.csv',[dict(zip(cols,x),positive=groups[y,x,1],negative=groups[y,x,0]) for x in basis]);(out/'cohort_qc.json').write_text(json.dumps(audit,indent=2))
manifest=dict(model_specification=spec,prior_sources=sources,cohort_audits=audits,new_outcomes_fitted=['Narrow6','BroadSX23'],Primary11_reused_without_refit=True,fit_implementation_sha256=hashlib.sha256((R/'work/firth_primary11/fit.R').read_bytes()).hexdigest(),covariate_or_definition_changes=False)
(O/'analysis_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2));print(json.dumps(audits,indent=2))
