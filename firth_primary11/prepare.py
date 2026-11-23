from pathlib import Path
import csv,gzip,json,hashlib
from collections import Counter
import numpy as np
R=Path(__file__).resolve().parents[2];O=R/'outputs/firth_primary11_main_v1_0';W=Path(__file__).parent
O.mkdir(parents=True,exist_ok=True)
sources=[R/'outputs/table1_missingness_v1_0/case_characteristics.csv.gz',R/'outputs/all_outcomes_v1_0/case_outcome_flags.csv.gz']
def read(p):
    with gzip.open(p,'rt',encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
hashes=[]
for p in sources:
    h=hashlib.sha256(p.read_bytes()).hexdigest()
    assert any(x.split()[0]==h and x.rstrip().endswith(p.name) for x in (p.parent/'checksums.sha256').read_text().splitlines())
    hashes.append(dict(path=str(p),sha256=h))
c,o=map(read,sources);oc={x['caseid']:x for x in o};assert len(oc)==len(o)==len(c)==174406
assert len({x['caseid'] for x in c})==len(c)
drugs=['dupilumab','tralokinumab','upadacitinib'];formal=dict(zip(drugs,[164515,3637,6254]));ex=Counter();ns=Counter();pos=Counter();rows=[];excluded=[];groups=Counter()
cols=['Intercept','drug_tralokinumab','drug_upadacitinib','sex_Male','sex_Missing','HCP_NonHCP','HCP_Missing','geo_NonUS','quarter_index']
for a in c:
    b=oc[a['caseid']]
    for key in ['drug','primaryid','caseversion','source_quarter','fda_dt','init_fda_dt']:assert a[key]==b[key]
    assert '20220401'<=a['fda_dt']<='20260630'
    q=(int(a['fda_dt'][:4])-2022)*4+(int(a['fda_dt'][4:6])-1)//3-1
    assert 0<=q<=16
    assert a['calendar_quarter']==f"{a['fda_dt'][:4]}Q{(int(a['fda_dt'][4:6])-1)//3+1}"
    y=int(b['Primary11']);assert y in [0,1]
    sex={'F':'Female','M':'Male','Missing / unknown':'Missing'}[a['sex']]
    hcp={'HCP':'HCP','Non-HCP':'Non-HCP','Missing / unknown':'Missing'}[a['hcp']]
    geo=a['us_reporter']
    row=dict(drug=a['drug'],caseid=a['caseid'],primaryid=a['primaryid'],caseversion=a['caseversion'],source_quarter=a['source_quarter'],fda_dt=a['fda_dt'],Primary11=y,sex=sex,HCP=hcp,geography=geo,quarter_index=q)
    if geo=='Missing / unknown':ex[a['drug']]+=1;excluded.append(row);continue
    assert geo in ['US','Non-US'];rows.append(row);ns[a['drug']]+=1;pos[a['drug']]+=y
    x=(1,int(a['drug']=='tralokinumab'),int(a['drug']=='upadacitinib'),int(sex=='Male'),int(sex=='Missing'),int(hcp=='Non-HCP'),int(hcp=='Missing'),int(geo=='Non-US'),q)
    groups[x,y]+=1
assert dict(ns)==dict(zip(drugs,[164513,3636,6192]));assert dict(ex)==dict(zip(drugs,[2,1,62]));assert list(pos.values())==[1068,19,21]
assert not any(a['Primary11'] for a in excluded)
def write(p,rs):
    op=gzip.open if str(p).endswith('.gz') else open
    with op(p,'wt',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rs[0]));w.writeheader();w.writerows(rs)
write(O/'analysis_cases.csv.gz',rows);write(O/'geography_exclusions.csv',excluded)
xs=sorted(set(x for x,y in groups));agg=[dict(zip(cols,x),positive=groups[x,1],negative=groups[x,0],N=groups[x,1]+groups[x,0]) for x in xs]
write(O/'grouped_design_counts.csv',agg)
X=np.array(xs,dtype=float);mixed=[x for x in xs if groups[x,0]>0 and groups[x,1]>0];rank=int(np.linalg.matrix_rank(X));mr=int(np.linalg.matrix_rank(mixed))
assert rank==mr==9
basis=[]
for x in mixed:
    if np.linalg.matrix_rank(basis+[x])>len(basis):basis.append(x)
    if len(basis)==9:break
# Exact integer determinant for a full-rank mixed-outcome basis.
from fractions import Fraction
B=[[Fraction(z) for z in x] for x in basis];det=Fraction(1)
for i in range(9):
    j=next(j for j in range(i,9) if B[j][i]);
    if j!=i:B[i],B[j]=B[j],B[i];det=-det
    pivot=B[i][i];det*=pivot
    for k in range(i+1,9):
        ratio=B[k][i]/pivot
        for j in range(i,9):B[k][j]-=ratio*B[i][j]
assert det!=0
write(O/'no_separation_certificate.csv',[dict(zip(cols,x),positive=groups[x,1],negative=groups[x,0]) for x in basis])
qc=dict(formal_N=formal,analysis_N=dict(ns),geography_exclusions=dict(ex),excluded_positive=0,analysis_total=sum(ns.values()),analysis_positive=sum(pos.values()),drug_positive=dict(pos),source_hashes=hashes,design_columns=cols,design_rank=rank,mixed_pattern_rank=mr,mixed_patterns=len(mixed),unique_design_patterns=len(xs),mixed_basis_exact_determinant=str(det),separation=False,separation_method='A mixed-outcome pattern forces x*b=0 for any weak separating direction. Nine linearly independent mixed-outcome patterns span the entire 9-column design, so b=0 is the only possible direction; complete and quasi-complete separation are both absent.',model='Primary11 ~ drug + sex + HCP + geography + quarter_index',reference_levels='dupilumab; Female; HCP; US',quarter_definition='(year(FDA_DT)-2022)*4 + quarter(FDA_DT)-2; 2022Q2=0 through 2026Q2=16',other_outcomes_fitted=False)
(O/'cohort_qc.json').write_text(json.dumps(qc,ensure_ascii=False,indent=2));print(json.dumps(qc,ensure_ascii=False,indent=2))
