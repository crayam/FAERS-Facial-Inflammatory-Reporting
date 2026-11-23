from pathlib import Path
import csv,gzip,json,hashlib,math
from collections import Counter,defaultdict
R=Path(__file__).resolve().parents[2];W=Path(__file__).parent;O=R/'outputs/prr_pt_composition_v1_0'
def read(p):
 op=gzip.open if str(p).endswith('.gz') else open
 with op(p,'rt',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def write(p,rows):
 with open(p,'w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
sources=[]
def source(rel):
 p=R/rel;h=hashlib.sha256(p.read_bytes()).hexdigest();lines=(p.parent/'checksums.sha256').read_text().splitlines();assert any(x.split()[0]==h and x.rstrip().endswith(p.name) for x in lines),rel
 sources.append(dict(path=rel,sha256=h));return p
fp=source('outputs/all_outcomes_v1_0/case_outcome_flags.csv.gz');pp=source('outputs/primary11_v1_0/cohort_case_PT.csv.gz');dp=source('outputs/all_outcomes_v1_0/phenotype_definition.json');cp=source('outputs/crude_comparative_v1_0/crude_comparative_results.csv');rp=source('outputs/primary11_v1_0/primary11_raw_reac_evidence.csv');pcp=source('outputs/all_outcomes_v1_0/dictionary_pt_counts.csv')
flags=read(fp);assert len(flags)==len({a['caseid'] for a in flags})==174406;byid={a['caseid']:a for a in flags};drugs=['dupilumab','tralokinumab','upadacitinib'];ys=['Primary11','Narrow6','BroadSX23'];defs=json.loads(dp.read_text());pts=defs['definitions']['Primary11'];assert len(pts)==len(set(pts))==11
N=Counter(a['drug'] for a in flags);assert [N[d] for d in drugs]==[164515,3637,6254]
counts={y:Counter(a['drug'] for a in flags if a[y]=='1') for y in ys};assert [[counts[y][d] for d in drugs] for y in ys]==[[1068,19,21],[754,13,14],[2347,51,30]]
casepts=defaultdict(set);pairs=set();allpairs=set();reconstructed={y:set() for y in ys}
for a in read(pp):
 b=byid[a['caseid']]
 for k in ['drug','primaryid','source_quarter']:assert a[k]==b[k]
 key=(a['drug'],a['caseid'],a['pt']);assert key not in allpairs;allpairs.add(key)
 for y in ys:
  if a['pt'] in defs['definitions'][y]:reconstructed[y].add(a['caseid'])
 if a['pt'] in pts:casepts[a['caseid']].add(a['pt']);pairs.add(key)
for y in ys:assert reconstructed[y]=={a['caseid'] for a in flags if a[y]=='1'}
for a in flags:
 p=casepts.get(a['caseid'],set());stored=set(a['Primary11_matched_pts'].split('; ')) if a['Primary11_matched_pts'] else set();assert p==stored and len(p)==int(a['Primary11_matched_pt_count']) and bool(p)==bool(int(a['Primary11']))
rawpairs=set();raw_rows=0
for a in read(rp):
 b=byid[a['caseid']];assert a['primaryid']==b['primaryid'] and a['source_quarter']==b['source_quarter'] and a['pt'] in pts
 rawpairs.add((b['drug'],a['caseid'],a['pt']));raw_rows+=1
assert rawpairs==pairs
priorpts={a['pt']:a for a in read(pcp) if a['outcome']=='Primary11'};composition=[];evidence=[]
for pt in pts:
 for drug in drugs:
  ids={cid for d,cid,p in pairs if d==drug and p==pt};n=len(ids);assert n==int(priorpts[pt][drug]);composition.append(dict(PT=pt,drug=drug,unique_CASEID_n=n,pct_Primary11_union=n/counts['Primary11'][drug],pct_formal_cohort=n/N[drug],Primary11_union_n=counts['Primary11'][drug],formal_N=N[drug]))
  for cid in sorted(ids):evidence.append(dict(drug=drug,CASEID=cid,PRIMARYID=byid[cid]['primaryid'],source_quarter=byid[cid]['source_quarter'],PT=pt))
overlap=[];hist={};uniontotals=[]
for drug in drugs:
 h=Counter(len(casepts[a['caseid']]) for a in flags if a['drug']==drug and a['Primary11']=='1');hist[drug]=dict(h);U=counts['Primary11'][drug]
 for cat,n in [('Exactly 1 PT',h[1]),('Exactly 2 PTs',h[2]),('Exactly 3 PTs',h[3]),('≥4 PTs',sum(v for k,v in h.items() if k>=4)),('Overall ≥2 PTs',sum(v for k,v in h.items() if k>=2))]:overlap.append(dict(drug=drug,category=cat,n=n,pct_Primary11_union=n/U,Primary11_union_n=U,mutually_exclusive=cat!='Overall ≥2 PTs'))
 assert sum(h.values())==U
 uniontotals.append(dict(drug=drug,formal_N=N[drug],Primary11_union_n=U,sum_PT_n=sum(k*v for k,v in h.items()),multi_PT_n=sum(v for k,v in h.items() if k>=2),max_PTs=max(h)))
prior={(a['outcome'],a['comparator']):a for a in read(cp)};results=[];z=1.9599639845400534
for y in ys:
 for comp in drugs[1:]:
  a=counts[y]['dupilumab'];b=N['dupilumab']-a;c=counts[y][comp];dd=N[comp]-c;prev=prior[y,comp]
  for k,v in zip(['a_index_positive','b_index_negative','c_comparator_positive','d_comparator_negative'],[a,b,c,dd]):assert int(prev[k])==v
  zero=min(a,b,c,dd)==0;assert not zero
  prr=(a/(a+b))/(c/(c+dd));se=math.sqrt(1/a-1/(a+b)+1/c-1/(c+dd));lo=math.exp(math.log(prr)-z*se);hi=math.exp(math.log(prr)+z*se)
  # Independently check algebraically equivalent binomial variance and ratio.
  p1=a/(a+b);p2=c/(c+dd);var=(1-p1)/((a+b)*p1)+(1-p2)/((c+dd)*p2)
  assert math.isclose(var,se*se,rel_tol=1e-14) and math.isclose(a*(c+dd)/((a+b)*c),prr,rel_tol=1e-14)
  ror=a*dd/(b*c);ror_se=math.sqrt(1/a+1/b+1/c+1/dd);rrlo=math.exp(math.log(ror)-z*ror_se);rrhi=math.exp(math.log(ror)+z*ror_se)
  for v,k in [(ror,'cROR'),(rrlo,'CI95_lower'),(rrhi,'CI95_upper')]:assert math.isclose(v,float(prev[k]),rel_tol=1e-13)
  results.append(dict(outcome=y,index_drug='dupilumab',comparator=comp,a=a,b=b,c=c,d=dd,N_index=N['dupilumab'],N_comparator=N[comp],PRR=prr,PRR_lower=lo,PRR_upper=hi,SE_log_PRR=se,crude_ROR=float(prev['cROR']),ROR_lower=float(prev['CI95_lower']),ROR_upper=float(prev['CI95_upper']),zero_cell=zero,continuity_correction=0))
qc=[]
def check(item,drug,expected,observed,details=''):
 assert expected==observed,(item,drug,expected,observed);qc.append(dict(check=item,drug=drug,expected=expected,observed=observed,status='PASS',details=details))
for drug,target in zip(drugs,[164515,3637,6254]):check('Formal denominator',drug,target,N[drug],'Matches prior formal crude ROR population; no geography exclusion')
for y in ys:
 for drug in drugs:check(y+' union count',drug,counts[y][drug],sum(byid[cid]['drug']==drug for cid in reconstructed[y]),'Rebuilt from retained CASEID–PT set, matched to frozen flags and prior crude counts')
for drug in drugs:
 p=[a for a in evidence if a['drug']==drug];check('Duplicate CASEID within each PT',drug,0,len(p)-len({(a['PT'],a['CASEID']) for a in p}),'Distinct PT-case pairs; no duplicate within any of 11 PTs')
 h=hist[drug];check('Exclusive overlap categories sum',drug,counts['Primary11'][drug],sum(a['n'] for a in overlap if a['drug']==drug and a['mutually_exclusive']),'Exactly 1 + exactly 2 + exactly 3 + ≥4; overall ≥2 is not added')
 check('Overall ≥2 matches union minus exactly 1',drug,counts['Primary11'][drug]-h.get(1,0),sum(v for k,v in h.items() if k>=2))
 check('PT count sum = multiplicity-weighted cases',drug,sum(a['unique_CASEID_n'] for a in composition if a['drug']==drug),sum(k*v for k,v in h.items()),'PT sum is not the union denominator')
check('Formal CASEID uniqueness','All',174406,len(byid));check('Raw retained REAC vs frozen PT pairs','All',len(pairs),len(rawpairs),'Exact set equality checked, including PRIMARYID and source quarter')
check('Crude 2×2 and ROR reproduce prior','All',6,len(results),'All 24 cells and six ROR/CI triples match prior output')
check('Zero cells in six PRR tables','All',0,sum(a['zero_cell'] for a in results));check('Continuity correction','All',0,sum(a['continuity_correction'] for a in results));check('Frozen source checksums','All',len(sources),len(sources),'Every source hash matched its existing checksums.sha256')
manifest=dict(sources=sources,definitions=defs['definitions'],formal_N=dict(N),union_counts={y:dict(counts[y]) for y in ys},PT_pair_count=len(pairs),raw_evidence_rows=raw_rows,raw_duplicate_rows=raw_rows-len(rawpairs),scope='PRR supplementary metric and Primary11 PT composition only',PRR_formula='[a/(a+b)]/[c/(c+d)]',CI_formula='exp(log(PRR) ± 1.9599639845400534*sqrt(1/a−1/(a+b)+1/c−1/(c+d)))',method_reference='https://search.r-project.org/CRAN/refmans/pvda/html/prr.html',no_other_analyses=True)
for name,rr in [('prr_results',results),('pt_composition',composition),('pt_overlap',overlap),('primary11_case_pt_evidence',evidence),('qc',qc)]:write(O/(name+'.csv'),rr)
(O/'analysis_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2));(W/'workbook_data.json').write_text(json.dumps(dict(manifest=manifest,prr=results,composition=composition,overlap=overlap,totals=uniontotals,qc=qc),ensure_ascii=False,indent=2))
print(json.dumps(dict(results=results,totals=uniontotals,histogram=hist,qc=len(qc),raw_duplicates=raw_rows-len(rawpairs)),ensure_ascii=False,indent=2))
