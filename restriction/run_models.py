from pathlib import Path
import json,subprocess,time
R=Path(__file__).resolve().parents[2];O=R/'outputs/restriction_sensitivity_v1_0';W=Path(__file__).parent
m=json.loads((O/'analysis_manifest.json').read_text())
for a in m['new_models']:
 p=O/'models'/a['model_id'];start=time.time()
 with open(p/'fit.log','w') as log:
  try:
   z=subprocess.run(['Rscript',str(W/'fit.R'),str(p)],stdout=log,stderr=subprocess.STDOUT,timeout=600)
   status='PASS' if z.returncode==0 and (p/'fit_results.json').exists() else 'FAILED'
   reason='' if status=='PASS' else (p/'fit.log').read_text()[-2000:]
  except subprocess.TimeoutExpired:status='FAILED';reason='600-second computation limit reached; design not changed.'
 (p/'run_status.json').write_text(json.dumps(dict(status=status,reason=reason,elapsed_seconds=time.time()-start),indent=2))
 print(a['model_id'],status,round(time.time()-start,2),flush=True)
