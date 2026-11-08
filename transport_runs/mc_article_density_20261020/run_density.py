from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import json,os,subprocess,time
R=Path(__file__).resolve().parent
jobs=[dict(n_cm2=n,seed=5102000+100*i+j,N=131072,duration_ps=2) for i,n in enumerate([3.2e12,6.4e12,1e13,1.5e13,2e13]) for j in range(12)]
(R/'cases.json').write_text(json.dumps(jobs,indent=2)+'\n')
def run(j):
 out=R/'cases'/f"n{j['n_cm2']:.4g}_s{j['seed']}";out.mkdir(exist_ok=False)
 cmd=[str(R/'optical_response'),str(j['n_cm2']),str(j['N']),str(R/'modes_7_11.txt')]
 p=subprocess.run(cmd,capture_output=True,text=True,env=dict(os.environ,OPTICAL_RESPONSE_SEED=str(j['seed'])))
 (out/'stdout.txt').write_text(p.stdout);(out/'stderr.txt').write_text(p.stderr)
 (out/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
 assert p.returncode==0,p.stderr
 d=json.loads(p.stdout);d['seed']=j['seed'];(out/'result.json').write_text(json.dumps(d,indent=2)+'\n')
 assert d['max_pair_balance_relative_residual']<1e-10 and d['max_energy_error_eV']<1e-10
 assert abs(sum(d['residence_ps'])/(d['N']*2)-1)<1e-10
 assert all(e>0 for e in d['events'])
 return dict(n_cm2=j['n_cm2'],seed=j['seed'],elapsed_s=d['elapsed_s'])
start=time.monotonic();results=[]
with ThreadPoolExecutor(max_workers=2) as pool:
 for f in as_completed([pool.submit(run,j) for j in jobs]):
  results.append(f.result());(R/'run_results.json').write_text(json.dumps(results,indent=2)+'\n')
  if len(results)%12==0:print(json.dumps(dict(completed=len(results),total=len(jobs),wall_s=time.monotonic()-start)),flush=True)
(R/'run_validation.json').write_text(json.dumps(dict(status='completed',new_cases=len(results),reused_cases=12,N_per_set=131072,sets_per_density=12,duration_ps=2,new_trajectories=len(results)*131072,wall_s=time.monotonic()-start),indent=2)+'\n')
print('Completed density series',flush=True)
