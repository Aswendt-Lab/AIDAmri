"""Launch the validation guidance-first interface detached, then open browser."""
import argparse,json,subprocess,sys,time,urllib.request,os,webbrowser
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def ready(url,dataset):
 try:
  with urllib.request.urlopen(url+'/api/cases',timeout=1) as response:data=json.load(response)
  return data.get('api_version')==4 and data.get('upload_supported') and data.get('demo_supported') and data.get('dataset')==str(dataset) and bool(data.get('cases'))
 except Exception:return False

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--port',type=int,default=8770);parser.add_argument('--dataset',type=Path,default=ROOT/'results/unseen_validation_20261125');parser.add_argument('--no-open',action='store_true');args=parser.parse_args();url=f'http://127.0.0.1:{args.port}';dataset=args.dataset.expanduser().resolve()
 if not (dataset/'manifest.csv').is_file():raise SystemExit('Dataset manifest missing: '+str(dataset))
 if not ready(url,dataset):
  runtime=ROOT/'.gui_runtime';runtime.mkdir(exist_ok=True)
  with (runtime/f'validation_{args.port}.log').open('a') as log:
   process=subprocess.Popen([sys.executable,'-u',str(ROOT/'scripts/validation_guidance_gui.py'),'--port',str(args.port)],cwd=ROOT,env={**os.environ,'LESION_GUIDANCE_DATASET':str(dataset)},stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,close_fds=True)
  (runtime/f'validation_{args.port}.pid').write_text(str(process.pid))
  for _ in range(120):
   if ready(url,dataset):break
   if process.poll() is not None:raise RuntimeError('Server stopped; inspect validation log')
   time.sleep(.25)
  else:raise RuntimeError('Server did not become ready; inspect validation log')
 print('Validation guidance interface: '+url,flush=True)
 if not args.no_open:webbrowser.open(url)
