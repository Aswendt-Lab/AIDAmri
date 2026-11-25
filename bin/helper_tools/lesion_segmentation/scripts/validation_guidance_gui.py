"""Guidance-first validation review, with explicit continuous-region selection."""
import argparse,base64,csv,hashlib,json,os,uuid,shutil,tempfile,re
from email.parser import BytesParser
from email import policy
from pathlib import Path
from datetime import datetime,timezone
from functools import lru_cache
from http.server import BaseHTTPRequestHandler,HTTPServer
from urllib.parse import urlsplit,parse_qs
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mplconfig'))
import nibabel as nib,numpy as np
from guidance_common import metric,sphere_mask,write_csv
from slice_restriction import continuous_seeded_region

BASE=Path(os.environ.get('LESION_GUIDANCE_DATASET',str(ROOT/'results/unseen_validation_20261125'))).expanduser().resolve()
OUT=BASE/'guided_continuity'

def resolve_file(value):
 path=Path(value).expanduser()
 return path if path.is_absolute() else BASE/path

def catalog():
 rows=[]
 for row in csv.DictReader((BASE/'manifest.csv').open()):
  if row['status']!='ready':continue
  for field in ['image','mask']:
   if row.get(field):row[field]=str(resolve_file(row[field]).resolve())
  rows.append(row)
 for path in sorted((BASE/'imported_scans').glob('*/case.json')):
  row=json.loads(path.read_text())
  for field in ['image','mask']:
   if row.get(field):row[field]=str(resolve_file(row[field]).resolve())
  rows.append(row)
 return {r['case_id']:r for r in rows}

def case_listing():
 return dict(api_version=4,dataset=str(BASE),upload_supported=True,demo_supported=True,cases=[dict(id=k,group=r['group'],name=r.get('display_name',r['case']),novelty=r.get('novelty','not_audited'),identity_uncertain=r.get('identity_uncertain','False')) for k,r in catalog().items()])

def demo_config():
 path=BASE/'demo.json'
 if not path.is_file():return dict(available=False,reason='The guided demo is available in the bundled example dataset.')
 config=json.loads(path.read_text());rows=catalog();row=rows.get(config['case'])
 if row is None:return dict(available=False,reason='The example scan is not available in this dataset.')
 folder=session_folder(config['session_id']);session=json.loads((folder/'session.json').read_text())
 if session['case']!=config['case'] or session['seeds']!=config['seeds'] or session['avoidance']!=config['avoidance']:raise ValueError('Demo guidance does not match its recorded source session')
 required=[Path(row['image']),Path(row['mask'])]+[BASE/'predictions'/row['group']/row['case']/v/'lesion_mask.nii.gz' for v in ['native','mirrored']]+[folder/(v+'_lesion_mask.nii.gz') for v in ['native','mirrored']]
 if not all(p.is_file() for p in required):return dict(available=False,reason='A required demo image or saved mask is missing.')
 return dict(**config,available=True,expected_metrics=list(csv.DictReader((folder/'metrics.csv').open())))

def demo_result():
 config=demo_config()
 if not config['available']:raise ValueError(config['reason'])
 result=preview(dict(case=config['case'],slice=config['slice'],seeds=config['seeds'],avoidance=config['avoidance'],session=config['session_id']))
 result.update(session=config['session_id'],demo_replay=True)
 return result

def import_upload(content_type,raw):
 """Validate uploads in staging, then atomically add a persistent case record."""
 if not content_type.startswith('multipart/form-data') or '\n' in content_type or '\r' in content_type:raise ValueError('Expected a file upload')
 message=BytesParser(policy=policy.default).parsebytes(('Content-Type: '+content_type+'\r\nMIME-Version: 1.0\r\n\r\n').encode()+raw)
 if not message.is_multipart() or message.defects:raise ValueError('Invalid file upload')
 files={}
 for part in message.iter_parts():
  field=part.get_param('name',header='content-disposition');filename=part.get_filename()
  if field not in ['image','mask']:raise ValueError('Unexpected upload field')
  if not filename:continue
  if field in files:raise ValueError('Select one T2w image and at most one reference mask')
  name=Path(filename.replace('\\','/')).name
  if not name.lower().endswith(('.nii','.nii.gz')):raise ValueError('Select NIfTI files ending in .nii or .nii.gz')
  payload=part.get_payload(decode=True)
  if not payload or len(payload)>128*1024**2:raise ValueError('Each NIfTI file must be nonempty and no larger than 128 MB')
  files[field]=(name,payload)
 if 'image' not in files:raise ValueError('Choose a T2w image')
 root=BASE/'imported_scans';root.mkdir(exist_ok=True);token=uuid.uuid4().hex
 stage=Path(tempfile.mkdtemp(prefix='.upload_',dir=root));dest=root/token
 try:
  objects={};arrays={};paths={}
  for field,(name,payload) in files.items():
   suffix='.nii.gz' if name.lower().endswith('.gz') else '.nii';path=stage/(field+suffix);path.write_bytes(payload);paths[field]=path.name
   try:
    ni=nib.load(path)
    if len(ni.shape)!=3 or any(n<=0 for n in ni.shape) or np.prod(ni.shape,dtype=np.int64)>32000000:raise ValueError('Expected a 3D volume with at most 32 million voxels')
    if ni.header.get_xyzt_units()[0]!='mm' or not np.isfinite(ni.affine).all() or abs(np.linalg.det(ni.affine[:3,:3]))<1e-12:raise ValueError('NIfTI requires a valid affine and millimetre spatial units')
    array=ni.get_fdata(dtype=np.float32)
    if not np.isfinite(array).all() or (field=='mask' and np.any(array<0)):raise ValueError('Image/mask contains invalid values')
   except Exception as exc:raise ValueError(f'Could not load {field}: {exc}') from exc
   objects[field]=ni;arrays[field]=array
  if 'mask' in objects and (objects['mask'].shape!=objects['image'].shape or not np.allclose(objects['mask'].affine,objects['image'].affine,rtol=0,atol=1e-5)):raise ValueError('The manual mask must match the T2w image shape and affine. No files were added; choose the matching mask or omit it.')
  name=files['image'][0];display=re.sub(r'(?i)\.nii(?:\.gz)?$','',name);safe=re.sub(r'[^A-Za-z0-9_.-]','_',display).strip('.')[:100] or 'scan'
  row=dict(case_id='uploaded_'+token,group='New',case=safe+'_'+token[:8],display_name=display,animal_id='not_verified',image=str((dest/paths['image']).relative_to(BASE)),mask=str((dest/paths['mask']).relative_to(BASE)) if 'mask' in files else '',source_variant='user_selected_T2w',status='ready',novelty='user_loaded_not_overlap_audited',identity_uncertain='True',image_sha256=hashlib.sha256(np.ascontiguousarray(arrays['image']).tobytes()).hexdigest(),mask_sha256=hashlib.sha256(np.ascontiguousarray(arrays['mask']).tobytes()).hexdigest() if 'mask' in files else '',original_image_filename=name,original_mask_filename=files['mask'][0] if 'mask' in files else '')
  (stage/'case.json').write_text(json.dumps(row,indent=2));stage.rename(dest)
  return dict(case=row['case_id'],has_reference='mask' in files,saved_path=str(dest),**case_listing())
 finally:
  if stage.exists():shutil.rmtree(stage)

def encoded(a):return base64.b64encode(np.ascontiguousarray(a.T,dtype=np.uint8).tobytes()).decode()

@lru_cache(maxsize=2)
def load(key):
 rows=catalog()
 if key not in rows:raise ValueError('Unknown dataset')
 row=rows[key];ni=nib.load(row['image']);x=ni.get_fdata(dtype=np.float32);ref=nib.load(row['mask']) if row.get('mask') else None
 if x.ndim!=3 or not np.isfinite(x).all() or ni.header.get_xyzt_units()[0]!='mm':raise ValueError('Expected a finite 3D MRI with millimetre units')
 if ref is not None and (ref.shape!=ni.shape or not np.allclose(ref.affine,ni.affine,rtol=0,atol=1e-5)):raise ValueError('Reference geometry mismatch')
 vmax=np.percentile(x[x>0],99.5) if np.any(x>0) else 1
 return dict(row=row,ni=ni,gray=(255*np.clip(x/max(vmax,1e-12),0,1)).astype(np.uint8),ref=ref.get_fdata()>0 if ref is not None else None,voxel=abs(np.linalg.det(ni.affine[:3,:3])))

def marks(body):
 data=load(body.get('case'));shape=data['ni'].shape
 def point(p):
  if not isinstance(p,list) or len(p)!=3 or any(type(v) is not int or v<0 or v>=shape[k] for k,v in enumerate(p)):raise ValueError('Invalid native voxel coordinate')
  return p
 seeds=body.get('seeds',[]);avoid=body.get('avoidance',[])
 if not isinstance(seeds,list) or not isinstance(avoid,list) or len(seeds)>200 or len(avoid)>200:raise ValueError('Provide up to 200 positive seeds and 200 avoidance ROIs')
 seeds=[point(p) for p in seeds];regions=[]
 for r in avoid:
  if not isinstance(r,dict):raise ValueError('Invalid avoidance ROI')
  radius=r.get('radius_mm')
  if type(radius) not in (int,float) or not np.isfinite(radius) or not .1<=radius<=5:raise ValueError('Radius must be 0.1–5 mm')
  regions.append(dict(point=point(r.get('point')),radius_mm=float(radius)))
 exclusion=sphere_mask(shape,data['ni'].header.get_zooms()[:3],regions)
 if any(exclusion[tuple(p)] for p in seeds):raise ValueError('A positive seed lies inside an avoidance ROI; move or remove the conflicting mark')
 return data,seeds,regions,exclusion

def preview(body):
 data,seeds,regions,exclusion=marks(body);z=body.get('slice',data['ni'].shape[2]//2)
 if type(z) is not int or not 0<=z<data['ni'].shape[2]:raise ValueError('Invalid slice')
 result=dict(case=body['case'],slice=z,shape=list(data['ni'].shape),spacing=list(map(float,data['ni'].header.get_zooms()[:3])),gray=encoded(data['gray'][:,:,z]),reference=encoded(data['ref'][:,:,z]) if data['ref'] is not None else encoded(np.zeros(data['ni'].shape[:2],bool)),has_reference=data['ref'] is not None,exclusion=encoded(exclusion[:,:,z]),versions=None)
 token=body.get('session')
 if token:
  folder=session_folder(token,body.get('demo_run',False));session=json.loads((folder/'session.json').read_text())
  if session['case']!=body['case'] or session['seeds']!=seeds or session['avoidance']!=regions:raise ValueError('Guidance changed; run guided segmentation again')
  result.update(versions={},metrics=list(csv.DictReader((folder/'metrics.csv').open())),warnings=session['warnings'],saved_path=str(folder))
  for version in ['native','mirrored']:
   original=nib.load(BASE/'predictions'/data['row']['group']/data['row']['case']/version/'lesion_mask.nii.gz').get_fdata()>0
   selected=nib.load(folder/(version+'_lesion_mask.nii.gz')).get_fdata()>0
   result['versions'][version]=dict(prediction=encoded(original[:,:,z]),filtered=encoded(selected[:,:,z]),seed_hits=session['versions'][version]['seed_hits'])
 return result

def session_folder(token,demo=False):
 if not isinstance(token,str) or len(token)>80 or any(c not in '0123456789TZ_abcdef' for c in token):raise ValueError('Invalid session ID')
 path=OUT/('demo_sessions' if demo else 'manual_sessions')/token
 if not (path/'session.json').is_file():raise ValueError('Unknown saved session')
 return path

def save_array(array,ni,path):
 hdr=ni.header.copy();hdr.set_data_dtype(np.uint8);obj=nib.Nifti1Image(array.astype(np.uint8),ni.affine,hdr)
 obj.set_qform(ni.get_qform(),int(ni.header['qform_code']));obj.set_sform(ni.get_sform(),int(ni.header['sform_code']));nib.save(obj,path)

def run(body):
 data,seeds,regions,exclusion=marks(body)
 demo=body.get('demo_run',False)
 if type(demo) is not bool:raise ValueError('Invalid demo mode')
 if demo:
  config=demo_config()
  if not config.get('available') or body['case']!=config['case']:raise ValueError('Demo practice requires the bundled tutorial scan')
 if not seeds and not regions:raise ValueError('Place a lesion seed and/or avoidance ROI before running')
 protocol=json.loads((BASE/'protocol.json').read_text());checkpoint=resolve_file(protocol['checkpoint'])
 if hashlib.sha256(checkpoint.read_bytes()).hexdigest()!=protocol['checkpoint_sha256']:raise ValueError('Frozen checkpoint identity changed')
 token=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]
 folder=OUT/('demo_sessions' if demo else 'manual_sessions')/token;folder.mkdir(parents=True)
 guidance=dict(case=body['case'],seeds=seeds,avoidance=regions,source_session=body.get('guidance_source_session'),guidance_saved_before_refinement=True)
 (folder/'guidance.json').write_text(json.dumps(guidance,indent=2))
 prediction_paths=[BASE/'predictions'/data['row']['group']/data['row']['case']/v/'lesion_mask.nii.gz' for v in ['native','mirrored']]
 reused=all(p.is_file() for p in prediction_paths)
 if not reused:
  from guidance_inference import ensure_predictions
  ensure_predictions(BASE,data['row'],checkpoint,protocol['checkpoint_sha256'])
 versions={};rows=[];warnings=[];prediction_hashes={}
 for version in ['native','mirrored']:
  path=BASE/'predictions'/data['row']['group']/data['row']['case']/version/'lesion_mask.nii.gz';ni=nib.load(path)
  if ni.shape!=data['ni'].shape or not np.allclose(ni.affine,data['ni'].affine,rtol=0,atol=1e-6):raise ValueError('Prediction geometry mismatch')
  prediction_hashes[version]=hashlib.sha256(path.read_bytes()).hexdigest()
  raw=ni.get_fdata()>0;selected,info=continuous_seeded_region(raw&~exclusion,seeds)
  versions[version]=info;warnings.extend(version+': '+w for w in info['warnings'])
  save_array(selected,data['ni'],folder/(version+'_lesion_mask.nii.gz'))
  for stage,pred in [('baseline',raw),('guided_continuity',selected)]:rows.append(dict(case=body['case'],group=data['row']['group'],scan=data['row']['case'],version=version,stage=stage,**metric(pred,data['ref'],data['voxel'])))
 save_array(exclusion,data['ni'],folder/'avoidance_mask.nii.gz');write_csv(folder/'metrics.csv',rows)
 session=dict(**guidance,demo_practice=demo,model='frozen clean nnU-Net v2',checkpoint_sha256=protocol['checkpoint_sha256'],baseline_prediction_sha256=prediction_hashes,algorithm='continuous_seeded_region_v1',algorithm_sha256=hashlib.sha256((ROOT/'scripts/slice_restriction.py').read_bytes()).hexdigest(),baseline_inference_reused=reused,network_adapted=False,reference_used_for_selection=False,source=data['row'],versions=versions,warnings=warnings,note='Tutorial practice, excluded from human review summaries.' if demo else 'Post-validation manual refinement. One continuous 3D track, at most one 8-connected region per native slice, no jumps across blank slices. Cannot add network-missed tissue.')
 (folder/'session.json').write_text(json.dumps(session,indent=2));body={**body,'session':token};result=preview(body);result['session']=token
 return result

def sessions():
 result=[]
 for p in (OUT/'manual_sessions').glob('*/session.json'):
  s=json.loads(p.read_text());result.append(dict(id=p.parent.name,case=s['case'],seeds=s['seeds'],avoidance=s['avoidance']))
 return sorted(result,key=lambda s:s['id'],reverse=True)

class Handler(BaseHTTPRequestHandler):
 def respond(self,status,payload,kind='application/json'):
  raw=payload if isinstance(payload,bytes) else json.dumps(payload).encode();self.send_response(status);self.send_header('Content-Type',kind);self.send_header('Content-Length',str(len(raw)));self.send_header('Cache-Control','no-store');self.end_headers();self.wfile.write(raw)
 def do_GET(self):
  if self.path=='/':self.respond(200,(ROOT/'gui/validation_guidance.html').read_bytes(),'text/html; charset=utf-8')
  elif self.path=='/api/cases':self.respond(200,case_listing())
  elif self.path=='/api/sessions':self.respond(200,dict(sessions=sessions()))
  elif self.path=='/api/demo':self.respond(200,demo_config())
  elif self.path=='/api/demo/result':self.respond(200,demo_result())
  else:self.respond(404,dict(error='Not found'))
 def do_POST(self):
  if self.headers.get('Origin') not in (None,'http://'+self.headers.get('Host','')):self.respond(403,dict(error='Use the local application origin'));return
  try:
   size=int(self.headers.get('Content-Length',0))
   if self.path=='/api/import':
    if not 0<size<=257*1024**2:raise ValueError('Upload is too large (maximum 128 MB per file)')
    result=import_upload(self.headers.get('Content-Type',''),self.rfile.read(size));self.respond(200,result);return
   if self.headers.get('Content-Type')!='application/json' or not 0<size<=200000:raise ValueError('Invalid JSON request')
   body=json.loads(self.rfile.read(size))
   if not isinstance(body,dict):raise ValueError('Expected JSON object')
   if self.path=='/api/preview':result=preview(body)
   elif self.path=='/api/run':result=run(body)
   else:self.respond(404,dict(error='Not found'));return
   self.respond(200,result)
  except (ValueError,TypeError,KeyError,ImportError) as exc:self.respond(400,dict(error=str(exc)))
  except Exception as exc:print(repr(exc),flush=True);self.respond(500,dict(error='Request failed; see server log'))
 def log_message(self,*args):pass

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--port',type=int,default=8770);args=parser.parse_args()
 server=HTTPServer(('127.0.0.1',args.port),Handler);print(f'Guidance-first validation interface: http://127.0.0.1:{args.port}',flush=True);server.serve_forever()
