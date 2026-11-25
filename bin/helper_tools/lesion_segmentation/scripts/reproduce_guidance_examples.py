"""Verify portable examples and replay every saved human-guidance session."""
import argparse,csv,hashlib,json,os,tempfile
from pathlib import Path
import nibabel as nib,numpy as np
from scipy.ndimage import label
from guidance_common import metric,sphere_mask
from slice_restriction import continuous_seeded_region
ROOT=Path(__file__).resolve().parents[1]

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--network',action='store_true',help='Also rerun the frozen network on every example (requires inference dependencies)');args=parser.parse_args()
 release_files=json.loads((ROOT/'SHA256SUMS.json').read_text())
 for name,digest in release_files.items():
  assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,'Changed or missing release file: '+name
 base=ROOT/'examples';sources={r['case_id']:r for r in csv.DictReader((base/'manifest.csv').open())};count=0
 for file in sorted((base/'guided_continuity/manual_sessions').glob('*/session.json')):
  if str(file.relative_to(ROOT)) not in release_files:continue
  session=json.loads(file.read_text());r=sources[session['case']];ni=nib.load(base/r['image']);ref=nib.load(base/r['mask']).get_fdata()>0
  exclusion=sphere_mask(ni.shape,ni.header.get_zooms()[:3],session['avoidance'])
  rows=list(csv.DictReader((file.parent/'metrics.csv').open()))
  for version in ['native','mirrored']:
   raw=nib.load(base/'predictions'/r['group']/r['case']/version/'lesion_mask.nii.gz').get_fdata()>0
   selected,info=continuous_seeded_region(raw&~exclusion,session['seeds']);saved=nib.load(file.parent/(version+'_lesion_mask.nii.gz'))
   assert np.array_equal(selected,saved.get_fdata()>0) and np.allclose(saved.affine,ni.affine,rtol=0,atol=1e-6)
   assert not np.any(selected&exclusion) and label(selected,np.ones((3,3,3)))[1]<=1
   for z in range(selected.shape[2]):assert label(selected[:,:,z],np.ones((3,3)))[1]<=1
   assert info['seed_hits']==session['versions'][version]['seed_hits']
   for stage,pred in [('baseline',raw),('guided_continuity',selected)]:
    fresh=metric(pred,ref,abs(np.linalg.det(ni.affine[:3,:3])));old=next(row for row in rows if row['version']==version and row['stage']==stage)
    for key,value in fresh.items():assert (value is None and old[key]=='') or (value is not None and np.isclose(value,float(old[key]),rtol=0,atol=1e-10)),(file,key)
   count+=1
 print(f'PASS: release file checksums, {len(sources)} example scans, {count} saved masks, all guidance, geometry, continuity and quantitative metrics',flush=True)
 if args.network:
  from guidance_inference import ensure_predictions
  protocol=json.loads((base/'protocol.json').read_text());checkpoint=(base/protocol['checkpoint']).resolve()
  scores=[]
  with tempfile.TemporaryDirectory() as folder:
   target=Path(folder)
   for key,r in sources.items():
    r={**r,'image':str(base/r['image'])};ensure_predictions(target,r,checkpoint,protocol['checkpoint_sha256'])
    for version in ['native','mirrored']:
     relative=Path('predictions')/r['group']/r['case']/version/'lesion_mask.nii.gz'
     a=nib.load(base/relative).get_fdata()>0;b=nib.load(target/relative).get_fdata()>0
     scores.append(dict(case=key,version=version,differing_voxels=int(np.count_nonzero(a!=b)),dice=metric(a,b,1)['dice']))
  (ROOT/'network_reproduction.json').write_text(json.dumps(scores,indent=2));print('Network rerun finished. See network_reproduction.json; small floating-point differences across devices may change boundary voxels.',flush=True)

if __name__=='__main__':main()
