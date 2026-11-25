"""Create a portable guidance dataset from skull-stripped T2 NIfTI images."""
import argparse,csv,hashlib,json,re,shutil
from pathlib import Path
import nibabel as nib,numpy as np
ROOT=Path(__file__).resolve().parents[1]
SUFFIXES=['BiasBet','BiasMicoThresBet','ThresBet','Bet']
FROZEN_SHA='f9bc7637b7aabf42ee3378e78d98b09f57bfd3f20a590d128d1e52f62552d33f'

def create(output,checkpoint,pairs):
    output=Path(output).expanduser().resolve();checkpoint=Path(checkpoint).expanduser().resolve()
    if output.exists():raise ValueError('Output folder already exists; choose a new folder to preserve existing work')
    if hashlib.sha256(checkpoint.read_bytes()).hexdigest()!=FROZEN_SHA:raise ValueError('Expected the bundled frozen clean-v2 checkpoint')
    rows=[]
    for index,(image,mask,name,variant) in enumerate(pairs):
        ni=nib.load(image);x=ni.get_fdata(dtype=np.float32)
        if x.ndim!=3 or not np.isfinite(x).all() or ni.header.get_xyzt_units()[0]!='mm':raise ValueError('Expected finite 3D T2 MRI with millimetre units: '+str(image))
        if mask:
            nm=nib.load(mask);ref=nm.get_fdata(dtype=np.float32)
            if nm.shape!=ni.shape or not np.allclose(nm.affine,ni.affine,rtol=0,atol=1e-5) or not np.isfinite(ref).all() or np.any(ref<0):raise ValueError('Reference mask geometry/data mismatch: '+str(mask))
        safe=re.sub(r'[^A-Za-z0-9_.-]','_',name).strip('.') or 'scan'
        if any(r['case']==safe for r in rows):safe+=f'_{index:03d}'
        rows.append(dict(case_id=f'scan_{index:03d}',group='New',case=safe,animal_id='not_verified',image=f'inputs/scan_{index:03d}/image.nii.gz',mask=f'inputs/scan_{index:03d}/reference.nii.gz' if mask else '',source_variant=variant,status='ready',novelty='user_imported_not_overlap_audited',identity_uncertain='True',image_sha256=hashlib.sha256(np.ascontiguousarray(x).tobytes()).hexdigest(),mask_sha256=hashlib.sha256(np.ascontiguousarray(ref).tobytes()).hexdigest() if mask else ''))
    if not rows:raise ValueError('No matching images found; use --image for a single arbitrary-named T2 NIfTI file')
    model=checkpoint.parent.parent
    for name in ['plans.json','dataset.json']:
        if not (model/name).is_file():raise ValueError('Model metadata missing: '+name)
    output.mkdir(parents=True)
    target=output/'model';(target/'fold_0').mkdir(parents=True)
    for name in ['plans.json','dataset.json']:shutil.copy2(model/name,target/name)
    shutil.copy2(checkpoint,target/'fold_0/checkpoint_best.pth')
    for row,(image,mask,_,_) in zip(rows,pairs):
        folder=output/row['image'];folder.parent.mkdir(parents=True)
        # Preserve bytes and compression: uncompressed input is written as genuine gzip.
        if str(image).endswith('.gz'):shutil.copy2(image,folder)
        else:nib.save(nib.load(image),folder)
        if mask:
            dest=output/row['mask']
            if str(mask).endswith('.gz'):shutil.copy2(mask,dest)
            else:nib.save(nib.load(mask),dest)
    with (output/'manifest.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (output/'protocol.json').write_text(json.dumps(dict(checkpoint='model/fold_0/checkpoint_best.pth',checkpoint_sha256=FROZEN_SHA,model='frozen clean nnU-Net 2D v2',native_input_is_primary=True,compute_both_versions=True,training=False,orientation_note='Preserves native axes; axis 2 must be the anatomical slice stack as in examples. No automatic orientation harmonization.'),indent=2))
    print(f'Created {len(rows)} scan(s): {output}\nNo detection has run. Launch this dataset, place guidance, then click Run.')
    return rows

def main():
    parser=argparse.ArgumentParser(description=__doc__);source=parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--image',type=Path,help='One skull-stripped T2 NIfTI file, any filename');source.add_argument('--input-folder',type=Path,help='Recursively pair *BiasBet / *Bet with optional *Stroke_mask')
    parser.add_argument('--mask',type=Path,help='Optional binary/manual reference for --image');parser.add_argument('--name',help='Optional case name for --image');parser.add_argument('--output',type=Path,required=True)
    default=ROOT/'model/fold_0/checkpoint_best.pth'
    if not default.exists():default=ROOT/'training/pilot_v2_clean/nnunet/results/Dataset502_MouseStrokeClean/nnUNetTrainer__nnUNetPlansPilot__2d/fold_0/checkpoint_best.pth'
    parser.add_argument('--checkpoint',type=Path,default=default);args=parser.parse_args()
    if args.mask and not args.image:parser.error('--mask is used only with --image')
    pairs=[]
    if args.image:pairs=[(args.image,args.mask,args.name or args.image.name.removesuffix('.gz').removesuffix('.nii'),'explicit_user_input')]
    else:
        chosen={}
        for suffix in SUFFIXES:
            for image in sorted(args.input_folder.rglob('*'+suffix+'.nii.gz')):
                prefix=image.name.removesuffix(suffix+'.nii.gz');key=(image.parent,prefix)
                if key in chosen:continue
                # Longest suffix is already handled before its Bet suffix.
                if any(image.name.endswith(other+'.nii.gz') for other in SUFFIXES[:SUFFIXES.index(suffix)]):continue
                mask=image.with_name(prefix+'Stroke_mask.nii.gz');chosen[key]=(image,mask if mask.is_file() else None,prefix,suffix)
        pairs=list(chosen.values())
    create(args.output,args.checkpoint,pairs)

if __name__=='__main__':main()
