"""Lazy frozen-checkpoint inference; called only after guidance is saved."""
import hashlib,json,os
from pathlib import Path
import nibabel as nib,numpy as np

def ensure_predictions(base,row,checkpoint,expected_sha):
    try:
        import torch
    except ImportError as exc:
        raise ImportError('New-scan detection dependencies are missing. Double-click Setup.command in the application folder, then retry. Your guidance has been saved.') from exc
    if hashlib.sha256(checkpoint.read_bytes()).hexdigest()!=expected_sha:
        raise ValueError('Checkpoint checksum mismatch')
    work=base/'inference_runtime';work.mkdir(exist_ok=True)
    for env,sub in [('nnUNet_raw','raw'),('nnUNet_preprocessed','preprocessed'),('nnUNet_results','results')]:
        folder=work/sub;folder.mkdir(exist_ok=True);os.environ[env]=str(folder)
    os.environ['nnUNet_compile']='false'
    try:
        from nnunetv2.inference.predict_from_raw_data import nnUNetPredictor
    except ImportError as exc:
        raise ImportError('New-scan detection dependencies are missing. Double-click Setup.command in the application folder, then retry. Your guidance has been saved.') from exc
    inputs=work/'inputs';outputs=work/'outputs';inputs.mkdir(exist_ok=True);outputs.mkdir(exist_ok=True)
    source=nib.load(row['image']);x=source.get_fdata(dtype=np.float32);identifiers=[]
    for version in ['native','mirrored']:
        identifier=row['case_id']+'_'+version;identifiers.append(identifier)
        array=np.flip(x,0).copy() if version=='mirrored' else x
        nib.save(nib.Nifti1Image(array,source.affine,source.header),inputs/(identifier+'_0000.nii.gz'))
    torch.set_num_threads(4)
    device=torch.device('cuda' if torch.cuda.is_available() else 'mps' if torch.backends.mps.is_available() else 'cpu')
    predictor=nnUNetPredictor(tile_step_size=.5,use_gaussian=True,use_mirroring=True,perform_everything_on_device=False,device=device,verbose=False,verbose_preprocessing=False,allow_tqdm=False)
    predictor.initialize_from_trained_model_folder(str(checkpoint.parent.parent),use_folds=(0,),checkpoint_name=checkpoint.name)
    predictor.predict_from_files([[str(inputs/(i+'_0000.nii.gz'))] for i in identifiers],[str(outputs/i) for i in identifiers],save_probabilities=False,overwrite=True,num_processes_preprocessing=1,num_processes_segmentation_export=1)
    for version,identifier in zip(['native','mirrored'],identifiers):
        pred=nib.load(outputs/(identifier+'.nii.gz'))
        if pred.shape!=source.shape or not np.allclose(pred.affine,source.affine,rtol=0,atol=1e-6):raise ValueError('Inference output geometry mismatch')
        mask=pred.get_fdata()>0
        if version=='mirrored':mask=np.flip(mask,0).copy()
        dest=base/'predictions'/row['group']/row['case']/version/'lesion_mask.nii.gz';dest.parent.mkdir(parents=True,exist_ok=True)
        hdr=source.header.copy();hdr.set_data_dtype(np.uint8);obj=nib.Nifti1Image(mask.astype(np.uint8),source.affine,hdr)
        obj.set_qform(source.get_qform(),int(source.header['qform_code']));obj.set_sform(source.get_sform(),int(source.header['sform_code']));nib.save(obj,dest)
    (outputs/(row['case_id']+'_provenance.json')).write_text(json.dumps(dict(checkpoint_sha256=expected_sha,device=str(device),guidance_not_network_input=True,reference_used_for_inference=False),indent=2))
