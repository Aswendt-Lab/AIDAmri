"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne



"""
from __future__ import print_function

import os
import sys
import numpy as np
import nibabel as nib
import nibabel.nifti1 as nii
import pv_parseBruker_md_np as pB
import P2_IDLt2_mapping as mapT2

class Bruker2Nifti:
    def __init__(self, study, expno, procno, rawfolder, procfolder, ftype='NIFTI_GZ'):
        self.study = study
        self.expno = str(expno)
        self.procno = str(procno)
        self.rawfolder = rawfolder
        self.procfolder = procfolder
        self.ftype = ftype

    def read_2dseq(self, map_raw=False, pv6=False, sc=1.0):
        study = self.study
        expno = self.expno
        procno = self.procno
        rawfolder = self.rawfolder

        self.acqp = pB.parsePV(os.path.join(rawfolder, study, expno, 'acqp'))
        self.method = pB.parsePV(os.path.join(rawfolder, study, expno, 'method'))
        self.subject = pB.parsePV(os.path.join(rawfolder, study, 'subject'))
        # get header information
        datadir = os.path.join(rawfolder, study, expno, 'pdata', procno)
        #self.d3proc = pB.parsePV(os.path.join(datadir, 'd3proc'))   # removed for PV6
        self.visu_pars = pB.parsePV(os.path.join(datadir, 'visu_pars'))
        hdr = pB.getNiftiHeader(self.visu_pars, sc=sc)
        #print("hdr:", hdr)

        if hdr is None or not isinstance(hdr[12], str):
            return

        # read '2dseq' file
        f_id = open(os.path.join(datadir, '2dseq'), 'rb')
        data = np.fromfile(f_id, dtype=np.dtype(hdr[12])).reshape(hdr[1], hdr[2], hdr[3], hdr[4], order='F')
        f_id.close()

        # map to raw data range (PV6)
        if map_raw:
            visu_core_data_slope = np.array(map(float, self.visu_pars['VisuCoreDataSlope'].split()), dtype=np.float32)
            visu_core_data_offs = np.array(map(float, self.visu_pars['VisuCoreDataOffs'].split()), dtype=np.float32)
            visu_core_data_shape = list(data.shape)
            visu_core_data_shape[:2] = (1, 1)
            if pv6:
                data = data / visu_core_data_slope.reshape(visu_core_data_shape)
            else:
                data = data * visu_core_data_slope.reshape(visu_core_data_shape)
            data = data + visu_core_data_offs.reshape(visu_core_data_shape)

        # NIfTI image
        nim = nii.Nifti1Image(data, None)

        # NIfTI header
        #header = nim.header
        header = nim.get_header()
        #print("header:"); print(header)
        header['pixdim'] = [0.0, hdr[5], hdr[6], hdr[7], hdr[8], 0.0, 0.0, 0.0]
        #nim.setXYZUnit('mm')
        header.set_xyzt_units(xyz='mm', t=None)
        #nim.header = header
        #header = nim.get_header()
        #print("header:"); print(header)

        # write header in xml structure
        #xml = pB.getXML(datadir + "/")
        xml = pB.getXML(os.path.join(datadir, 'visu_pars'))
        #print("xml:"); print(xml)

        # add protocol information (method, acqp, visu_pars, d3proc) to Nifti's header extensions
        #nim.extensions += ('comment', xml)
        #extension = nii.Nifti1Extension('comment', xml)

        self.hdr = hdr
        self.nim = nim
        self.xml = xml

    def save_nifti(self, subfolder=''):


        procfolder = os.path.join(self.procfolder, self.study)
        if not os.path.isdir(procfolder):
            os.mkdir(procfolder)

        if "Localizer" in self.acqp['ACQ_protocol_name']:
            procfolder = os.path.join(self.procfolder, self.study, "Localizer")
        elif "DTI" in self.acqp['ACQ_protocol_name'] or "Diffusion" in self.acqp['ACQ_protocol_name']:
            procfolder = os.path.join(self.procfolder, self.study, "DTI")
        elif "fMRI" in self.acqp['ACQ_protocol_name']:
            procfolder = os.path.join(self.procfolder, self.study, "fMRI")
        elif "Turbo" in self.acqp['ACQ_protocol_name']:
            procfolder = os.path.join(self.procfolder, self.study, "T2w")
        elif "MSME" in self.acqp['ACQ_protocol_name']:
            procfolder = os.path.join(self.procfolder, self.study, "T2map")
        else:
            procfolder = os.path.join(self.procfolder, self.study, subfolder)



        if not os.path.isdir(procfolder):
            os.mkdir(procfolder)

        if self.ftype   == 'NIFTI_GZ': ext = 'nii.gz'
        elif self.ftype == 'NIFTI':    ext = 'nii'
        elif self.ftype == 'ANALYZE':  ext = 'img'
        else: ext = 'nii.gz'

        fname = '.'.join([self.study, self.expno, self.procno, ext])

        # write Nifti file

        print(os.path.join(procfolder, fname))
        if not hasattr(self, 'nim'):
            return
        nib.save(self.nim, os.path.join(procfolder, fname))

        return os.path.join(procfolder, fname)

    def save_table(self, subfolder=''):
        procfolder = os.path.join(self.procfolder, self.study)
        if not os.path.isdir(procfolder):
            os.mkdir(procfolder)

        procfolder = os.path.join(self.procfolder, self.study, subfolder)
        if not os.path.isdir(procfolder):
            os.mkdir(procfolder)

        #dw_bval_each = float(self.method['PVM_DwBvalEach'])
        if 'PVM_DwEffBval' in self.method:
            dw_eff_bval = np.array(list(map(float, self.method['PVM_DwEffBval'].split())), dtype=np.float32)
        #print("dw_bval_each:", dw_bval_each)
        #print("dw_eff_bval:"); print(dw_eff_bval)

        if 'PVM_DwAoImages' in self.method:
            dw_ao_images = int(self.method['PVM_DwAoImages'])

        if 'PVM_DwNDiffDir' in self.method:
            dw_n_diff_dir = int(self.method['PVM_DwNDiffDir'])
        #print("dw_ao_images:", dw_ao_images)
        #print("dw_n_diff_dir:", dw_n_diff_dir)

            if 'PVM_DwDir' in self.method:
                dw_dir = np.array(list(map(float, self.method['PVM_DwDir'].split())), dtype=np.float32)
                dw_dir = dw_dir.reshape((dw_n_diff_dir, 3))

                nd = dw_ao_images + dw_n_diff_dir
                bvals = np.zeros(nd, dtype=np.float32)
                dwdir = np.zeros((nd, 3), dtype=np.float32)

                bvals[dw_ao_images:] = dw_eff_bval[dw_ao_images:]
                dwdir[dw_ao_images:] = dw_dir

                fname = '.'.join([self.study, self.expno, self.procno, 'btable', 'txt'])
                print(os.path.join(procfolder, fname))

                # Open btable file to write binary (windows format)

                #fid = open(os.path.join(procfolder, fname), 'wb') - py 2.6
                fid = open(os.path.join(procfolder, fname),mode='w',buffering=-1)

                for i in range(nd):
                    fid.write("%.4f" % (bvals[i],) + " %.8f %.8f %.8f" % tuple(dwdir[i]))
                    #print("%.4f" % (bvals[i],) + " %.8f %.8f %.8f" % tuple(dwdir[i]), end="\r\n", file=fid) - py 2.6

                # Close file
                fid.close()

                fname = '.'.join([self.study, self.expno, self.procno, 'bvals', 'txt'])
                print(os.path.join(procfolder, fname))

                # Open bvals file to write binary (unix format)
                fid = open(os.path.join(procfolder, fname), mode='w', buffering=-1)
                #fid = open(os.path.join(procfolder, fname), 'wb') - py 2.6


                fid.write(" ".join("%.4f" % (bvals[i],) for i in range(nd)))
                #print(" ".join("%.4f" % (bvals[i],) for i in range(nd)), end=chr(10), file=fid) - py 2.6

                # Close bvals file
                fid.close()

                fname = '.'.join([self.study, self.expno, self.procno, 'bvecs', 'txt'])
                print(os.path.join(procfolder, fname))

                # Open bvecs file to write binary (unix format)
                fid = open(os.path.join(procfolder, fname), mode='w', buffering=-1)
                #fid = open(os.path.join(procfolder, fname), 'wb') - py 2.6

                for k in range(3):
                    fid.write(" ".join("%.8f" % (dwdir[i,k],) for i in range(nd)))
                    #print(" ".join("%.8f" % (dwdir[i,k],) for i in range(nd)), end=chr(10), file=fid)- py 2.6

                # Close bvecs file
                fid.close()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Convert ParaVision to NIfTI')

    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument('-i','--input_folder', help='raw data folder')
    # parser.add_argument('-o','--output_folder', help='output data folder')
    # parser.add_argument('study', help='study name')
    # parser.add_argument('expno', help='experiment number')
    # parser.add_argument('procno', help='processed (reconstructed) images number')
    parser.add_argument('-f','--model',
                        help='T2_2p  (default)  : Two   parameter T2 decay S(t) = S0 * exp(-t/T2)\n'
                             'T2_3p             : Three parameter T2 decay S(t) = S0 * exp(-t/T2) + C'
                        , nargs='?', const='T2_2p', type=str, default='T2_2p')
    parser.add_argument('-u','--upLim', help='upper limit of TE - default: 100', nargs='?', const=100, type=int, default=100)
    parser.add_argument('-s','--snrLim', help='upper limit of SNR - default: 1.5', nargs='?', const=1.5, type=float,
                        default=1.5)
    parser.add_argument('-k','--snrMethod', help='Brummer ,Chang, Sijbers', nargs='?', const='Brummer', type=str,
                        default='Brummer')
    parser.add_argument('-m', '--map_raw', action='store_true', help='get the real values')
    parser.add_argument('-p', '--pv6', action='store_true', help='ParaVision 6')
    parser.add_argument('-t', '--table', action='store_true', help='save b-values and diffusion directions')
    args = parser.parse_args()

    input_folder = None
    # raw data folder
    if args.input_folder is not None:
        input_folder = args.input_folder
    if not os.path.isdir(input_folder):
        sys.exit("Error: '%s' is not an existing directory." % (input_folder,))



    list = os.listdir(input_folder)
    listOfScans = [s for s in list if s.isdigit()]

    if len(listOfScans) is 0:
        sys.exit("Error: '%s' contains no numbered scans." % (input_folder,))

    print('Start to process '+str(len(listOfScans))+' scans...')
    procno ='1'
    study=input_folder.split('/')[len(input_folder.split('/'))-1]
    print(study)
    img = []
    for expno in np.sort(listOfScans):
        path = os.path.join(input_folder, expno, 'pdata', procno)
        if not os.path.isdir(path):
            sys.exit("Error: '%s' is not an existing directory." % (path,))

        if os.path.exists(os.path.join(path,'2dseq')):

            img = Bruker2Nifti(study, expno, procno, os.path.split(input_folder)[0], input_folder, ftype='NIFTI_GZ')
            img.read_2dseq(map_raw=args.map_raw, pv6=args.pv6)
            resPath = img.save_nifti()

            if 'VisuAcqEchoTime' in img.visu_pars:

                echoTime = img.visu_pars['VisuAcqEchoTime']
                echoTime = np.fromstring(echoTime, dtype=float, sep=' ')
                if len(echoTime) > 3:
                    mapT2.getT2mapping(resPath,args.model,args.upLim,args.snrLim,args.snrMethod,echoTime)
    pathlog = os.path.dirname(os.path.dirname(resPath))
    pathlog = os.path.join(pathlog, 'data.log')

    logfile = open(pathlog, 'w')
    logfile.write(img.subject['coilname'])
    logfile.close()
