import nibabel as nii
import os
import argparse



def correct_orientation(qform,sform, img):
    # overwrite img with correct orienation
    old_img = nii.load(img)
    imgTemp = old_img.dataobj.get_unscaled()

    old_img.header.set_qform(qform)
    old_img.header.set_sform(sform)

    new_img = nii.Nifti1Image(imgTemp, None, old_img.header)
    nii.save(new_img, img)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='This script automates the adjustment of the orienation of nifti files in order to get a homogenous orienation amongst all images')
	parser.add_argument('-i', '--input', required=True, type=str)
	parser.add_argument('-t', '--template', required=True,type=str)

	args = parser.parse_args()


	if args.input is not None:
		input_img = args.input

	template_img = args.template
	data = nii.load(template_img)
	header = data.header

	sform = header.get_sform()
	qform = header.get_qform()

	correct_orientation(qform, sform, input_img)


