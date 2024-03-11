import os
import glob
import argparse

def main(inputPath):
    SearchPath = os.path.join(inputPath, "**","anat", "*Stroke_mask.nii.gz")
    List_of_Stroke_rois = glob.glob(SearchPath, recursive=True)
    print(List_of_Stroke_rois)
    for ss in List_of_Stroke_rois:
        tempSplit = ss.split(os.sep)
        print(tempSplit)
        modality = tempSplit[-2]
        timepoint = tempSplit[-3]
        Subject = tempSplit[-4]
        IndicdencPath = glob.glob(os.path.join(os.path.dirname(ss),"*IncidenceData.nii.gz"))[0]
        TransMatInv = glob.glob(os.path.join(os.path.dirname(ss),"*MatrixInv.txt"))[0]
        OutputStrokeIncidence = os.path.join(os.path.dirname(ss),Subject + "_"+timepoint+"_"+"StrokeM_IncidenceSpace.nii.gz")
        SubjectPath = os.path.join(inputPath, Subject)
        List_of_all_tp = glob.glob(os.path.join(SubjectPath, "*"))
        
        command1 = f"reg_resample -ref {IndicdencPath} -flo {ss} -inter 0 -trans {TransMatInv} -res {OutputStrokeIncidence}"
        os.system(command1)
        
        for tp in List_of_all_tp:
            if tp != timepoint:
                anat_path_for_tp_Affine = os.path.join(tp, "anat", "*MatrixAff.txt")
                anat_path_for_tp_Bspline = os.path.join(tp, "anat", "*MatrixBspline.nii")
                anat_path_for_tp_BetFile = os.path.join(tp, "anat", "*BiasBet.nii.gz")
                
                MatrixAff = glob.glob(anat_path_for_tp_Affine)[0] 
                MatrixBspline = glob.glob(anat_path_for_tp_Bspline)[0]
                BetFile = glob.glob(anat_path_for_tp_BetFile)[0]
                OutputStroke = os.path.join(tp, "anat", Subject + "_"+ timepoint + "_" + "Stroke_mask.nii.gz")
                CheckIfStroke = glob.glob(os.path.join(tp, "anat", "*Stroke_mask.nii.gz"))
                
                if not CheckIfStroke:
                
                    command2 = f"reg_resample -ref {BetFile} -flo {OutputStrokeIncidence} -inter 0 -trans {MatrixBspline} -res {OutputStroke}"
                    os.system(command2)
                    #command = f"reg_resample -ref {BetFile} -flo {ss} -inter 0 -trans {MatrixBspline} -res {OutputStroke}"
                else:
                    continue
                
                print("done")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process stroke mask files.")
    parser.add_argument("-i", "--input", type=str, help="Input path", required=True)
    args = parser.parse_args()

    main(args.input)
