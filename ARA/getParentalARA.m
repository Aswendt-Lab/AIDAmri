%% Method generate a Atlas with all parental Regions of the given xml_file
% getParentalARA('/Volumes/AG_Aswendt_Share/Scratch/Asw_fMRI2AllenBrain_Data/annotation_50CHANGEDanno_label_IDs_valid.xlsx','/Users/pallastn/ImageProcessing/MatLabEvaluation/Tools/ARA/annotation.nii')
function getParentalARA(xml_file,atlasNii_file)
addpath('./AllenBrainAPI-master/');
labelsStrArray = char(readXML_Lables(xml_file));
atlasData = load_nii(atlasNii_file);
parentalAtlasVolume = zeros(size(atlasData.img));
for label_idx = 1:length(labelsStrArray)
    disp(labelsStrArray(label_idx,:));
    childTable = getAllenStructureList('childrenOf',labelsStrArray(label_idx,:));
    if isempty(childTable)
        continue
    end
    childIDs = childTable.id;
   
    parentalID = name2structureID(labelsStrArray(label_idx,:));
    for child_idx = 1:length(childIDs)
        parentalAtlasVolume(atlasData.img==childIDs(child_idx)) = parentalID;
    end
    
end
% change large Values
parentalAtlasVolume(atlasData.img==182305689)= 1098;
%splitDataset


fileName = strsplit(atlasNii_file,'.');
output_file = [fileName{1} '_parent.' fileName{2}];
parentalAtlasVolume = flip(parentalAtlasVolume,3);
atlasData.img = parentalAtlasVolume;
save_nii(atlasData,output_file);
end

