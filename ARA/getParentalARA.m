%% Method generate a atlas with all parental regions of the given xlsx file
% getParentalARA('./annotation_label_IDs_valid.xlsx','./annotation/annotation.nii.gz')
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

% change large annotation value of Primary somatosensory area, unassigned   
parentalAtlasVolume(parentalAtlasVolume==182305689)= 1098;

file=dir(atlasNii_file);
fileName = strsplit(string(file.name),'.');
output_file = [fileName{1} '_parent.' fileName{2} '.' fileName{3}];
atlasData.img = parentalAtlasVolume;
save_nii(atlasData,output_file);
end

