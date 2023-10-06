function annotation50CHANGEDannolabelIDs = readXML_Lables(xml_file)
%% Import the data
[~, ~, annotation50CHANGEDannolabelIDs] = xlsread(xml_file);
annotation50CHANGEDannolabelIDs = annotation50CHANGEDannolabelIDs(2:end,5);

annotation50CHANGEDannolabelIDs = string(annotation50CHANGEDannolabelIDs);
annotation50CHANGEDannolabelIDs(ismissing(annotation50CHANGEDannolabelIDs)) = '';
end