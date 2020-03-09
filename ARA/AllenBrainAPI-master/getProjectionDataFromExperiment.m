function result = getProjectionDataFromExperiment(expID)
% get projection data from Allen experiment ID
%
% function data = getProjectionDataFromExperiment(expID)
%
% Purpose
% Get projection information from ARA sample brain(s) given one or more experiment IDs.
% These can be searched for using findAllenExperiments
%
%
% Inputs
% 'expID' - a scalar or vector defining one or more experiment IDs to extract using the ARA API.
%
%
% Outputs
% data - the connectivity data in a cell array [of length equal to length(expID)] of structures.
%        The output structure as the following fields
%    hemisphere_id
%    id  .  .  .  .  .  .  .  .  . experiment ID
%    is_injection .  .  .  .  .  . the labeling here part of the injection volume?
%    max_voxel_density
%    max_voxel_x
%    max_voxel_y
%    max_voxel_z
%    normalized_projection_volume
%    projection_density .  .  .  . detected signal volume / structure volume  sum detected pixels / sum all pixels in a grid 
%    projection_intensity. .  .  . sum detected pixel intensity / sum all detected pixels
%    projection_energy. .  .  .  . sum detected pixel intensity / structure volume
%    projection_volume
%    section_data_set_id
%    structure_id .  .  .  .  .  . the ID of the brain area with which these data are associated.
%    sum_pixel_intensity
%    sum_pixels
%    sum_projection_pixel_intensity
%    sum_projection_pixels
%    volume .  .  .  .  .  .  .  . sum of detected signal volume in mm^3 (I think)
%
% 
% Rob Campbell - Basel 2015
%
%
% Also see: 
%  findAllenExperiments
%
% requires JSONlab from the


if ~exist('loadjson')
   disp('Please install JSONlab from the FEX') 
   return
end


if ~isnumeric(expID)
    error('expID should be numeric')
end


%Build the URL 
url = 'http://connectivity.brain-map.org/api/v2/data/ProjectionStructureUnionize/query.json?criteria=[section_data_set_id$eq%d]&num_rows=all';


%Get projection data (this can be slow)
for ii=1:length(expID)
    if length(expID)>1
        fprintf('%d/%d. Getting data for experiment ID %d\n',ii,length(expID),expID(ii))
    end
    try 
        page{ii} = urlread(sprintf(url,expID(ii)));
    catch
        fprintf('Failed to get data for ID %d\n', expID(ii))
        page{ii} = [];
    end
end


%parse the JSON data into a cell array of structures that the function returns
n=1;
for ii=1:length(page)
    if isempty(page{ii}), continue, end        

    tmp = loadjson(page{ii});

    if tmp.success %skip any failures 
        %Convert to an array of structures because cell arrays are annoying
        for jj=1:length(tmp.msg)
            result{n}(jj)=tmp.msg{jj};
        end
    end

    n=n+1;
    
end

