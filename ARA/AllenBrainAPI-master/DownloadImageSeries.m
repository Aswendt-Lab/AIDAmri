function DownloadImageSeries(outdir, expid, varargin)
% download Allen sample brain using the Allen API
%
% function DownloadImageSeries(outdir, expid, varargin)
%
% 
% Inputs [required]
% outdir - where to put the JPEGs
% expid  - [numerical scalar] experiment ID assigned by Allen
%
% Inputs [optional]
% 'downsample' sets how many times the image will be downsampled and scaled down, 
%   e.g. downsample=3 means the image will be 1/2^3 = 1/8 of original size. 
% 'range' specifies the range of 16 bit RGB values that will be mapped onto 8 bit
%
%
% Example
% - Pull in a nicely downsampled version of experiment ID 479701339 into the current directory:
% >> DownloadImageSeries(pwd,479701339,'downsample',4)
%
%
% PZ


params = inputParser;
params.addParamValue('range', '0,2500,0,2500,0,4095', @ischar);
params.addParamValue('downsample', 2, @isnumeric);
params.parse(varargin{:});

% download XML data for experiment
exp_url = ['http://api.brain-map.org/api/v2/data/SectionDataSet/' ...
    num2str(expid) '.xml?include=section_images'];
disp(['Accessing ',exp_url]);
doc = xmlread(exp_url);

% select the image IDs for all images
list = doc.getElementsByTagName('section-image');

nImages = list.getLength;

% cycle through images, zero-indexed
for ind = 0:nImages-1
    thisImage = list.item(ind);
    imageid = thisImage.getElementsByTagName('id');
    % number of the brain section (they aren't in order for some reason)
    sectionnum = thisImage.getElementsByTagName('section-number');
    
    query = ['http://api.brain-map.org/api/v2/section_image_download/' ...
        char(imageid.item(0).getFirstChild.getData) ...
        '?range=' params.Results.range...
        '&downsample=' num2str(params.Results.downsample)];
    disp([ num2str(ind+1) ' of ' num2str(nImages) ': ' query]);
    
    try
        urlwrite(query,...
            [ outdir filesep ...
            char(sectionnum.item(0).getFirstChild.getData) '.jpg']);
    catch
        disp(['Error getting image from ',query]);
    end
end
