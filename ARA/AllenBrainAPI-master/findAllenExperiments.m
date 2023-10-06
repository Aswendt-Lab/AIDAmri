function varargout = findAllenExperiments(varargin)
% find all Allen experiments defined parameters
%
% function [IDs,json]=findAllenExperiments('param1','val1','param2','val2',...)
%
% 
% Inputs
% 'injection' - search for experiments with injections in this location. not searched for by default.
%               this is case senstive. So 'VISp' searches for V1 but 'visp' produces an error.
% 'line' - search for experiments on this transgenic line. not searched for by default. use '0' for wild-type, 
% 'primary' - true/false. true by default, search for injections where 'injection' was the primary
%             injection structure. if false it it will search for cases where 'injection' was not the primary
%             structure.
%
%
% Outputs
% IDs - a vector of numbers corresponding to experiment IDs
% json - cell array of structures containing all data pulled out of the JSON returned by the Allen API.
%
%
% Examples
% a) Return all experiments on wild type animals:
% findAllenExperiments('line','0'); 
%
% b) Return all experiments with injections in V1
% findAllenExperiments('injection','VISp'); 
% or:
% findAllenExperiments('injection','385'); 
%
% c) Return all experiments with injections in cortex
% findAllenExperiments('injection','Isocortex'); 
% 
% Rob Campbell - Basel 2015
%
% requires JSONlab from the FEX
%
% See Also:
% getInjectionIDfromExperiment


if ~exist('loadjson')
   disp('Please install JSONlab from the FEX') 
   return
end

%Handle input arguments
params = inputParser;
params.CaseSensitive = false;
params.addParamValue('injection','',@ischar)
params.addParamValue('line','',@ischar)
params.addParamValue('primary',true,@islogical)

params.parse(varargin{:});


%Build the URL 
%this is the base URL
url = 'http://api.brain-map.org/api/v2/data/query.json?criteria=service::mouse_connectivity_injection_structure';

%now we extend it according to what is being searched
if ~isempty(params.Results.injection)
    url = [url,'[injection_structures$eq',params.Results.injection,']'];
end

if ~isempty(params.Results.line)
    url = [url,'[transgenic_lines$eq',params.Results.line,']'];
end

if params.Results.primary
    primary='true';
else
    primary='false';
end
url = [url,'[primary_structure_only$eq',primary,']'];



%Get data from Allen
page=urlread(url);
result=loadjson(page);
if ~result.success
    fprintf('Query failed!\n%s\nAt URL: %s\n\n',result.msg,url);
end


%Return data
IDs = ones(1,length(result.msg));
for ii=1:length(IDs)
    IDs(ii) = result.msg{ii}.id;
end
fprintf('Found %d experiments\n',length(IDs))



if nargout>0
    varargout{1}=IDs;
end
if nargout>1
    varargout{2}=result.msg;
end


