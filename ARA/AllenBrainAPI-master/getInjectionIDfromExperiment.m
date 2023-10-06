function [IDs,names] = getInjectionIDfromExperiment(expIDs)
% Download structure ID of the primary injection structure from an Allen experiment
% 
% function IDs = getInjectionIDfromExperiment(expIDs)
% 
% Purpose
% Make an API query that downloads the structure id of the primary injection 
% structure of each experiment in the list expIDs. 
%
% 
% Inputs
% expIDs -  a vector (list) of integers corresponding to Allen experiment IDs
%
%
% Outputs
% IDs - a list of brain area index values associated with each experiment
% names - the names associated with the IDs
%
%
% Notes
% Based on R example at http://api.brain-map.org/examples/doc/thalamus/thalamus.R.html
% 
%
% Rob Campbell - Basel 2015
%
% See Also:
% findAllenExperiments




csvQueryURl = 'http://api.brain-map.org/api/v2/data/query.csv';

% This uses a tabular query to select a few columns of interest for the results. 
% The first criteria chooses seven experiments by their ids. The specimen and 
% injection models are included in the criteria so that the subsequent tabular 
% query can access their columns. Each experiment can have multiple injections 
% but only one primary injection structure, so the 'distinct' clause limits the 
% results to one data set and primary injection structure per row.
% two percent sizes instead of one because we need to escape them for sprintf to work.
% TODO: perhaps make this more flexible or create a function that allows more elaborate API queries.
queryURL = '?criteria=model::SectionDataSet,rma::criteria,[id$in%s],specimen%%28injections%%29,rma::options[tabular$eq%%27distinct%%20data_sets.id%%20as%%20section_data_set_id,injections.primary_injection_structure_id%%27]';


%Build a comma-separated text string of experiment IDs
csvID = sprintf('%d,',expIDs);
csvID(end)=[];
finalURL = sprintf([csvQueryURl,queryURL],csvID);



result = strsplit(urlread(finalURL),'\n');

IDs = [];
n=1;

for ii = 2:length(result)
    if length(result{ii})==0
        continue
    end

    tmp=strsplit(result{ii},',');

    IDs(n) = str2num(tmp{2});
    n=n+1;
end


if nargout>1
    names = structureID2name(IDs);
end