function [IDs,ARA_LIST]=acronym2structureID(acronyms,ARA_LIST,quiet)
% Convert a list of ARA (Allen Reference Atlas) acronyms to a vector of structure IDs
%
% function [IDs,ARA_LIST]=acronym2structureID(acronyms,ARA_LIST,quiet)
%
% Purpose
% Each Allen Reference Atlas (ARA) brain area is associated with a unique
% number (structure ID), a long name and an acronym. This function converts 
% the acronym to an area structure ID number.
% 
%
% Inputs
% acronyms - a cell array (if a list) of brain area acronyms or a single string. 
% ARA_LIST - [optional] the first output of getAllenStructureList
% quiet - [optional, false by default] if true, do not print any warning messages
%
%
% Outputs
% IDs - a vector of brain area structure IDs (if more than one acronym was provided)
% ARA_LIST - the CSV data from getAllenStructureList in the form of a cell array
%
%
% Examples
%
% >> acronym2structureID({'VISp','VISal'}) 
%
% ans =
%
%  1x2 int32 row vector
%
%   385   402
%
%
%
% Rob Campbell
%
% See also:
% getAllenStructureList, structureID2name


if isstr(acronyms)
    acronyms={acronyms};
end

if nargin<2 || isempty(ARA_LIST)
    ARA_LIST = getAllenStructureList;
end

if nargin<3
    quiet = false;
end

%loop through and find all the acronyms
for ii=1:length(acronyms)

    f=strmatch(lower(acronyms{ii}),lower(ARA_LIST.acronym),'exact');
    if isempty(f)
        if ~quiet
            fprintf('%s finds no acronym %s in the atlas\n',mfilename, acronyms{ii})
        end
        continue
    end
    if length(f)>1
        if ~quiet
            error('found more than one ID index')
        end
    end
    IDs(ii) = ARA_LIST.id(f);
end


