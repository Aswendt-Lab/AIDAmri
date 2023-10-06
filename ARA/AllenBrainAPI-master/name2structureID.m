function [IDs,ARA_LIST]=name2structureID(names,ARA_LIST,quiet)
% Convert a list of ARA (Allen Reference Atlas) area names to a vector of structure IDs
%
% function [IDs,ARA_LIST]=name2structureID(names,ARA_LIST,quiet)
%
% Purpose
% Each Allen Reference Atlas (ARA) brain area is associated with a unique
% number (structure ID), a long name and an acronym. This function converts 
% the acronym to an area structure ID number.
% 
%
% Inputs
% names - a cell array (if a list) of brain area names or a single string. 
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
% >> name2structureID('Paraflocculus')
%
%ans =
%
%  int32
%
%   1041
%
%
%
% Rob Campbell
%
% See also:
% getAllenStructureList, structureID2name


if isstr(names)
    names={names};
end

if nargin<2 || isempty(ARA_LIST)
    ARA_LIST = getAllenStructureList;
end

if nargin<3
    quiet = false;
end

%loop through and find all the names
for ii=1:length(names)

    f=strmatch(lower(names{ii}),lower(ARA_LIST.name),'exact');
    if isempty(f)
        if ~quiet
            fprintf('%s finds no name %s in the atlas\n',mfilename, names{ii})
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


