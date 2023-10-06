function [names,acronyms,ARA_LIST]=structureID2name(structIDs,ARA_LIST,quiet)
% convert a list of ARA (Allen Reference Atlas) structure IDs to a cell array of names
%
% function [names,acronyms,ARA_LIST]=structureID2name(structIDs,ARA_LIST,quiet)
%
% Purpose
% Each Allen Reference Atlas (ARA) brain area is associated with a unique
% number (structure ID). This function converts the ID to an area name.
% 
%
% Inputs
% structIDs - a vector (list) of integers corresponding to brain structure ids 
% ARA_LIST - [optional] the first output of getAllenStructureList
% quiet - [optional, false by default] if true, do not print any warning messages
%
%
% Outputs
% names - a list of brain area names in a cell array (if there is more than one name)
% acronyms - a list of brain area abbreviations in a cell array (if there is more than one)
% ARA_LIST - the CSV data from getAllenStructureList in the form of a cell array
%
%
% Examples
%
% >> structureID2name(644) 
% ans = 
%    'Somatomotor areas, Layer 6a'
%
% >> structureID2name([60,33])
% ans = 
%    'Entorhinal area, lateral part, layer 6b'    'Primary visual area, layer 6a'
%
% >> structureID2name([60,33],[],true) %for quiet operation
% ans = 
%    'Entorhinal area, lateral part, layer 6b'    'Primary visual area, layer 6a'
%
%
% Rob Campbell
%
% See also:
% getAllenStructureList, acronym2structureID


if nargin<2 || isempty(ARA_LIST)
    ARA_LIST = getAllenStructureList;
end

if nargin<3
    quiet = false;
end

%loop through and find all the IDs
names={};

for ii=1:length(structIDs)
    if structIDs(ii)==0
        if ~quiet
            names{ii}='Out of brain';
        end
        continue
    end

    f=find(ARA_LIST.id == structIDs(ii));
    if isempty(f)
        if ~quiet
            fprintf('%s finds no name for ARA ID %d\n',mfilename,structIDs(ii))
        end
        continue
    end
    if length(f)>1
        if ~quiet
            error('found more than one ID index')
        end
    end
    names{ii} = ARA_LIST.name{f};
    acronyms{ii} = ARA_LIST.acronym{f};
end

if ~isempty(names) & length(names)==1
    names=names{1};
end

