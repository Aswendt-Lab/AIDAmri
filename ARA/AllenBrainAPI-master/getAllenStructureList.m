function [ARA_table,tableRowInds] = getAllenStructureList(varargin)
% Download the list of adult mouse structures from the Allen API. 
%
% function ARA_table = getAllenStructureList('param1',val1,...)
%
%
% Purpose
% Make an API query to read in the Allen Reference Atlas (ARA) brain area
% list. All areas and data are read. Data are cached in the system's temp
% directory and re-read from here if possible to improve speed. 
% 
%
% Inputs (all optional param/val pair)
% 'downloadAgain' -  [false]  If true, the function wipes cached data and 
%                    re-reads. zero by default. 
% 'ancestorsOf'    -  [empty] Returns only those areas that are ancestors of
%                    of the named area. You may supply a string, numeric scalar, 
%                    or a cell array that combines these to return a table that 
%                    contains the acestors of multiple areas. If the ID or 
%                    structure name can not be found, the function returns an 
%                    empty array and displays a warning on-screen.
% 'childrenOf'     -  [empty] Returns only those areas that are children of
%                    of the named area. As above, you may supply a string, numeric 
%                    scalar, or a cell array that combines these.
% 'excludeReferenceArea' - [false] if true, the areas supplied by the childrenOf
%                           and ancestorOf arguments are removed from the 
%                           outputs.
%
%
% Outputs
% ARA_table - a table containing the imported data. 
% tableRowInds - The rows in the original table that are present in ARA_table.
%                If the user has selected a subset of data using the childrenOf
%                of or ancestorsOf arguments, then tableRowInds is:
%                ARA_table = ORIGINAL_TABLE(tableRowInds,:);
%
%
% Examples
%
%  a) Basic usage
%  S=getAllenStructureList;
%
%  b) Returning subsets of the data
%  S=getAllenStructureList('ancestorsOf',{'Posterior auditory area, layer 1',1017})
%  S=getAllenStructureList('ancestorsOf','Posterior auditory area, layer 1')
%  S=getAllenStructureList('childrenOf','Cerebellum')
%  S=getAllenStructureList('childrenOf','Cerebellum','ancestorsOf','Posterior auditory area, layer 1')
%
%  c) Remove the cerebellum
%  S=getAllenStructureList;
%  [~,ind]=getAllenStructureList('childrenOf','Cerebellum');
%  S(ind,:)=[];
%
%  d) Remove the cerebellum's children but keep the cerebellum
%  S=getAllenStructureList;
%  [~,ind]=getAllenStructureList('childrenOf','Cerebellum', 'excludeReferenceArea',true);
%  S(ind,:)=[];  
%
%
%
% Rob Campbell - Basel 2015



if nargin==1
    %Alert anyone who might be using the old scheme
    fprintf('\n\n')
    help mfilename
    error('You supplied only one input argument')
    return
end

% Parse input aruments
params = inputParser;
params.CaseSensitive = false;
params.addParamValue('downloadAgain', false, @(x) islogical(x) || x==0 || x==1);
params.addParamValue('ancestorsOf', {}, @(x) ischar(x) || isnumeric(x) || iscell(x))
params.addParamValue('childrenOf', {}, @(x) ischar(x) || isnumeric(x) || iscell(x))
params.addParamValue('excludeReferenceArea', false, @(x) islogical(x) || x==0 || x==1);
params.parse(varargin{:})

downloadAgain = params.Results.downloadAgain;

%Ensure that ancestorsOf and chilrenOf are cell arrays of IDs or names in order to simplify later code
ancestorsOf = checkFilteringInput(params.Results.ancestorsOf);
childrenOf = checkFilteringInput(params.Results.childrenOf);
excludeReferenceArea = params.Results.excludeReferenceArea;


%Cached files will be stored here
cachedCSV = fullfile(tempdir,sprintf('%s_CACHED.csv',mfilename));
cachedMAT = fullfile(tempdir,sprintf('%s_CACHED.mat',mfilename));



if ~exist(cachedMAT,'file') || downloadAgain
    % The data are to be re-read or we couldn't find any cached data

    % The adult mouse structure graph has an id of 1.       
    url='http://api.brain-map.org/api/v2/data/Structure/query.csv?criteria=[graph_id$eq1]&num_rows=all';
    [~,status] = urlwrite(url,cachedCSV);
    if ~status
        error('Failed to get CSV file from URL %s', url)
    end

    fid = fopen(cachedCSV);
    if fid<0
        error('Failed to open CSV file at %s\n', cachedCSV)
    end

    col_names = strsplit(fgetl(fid),','); %The names of the columns in the main cell array

    %Loop through and read each data row
    readParams={'%d%d%q%q%d%d%d%d%d%d%d%d%s%s%s%s%s%d%d%d%s\n','delimiter',','};
    ARA_table=textscan(fid,readParams{:});
    fclose(fid);

    ARA_table=readtable(cachedCSV,'format',readParams{:});

    %cache to disk in temporary location
    save(cachedMAT,'ARA_table')

else
    %If the data have been imported before we can just return them
    load(cachedMAT)
end





%Filter the structure list if needed
if isempty(ancestorsOf) && ~isempty(childrenOf)
    [ARA_table,tableRowInds] = returnChildrenOnly(ARA_table,childrenOf,excludeReferenceArea);
elseif ~isempty(ancestorsOf) && isempty(childrenOf)
    [ARA_table,tableRowInds] = returnAncestorsOnly(ARA_table,ancestorsOf,excludeReferenceArea);
elseif ~isempty(ancestorsOf) && ~isempty(childrenOf)
    [ARA_tableC,tableRowIndsC] = returnChildrenOnly(ARA_table,childrenOf,excludeReferenceArea);
    [ARA_tableA,tableRowIndsA] = returnAncestorsOnly(ARA_table,ancestorsOf,excludeReferenceArea);
    ARA_table = unique([ARA_tableC;ARA_tableA]);
    tableRowInds = unique([tableRowIndsC;tableRowIndsA]);
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fInput = checkFilteringInput(fInput)
    %Ensure that ancestorsOf or childrenOf are a suitable cell array 
    if iscell(fInput)
        %Do nothing
    elseif isvector(fInput) && ischar(fInput)
        fInput={fInput};
    elseif isnumeric(fInput)
        fInput=num2cell(fInput);
    else
        fprintf('\n *** %s - Unknown value of input variable \n', mfilename)
        help mfilename
        return
    end



function [returnedTable,tableRowInds] = returnAncestorsOnly(ARA_table,ancestorsOf,excludeReferenceArea)
    % If the user asked for only the ancestors of an area, we search for these here and
    % return an empty array with an on-screen warning if nothing could be found. 
    % 
    % Outputs
    % returnedTable - the table containing only the ancestors of the selected area
    % tableRowInds - these are: returnedTable = ARA_table(tableRowInds,:)

    if isempty(ancestorsOf)
        returnedTable=ARA_table;
        tableRowInds=[];
        return
    end

    childRows=[]; %Rows of the table for which we will find ancestors
    for ii=1:length(ancestorsOf)
        tChild=ancestorsOf{ii}; %This child for which we will look for parents
        if isnumeric(tChild)
            childRows=[childRows;find(ARA_table.id==tChild)];
        elseif ischar(tChild)
            childRows=[childRows;strmatch(tChild,ARA_table.name)];
        end
    end

    %Loop through childRows and collect the table rows of all ancestors
    ancestors=[];
    for ii=1:length(childRows)
        grandpa = ARA_table.structure_id_path(childRows(ii));
        grandpa = strsplit(grandpa{1},'/'); %produce a cell array of character arrays that are area index values
        grandpa = cell2mat(cellfun(@str2num,grandpa,'UniformOutput',false));
        ancestors = [ancestors, grandpa];
    end

    if excludeReferenceArea
        %Don't keep if this is the root area whose children we are looking for if the user asked to discard this
        for ii = 1:length(childRows)
            thisID = ARA_table.id(childRows(ii));
            ancestors(ancestors==thisID)=[];
        end
    end

    if isempty(childRows) || isempty(ancestors)
        fprintf('\n\n *** NO ANCESTORS FOUND. RETURNING EMPTY ARRAY ***\n\n')
        returnedTable=[];
        tableRowInds=[];
        return
    end

    tableRowInds = unique(ancestors);
    for ii=1:length(ancestors)
        tableRowInds(ii)=find(ARA_table.id==ancestors(ii));
    end

    returnedTable = ARA_table(tableRowInds,:); %filter it





function [returnedTable,tableRowInds] = returnChildrenOnly(ARA_table,childrenOf,excludeReferenceArea)
    % If the user asked for only the children of an area, we search for these here and
    % return an empty array with an on-screen warning if nothing could be found. 
    %
    % Outputs
    % returnedTable - the table containing only the children of the selected area
    % tableRowInds - these are: returnedTable = ARA_table(tableRowInds,:)

    if isempty(childrenOf)
        returnedTable=ARA_table;
        tableRowInds=[];
        return
    end

    childRows=[]; %Rows of the table for which we will find children
    for ii=1:length(childrenOf)
        tChild=childrenOf{ii}; %This child for which we will look for parents
        if isnumeric(tChild)
            childRows=[childRows;find(ARA_table.id==tChild)];
        elseif ischar(tChild)
            childRows=[childRows;strmatch(tChild,ARA_table.name)];
        end
    end

    %Get the index values associated with these rows
    childRows = unique(childRows);
    ind = zeros(size(childRows));
    for ii=1:length(childRows)
        ind(ii)=ARA_table.id(childRows(ii));
    end

    % Now we will loop through the whole table and look for rows that list each of these
    % values in their structure_id_path
    tableRowInds = [];
    for thisInd = 1:length(ind)
        for thisRow = 1:height(ARA_table)

            sID = strsplit(ARA_table.structure_id_path{thisRow},'/');
            sID = cell2mat(cellfun(@str2num,sID,'UniformOutput',false));
            f=find(sID==ind(thisInd));
            if ~isempty(f)
                %Don't keep if this is the root area whose children we are looking for if the user asked to discard this
                if excludeReferenceArea && max(f)==length(sID) %this works because it returns true when the target area is the last in the list
                    continue
                end
                tableRowInds(end+1)=thisRow;
            end

        end    
    end

    if isempty(tableRowInds)
        fprintf('\n\n *** NO CHILDREN FOUND. RETURNING EMPTY ARRAY ***\n\n')
        returnedTable=[];
        return
    end

    returnedTable = ARA_table(tableRowInds,:); %filter it

