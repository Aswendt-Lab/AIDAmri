function thalamus
% make a projection density plot
%
% This is the MATLAB version of the R example found at
% http://api.brain-map.org/examples/doc/thalamus/thalamus.R.html
% 
% We have some more robust MATLAB functions that encapsulate 
% some operations that were performed in-line in the R example. 
% These standalone functions are used where appropriate here. 
%
%
% Rob Campbell - Basel 2015


structures = getAllenStructureList; %read structure data as a table

data_sets = url2table('http://api.brain-map.org/api/v2/data/query.csv?criteria=model::SectionDataSet,rma::criteria,[id$in100141219,112423392,127084296,127866392,139426984,146858006,112424813],specimen%28injections%29,rma::options[tabular$eq%27distinct%20data_sets.id%20as%20section_data_set_id,injections.primary_injection_structure_id%27]');

%Here we use arrayfun instead of the R sapply function
data_sets.graph_order = arrayfun(@(x) structures.graph_order(structures.id == x), data_sets.primary_injection_structure_id);
data_sets.acronym 	  = arrayfun(@(x) structures.acronym(structures.id == x), data_sets.primary_injection_structure_id);

data_sets = sortrows(data_sets,'graph_order');


unionizes = url2table('http://api.brain-map.org/api/v2/data/query.csv?criteria=model::ProjectionStructureUnionize,rma::criteria,section_data_set[id$in100141219,112423392,127084296,127866392,139426984,146858006,112424813],structure[acronym$in''VPL'',''VPM'',''PO'',''VAL'',''PF'',''VM'',''CM'',''RH'',''MD'',''PVT'',''RE'',''AM'',''AV'',''AD'',''LD'',''LP'',''LGv'',''LGd'',''MG''],rma::include,structure,rma::options[num_rows$eqall]');



%This is much more long-winded than the R way: 
%m = xtabs(projection_volume ~ section_data_set_id + structure_id, unionizes)
setID = unique(unionizes.section_data_set_id);
strID = unique(unionizes.structure_id);
for ii=1:length(setID)
	for jj=1:length(strID)
		f=find(unionizes.section_data_set_id == setID(ii) & unionizes.structure_id == strID(jj));
		m(ii,jj) = sum(unionizes.projection_volume(f));
	end
end


%We now get the row and column names rather than use the ID numbers
row_acronyms = arrayfun(@(x) data_sets.acronym(data_sets.section_data_set_id == x), setID);
col_acronyms = arrayfun(@(x) structures.acronym(structures.id == x), strID);



%Sort the rows and columns by sorted graph_order
[~,row_order] = sort(arrayfun(@(x) data_sets.graph_order(data_sets.section_data_set_id == x), setID));
[~,col_order] = sort(arrayfun(@(x) structures.graph_order(structures.id == x), strID));

om = m(row_order,col_order);
row_acronyms = row_acronyms(row_order);
col_acronyms = col_acronyms(col_order);


%Make the plot
clf
imagesc(om)
set(gca,'YTick',1:size(om,1), 'YTickLabel',row_acronyms, 'XTick',1:size(om,2), 'XTickLabel',col_acronyms)
colormap hot
axis equal tight

xlabel('target structure')
ylabel('primary injection structure')
title('Cortico-thalamic projection')




%----------------------------------------------
function thisTable=url2table(url)
	%Read a URL into a table
	thalamusTMP=fullfile(tempdir,'thalamus.csv');
	urlwrite(url,thalamusTMP);
    thisTable=readtable(thalamusTMP);
    delete(thalamusTMP)