function AutomizeSimsSort(pathnames)
outputs=cell(numel(pathnames),1);
parfor i=1:numel(pathnames)
    outputs{i,1}=IDGroupingT(InputFix(loadgdf([pathnames{i,1},'.gdf'])));
end
for i=1:numel(pathnames)
    output=outputs{i,1};
    save(strcmp(pathnames(i), '.mat'), 'output','-v7.3');
end