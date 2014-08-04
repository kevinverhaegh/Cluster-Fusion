function sortrun(pathname)
input=loadgdf([pathname, '.gdf']);
input=InputFix(input);
output=IDGroupingT(input);
save([pathname, '.mat'], 'output', '-v7.3');