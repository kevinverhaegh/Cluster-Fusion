function sortrun(pathname, OS)
input=loadgdf([pathname, '.gdf']);
input=InputFix(input);
output=IDGroupingT(input);
save([pathname, '.mat'], 'output', '-v7.3');