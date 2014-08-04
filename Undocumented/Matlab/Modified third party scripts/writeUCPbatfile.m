function [] = writeUCPbatfile(directosave,N, IonModel)

filename = 'UCP.bat';
fid = fopen([directosave,filename],'wt');

for i = 1:N
fprintf(fid,'mr -o UCP%03.f.gdf scan%03.0f.mr gpt UCP.in\n',i,i); 
fprintf(fid,'\n\n');
end

if IonModel == 3
    fprintf(fid, 'asci2gdf –o LandauIonization.gdf LandauIonization.txt\n\n');
end

fclose(fid);