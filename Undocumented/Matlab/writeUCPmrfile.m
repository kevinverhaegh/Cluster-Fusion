function [] = writeUCPmrfile(directosave,mr,i)

filename = sprintf('scan%03.0f.mr',i);
fid = fopen([directosave,filename],'wt');
fprintf(fid,'nps %d\n',mr.nps);
fprintf(fid,'dens %0.3e\n',mr.dens);
fprintf(fid,'rseed %0.3e\n',mr.rseed);
fprintf(fid,'eps %d\n',mr.eps);
fprintf(fid,'wd %d\n',mr.wd);
fprintf(fid,'acc %0.3e\n',mr.acc);
fprintf(fid,'theta %0.3e\n\n',mr.theta);
fprintf(fid,'verbose %0.3e\n',mr.verbose);
fprintf(fid,'Eioncase %0.3e\n',mr.Eioncase);
fprintf(fid,'tion %0.3e\n',mr.tion);
fprintf(fid,'tau %0.3e\n', mr.tau);
fprintf(fid,'E0 %0.3e\n',mr.E0);
fprintf(fid,'IoniTim %0.3e\n',mr.IoniTim);
fprintf(fid,'tion %0.3e\n',mr.tion);
fprintf(fid,'Ez %0.3e\n',mr.Ez);
fprintf(fid,'thalf %0.3e\n',mr.thalf);
fprintf(fid,'IonModel %0.3e\n',mr.IonModel);
fprintf(fid,'SpaceCharge %0.3e\n',mr.SpaceCharge);

fclose(fid);