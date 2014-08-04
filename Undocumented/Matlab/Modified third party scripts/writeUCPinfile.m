function [] = writeUCPinfile(directosave,cut_off,snapshot,t1,t2,step1,step2,whichmodel, SpaceCharge, IonModel, Eioncase)

filename = 'UCP.in';
fid = fopen([directosave,filename],'wt');

thalf = t2/2;

fprintf(fid,'#include <stdio.h> \n');
fprintf(fid,'#include <math.h> \n\n\n');
fprintf(fid,'auv=2.187e7; \n');
fprintf(fid,'auE=5.14220652e11; \n');
fprintf(fid,'aut=2.418884326505e-17; \n');
fprintf(fid,'auEn=4.35974417e-18; \n');
fprintf(fid,'auL=5.2917720859e-11; \n\n\n');
fprintf(fid,'setparticles("electrons",nps,me,qe,nps*qe);\n');
if SpaceCharge ~= 1
fprintf(fid,'setparticles("ions",nps,2*mp,-qe,-nps*qe);\n\n\n');
end
fprintf(fid,'rb = ((3*nps)/(4*pi*dens))^(1/3);\n');
fprintf(fid,'a = (3/(4*pi*2*dens))^(1/3);\n');
fprintf(fid,'wp = sqrt((dens*qe*qe)/(me*eps0));\n');
fprintf(fid,'wm = wp/sqrt(3);\n\n\n');

fprintf(fid,'if (rseed==-1) {\n');
fprintf(fid,'\t randomize();\n');
fprintf(fid,'\t }\n');
fprintf(fid,'\t else{\n');
fprintf(fid,'\t randomize(rseed);\n');
fprintf(fid,'\t } \n\n\n');

if SpaceCharge ~= 1
fprintf(fid,'setellipse("ions",rb,rb,rb);\n');
end
fprintf(fid,'setellipse("electrons",rb,rb,rb);\n');

if IonModel == 0 || IonModel == 1
fprintf(fid,'settdist("electrons","u",tion+0.5*IoniTim,IoniTim);\n');
if SpaceCharge ~= 1
fprintf(fid,'settdist("ions","u",tion+0.5*IoniTim,IoniTim);\n');
end
end

if IonModel == 2
fprintf(fid,'settdist("electrons",“F”,“LandauIonization.gdf”,”t”,”Pt”,1,0);\n');
if SpaceCharge ~= 1
fprintf(fid,'settdist("ions",“F”,“LandauIonization.gdf”,”t”,”Pt”,1,0);\n');
end
end
    

fprintf(fid,'kb = 1.38065e-23;\n');
if ((Eioncase~=0) && (Eioncase~=1))
fprintf(fid,'\t setGBxdist("electrons","g",0,sqrt(kb*(Eioncase/e)/me)/c,3,3);\n');
fprintf(fid,'\t setGBydist("electrons","g",0,sqrt(kb*(Eioncase/e)/me)/c,3,3);\n');
fprintf(fid,'\t setGBzdist("electrons","g",0,sqrt(kb*(Eioncase/e)/me)/c,3,3);\n');
end
if (Eioncase == 1)
fprintf(fid,'\t setGBxdist("electrons","g",0,sqrt((3*(Ez/auEn)^4.5)/(128*sqrt(2.0)*(omega/aut)^2))*auv/c,3,3);\n');
fprintf(fid,'\t setGBydist("electrons","g",0,sqrt((Ez/auEn)^1.5 / (4*sqrt(2.0)))*auv/c,3,3);\n');
fprintf(fid,'\t setGBzdist("electrons","g",0,sqrt((Ez/auEn)^1.5 / (4*sqrt(2.0)))*auv/c,3,3);\n');
end

if cut_off == 1
fprintf(fid,'removeR("wcs","I",ZL);\n\n\n');
end

fprintf(fid,'R = eps*a;\n');

fprintf(fid,'setrmacrodist("electrons","u",R,0);\n');
if SpaceCharge ~= 1
fprintf(fid,'setrmacrodist("ions","u",R,0);\n\n\n');
end

fprintf(fid,'FemtoLinPlane("wcs","I",wd,tau,E0,thalf); \n');

fprintf(fid,'outputvalue("wd",wd);\n');
fprintf(fid,'outputvalue("E0",E0);\n');
fprintf(fid,'outputvalue("R",R);\n\n\n');

if SpaceCharge ~= 0
if whichmodel==0
    fprintf(fid,'spacecharge3Dtree(theta);\n');
elseif whichmodel==1
    fprintf(fid,'spacecharge3D();\n\n\n');
end
end

fprintf(fid,'accuracy(acc);\n\n\n');

fprintf(fid,'if (verbose){\n');
fprintf(fid,'\t pp("rb=",rb);\n');
fprintf(fid,'\t pp("nps=",nps);\n');
fprintf(fid,'\t pp("a=",a);\n');
fprintf(fid,'\t pp("rmacrodist:",R);\n');
fprintf(fid,'\t pp("wp=",wp);\n');
fprintf(fid,'\t pp("wd=",wd);\n');
fprintf(fid,'\t pp("E0=",E0);\n');
fprintf(fid,'\t pp("E0ioncase=",Eioncase);\n');
fprintf(fid,'\t pp("tion=",tion);\n');
fprintf(fid,'}\n\n\n');

if snapshot == 1
    fprintf(fid,'snapshot(0,%0.3e,%0.3e);\n',t1,step1);
    fprintf(fid,'snapshot(%0.3e,%0.3e,%0.3e);\n',t1,t2,step2);
    
elseif snapshot == 0
    
    fprintf(fid,'tout(0,%0.3e,%0.3e);\n',t1,step1);
    fprintf(fid,'tout(%0.3e,%0.3e,%0.3e);\n',t1,t2,step2);
    
    
end
fclose(fid);



