function LandauIonization(Field, E0, wd, tau, directosave, thalf,Ez)

%constants
e0=1.602e-19;
mu0=4*pi*10^-7;
auE=5.14220652e11;
aut=2.418884326505e-17;
auEn=4.35974417e-18;
auL=5.2917720859e-11;
auv=2.1876912633e6;

%Assumptions: No explosion model (e.q. field is the field delivered by the
%laser, full penetration, etc.) + Z dependence field negligable, so the
%ionization only depends on time. 

Ez=Ez/auEn;

F0 = E0/auE;

tvec = 0:0.01*(2*pi/wd):thalf;

if Field == 1
F=F0 * exp(-(1/(2*(tau*tau))).*(tvec- thalf).*(tvec - thalf)).*abs(cos(wd.*(tvec-thalf))) ;
end

dwdtLandau=0.*tvec;
kappa = sqrt(2*Ez);
for i=1:numel(tvec)
    dwdtLandau(i)= (4*(kappa^5)/F(i)) * exp(-(2*(kappa^3))/(3*F(i))) * ((0.01*(2*pi/wd))/aut);
    if sum(dwdtLandau) >= 1
        maxi=i;
        break
    end
end

fid = fopen([directosave, 'LandauIonization.txt'], 'w');
fprintf(fid,'t            Pt\n');
for i=1:maxi
fprintf(fid, '%E %E\n', tvec(i), dwdtLandau(i));
end
fclose(fid);