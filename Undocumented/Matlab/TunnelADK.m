function info=TunnelADK(I,tFWHM, omega)
    
%Calculates Keldysh parameter
%For hydrogen?????? Some sources say FBSI = Ez/4 is wrong for hydrogen and
%it should be Ec = (sqrt(2) - 1) Ez^(3/2), instead of Ec=Ez^2 / 4, how
%Keldysh then ????????

%constants
e0=1.602e-19;
mu0=4*pi*10^-7;
auE=5.14220652e11;
aut=2.418884326505e-17;
auEn=4.35974417e-18;
auL=5.2917720859e-11;
auv=2.1876912633e6;

%assume deuterium
%Ez=15.4667;
Ez=13.5984;
Ez=Ez*PhysConst.e / auEn ;

tau=tFWHM/(2*sqrt(log(2)));
tt=-3*tFWHM:tFWHM/10000:3*tFWHM;
F0=sqrt((2*I)/(PhysConst.c*PhysConst.epsilon0))/auE;
Ftri=F0*(exp(-(tt./tau).^2)).*abs(cos(omega*tt));
[~,imax]=min(abs(Ftri-1));
t=-3*tFWHM:tFWHM/10000:3*tFWHM;
env.A = F0;
env.omega=omega;
env.tau = tau;
env.type='gaussian';
Fs = omega*enveta(omega.*t, env);
F = abs(omega*enveta(omega.*t, env).*cos(omega.*t));
omega = omega*aut;

%PPT parameters
kappa = sqrt(2*Ez);
Fn = Fs/(kappa^3);
nstar = 1/kappa;
Keldysh = ((kappa*omega)./Fs);

for i=1:numel(t)
    if F(i)>((Ez^2)/4)
        tBethe=t(i);
        iBethe = i;
        break
    end
end

for i=1:numel(t)
    if F(i) > (sqrt(2) - 1)*(Ez^(3/2))
        disp('Barrier suppression at:')
        disp(t(i))
        disp('s')
        tBSI=t(i);
        iBSI=i;
        break
    end
end

dwdtADK=0.*t;
dwdtADK2=0.*t;
dwdtLandau=0.*t;
dwdtPPT=0.*t;
dwdtBauer = 0.*t;
WcumADK = 0.*t;
WcumADK2=0.*t;
WcumLandau = 0.*t;
WcumPPT=0.*t;
ib=0;
for i=1:numel(t)
    if Keldysh(i) < 1
        if ib==0
            ib=i;
        end
    dwdtADK(i) = ADK(0,0, Ez, 1, F(i));
    dwdtADK2(i) = ADK2(0,0, Ez, 1, F(i));
    dwdtLandau(i)= (4*(kappa^5)/Fs(i)) * exp(-(2*(kappa^3))/(3*Fs(i)));
    dwdtBauer(i) = 2.4*(Fs(i)^2);
    dwdtPPT(i) = (kappa^2)*(((2^((2*nstar) - 2))/(nstar*factorialni(nstar)*factorialni(nstar-1))))*sqrt((3*Fn(i))/pi)*(2^(2*nstar))*(Fn(i)^(1-(2*nstar)))*exp(-(2/(3*Fn(i)))*(1-(Keldysh(i)^2)/10));
    WcumADK(i) = sum(dwdtADK)*((tFWHM/10000)/aut);
    WcumADK2(i) = sum(dwdtADK2)*((tFWHM/10000)/aut);
    WcumLandau(i) = sum(dwdtLandau)*((tFWHM/10000)/aut);
    WcumPPT(i) = sum(dwdtPPT)*((tFWHM/10000)/aut);
    else if i~=1
        WcumADK(i)=WcumADK(i-1);
        WcumADK2(i)=WcumADK2(i-1);
        WcumLandau(i)=WcumLandau(i-1);
        WcumPPT(i)=WcumPPT(i-1);
        end
    end
end

%check where F>0.05
%rADK = findfirstroot((dwdtADK(ib:1:numel(t)-1)-dwdtBauer(ib:1:numel(t)-1)));
%rLandau = findfirstroot((dwdtLandau(ib:1:numel(t)-1)-dwdtBauer(ib:1:numel(t)-1)));
%rPPT = findfirstroot((dwdtPPT(ib:1:numel(t)-1)-dwdtBauer(ib:1:numel(t)-1)));
%disp('ADK')
%disp(t(rADK+ib))
%disp('Field')
%disp(Fs((rADK+ib)))
%disp('Landau')
%disp(t(rLandau+ib))
%disp('Field')
%disp(F((rLandau+ib)))
%disp('PPT')
%disp(t(rPPT+ib))
%disp('Field')
%disp(Fs((rPPT+ib)))
%dwdtADKBauer=dwdtADK;
%dwdtLandauBauer=dwdtLandau;
%dwdtPPTBauer = dwdtPPT;
%WcumADKBauer = WcumADK;
%WcumPPTBauer = WcumPPT;
%WcumLandauBauer = WcumLandau;
%for i=1:numel(t)
%    if i>(rADK+ib)
%        dwdtADKBauer(i) = dwdtBauer(i);
%        WcumADKBauer(i)=sum(dwdtADKBauer(1:i))*((tFWHM/10000)/aut);
%    end
%    if i>(rPPT+ib)
%        dwdtPPTBauer(i) = dwdtBauer(i);
%        WcumPPTBauer(i)=sum(dwdtPPTBauer(1:i))*((tFWHM/10000)/aut);
%    end
%    if i>(rLandau+ib)
%        dwdtLandauBauer(i) = dwdtBauer(i);
%        WcumLandauBauer(i)=sum(dwdtADKBauer(1:i))*((tFWHM/10000)/aut);
%    end
%end


subplot(4,1,1)
plot(t,F(1:numel(t)),t,Fs(1:numel(t)))
%hold
axis([t(ib), max(t)*0.95, 0, 1.2*max(Fs(1:numel(t)))])
try
line([tBSI, tBSI], [0, 1.2*max(Fs(1:numel(t)))])
line([tBethe, tBethe], [0, 1.2*max(Fs(1:numel(t)))])
catch
end
subplot(4,1,2)
plot(t,(dwdtADK), t, (dwdtLandau), t, (dwdtPPT), t, dwdtADK2)
%plot(t,log(dwdtADK), t, log(dwdtLandau))
ylabel('dW/dt')
xlabel('Time')
legend('ADK', 'Landau', 'PPT','ADK2')
axis([t(ib), max(t)*0.95, 0, 1])
try
line([tBSI, tBSI], [-100, 0])
line([tBethe, tBethe], [-100, 0])
catch
end
%legend('ADK', 'Landau')
%hold
subplot(4,1,3)
plot(t, log(WcumADK), t, log(WcumLandau), t, log(WcumPPT),t,log(WcumADK2))
%plot(t, log(WcumADK), t, log(WcumLandau))
ylabel('Log tunnel ionization probability')
xlabel('Time')
legend('ADK', 'Landau', 'PPT','ADK2')
axis([t(ib), max(t)*0.95, -11, 1])
try
line([tBSI, tBSI], [-11, 0])
line([tBethe, tBethe], [-11, 0])
line([t(ib), tBethe], [-11, -11])
line([tBethe, max(t)*0.95], [0, 0])
line([t(ib), tBSI], [-11, -11])
line([tBSI, max(t)*0.95], [0, 0])
line([t(ib), max(t)*0.95], [0, 0])
catch
end
%legend('ADK', 'Landau')
%hold
subplot(4,1,4)
plot(t, (WcumADK), t, (WcumLandau), t, (WcumPPT),t,WcumADK2)
%plot(t, log(WcumADK), t, log(WcumLandau))
ylabel('Log tunnel ionization probability')
xlabel('Time')
legend('ADK', 'Landau', 'PPT','ADK2')
axis([t(ib), max(t)*0.95, 0, 1.1])
try
line([tBSI, tBSI], [0, 1])
line([tBethe, tBethe], [0, 1])
line([t(ib), tBethe], [0, 0])
line([tBethe, max(t)*0.95], [1, 1])
line([t(ib), tBSI], [0, 0])
line([tBSI, max(t)*0.95], [1, 1])
line([t(ib), max(t)*0.95], [1, 1])
catch
end
%legend('ADK', 'Landau')
%hold
try
par = -2*sqrt((3*omega)/Keldysh(iBSI)^3):2*sqrt((3*omega)/Keldysh(iBSI)^3)/1000:2*sqrt((3*omega)/Keldysh(iBSI)^3);
pen = 0.5*((0:2*sqrt((3*omega)/Keldysh(iBSI)^3)/1000:2*sqrt((3*omega)/Keldysh(iBSI)^3)).^2);
fppar=(((sqrt((4.*(Keldysh.^3))./(3*pi*omega))')* (zeros(1,numel(par))+1))).*exp(-(((Keldysh.^3)/(3*omega))')*((par.^2)));
fpperp=((sqrt((4*pi*kappa)./Fs)')*(1+zeros(1, numel(par)))).*exp(-((kappa./Fs)')*(par.^2));
fpen = ((((2.*(Keldysh.^3))./(3*omega))')*(1+zeros(1, numel(pen)))).*exp(-(((2*Keldysh.^3)/(3*omega))')*(pen));
figure
subplot(3,1,1)
plot(par, fppar(iBSI,:))
subplot(3,1,2)
plot(par,fpperp(iBSI,:))
subplot(3,1,3)
plot(pen, fpen(iBSI,:))
disp('Values distribution')
disp('Parallel momentum')
disp(sqrt((3*omega)/(Keldysh(iBSI)^3)))
disp('a.u.')
disp('Parallel velocity')
disp(sqrt((3*omega)/(Keldysh(iBSI)^3))*auv)
disp('m/s')
disp('Perpendicular momentum')
disp(sqrt(Fs(iBSI)/kappa))
disp('a.u.')
disp('Perpendicular velocity')
disp(sqrt(Fs(iBSI)/kappa)*auv)
disp('m/s')
disp('Angle')
disp(sqrt(Fs(iBSI)/kappa) / sqrt((3*omega)/(Keldysh(iBSI)^3)))
disp('Energy')
disp(Keldysh(iBSI)^3 / (3*omega))
disp('a.u.')
disp((Keldysh(iBSI)^3 / (3*omega)) * auEn / e0)
disp('eV')
disp('Quiver velocity')
disp(Fs(iBSI)/omega)
disp('a.u.')
disp(Fs(iBSI)/omega * auv)
disp('m/s')
catch
end
%info = struct('t', t, 'dwdtADK', dwdtADK, 'WcumADK', WcumADK, 'dwdtLandau', dwdtLandau, 'WcumLandau', WcumLandau, 'dwdtPPT', dwdtPPT, 'WcumPPT', WcumPPT, 'tnotadj', t, 'Keldysh', Keldysh, 'F',F, 'fppar',fppar,'fpperp',fpperp, 'par',par,'pen',pen);
info = struct('F',F','t', t', 'dwdtADK', dwdtADK, 'WcumADK', WcumADK', 'dwdtLandau', dwdtLandau, 'WcumLandau', WcumLandau');

