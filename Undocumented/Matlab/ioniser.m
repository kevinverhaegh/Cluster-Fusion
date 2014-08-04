function info=ioniser(N,n,I,tFWHM,Ez,valence)

%Ignition model Krainov, with added collisional ionization and tunnel
%ionization

%constants
mu0=4*pi*10^-7;
e0=1.602e-19;
auE=5.14220652e11;
aut=2.418884326505e-17;
auEn=4.35974417e-18;
auL=5.2917720859e-11;


%timevector
t=-2*tFWHM:tFWHM/100:2*tFWHM;

%other quantities/transform to atomic untis
F0=sqrt(I*sqrt(mu0/e0))/auE;
tau=tFWHM/(2*sqrt(log(2)));
F = F0.*exp(-(t./tau).^2);
Ez=Ez./auEn;
Zb=1:1:numel(Ez);
R0 = ((3*N/n)/(4*pi))^(1/3)/auL;
Rmax = 10*R0;
rvec = R0:0.1*R0:Rmax;
tvec = t;

omega_mie = sqrt(4*pi*n/3)*aut;

%create table alpha(R, Z)
alphatab=alphasolver(N, Zb, rvec, Ez);  
%create table ionization tims Z(t,R)
Ztab=tsolver(N,t,rvec,F, alphatab, Zb);

%save settings for differential equation solving
save('ionisercalc.mat', 'rvec', 'alphatab', 'Ztab',  'tvec', 'Zb', 'N');
%use ODE45 to solve DV
[t,R] = ode45('eqsys',[min(t), max(t)],[R0, 0]);

%calculate parameters of itnerest
Z=0.*t;
R=R(:,1);
alpha=0.*t;
Q=0.*t;
WcumTun=0.*t;
Fmax=0.*t;
Ne=0.*t;
sigma=0.*t;
WcumCol=0.*t;
for i=1:numel(t)
    [~,dimt]=min(abs(t(i)-tvec));
    [~,dimR]=min(abs(R(i)-rvec));
    %Z(t)
    Z(i)=Ztab(dimt, dimR);
    %alpha(t)
    alpha(i)=alphatab(dimR, Z(i));
    %Q(t)
    Q(i)=((N*Z(i))/4)*(2-3*cos(alpha(i)) + cos(alpha(i))^3);
    %Maximum field at B
    Fmax(i)=F(dimt)+ ((N*Z(i)) / (4*R(i)^2)).*(2-3.*cos(alpha(i)) + cos(alpha(i)).^3 ).*(3 - sin(alpha(i)./2))*(sin(alpha(i)./2)^2);
    %Calculate cummulative tunneling chance (ADK) during bethe ionization
    %(assumed tunneling over the same energy levels, can be altered to
    %achieve higher Zmax. 
    %Calculate collisional crossection (Lotz)
    sigma(i) = 2.17*valence*log((F(i)^2 / (4*(omega_mie)^2))/Ez(Z(i)))/((F(i)^2 / (4*(omega_mie)^2))*Ez(Z(i)));
    Ne(i)=(N*Z(i) - Q(i))/(4*pi*(R(i)^3)/3);
    if i~=numel(t) 
        WcumCol(i) =sum(WcumCol)+ Ne(i)*sigma(i)*sqrt(2*(F(i)^2 / (4*(omega_mie)^2)))*(abs(t(i)-t(i+1))*aut);
        WcumTun(i) = sum(WcumTun) + (abs(t(i)-t(i+1))*aut)*((exp(1)/pi)^(3/2))*(sqrt(3)*(2*Ez(Z(i)))^(9/4))/(Z(i)^(5/2))*((16*exp(1)*(Ez(Z(i))^2))/(Z(i)*Fmax(i)))^(2*Z(i) / sqrt(2*Ez(Z(i)))-3/2)*exp(-(2*((2*Ez(Z(i)))^(3/2)))/(3*Fmax(i)));
    end
end

t_esc1 = 2*sqrt(2*(R0^3)/(N*max(Z)));

disp('Estimated escape time', t_esc1);
disp('Cumulated colision ionization chance during discharge', WcumCol(numel(WcumCol)));
disp('Cumulated tunneling ionization chance', WcumTun(numel(WcumTun)));
disp('Final ionization degree', Max(Q)/(N*Max(Z)));
disp('Maximum ionization', Max(Z));

disp('Plotting result')

 figure(1)
  clf;

  subplot(2,3,1)
  plot(t,Z, 'b')
  xlabel('Time')
  ylabel('Z [inner-barrier ionization]')
  
  subplot(2,3,2)
  plot(t,Q, 'b')
  xlabel('Time')
  ylabel('Charge (a.u.) [outer-barrier ionization]')

  subplot(2,3,3)
  plot(t,R, 'b')
  xlabel('Time')
  ylabel('Radius (a.u.) [Coulomb expansion]')

  subplot(2,3,4)
  plot(t,alpha, 'b')
  xlabel('Time')
  ylabel('Alpha')

  subplot(2,3,5)
  plot(t,WcumCol, 'b', t,WcumTun, 'r')
  xlabel('Time')
  ylabel('Below-barrier ionization chance')
  legend('Collisional ionization', 'Tunneling ionization')

  subplot(2,3,6)
  plot(t, Fmax, 'b', t, F, 'r')
  xlabel('Time')
  ylabel('Electric field [a.u.]')
  legend('Boosted ignition field', 'Laser field')

info.input=struct('StartDens', n, 'Particles', N, 'tFWHM', tFWHM, 'LaserField', F, 'R0', R0, 'Ez', Ez);
info.output=struct('Q', Q, 'Z', Z, 'time', t, 'R', R, 'ElectronDens', Ne, 'MaxField', Fmax, 'alpha', alpha, 'pTunCum', WcumTun, 'pColCum', WcumCol, 'EscapeTime', t_esc1, 'FinalIon', Max(Q)/(N*Max(Z)), 'ExpectedTunnelIonization', WcumTun(numel(WcumTun))*N, 'ExpectedColIonization', Wcumcol(numel(WcumCol))*N, 'Dens', N/((4*pi*R.^3) / 3));
