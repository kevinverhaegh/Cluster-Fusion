function [t, Keldysh]=keldysh(I,omega,tFWHM)
    
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

%assume deuterium
Ez=15.5;

t=-2*tFWHM:tFWHM/100:2*tFWHM;
tau=tFWHM/(2*sqrt(log(2)));
F0=sqrt(I*sqrt(mu0/e0))/auE;
F=F0*(exp(-(t./tau).^2));
Ez=Ez*e0/auEn;
omega=omega*aut;
%Calculate Keldysh parameter
Keldysh = omega*sqrt((2*Ez)./F);
t=(t+2*tFWHM)*1e15;
for i=1:numel(t)
    if F(i)>(Ez^2/4)
        disp('Barrier suppression at:')
        disp(t(i))
        disp('fs')
        tBSI=t(i);
        break
    end
end

subplot(2,1,1)
plot(t,F)
axis([0, max(t), min(F), 1.2*max(F)])

subplot(2,1,2)
plot(t,Keldysh)
ylabel('Keldysh')
xlabel('Time')
line([tBSI, tBSI], [0, 2])
line([0, max(t)], [1, 1])
axis([0, max(t), 0, 2])

%legend('Estimate','Upper','Lower','One')



