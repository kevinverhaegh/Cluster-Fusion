function CalcZt(data, flatIB, flatIE, Spike)

I=1e-1; %current is 100 mA
Ierr=5e-3; %assume 5% error in the current
dt = 1/200000000;
omega=2*pi*100e3;

R=zeros(numel(flatIB),1);
Rerr=zeros(numel(flatIB),1);

for i=1:numel(flatIB)
Vf=abs(mean(data(flatIB(i):flatIE(i)))); %Voltage flat top
Vferr=std(data(flatIB(i):flatIE(i))); %Error in the voltage flat top
R(i) = Vf/I;
Rerr(i) = R(i)*sqrt((Vferr/Vf)^2 + (Ierr/I)^2);
end

Rval = sum(R./(Rerr.^2))/sum(1./(Rerr.^2));
Rvalerr = 1/sqrt(sum(1./(Rerr.^2)));

disp('Resistance is (Ohm):')
disp(Rval)
disp('With an error of (Ohm):')
disp(Rvalerr)

%Assume slope between flatIB(i) and flatIE(i+1)

M=zeros(numel(flatIB)-1,1);
Merr=zeros(numel(flatIB)-1,1);


for i=2:numel(flatIB)
    [p,s]=polyfit(dt*(flatIE(i-1):1:Spike(i)), data(flatIE(i-1):Spike(i)), 1);
    perr = sqrt(diag(inv(s.R)*inv(s.R')).*s.normr.^2./s.df);
    M(i-1)=abs(p(1))/I;
    Merr(i-1)=M(i-1)*sqrt((perr(1)/p(1))^2 + (Ierr/I)^2);
end

Mval = sum(M./(Merr.^2))/sum(1./(Merr.^2));
Mvalerr = 1/sqrt(sum(1./(Merr.^2)));

disp('Mutual inductance is (H):')
disp(Mval)
disp('With an error of (H):')
disp(Mvalerr)

disp('Transfer impedance (Ohm):')
Zt = Rval + 1i*omega*Mval;
disp(Zt)
disp('With an error of (Ohm):')
Zterr = Rvalerr + 1i*omega*Mvalerr;
disp(Zterr)
