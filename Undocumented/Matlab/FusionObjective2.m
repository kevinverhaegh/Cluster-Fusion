function optim = FusionObjective2(x)

n0=x(1)*1e29;
r0=x(2)*1e-7;

r=1/100:1/100:1;
R0=[r0, r0]';
Z=[1,1]';
M=[2,3]' .*PhysConst.AMU ;

r0t=r0.*r';
rho=n0 + zeros(2,numel(r))';

N=zeros(2,1);

for i=1:2
N(i)=trapz(r0t, rho(:,i).*4.*pi.*r0t.^2);
end

densavg=sum(Z.*N)./(4/3 * pi .* R0.^3); 

t0 = sqrt(3.*M.*PhysConst.epsilon0./(densavg .* Z * PhysConst.e^2));

t=0:t0/200:3*t0;

ReacCross=[0,2;0,0];

omegaL = (2*pi*PhysConst.c)/(700e-9)
omega0 = sqrt((PhysConst.e^2 * n0 * sum(Z))/(3*PhysConst.epsilon0.*PhysConst.me)); %eigenfrequency oscillator
E0 = (2*R0(1)*PhysConst.me)/PhysConst.e * abs(omegaL^2 - omega0^2); %mimimum required electric field
I0 = 0.5 * PhysConst.epsilon0 * PhysConst.c * E0^2

output=ShockShellModelGeneralEuler(rho, R0, Z, M, t, r0t, ReacCross);

optim=(r0/1e-7)^3 / output.FusY;

dlmwrite('FusionObjective2Results',horzcat(n0,r0,output.FusY,optim),'delimiter','\t','-append')


