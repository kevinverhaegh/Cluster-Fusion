function optim = FusionObjective2L(x)

n0=1e29;
r0=1e-7;

r1D=x(1)

r=1/100:1/100:1;
R0=[r0, r0]';
Z=[1,1]';
M=[2,3]' .*PhysConst.AMU ;

r0t=r0.*r';
rho(:,1)=(2*(n0 + zeros(1,numel(r))).*(r<=r1D))';
rho(:,2)=(2*(n0 + zeros(1,numel(r))).*(r>r1D))';

N=zeros(2,1);

for i=1:2
N(i)=trapz(r0t, rho(:,i).*4.*pi.*r0t.^2);
end

densavg=sum(Z.*N)./(4/3 * pi .* R0.^3); 

t0 = sqrt(3.*M.*PhysConst.epsilon0./(densavg .* Z * PhysConst.e^2));

t=0:t0/200:3*t0;

ReacCross=[0,2;0,0];

output=ShockShellModelGeneralEuler(rho, R0, Z, M, t, r0t, ReacCross);

optim=(r0/1e-7)^3 / output.FusY;

omegaL = (2*pi*PhysConst.c)/(700e-9)
omega0 = sqrt((PhysConst.e^2 * (2*n0))/(3*PhysConst.epsilon0.*PhysConst.me))
E0 = (2*r0*PhysConst.me)/PhysConst.e * abs(omegaL^2 - omega0^2); %mimimum required electric field
I0 = 0.5 * PhysConst.epsilon0 * PhysConst.c * E0^2

dlmwrite('FusionObjective2Results',horzcat(n0,r0,output.FusY,optim),'delimiter','\t','-append')


