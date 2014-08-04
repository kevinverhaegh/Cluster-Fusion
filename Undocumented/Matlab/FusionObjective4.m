function optim = FusionObjective4(x)

n0Xe=x(1)*1e25;
r0Xe=x(2)*1e-7;
ZXe = x(3);
n0D = x(4)*1e25;
r0D = x(5)*1e-7;
n0T = x(6)*1e25;
r0T = x(7)*1e-7;

r=1/100:1/100:1;
R0=[r0Xe+r0D+r0T, r0Xe+r0D+r0T, r0Xe+r0D+r0T]';
Z=[1,1,ZXe]';
M=[2,3,131.3]' .*PhysConst.AMU ;

r0t=(r0Xe+r0D+r0T).*r';
rho=zeros(3,numel(r))';
rho(:,3)=n0Xe.*(r0t<r0Xe);
rho(:,2)=n0T.*(r0t>r0Xe+r0D);
rho(:,1)=n0D.*((r0t>r0Xe) & (r0t<(r0Xe+r0D)));

N=zeros(3,1);

for i=1:3
N(i)=trapz(r0t, rho(:,i).*4.*pi.*r0t.^2);
end

densavg=sum(Z.*N)./(4/3 * pi .* R0.^3); 

t0 = sqrt(3.*M.*PhysConst.epsilon0./(densavg .* Z * PhysConst.e^2));

t=0:t0/20:4*t0;

ReacCross=[0,0,0;0,1,2;0,0,0];

output=ShockShellModelGeneralEuler(rho, R0, Z, M, t, r0t, ReacCross);

optim=(R0(1)/1e-7)^3 / output.FusY;

dlmwrite('FusionObjective4Results',horzcat(n0Xe,r0Xe,ZXe,n0D,r0D,n0T,r0T,output.FusY,optim),'delimiter','\t','-append')


