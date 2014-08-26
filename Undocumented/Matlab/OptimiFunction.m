function m=OptimiFunction(x)
r0t=x(2)/20:x(2)/20:x(2);
rho=r0t.^-x(3);
A=trapz(r0t, rho.*.4.*pi.*r0t.^2);
rho=x(1)/A * r0t.^-x(3);
t0 = (1/(x(1)*PhysConst.e))*sqrt(2*PhysConst.AMU*(x(2)^3)/x(1));
t=0:t0/10:20*t0;
output=ShockShellModel(rho, x(2), 1, 2*PhysConst.AMU, t, r0t);
if output.FusY~=0
m=x(1)/output.FusY;
else
    m=1e90;
end