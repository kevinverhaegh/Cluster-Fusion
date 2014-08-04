function [Ex,Ey,Ez,Bx,By,Bz] = QueMorLaser(wd,tFWHM,I,w0,x,y,z,t)

b=0.0001:1/1000:0.9999;
k0= wd/PhysConst.c;
epsilon = 1/(k0*w0);
E0 = sqrt(I*sqrt(PhysConst.mu0/PhysConst.e));

Ex = zeros(numel(x),numel(y),numel(z),numel(t));
Ey = zeros(numel(x),numel(y),numel(z),numel(t));
Ez = zeros(numel(x),numel(y),numel(z),numel(t));
Bx = zeros(numel(x),numel(y),numel(z),numel(t));
By = zeros(numel(x),numel(y),numel(z),numel(t));
Bz = zeros(numel(x),numel(y),numel(z),numel(t));

tau=tFWHM/(2*sqrt(log(2)));

Env = exp(-(t.*t)/(tau*tau));

for i=1:numel(t)
    for j=1:numel(z)
    phib = wd.*t(i) - k0.*z(j).*sqrt(1-b.*b);
        for k=1:numel(x)
            for l=1:numel(y)
            r = sqrt(x(k)*x(k) + y(l)*y(l));
            I1 = sum(exp(-((b.*b)./(4.*epsilon.*epsilon))).*(1+sqrt(1-b.*b)).*sin(phib).*b.*besselj(0,(k0.*r.*b)).*(1./1000));
            I2 = sum(exp(-((b.*b)./(4.*epsilon.*epsilon))).*(sin(phib)./sqrt(1-b.*b)).*b.*b.*besselj(1,(k0.*r.*b)).*(1./1000));
            I3 = sum(exp(-((b.*b)./(4.*epsilon.*epsilon))).*(sin(phib)./sqrt(1-b.*b)).*b.*b.*b.*besselj(0,(k0.*r.*b)).*(1./1000));
            I4 = sum(exp(-((b.*b)./(4.*epsilon.*epsilon))).*cos(phib).*(1 + (1./sqrt(1-b.*b))).*b.*b.*besselj(1,(k0.*r.*b)).*(1./1000));
            
            Ex(k,l,j,i) = Env(i)*(E0/(4*epsilon*epsilon))*(I1 + ((x(k)*x(k) - y(l)*y(l))/(k0*r*r*r))*I3 + ((y(l)*y(l))/(r*r))*I3);
            Ey(k,l,j,i) = - Env(i)*(E0/(4*epsilon*epsilon))*((x(k)*y(l))/(k0*r*r*r))*(k0*r*I3 - 2*I2);
            Ez(k,l,j,i) = Env(i)*(E0/(4*epsilon*epsilon))*(x(k)/r)*I4;
            Bx(k,l,j,i) = Env(i)*Ey(k,l,j,i)/PhysConst.c;
            By(k,l,j,i) = Env(i)*(E0/(4*PhysConst.c*epsilon*epsilon))*(I1 + ((y(l)*y(l) - x(k)*x(k))/(k0*r*r*r))*I2 + ((x(k)*x(k))/(r*r))*I3);
            Bz(k,l,j,i) = Env(i)*(E0/(4*PhysConst.c*epsilon*epsilon))*(y(l)/r)*I4;
            
            end
        end
    end
    disp((i/numel(t))*100);
end
            