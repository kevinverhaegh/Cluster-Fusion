function pout = LaserMomentaModel(x, eta0)
%function pout = LaserMomentaModelCE(env, x, eta0, px0, py0, pz0)
A = 1.9410e+12;
tau = 3.0030e-14;
omega = 2.3546e+15;
envi0 = -(A/omega) * exp(-(eta0/omega)^2 / (2 * tau * tau));
envi = -(A/omega) * exp(-(x./omega).^2 ./ (2 * tau * tau));
%pout = [PhysConst.e*(envi*sin(x) - envi0*sin(eta0)) + px0, ...
    %pz0 + (0.5/(PhysConst.me*PhysConst.c))*((PhysConst.e^2)*((envi*sin(x))^2 -2*envi*sin(x)*envi0*sin(eta0) - (envi0*sin(eta0))^2)), ...
    %py0];
    
pout = PhysConst.e.*(envi.*sin(x) - envi0*sin(eta0));