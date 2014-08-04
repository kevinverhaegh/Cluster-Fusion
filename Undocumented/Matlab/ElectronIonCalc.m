function output = ElectronIonCalc(t, info)
%atomic unit parameters
e0=1.602e-19;
mu0=4*pi*10^-7;
auE=5.14220652e11;
aut=2.418884326505e-17;
auEn=4.35974417e-18;
auL=5.2917720859e-11;
auv=2.1876912633e6;

Ez=15.5; %assume deuterium
Z = info.Z;
N = info.N;
ionMod = info.ionMod;
env = info.env;
omega = env.omega;
thalf = info.thalf;
%ionization dynamics determine x0, y0, z0
%polarization x in
if strcmp(ionMod, 'Bethe')==1
    x0 = (2*Z /(Ez*PhysConst.e / auEn))*auL;
    y0 = 0;
    z0 = 0;
    vx0 = 0;
    vy0 = 0;
    vz0 = 0;
    eta0v = omega.*(t - z0/PhysConst.c);
    env0 = omega.*enveta(eta0v, env).*cos(eta0v);
    for i=1:numel(eta0v)
        if env0(i)> ((Ez*PhysConst.e/auEn)^2 / (4*Z))*auE
        break
        end
    end
    thelp = t(i:end);
    clear t
    t = thelp;
end

[tn, y] = ode23(@(t, y) IonElectronModel(t, y, env, vz0, vx0, omega.*(t(1) - z0/PhysConst.c), N), [min(t), max(t)], [x0; vx0; z0; vz0]);
output = struct('t', tn, 'Trajec', y);
