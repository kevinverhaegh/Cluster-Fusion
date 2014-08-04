function output =CompareHellical(env, delta, t, acc)

dum = size(t);
if dum(1) == 1
    t=t';
end

pred1 = ElectronField('GeneralTravel', 0,0,0,env,0,0,0,t);
dz = delta*(2*pi*PhysConst.c/env.omega);
pred2 = ElectronField('GeneralTravel', 0,0,dz,env,0,0,0,t);

frac = 0.7;

R = 1e-6;
omega = frac*PhysConst.c/R;

fracH = 0.7;
vH = fracH.*PhysConst.c;

pred1.x = 0.*t;
pred1.y = 0.*t;
pred1.z = vH.*t;
pred1.vz = vH.*t./t;
pred1.vy = 0.*t;
pred1.vx = 0.*t;

pred2.z = vH.*t;
pred2.y = R.*sin(omega.*t);
pred2.x = R.*cos(omega.*t);
pred2.vz = vH.*t./t;
pred2.vy = R.*omega.*cos(omega.*t);
pred2.vx = - R.*omega.*sin(omega.*t);

output = CompareRetard3(pred1, pred2, t);