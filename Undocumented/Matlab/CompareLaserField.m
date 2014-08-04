function output = CompareLaserField(env, delta, info, t)

dum = size(t);
if dum(1) == 1
    t=t';
end

pred1 = ElectronField('GeneralTravel', 0,0,0,env,0,0,0,t-info{1,1}.param.thalf);
dz = delta*(2*pi*PhysConst.c/env.omega);
pred2 = ElectronField('GeneralTravel', 0,0,dz,env,0,0,0,t-info{1,1}.param.thalf);

output=CompareFRet(pred1, pred2, acc, t);