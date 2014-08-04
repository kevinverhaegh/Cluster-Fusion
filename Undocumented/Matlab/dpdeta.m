function deriva = dpdeta(typ, etaval, env, eta0, vz0, vx0)
gamma0 = 1/sqrt(1 - ((vz0^2 + vx0^2)/(PhysConst.c^2)));
etat = [etaval-0.2, etaval+0.2];
feta = enveta(etaval, env);
feta0 = enveta(eta0, env);
dfdeta = diff(PhysConst.c .* enveta(etat, env));
if typ == 1
    %x component
    deriva = (PhysConst.e/PhysConst.c)*(dfdeta * sin(etaval) + feta * cos(etaval));
end
if typ ==2
    %z component
    deriva = (1/(gamma0.*PhysConst.me.*(PhysConst.c - vz0) ))*((((PhysConst.e/PhysConst.c)^2) * ...
        dfdeta * sin(etaval)^2 + feta * cos(etaval) * sin(etaval) - dfdeta * feta0 * sin(etaval) ...
        * sin(eta0) - feta * feta0 * cos(etaval) * sin(etaval)) + (PhysConst.e/PhysConst.c) * ...
        gamma0 * PhysConst.me * vx0 * (dfdeta * sin(etaval) + feta * cos(etaval)));
end