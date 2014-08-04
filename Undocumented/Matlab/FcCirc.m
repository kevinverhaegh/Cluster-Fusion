function output=FcCirc(R, frac)
K = - PhysConst.e^2 / (4*pi*PhysConst.epsilon0);
v = frac*PhysConst.c;
omega = frac*PhysConst.c/R;
t = 0:4*pi/(omega * 10^5):4*pi/omega;
tr = t - R/PhysConst.c;
er = [-cos(omega.*t)', - sin(omega.*t)'];

Fdir = K / R^2 .* er;

Fret = K / R^2 .* [ (v/PhysConst.c .* sin(omega.*tr) - (1 - (v/PhysConst.c)^2) .*cos(omega.*tr))', ( -(v/PhysConst.c .*cos(omega.*tr) ...
    + (1 - (v/PhysConst.c)^2) .* sin(omega.*tr)))'];

figure
plot(t, Fret(:,1), t, Fdir(:,1))
legend('Retarded force', 'Direct force')
xlabel('t')
ylabel('F')

figure
plot(t, Fret(:,2), t, Fdir(:,2))
legend('Retarded force', 'Direct force')
xlabel('t')
ylabel('Fy')

output = struct('t', t, 'Fdir', Fdir, 'Fret', Fret);