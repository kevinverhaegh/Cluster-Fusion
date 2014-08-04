function output=FcHeli(R, frac, fracH)

if (frac^2 + fracH^2 >= 1)
    disp('Velocity larger than light speed !!!!')
end

K =  PhysConst.e^2 / (4*pi*PhysConst.epsilon0);
v = frac*PhysConst.c;
omega = frac*PhysConst.c/R;
t = 0:4*pi/(omega * 10^5):4*pi/omega;
trCirc = t - R/PhysConst.c;
tr = t - R/(PhysConst.c * (sqrt(1 - fracH^2)));
er = [-cos(omega.*t)', - sin(omega.*t)'];

Fdir = K / R^2 .* er .*(1 - frac*fracH);

FretCirc = K / R^2 .* [ (v/PhysConst.c .* sin(omega.*trCirc) - (1 - (v/PhysConst.c)^2) .*cos(omega.*trCirc))', ( -(v/PhysConst.c .*cos(omega.*trCirc) ...
    + (1 - (v/PhysConst.c)^2) .* sin(omega.*trCirc)))'];
FretHeli =  K / R^2 .* [ (  frac .* sin(omega.*tr) - ((1 - (frac)^2 - fracH^2)/sqrt(1 - fracH^2))*cos(omega.*tr) )', (-frac.*cos(omega.*tr) ...
    - ((1 - (frac^2) - fracH^2)/(sqrt(1 - fracH^2)))*sin(omega.*tr))'];

figure
plot(t, FretCirc(:,1), t, Fdir(:,1), t, FretHeli(:,1))
legend('Retarded force circ', 'Direct force', 'Retarded hellical force')
xlabel('t')
ylabel('F')

figure
plot(t, FretCirc(:,2), t, Fdir(:,2), t, FretHeli(:,2))
legend('Retarded force', 'Direct force', 'Retarded Hellical force')
xlabel('t')
ylabel('Fy')

output = struct('t', t, 'Fdir', Fdir, 'FretC', FretCirc, 'FretH', FretHeli);