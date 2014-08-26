function optim = OptimCross(x)
load('DD_crossection.mat')
rho = x(1)*1e29;
R = x(2)*1e-7;
Q = 4/3 * pi * rho * PhysConst.e * R^3 ;
emax = (PhysConst.e*Q/(4*pi*PhysConst.epsilon0*(R)))/(1e3 * PhysConst.e) ; 
coulomb = (3/2 .* sqrt(EnDDsp)./(emax^(3/2))).*(EnDDsp<emax);
optimT = trapz(EnDDsp, coulomb.*DDsp.*(10^-28))/(rho^2 * R^2);
optim = 1/optimT ; 