function U = Ujellium(r, R, N)
rho = N/((4/3)*pi*(R^3));
BoolrsR = (r<R);
BoolrbR = abs(1-BoolrsR);

U = ((PhysConst.e^2 * rho)/PhysConst.epsilon0)*(((r.^2)/6 - (R^2)/2).*BoolrsR - ((R^3)./(3.*r)).*BoolrbR);