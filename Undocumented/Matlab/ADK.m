function w=ADK(l,m, Ez, Z, F)
nstar = Z/sqrt(2*Ez);
Cnstar = (((2*exp(1))/nstar)^nstar)*((2*pi*nstar)^(-1/2));
flm = ((2*l+1)*factorial(l+abs(m)))/((2^abs(m))*(factorial(abs(m)))*factorial(l-abs(m)));
w=(Cnstar^2)*flm*((Z^2)/(2*(nstar^2)))*sqrt((3*F*(nstar^3))/(pi*(Z^3)))*(((2*(Z^3))/(F*(nstar^3)))^(2*nstar-abs(m)-1))*exp(-(2*(Z^3))/(3*(nstar^3)*F));

